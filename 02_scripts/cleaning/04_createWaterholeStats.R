# cleanWaterholeStats.R
# Created 23 Feb 2023
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am(file.path('02_scripts', 'cleaning', '04_createWaterholeStats.R'))
source(here::here("02_scripts", "utilities.R"))
pacman::p_load(rgdal, units, lwgeom, googledrive)
drivepath="GEEData/khauStatsOutput"

# ******************************************************************************
#                                 STATIC STATS
# ******************************************************************************

stat.file <- drive_ls( path=drivepath, pattern="staticStatsKHAU.csv" )
wh.stats <- stat.file %>% 
  drive_read_string(encoding="UTF-8") %>%
  read.csv(text = .) %>% 
  dplyr::select(id, latitude_mean, longitude_mean, elevation_m_mean, 
                is_water_count) %>%
  rename_with(.fn=function(e) {gsub('_mean', '', e)}) %>%
  rename(area_px = is_water_count) %>%
  mutate(id = as.character(id),
         n_fill_wet = 0,
         n_fill_dry = 0,
         p_fill_wet = 0,
         p_fill_dry = 0,
         )

# ******************************************************************************
#                                 DYNAMIC STATS
# ******************************************************************************

files <- drive_ls( path=drivepath, pattern="dynamicStats.*.csv" )
LD <- nrow(files)
wetcount = 0
for (i in 1:LD) {
  message('file ', i, ' of ', LD)
  fname <- files[i,]
  df <- fname %>% 
    drive_read_string(encoding="UTF-8") %>%
    read.csv(text = .) %>% 
    mutate(id = as.character(labels_median),
           is_full = is_water_dyn_sum > 0,
           MONTH = month(date_begin),
           SEASON = ifelse(MONTH %in% 5:10, "dry", "wet")) %>%
    filter(id %in% wh.stats$id)
  inx <- match(df$id, wh.stats$id)
  
  # get season
  szn <- unique(df$SEASON)
  if (length(szn) == 2) warning('more than one season for ', i)
  else {
    if (szn == "wet") wetcount <- wetcount + 1
    cols <- which(grepl(szn,names(wh.stats)))
    
    # number of times there is some water
    wh.stats[inx, cols[1]] <- wh.stats[inx, cols[1]] + df$is_full
    
    # percent full
    p_full <- df$is_water_dyn_sum / wh.stats$area_px[inx]
    wh.stats[inx, cols[2]] <- wh.stats[inx, cols[1]] + p_full
  }
}

# Extra processing post-drive
wh.stats <- wh.stats %>% 
  mutate(
         # dry
         time_dry = n_fill_dry / (LD-wetcount),
         p_fill_dry = p_fill_dry / (LD-wetcount),
         px_dry = (p_fill_dry!=0) * area_px,
         rly_dry = p_fill_dry * px_dry,
         
         # wet
         time_wet = n_fill_wet / wetcount,
         p_fill_wet = p_fill_wet / wetcount,
         px_wet = (p_fill_wet!=0) * area_px,
         rly_wet = p_fill_wet * px_wet,
         
         # both
         keep = (n_fill_wet + n_fill_dry) >= 6
         ) %>%
  filter(keep) 

# save stats data
outfile = here(procpath, "waterholeStats.RData")
save(wh.stats, file=outfile)

# write to shapefile
wh.sf <- wh.stats %>% 
  st_as_sf(coords=c('longitude', 'latitude'), crs=st_crs(collar.df)) %>% 
  dplyr::select(id, px_dry, px_wet)
outfile = here(procpath, "waterholeStatsSF.shp")
st_write(obj=wh.sf, dsn=outfile)

# EOF