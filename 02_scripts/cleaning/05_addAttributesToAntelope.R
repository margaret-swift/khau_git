# addAttributesToAntelope.R
# Created 21 Feb 2023
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(lubridate, suncalc, units, exactextractr, googledrive) 
load(antfile)
setDataPaths('khaudum_geography', verbose=TRUE)
load(here(procpath, "geographicData.RData"))

# ******************************************************************************
#                            IS POINT WITHIN KHAUDUM?
# ******************************************************************************
pts <- st_filter( collar.df, khau ) %>% st_drop_geometry() 
collar.df$IN_BORDER <- collar.df$INX %in% pts$INX


# ******************************************************************************
#                            ADDING TIME OF DAY (FANCY)
# ******************************************************************************
calcTOD <- function(t, times) {
  # find the time of day compared to the calculated solar times
  # this isn't as efficient as it could be but... eh
  times <- times %>% dplyr::select( dawn, sunriseEnd, sunset, night )
  time.names <- c('DAWN', 'DAY', 'DUSK', "NIGHT")
  
  diffs <- unlist( lapply( times, function(e) difftime(e, t, units="mins") ) )
  position <- which( diffs < 0 )
  
  if (!length(position)) { inx <- 4 # if none are <0 it's before dawn!
  } else { inx <- max( position ) }
  TOD <- time.names[inx]
  return( TOD )
}

### THIS TAKES ABOUT 20 MINUTES~~~~~~~~~~~
data <- collar.df %>% rename(lat=LAT, lon=LON, date=DATE)
sun.times <- getSunlightTimes(data=data[1:10,])
TOD <- vector(mode="character", length=nrow(collar.df))
t0 <- proc.time()[['elapsed']]
for (i in 1:nrow(collar.df)) { # hate to do it in a loop but... sigh...
  if ((i %% 1000) == 0) {
    message( round(100 *i/nrow(collar.df)), "% complete")
    message("time elapsed: ", round((proc.time()[['elapsed']]-t0)/60), " minutes")
  }
  times <- sun.times[i,]
  t <- collar.df$DATETIME[i]
  TOD[i] <- calcTOD(t, times)
}

# treated a slightly different way
collar.df$TOD <- TOD
collar.df <- collar.df %>% 
  mutate(
    TOD_ACTIV = ifelse(TOD %in% c("DAWN", "DUSK"), "CREPUSCULAR", 
                    ifelse(TOD == "NIGHT", "NOCTURNAL", "DIURNAL")))

# ******************************************************************************
#                                  ADDING SEASON
# ******************************************************************************

dry.season = yday("01-05-2022"):yday("01-10-2022")
collar.df <- collar.df %>%
  mutate(SEASON = ifelse(yday(DATETIME) %in% dry.season, "DRY", "WET")) %>%
  arrange(ANIMAL_ID, DATETIME)


# ******************************************************************************
#                                 ADDING LANDCOVER
# ******************************************************************************

# Putting everything together form split LC's
LC.df <- data.frame(INX=collar.df$INX, 
                    lc1=raster::extract(lands.raster.1, collar.df),
                    lc2=raster::extract(lands.raster.2, collar.df)) %>% 
  mutate(LC_CLASS = ifelse(is.na(lc1), lc2, lc1))
LC.df$LC_CATEG<- lands.meta$class[ unlist( lapply(LC.df$LC, function(e) findLC(e)) ) ]
collar.df$LC_CLASS <- LC.df$LC_CLASS
collar.df$LC_CATEG <- factor( LC.df$LC_CATEG, levels=lands.meta$class )


# ******************************************************************************
#                                 ADDING TEMPERATURES
# ******************************************************************************
cdf <- collar.df.1h
setDataPaths('temperature')
load(here(procpath, 'era5estimates.rds'))
cols <- c('ERA5_TEMP_K', 'ERA5_TEMP_C', 'DAYTYPE')
inx <- match(cdf$INX, collar.df.1h$INX)
cdf[,cols] <- collar.df.1h[inx,cols] %>% nog()
cdf <- cdf %>% relocate(ERA5_TEMP_K, ERA5_TEMP_C, TEMPGROUP, DAYTYPE, .after=TEMP_DEG_C) %>% 
  dplyr::select(-ERA5_TEMP_INX)
collar.df.1h <- cdf
# mean <- mean(day.temps$MAXTEMP)
# sd <- sd(day.temps$MAXTEMP)
# day.types <- day.temps %>%
#   mutate(DAYTYPE = factor( ifelse(MAXTEMP < (mean-sd), "MILD",
#                                   ifelse(MAXTEMP < (mean+sd), "AVG", "HOT")),
#                            levels=c('MILD', 'AVG', 'HOT')))
# data <- data %>% 
#   mutate(DAYTYPE = day.types$DAYTYPE[ match(DATE, day.types$DATE) ])



# ******************************************************************************
#                       PULL EVI DATA FROM DRIVE AND MERGE
# ******************************************************************************

files <- drive_ls(path="GEEData/antelope_covariates", pattern=".csv$")
LF <- nrow(files)
drive.data <- NULL
for (i in 1:LF) {
  message(i, ' of ', LF)
  df <- files[i,] %>% 
    drive_read_string(encoding="UTF-8") %>%
    read.csv(text = .) %>% 
    rename(WATERDIST = waterDist)
  if (is.null(drive.data)) drive.data <- df
  else drive.data <- rbind(drive.data, df)
}

inx <- match(collar.df$INX, drive.data$INX)
collar.df$EVI <-  drive.data$EVI[inx]
collar.df$DIST_NAT <-  drive.data$WATERDIST[inx]
collar.df <- collar.df %>% arrange(INX)

# ******************************************************************************
#                     ADDING WATERHOLE AND BORDER DISTANCES
# ******************************************************************************

# DIST FROM NEAREST ARTIFICIAL WATERHOLE
water.artificial <- st_transform(water.artificial, st_crs(collar.df))
dists <- st_distance(collar.df, water.artificial)
mindists.art <- apply(dists, 1, min)

# DISTANCE FROM NEAREST RIVER
# TIMING: 2.5min
# dists.riv <- st_distance(collar.df, rivers)
# mindists.riv <- apply(dists.riv, 1, min)

# DISTANCE FROM FENCE
# dists.fence <- st_distance(collar.df, fence %>% st_cast("MULTILINESTRING"))
# mindists.fence <- apply(dists.fence, 1, min)

# SAVE TO DF
collar.df <- collar.df %>%
  mutate(
         # DIST_ART  = mindists.art,
         # DIST_NAT  = mindists.nat,
         # DIST_RIV  = mindists.riv,
         # DIST_WATER = pmin(DIST_ART, DIST_NAT, DIST_RIV),
         # DIST_FENCE = mindists.fence
  )


# ******************************************************************************
#                           REORDERING AND DROPPING COLS
# ******************************************************************************

# SAVE TO DF
collar.df <- collar.df %>%
  dplyr::select(-WATERDIST, -EVI_300m
                #-KPH, 
                # -DIR, 
                # -COV, 
                # -HDOP, 
                # -DISTANCE_M, 
                # -COUNT, 
                # -DIST_BORD
                ) %>%
  relocate(TOD, TOD_ACTIV, SEASON, .after=DATETIME) %>%
  relocate(ALT_M,
           DIST, 
           BEARING,
           AREA,
           LC_CLASS, 
           LC_CATEG,
           IN_BORDER,
           EVI,
           DIST_ART,
           DIST_NAT,
           DIST_RIV,
           DIST_WATER,
           DIST_FENCE,
           .after=LON )


# ******************************************************************************
#                                    FAKE DATASET
# ******************************************************************************

# create fake dataset
set.seed(999)
fake.df <- st_sample(khau, size=10000) %>% 
  st_as_sf() %>% st_transform(st_crs(collar.df))
fake.df <- cbind(fake.df, st_coordinates(fake.df)) %>% 
  rename(LAT=Y, LON=X)

# landcover class
fake.df$LC <- raster::extract(lands.raster, fake.df)
inx <- unlist( lapply(fake.df$LC, function(e) findLC(e)) )
fake.df$LC_CATEG <- factor( lands.meta$class[inx], levels=lands.meta$class )
fake.df <- fake.df %>% rename(LC_CLASS = LC)

# distances from water
water.artificial <- st_transform(water.artificial, st_crs(collar.df))
fake.df$DIST_ART <- apply(st_distance(fake.df, water.artificial), 1, min)
fake.df$DIST_FENCE <- apply(st_distance(fake.df, 
                                       khau %>% st_cast("MULTILINESTRING")), 1, min)
fake.df$DIST_RIV <- apply(st_distance(fake.df, rivers), 1, min)
fake.df$DIST_NAT <- apply(fake.df, 1, findDists)
fake.df$DIST_WATER = with(fake.df, pmin(DIST_ART, DIST_NAT, DIST_RIV))

# reorder
fake.df <- fake.df %>% 
  relocate(DIST_ART, DIST_NAT, DIST_RIV, DIST_WATER, DIST_FENCE, .after=LC_CATEG)


# ******************************************************************************
#                                     SAVING
# ******************************************************************************

setDataPaths('antelope')
collar.df.1h <- collar.df %>% filter(COLLAR_TYPE == "SAT+")
collar.df.5h <- collar.df %>% filter(COLLAR_TYPE == "SAT")
save(collar.df, collar.meta, #collar.df.resamp, 
     collar.df.1h, collar.df.5h,
     fake.df, file=antfile)


# need to save individual pieces as well!!!!!!
# outfile <- here::here(procpath, "collar_pieces.RData")
# save(
#   mindists.art, 
#   mindists.nat, 
#   mindists.riv, 
#   mindists.bord,
#   TOD, file=outfile
#   )
# 

