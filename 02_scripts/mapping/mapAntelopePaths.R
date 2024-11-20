# mapAntelopePaths.R
# Created 08 Dec 2022
# Margaret Swift <mes114@duke.edu>


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(units, lwgeom, gganimate, gifski, transformr, move, moveVis)
load(antfile)


# ******************************************************************************
#                              PLOTTING PATHS ETC
# ******************************************************************************

# ------------------------------------------------------------------------------
# Plotting a single individual colored by months since start (MONTHSSS)

plotIndiv <- function(id) {
  mss.lines <- collar.df %>%
    filter(!is.na(DATE), ANIMAL_ID == id) %>%
    mutate(MONTHYEAR = gsub("-[0-9]{2}$", "", DATE),
           MONTHSSS = interval("2013-08-01", DATE) %/% months(1) ) %>%
    group_by(ID, MONTHSSS, SEASON) %>% 
    summarise(do_union = FALSE) %>%
    st_cast("LINESTRING")
  bb <- st_bbox(mss.lines)
  p1 <- ggplot() + 
    geom_sf(data=khau, color='purple', linewidth=1, fill="transparent") +
    # geom_tile(data=lands, mapping=aes(x=x, y=y, fill=class)) +
    geom_sf(data=mss.lines, aes(color=SEASON), alpha=0.4, linewidth=1) +
    scale_color_brewer(palette="Dark2") +
    geom_tile(data=water.nat.df, mapping=aes(x=x, y=y), color='#21ebeb') +
    geom_sf(data=water.artificial, size=4, color='#eb1362') +
    coord_sf(xlim=c(bb$xmin, bb$xmax), ylim=c(bb$ymin, bb$ymax)) +
    plot.theme + ggtitle(id) + blank.theme +
    theme( panel.background = element_rect(fill = "black") )
  
  p2 <- khauPlot(color="purple") +  
    theme(
      panel.grid.major = element_line(color="#41444a"),
      panel.grid.minor = element_line(color='#41444a'),
      panel.background = element_rect(fill = "white"),
    ) +
    annotate("rect", alpha = .2, fill = "black", color='black',
             xmin = bb['xmin'], xmax = bb['xmax'], 
             ymin = bb['ymin'], ymax = bb['ymax'])
  p1 + p2 + 
    plot_layout(widths=c(4,1))
}
# 
# ids <- unique(collar.df$ANIMAL_ID)
# for (id in ids) {
#   plotIndiv(id)
#   fname = file.path(outpath, 'individual_maps_seasons', paste0(id, '.png'))
#   ggsave(filename = fname)
# }



# ------------------------------------------------------------------------------
# Plotting each individual on the Khaudum NP map
id.lines <- collar.df %>%
  mutate( TYPE = paste(SPECIES, SEX, sep="_")) %>%
  group_by(ANIMAL_ID, TYPE) %>%
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")

ggplot() + 
  geom_sf(data=khau, fill='#abed9f', color='#71de5d', linewidth=2) + 
  geom_sf(data=id.lines, aes(color=ANIMAL_ID, group=ANIMAL_ID),
          alpha=0.5, linewidth=1) + 
  facet_wrap(~TYPE, nrow=1) +
  guides(color="none")
  # scale_color_brewer(palette="Set3")
ggsave(filename=here::here(outpath, 'antelope_paths_mf_sp.png'))

# ------------------------------------------------------------------------------
# Creating an animated plot of an individual's movement

# first, static plot with waterholes
my.ids <- c("Og_1721")
indiv <- collar.df %>% filter(ANIMAL_ID %in% my.ids)
bb <- st_bbox(indiv)
ggplot() +
  geom_sf(data=indiv, size=2, mapping=aes(color=ANIMAL_ID), alpha=0.5) +
  geom_sf(data=water.artificial, size=5, color="black") +
  geom_tile(data=water.natural, mapping=aes(x=x, y=y), color="black") +
  coord_sf(xlim=c(bb$xmin, bb$xmax), ylim=c(bb$ymin, bb$ymax))

# convert to Move data type
# WHY IS POSIXCT NOT WORKING AHHHHHHH
indiv$DATETIME <- paste(indiv$DATE, indiv$TIME)
indiv$DATETIME <- parse_date_time(indiv$DATETIME, orders="%Y-%m-%d %H:%M:%S")

indiv.mv <- df2move(indiv, crs, "LON", "LAT", "DATETIME", "ANIMAL_ID")

# align to a uniform time scale
m <- align_move(indiv.mv, res = 5, unit = "hours")

# create spatial frames with a OpenStreetMap watercolour map
frames <- frames_spatial(m, path_colours = c("#eb1362"),
                         map_service = "carto", map_type = "light", alpha = 0.5) %>% 
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_northarrow() %>% 
  add_scalebar() %>% 
  add_timestamps(type = "label") %>% 
  add_progress()

frames[[1000]] # preview one of the frames, e.g. the 100th frame

# animate frames
filename <- here::here(outpath, "gifs", "SAT975_SAT977_full.gif")
animate_frames(frames, out_file = filename)
# ------------------------------------------------------------------------------




