# covariateStats.R
# Created 28 March 2023
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load('ggspatial')
load(antfile)

blank.hist.theme <- 
  theme( panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank() )

data <- collar.df.1h %>%
  st_drop_geometry() %>% 
  mutate(TEMPGROUP = factor(plyr::round_any(TEMP_DEG_C, 10)),
         EVIGROUP = round(EVI, 1),
         EVIGROUP = factor(ifelse(EVIGROUP >= 0.7, "0.7+", EVIGROUP)))

  
# ******************************************************************************
#                           TEMPERATURE DATA BY MONTH
# ******************************************************************************
data %>% 
  mutate(MONTH=factor(MONTH, levels=c(5:12, 1:4))) %>% 
  ggplot(aes(x=TEMP_DEG_C)) +
    geom_histogram(bins=100, aes(fill=TEMPGROUP, group=TEMPGROUP)) +
    facet_wrap(~MONTH, nrow=12, strip.position='left') + 
    scale_fill_brewer(palette="YlOrRd", direction=1, name="Temp (C)") +
    scale_x_continuous(breaks=seq(0, 60, by=5)) +
    xlab('temperature (C)') + 
    theme(strip.text=element_blank())
outfile=here::here('03_output', 'eda', 'temperatureByMonth.png')
ggsave(file=outfile)

# by hour and species
ggplot(data, aes(x=hour(DATETIME), y=TEMP_DEG_C)) +
  geom_point(aes(y=TEMP_DEG_C, color=TEMPGROUP)) +
  geom_smooth(color='black') +
  facet_wrap(~SPECIES+ MONTH, nrow=2) +
  plot.theme + 
  theme(strip.text=element_blank()) +
  scale_color_brewer(palette="RdBu", direction=-1, name="Temp (C)") +
  xlab('hour') + ylab('temperature (C)')
outfile=here::here('03_output', 'eda', 'temperatureBySpecies.png')
ggsave(file=outfile)


# ******************************************************************************
#                             EVI DATA BY MONTH
# ******************************************************************************

ggplot(data, aes(x=EVI)) +
  geom_histogram(bins=100, aes(fill=EVIGROUP, group=EVIGROUP)) +
  facet_wrap(~MONTH, nrow=12) + 
  scale_fill_brewer(palette="Greens", direction=1, name="EVI") +
  scale_x_continuous(breaks=seq(-1,1,by=0.1)) +
  xlab('Enhanced Vegetation Index') + ylab('') +
  blank.hist.theme 

outfile=file.path('03_output', 'eda', 'temperatureByMonth.png')
save(file=outfile)


# ******************************************************************************
#                             MAPPING WATER SOURCES
# ******************************************************************************

# Plot water stuff for chapter 2
khauPlot(fill='#EBFBA6') +
  geom_sf(data=khau, color='#BDCE76', fill="transparent", linewidth=1.5) +
  geom_sf(data=water.nat.df, color='darkgray', size=0.2) +
  geom_sf(data=fence, color='red', linewidth=1.5) +
  geom_sf(data=rivers, color='#63a5bf', linewidth=0.75) +
  geom_sf(data=rivers %>% filter(RIV_ORD <= 6), 
          color='darkblue', linewidth=1.5) +
  geom_sf(data=water.artificial, color='black', size=6) +
  geom_sf(data=water.artificial %>% filter(CODE==29), 
          color='white', size=3) +
  coord_sf(xlim=c(20.4, 21.2), ylim=c(-19.2, -18.3)) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.0, "in"), 
                         pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)
outfile=here::here('03_output', 'covariates', 'KhauWaters.png')
ggsave(file=outfile, width=9.21, height=10)

# plotting for chapter 3
# Plot water stuff
khauPlot(fill='#EBFBA6', color="transparent") +
  geom_sf(data=khau, color='#BDCE76', fill="transparent", linewidth=1, linetype="dashed") +
  geom_sf(data=water.nat.df, color='lightblue', size=0.05) +
  geom_sf(data=rivers %>% filter(RIV_ORD <= 6), color='#63a5bf', linewidth=0.75) +
  geom_sf(data=water.artificial, color='darkblue', size=4) +
  geom_sf(data=water.artificial %>% filter(CODE==29), color='white', size=2) +
  geom_sf(data=fence, color='red', linewidth=1.5) +
  geom_sf(data=cbpp, color='red', linewidth=1.5) +
  coord_sf(xlim=c(20.48, 21.4), ylim=c(-19.2, -18.3)) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.0, "in"), 
                         pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)
outfile=here::here('03_output', 'covariates', 'KhauFences.png')
ggsave(file=outfile, width=9.21, height=10)



# ******************************************************************************
# MAPPING LANDCOVER

# change to points for plotting purposes
lands.df <- lands.raster.1 %>%
  rasterToPoints(spatial = TRUE) %>%
  as.data.frame()
names(lands.df) <- c("lc.class", "x", "y")

# find landcover classes
inx <- unlist( lapply(lands.df$lc.class, function(e) findLC(e)) )
lands.df$class <- factor( lands.meta$class[inx], levels=lands.meta$class )
lands.df$category <- factor( lands.meta$category[inx], levels=lands.meta$category )
lands.df <- lands.df %>% dplyr::select(-lc.class)