# mapWaterholes.R
# Created 22 July 2022
# Author: Margaret Swift <mes114@duke.edu>

# Script reads the waterhole data sent by Piet Beytell & does some EDA with these.
# Waterholes are located in Khaudum national park, Namibia.

################################################################################
### SETUP ###
################################################################################
source(here::here("02_scripts", "utilities.R"))

# Load libraries
pacman::p_load(raster, tmap, leaflet, cartogram, nationalparkcolors)

# Load datasets
load('./zz_trash/mapsData.RData')

################################################################################
### FUNCTIONS ###
################################################################################


################################################################################
### PLOT PARAMS ###
################################################################################

# themes
my.pal = park_palette("Voyageurs")
kcol=my.pal[2]; bcol=my.pal[2]; zcol=my.pal[4];
my.plot.theme = theme(
  plot.background  = element_rect(fill = "transparent",colour = NA),
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.major = element_line(color='#e6e6e6'),
  panel.border = element_blank(),
  axis.ticks   = element_blank(),
  axis.text    = element_text(color='darkgray'),
  text = element_text(size=14),
  title = element_text(size=20),
)

# maps
my.base <- ggplot() + xlab('') + ylab('') + my.plot.theme
africa.map <- my.base + geom_sf(data=africa, fill='#e6e6e6')
khau.map  <- geom_sf(data=khau,  fill=kcol)
bwab.map <- geom_sf(data=bwab, fill=bcol)
kaza.map <- geom_sf(data=kaza, fill=zcol, color=zcol)

# Get shapefile for only countries within KAZA
kaza_countries_list <- c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
countries <- africa %>% 
  filter(name_long %in% kaza_countries_list)


################################################################################
### PLOT ###
################################################################################

country.lab <- data.frame(
  x = c( 19.4,       17,        27.5,     24.8,       27.55), 
  y = c(-15,        -22,       -15,      -22,        -18.8),
  val = c("ANGOLA", "NAMIBIA", "ZAMBIA", "BOTSWANA", "MOZAMBIQUE")
) %>% st_as_sf(coords=c('x', 'y')) %>% st_set_crs(4326)

kaza.map.final <- africa.map + 
  # KAZA
  geom_sf(data=kaza, fill=zcol, color=alpha(zcol, 0.4), alpha=0.8) +
  geom_sf_text(data=kaza, label="KAVANGO-ZAMBEZI\nTRANSFRONTIER AREA", 
               hjust=0.5, vjust=-1, fontface="italic", color='white') +

  # Country labels
  geom_sf_text(data=country.lab, aes(label=val),
               color="gray55", size=20/.pt, vjust=0, hjust=0) +
  coord_sf(xlim=c(17, 30), ylim=c(-12, -22))
ggsave('../03_output/kaza_map.png')

# add Khaudum NP
kaza.map.final + khau.map + coord_sf(xlim=c(17, 30), ylim=c(-12, -22)) +
  geom_sf_label(data=khau, label="Khaudum NP", hjust=1.1, vjust=0.6, size=18/.pt)
ggsave('../03_output/khaudum_map.png')

my.base + khau.map + 
  geom_sf(data=waterpts, color='black', size=5) +
  geom_sf(data=waterpts, mapping=aes(color=IDENTITY), size=4) +
  geom_sf_label(data=waterpts, mapping=aes(label=NAME), vjust=1.3, hjust=-0.05) +
  scale_color_manual('', values=c('black', 'white')) +
  ggtitle('KHAUDUM WATERPOINTS')
ggsave('../03_output/waterpoints.png')
  

my.base + 
  geom_sf(data=bwab, fill='transparent', color='red') +
  theme(
    plot.background  = element_rect(fill = "transparent",colour = NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_line(color='transparent'),
    panel.border = element_blank(),
    axis.ticks   = element_blank(),
    axis.text    = element_blank(),
    text = element_text(size=14),
    title = element_text(size=20),
  )
my.base + 
  geom_sf(data=khau, fill='transparent', color='red') +
  theme(
    plot.background  = element_rect(fill = "transparent",colour = NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_line(color='transparent'),
    panel.border = element_blank(),
    axis.ticks   = element_blank(),
    axis.text    = element_blank(),
    text = element_text(size=14),
    title = element_text(size=20),
  )
ggsave('../03_output/khau_outline.png')
