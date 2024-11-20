# geographicMaps.R
# Created 08 Dec 2021
# Author: Margaret Swift <mes114@duke.edu>

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SOURCE & SETUP

pacman::p_load(sf, raster, dplyr, spData, ggplot2, here)

# some useful lists of countries, species, and rivers of interest
countries <- c("Botswana","Namibia", "Zimbabwe",
               "South Africa", "Lesotho", "eSwatini", "Mozambique", 
               "Malawi")

africa <- world %>% 
    filter(continent == "Africa", !is.na(iso_a2))
s.africa <- africa %>% filter(name_long %in% countries)
region <- africa %>% filter(name_long %in% c('Namibia', 'Botswana'))
kaza = st_read('../01_data/kaza_boundary/kaza_boundary_2014.shp')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plotting parameters

# function to zoom to a shape with a user-provided buffer
zoomTo <- function(shape, buffer=NULL, buffer.x=NULL, buffer.y=NULL) {
  if (is.null(buffer.x)) buffer.x <- buffer
  if (is.null(buffer.y)) buffer.y <- buffer
  bb <- st_bbox(shape)
  bb['xmax'] <- bb['xmax'] + buffer.x
  bb['xmin'] <- bb['xmin'] - buffer.x
  bb['ymax'] <- bb['ymax'] + buffer.y
  bb['ymin'] <- bb['ymin'] - buffer.y
  coord_sf(xlim=bb[c('xmin', 'xmax')],
           ylim=bb[c('ymin', 'ymax')])
}

# function to add transparent background to plot
#   NOTE: figure must be saved as PNG to have transparent BG!
transparentBg <- function() {
  theme(
    axis.title.x=element_blank(),
    axis.text.x =element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y =element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background= element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )
}


# themes
my.pal = c('#5CC7B5', '#EA3365', "#ffaa00") #rev(saloon$seanchai)
kcol=my.pal[2]; bcol=my.pal[1]; zcol=my.pal[3];
my.plot.theme = theme(
    plot.background  = element_rect(fill = "transparent",colour = NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_line(color='#e6e6e6'),
    panel.border = element_blank(),
    axis.ticks   = element_blank(),
    axis.text    = element_text(color='darkgray'),
    text = element_text(size=14),
  )
blank.theme <- theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
# maps
my.base <- ggplot() + xlab('') + ylab('') + my.plot.theme
africa.map <- my.base + geom_sf(data=africa, fill='#e6e6e6')
s.africa.map <- my.base + geom_sf(data=s.africa, fill='#e6e6e6')
kaza.map <- geom_sf(data=kaza, fill=zcol, color=zcol)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#gray southern africa
ymin1=-15; ymax1=-35
xmin1=11; xmax1=41;
my.base + 
  geom_sf(data=africa, fill='black', color='black', linewidth=2) +
  geom_sf(data=s.africa, fill='#e6e6e6', color="#e6e6e6") +
  blank.theme + 
  theme(rect = element_rect(fill = "transparent"))
ggsave(here("03_output", "maps", "africa_map.png"), width=9.6, height=9.9)
  
my.base + 
  geom_sf(data=s.africa, fill='#e6e6e6', color='black', linewidth=1) +
  geom_rect(aes(xmin = 20.4, xmax = 21.4, ymin = -18.3, ymax = -19.2), color = "#15ab13", fill = '#15ab13', linewidth=2) +
  coord_sf(xlim = c(xmin1, xmax1)) + 
  theme(rect = element_rect(fill = "transparent")) +
  plot.theme
ggsave(here("03_output", "maps", "s_africa_map.png"), width=9.6, height=9.9)
  
my.base + geom_sf(data=s.africa, fill='#e6e6e6') + 
  geom_sf(data=kaza, fill="#15ab13", color="#15ab13", alpha=0.7) + 
  blank.theme + theme(rect = element_rect(fill = "transparent"))
ggsave("~/Desktop/southAfricaMap.png")


## KHAUDUM FIG 1
p = ggplot() + 
  geom_sf(data=khau, fill='#abed9f', color='#71de5d', linewidth=2) + 
  geom_sf(data=water.nat.df, color='#54dceb', size=0.2) +
  geom_sf(data=rivers %>% filter(BB_DIS_ORD<=6), color='#02a8ba', linewidth=2) +
  geom_sf(data=water.artificial, size=4, color='blue') +
  geom_sf(data=water.artificial %>% filter(CODE==29), size=2, color='white') +
  geom_sf(data=cbpp, color='darkgray', linewidth=1.7) +
  geom_sf(data=fence, color='black', linewidth=2) +
  coord_sf(xlim=c(20.5, 21.4), ylim=c(-19.2, -18.3)) +
  plot.theme
outdir = here('03_output')
ggsave(p, filename=file.path(outdir, 'khau.png'),
       width=10, height=10, bg='transparent')

## KHAUDUM SIMPLE
ggplot() + 
  geom_sf(data=khau, fill='#abed9f', color='#71de5d', linewidth=2) + 
  geom_sf(data=rivers %>% filter(BB_DIS_ORD<=6), color='blue', linewidth=2) +
  geom_sf(data=fence, color='black', linewidth=2) +
  geom_sf(data=cbpp, color='black', linewidth=1) +
  transparentBg() + 
  zoomTo(khau, buffer=0.04)
ggsave(filename=file.path(outdir, 'khau.png'),
       width=10, height=10, bg='transparent')
