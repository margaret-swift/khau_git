# mapPlots.R
# Created 5 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

# tell R where you are & load utilities
# here::i_am('02_scripts/2_eda/mapPlots.R')
# source(here::here('02_scripts','utilities.R')) #this is for other custom funcs

# load the "world" dataset from spData
library(spData)

# pull geographic data (KAZA boundary, fences, roads, etc.)
setDataPaths('khaudum_geography')
load(here(procpath, 'geographicData.RData'))

# where to save
outdir = "./03_output/maps"


# ******************************************************************************
#                                 FUNCTIONS & PARAMETERS
# ******************************************************************************

# function to zoom to a shape with a user-provided buffer
zoomTo <- function(shape, buffer=NULL, buffer.x=NULL, buffer.y=NULL) {
  if (is.null(buffer.x)) buffer.x <- buffer
  if (is.null(buffer.y)) buffer.y <- buffer
  bb <- st_bbox(shape)
  bb[c('xmin', 'ymin')] <- bb[c('xmin', 'ymin')] - buffer.x
  bb[c('xmax', 'ymax')] <- bb[c('xmax', 'ymax')] + buffer.y
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

## GEOGRAPHIC OBJECTS
sadc_list <- c("south africa", "lesotho", "eswatini", "angola",
               'mozambique', 'malawi', 'zambia', 'zimbabwe', 
               'botswana', 'namibia', 'tanzania', 'mauritius', 
               'democratic republic of the congo' , 'madagascar')
kaza_members <- c('angola', 'zambia', 'zimbabwe', 'botswana', 'namibia')
africa <- world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  st_transform(crs=st_crs(fences))
sadc <- africa %>% filter(tolower(name_long) %in% sadc_list)
sadc_kaza <- africa %>% filter(tolower(name_long) %in% kaza_members)
bots <- africa %>% filter(name_long=="Botswana")

# ******************************************************************************
#                                    PLOTTING
# ******************************************************************************

## KHAUDUM
ggplot() + 
  geom_sf(data=khau, fill='#abed9f', color='#71de5d', linewidth=2) + 
  geom_sf(data=rivers %>% filter(BB_DIS_ORD<=6), color='blue', linewidth=2) +
  transparentBg() + 
  zoomTo(khau, buffer=0)
ggsave(filename=file.path(outdir, 'khau.png'),
       width=10, height=10, bg='transparent')
