# cleanMaps.R
# Created 08 Dec 2022
# Margaret Swift <mes114@duke.edu>

# Code to clean and pare down protected areas datafiles and waterpoints shapes

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(rgdal, units, lwgeom, geojsonR, geojsonsf, raster,nominatimlite)
load(antfile)
setDataPaths('khaudum_geography')
files <- list.files( rawpath )

# ******************************************************************************
#                                   FENCELINE
# ******************************************************************************

# Only the eastern border with Botswana is fenced in Khaudum
nb <- geodata::gadm("Namibia", level=0, path='tmp') %>% 
  st_as_sf()
fence <- nb %>% 
  st_as_sf() %>% 
  st_crop(c(xmin = 20.8, ymin = -18.3, xmax = 21.5, ymax = -19.5)) %>% 
  nngeo::st_segments()
fence <- fence[2:7,] %>% st_union() %>% st_cast("MULTILINESTRING")

# Ikoga and Samochima fences
cbpp <- st_read(file.path(rawpath, "botswana veterinary fences.kml")) %>% 
  filter(grepl("Cordon Fence", Name)) %>% 
  st_transform(st_crs(khau)) %>% 
  st_zm(drop = TRUE, what = "ZM")


fences <- rbind(as.data.frame(cbpp) %>% dplyr::select(-Description), 
                data.frame(Name="Border Cordon Fence", geometry=fence))
st_write()


# ******************************************************************************
#                                       WATER
# ******************************************************************************

# -----------------------------------------------------------------
# ARTIFICIAL WATERPOINTS
# -----------------------------------------------------------------
file <- files[grepl("water", files)][1]
shp <- here(rawpath, file, "All Khaudum Waters_Edit.shp")
water.artificial <- st_read(shp) %>% 
  st_as_sf(coords=c("X_COORD", "Y_COORD"), crs=crs)

# -----------------------------------------------------------------
# RIVERS
# -----------------------------------------------------------------
# file <- files[grepl("water", files)][1]
# fname <- file.path(rawpath, file, 'KhaudumRivers.geojson')
# rivers <- geojson_sf(fname, expand_geometries=TRUE)
  
# -----------------------------------------------------------------
# NATURAL WATERPOINTS
# -----------------------------------------------------------------
load(here(procpath, "waterholeStats.RData"))

# filter out unreliable waterholes (3 or fewer time steps of fill)
water.nat.df <- wh.stats %>% 
  st_as_sf(coords=c('longitude', 'latitude'), remove=FALSE, 
           crs=st_crs(khau)) %>% 
  filter(n_fill >= 3)


# ******************************************************************************
#                                   BORDERS
# ******************************************************************************

# -----------------------------------------------------------------
# KHAUDUM
file <- files[grepl("WDPA", files)][1]
base <- gsub("_0.*$", "", file)
shp <- here(rawpath, file, paste0(base, '-polygons.shp'))
khau <- st_read(shp) %>% 
  filter( NAME == "Khaudum" ) %>%
  select(NAME, ORIG_NAME, REP_AREA, GIS_AREA, STATUS_YR, geometry)
write_sf(khau, file.path(rawpath, "khau.shp"))

# -----------------------------------------------------------------
# KAZA
# -----------------------------------------------------------------
file <- files[grepl("kaza", files)][1]
shp <- here(rawpath, file, "kaza_boundary_2014.shp")
kaza <- st_read(shp)

# ******************************************************************************
#                                   LANDCOVER
# ******************************************************************************

path <- here(rawpath, files[grepl("landcover", files)])
land.files <- list.files(path)

# read LC categories metadata
meta.file <- file.path(path, land.files[grepl('csv', land.files)])
lands.meta <- read.csv(meta.file)

# read shapefile
shp.file <- file.path(path, land.files[grepl('tif', land.files)])
lc.r <- raster(shp.file)

# crop LC to a reprojected khaudum map, then reproject smaller raster to khau crs

bb1 <- bb2 <- st_bbox(collar.df) 
bb1 [['ymin']] <- bb2 [['ymax']] <- -18.59715

bb1 <- bb1 %>% 
  nominatimlite::bbox_to_poly() %>% 
  st_transform(crs(lc.r)) %>% 
  st_as_sf()

bb2 <- bb2 %>% 
  nominatimlite::bbox_to_poly() %>% 
  st_transform(crs(lc.r)) %>% 
  st_as_sf()

t0 <- proc.time()[['elapsed']]
lands.raster.1 <- crop(lc.r, bb1) %>%
  projectRaster(crs=crs(khau), method="ngb") #### THIS TAKES FOREVER FYI ! ####
lands.raster.2 <- crop(lc.r, bb2)%>%
  projectRaster(crs=crs(khau), method="ngb")
lands.raster <- stack(lands.raster.1, lands.raster.2)
t <- proc.time()[['elapsed']] - t0
message(round(t/60, 2), " minutes elapsed")

# ******************************************************************************
#                                   ELE CROSSINGS
# ******************************************************************************

setDataPaths("ele")
f.bots <- list.files(rawpath, pattern="bots.*.shp", full.names=TRUE)
f.namib <- list.files(rawpath, pattern="nam.*.shp", full.names=TRUE)

ele.bots <- st_read(f.bots) %>% 
  rename_all(.fun=c(toupper)) %>% 
  mutate(ORIGIN="BOTSWANA") %>% 
  dplyr::select(ID, LON, LAT, DATE, SEX)
ele.namib <- st_read(f.namib) %>% 
  rename_all(.fun=c(toupper)) %>% 
  mutate(ORIGIN="NAMIBIA") %>% 
  dplyr::select(ID, LON, LAT, DATE, SEX)
ele <- rbind(ele.bots, ele.namib) %>% 
  mutate(DATE=as.Date(DATE, "%d/%m/%Y"),
         SEX=toupper(SEX)) %>% 
  filter(LON < 21.2) %>% 
  st_transform(crs=st_crs(khau))


# ******************************************************************************
#                                   SAVE ALL
# ******************************************************************************

save(khau, kaza, cbpp, fence, ele,
     water.artificial, water.nat.df, rivers,
     lands.raster.1, lands.raster.2, lands.meta,
     file=here(procpath, "geographicData.RData"))

# EOF
