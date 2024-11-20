# createCollarDataset.R
# Created 08 Dec 2022
# Margaret Swift <mes114@duke.edu>

# Code to run EDA on antelope collar data from Piet Beytell

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(units, lwgeom, openxlsx, readxl, measurements,
               amt, lubridate)
setDataPaths('antelope')

# loading metadata from Piet
load(here::here(procpath, "Collar_Metadata.RData"))

# loading collar data from files
files <- list.files( rawpath, full.names=TRUE)
files <- files[grepl("\\/AWT", files)]
for (i in 1:length(files)) {
  print(i)
  f <- files[i]
  d <- read_excel( f, skip=2 ) %>% fixNames()
  
  # get metadata info about the animal
  id <- gsub("[A-Z]", "", unique(d$ID))
  meta.row <- collar.meta %>% filter(COLLAR_ID == id)
  d <- d %>% mutate(ANIMAL_ID=meta.row$ANIMAL_ID,
                    SEX=meta.row$SEX, 
                    SPECIES=meta.row$ANIMAL,
                    AGE=meta.row$LENGTH_AGE,
                    AREA=meta.row$AREA_DEPLOYED,
                    NOTES=meta.row$REMARKS,
                    IS_SOLITARY=meta.row$IS_SOLITARY,
                    HERD_SIZE=meta.row$HERD_SIZE,
                    DATE_START=meta.row$DATE_START,
                    DATE_END=meta.row$DATE_END,
                    ALT_M=as.numeric(ALT_M),
                    COV=as.numeric(COV))
  if (i == 1) { data <- d
  } else { data <- bind_rows(data, d)}
}


# ******************************************************************************
#                          CLEAN ANTELOPE COLLAR DATA
# ******************************************************************************
# only use data that is within the bounds (not in transit to site)
inx.1 <- data$LOCAL_DATE <= data$DATE_START
inx.2 <- data$LOCAL_DATE >= data$DATE_END
keep <- inx.1 + inx.2 == 0
data <- data[keep,]

# A few data points are missing a leading 1 on Lat for some reason
inx <- which(data$LAT > -15)
data$LAT[inx] <- data$LAT[inx] - 10

# some data is way way off; remove
inx.rm <- which(data$LAT < -19.4)
data <- data[-inx.rm,]

# And now for the big cleanup
collar.df <- data %>%
  
  # clean up data names
  rename(KPH = TRUE_SPEED_KMH) %>%
  
  # make spatial
  st_as_sf(coords=c("LON", "LAT"), remove=FALSE, crs=crs) %>%
  
  mutate(# Make date/time into proper formats
         DATE = as.Date(data$LOCAL_DATE),
         TIME = gsub("^.* ", "", as.character(data$LOCAL_TIME)),
         
         # get collar type
         COLLAR_TYPE = collar.meta$TYPE[match(gsub("SAT", "", ID), 
                                       collar.meta$COLLAR_ID)],
         
         # Get distance and bearing from last point unless it's 1st of a new ID
         FIRST = ifelse( !is.na(lag(ID)), ID != lag(ID), TRUE),
         DIST = ifelse(FIRST, NA, 
                       st_distance( geometry, lag(geometry), by_element = TRUE)),
         DIST = set_units(DIST, "meters"), 
         BEARING = ifelse(FIRST, NA, c(st_geod_azimuth(.), set_units(NA, "radians"))),
         BEARING = ifelse(BEARING < 0, BEARING + 2*pi, BEARING),
         BEARING = set_units(BEARING, "radians"),
         BEARING = set_units(BEARING, "degrees") ) %>%
  rename(COLLAR_ID = ID) %>% 
  
  # remove unimportant columns an rearrange
  dplyr::select(-c(EXT_TEMP_DEG_C, ACTIVITY, DATE_START, DATE_END,
                   GMT_DATE, GMT_TIME, LOCAL_DATE, LOCAL_TIME)) %>% 
  relocate(COLLAR_TYPE, ANIMAL_ID, SEX, SPECIES, AGE, AREA, .after=COLLAR_ID) %>%
  relocate(DIST, BEARING, .after=LON)

# FIX DATETIME
collar.df$DATETIME = lubridate::as_datetime(paste(collar.df$DATE, collar.df$TIME))
collar.df <- collar.df %>% 
  mutate(MONTH = factor(month(DATETIME)),
         HOUR = factor(hour(DATETIME))) %>% 
  relocate(DATETIME, DATE, MONTH, TIME, HOUR, .after=AREA)

# Remove duplicates
inx.rm <- which(duplicated(collar.df[,c('ANIMAL_ID', 'DATETIME')]))
collar.df <- collar.df[-inx.rm,]

# REMOVE TEMPERATURE OUTLIERS
temps <- collar.df %>% 
  st_drop_geometry() %>% 
  group_by(MONTH) %>% 
  summarize(IQR = IQR(TEMP_DEG_C),
            min5 = quantile(TEMP_DEG_C, 0.05)-IQR,
            max5 = quantile(TEMP_DEG_C, 0.95)+IQR,
            meanT = mean(TEMP_DEG_C)) %>% 
  as.data.frame()
rownames(temps) <- temps$MONTH
collar.df <- collar.df %>% 
  mutate(TOOCOLD = TEMP_DEG_C < temps[MONTH, 'min5'],
         TOOHOT  = TEMP_DEG_C > temps[MONTH, 'max5'],
         TEMP_DEG_C = ifelse(TOOHOT + TOOCOLD == 0, TEMP_DEG_C, NA)) %>% 
  dplyr::select(-TOOHOT, -TOOCOLD)

# add index for later processing
collar.df$INX <- 1:nrow(collar.df)
collar.df <- collar.df %>% relocate(INX, .before=COLLAR_ID)


# ******************************************************************************
#                                   SAVING
# ******************************************************************************

# save
save(collar.df, collar.meta, collar.df.resamp, fake.df, file=antfile)

# write shapefile for GEE stuff
# outfile = file.path(rawpath, 'antelope_collar.shp')
# st_write(collar.df, dsn=outfile)

# ******************************************************************************
#                                   EOF
# ******************************************************************************


# PLOTTING TO FIND OUTLIERS IN TEMPERATURE
# collar.df %>% 
#   st_drop_geometry() %>% 
#   mutate(MONTH = factor(month(DATETIME))) %>% 
#   filter(ANIMAL_ID %in% ids[26:30]) %>% 
#   mutate(toocold = TEMP_DEG_C < temps[MONTH, 'min5'],
#          toohot  = TEMP_DEG_C > temps[MONTH, 'max5'],
#          extremes = toohot+toocold >= 1) %>% 
#   # plotting now
#   ggplot(aes(x=DATETIME, y=TEMP_DEG_C, color=extremes)) + 
#   geom_point(alpha=0.5, size=4) + 
#   guides(color="none") +
#   facet_wrap(~ANIMAL_ID, nrow=5)
# 

