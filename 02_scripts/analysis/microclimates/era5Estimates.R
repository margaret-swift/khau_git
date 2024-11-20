
# ******************************************************************************
#             TEMPERATURE DATA -- LINKING ERA5 DATA WITH POINTS
# ******************************************************************************

setDataPaths('temperature')
library(raster)
fname <- list.files(rawpath, pattern='.grib', full.names=TRUE)
grib <- brick(fname)
times <- seq(ymd_hm("2015-01-01 00:00"), ymd_hm("2017-12-31 23:00"), by="hour")


# First, get the indices of raster layers that correspond with points
# t0 <- proc.time()[['elapsed']]
# for (i in 1:nrow(collar.df.1h)) {
#   if (i %% 1000 == 0) print(i)
#   pt <- collar.df.1h[i,]
#   inx <- which.min(abs(times - pt$DATETIME))
#   collar.df.1h$ERA5_TEMP_INX[i] <- inx
# }
# t1 <- proc.time()[['elapsed']]
# t1 - t0


# now that I've got the indices i can read each raster only once to save time
# indices <- unique(collar.df.1h$ERA5_TEMP_INX) #that's a lot fewer calls to extract()!
inx = which(is.na(collar.df.1h$ERA5_TEMP_K))
indices = sort( unique( collar.df.1h$ERA5_TEMP_INX[inx] ) )

t0 <- proc.time()[['elapsed']]
# 10% completed in one hour... :( 
for (i in 1:length(indices)) {
  if (i %% 500 == 0) {
    print(i)
    save(collar.df.1h, file=here(procpath, paste0('era5estimates', i, '.rds')))
  }
  inx <- indices[i]
  rows <- which(collar.df.1h$ERA5_TEMP_INX == inx)
  pts <- collar.df.1h[rows,]
  r <- grib[[inx]]
  temps <- extract(r, pts)
  collar.df.1h$ERA5_TEMP_K[rows] <- temps
}
t1 <- proc.time()[['elapsed']]
t1 - t0
collar.df.1h <- collar.df.1h %>% 
  mutate(TEMPGROUP = factor(plyr::round_any(TEMP_DEG_C, 10)),
         SPECIES = ifelse(SPECIES=="Roan", "Roan", "Gemsbok"),
         ERA5_TEMP_C = ERA5_TEMP_K-273.15) 

# find daily maximums
day.temps <- collar.df.1h %>% nog() %>% 
  group_by(DATE) %>%
  summarize(MAXTEMP = max(ERA5_TEMP_C, na.rm=TRUE))

# bin days from mild to hot based on day temperatures
mean <- mean(day.temps$MAXTEMP)
sd <- sd(day.temps$MAXTEMP)
day.types <- day.temps %>%
  mutate(DAYTYPE = factor( ifelse(MAXTEMP < (mean-sd), "MILD",
                                  ifelse(MAXTEMP < (mean+sd), "AVG", "HOT")),
                           levels=c('MILD', 'AVG', 'HOT')))
types <- day.types$DAYTYPE[ match(collar.df.1h$DATE, day.types$DATE) ]
collar.df.1h$DAYTYPE <- types

### STS ###
save(collar.df.1h, file=here(procpath, 'era5estimates.rds'))