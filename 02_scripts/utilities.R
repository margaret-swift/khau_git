# utilities.R
# Created 08 Dec 2022
# Margaret Swift <mes114@duke.edu>

# Utilities file for KAZA antelope project

# ******************************************************************************
#                    LIBRARY LOADING & DEFINE GLOBAL OBJECTS
# ******************************************************************************
message("\nLoading all utility functions and parameters from utilities.R\n")

message('Loading base packages...')
pacman::p_load(tidyverse, patchwork, reshape2, 
               sf, lubridate, here)
# library(rgdal)
message('   ...Packages loaded.')

message("Setting base objects...")
datapath <- here::here("01_data")
scriptpath <- here::here("02_scripts")
outpath <- here::here("03_output")
# crs <- CRS("+proj=longlat +zone=34S +datum=WGS84")
message("   ...Base objects loaded.")

# Plot Params
message("Setting plot params...")
no.y.theme <- theme(axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.title.y=element_blank())
big.theme <- theme( text = element_text(size=18), 
                     title = element_text(size=20))
plot.theme <- theme( text = element_text(size=18), 
                     title = element_text(size=20),
                     axis.text.x = element_text(angle=45, vjust=0.5))
blank.theme <- theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
scalebar <- ggspatial::annotation_scale( pad_x = unit(0.05, "in"), 
                                         pad_y = unit(0.05, "in"),
                                         text_face="bold") 
message("   ...Plot params loaded.")

# ******************************************************************************
#                             BASE UTILITY FUNCTIONS
# ******************************************************************************
message("Loading functions...")

as.c <- as.character
as.n <- as.numeric
as.v <- as.vector
nog <- st_drop_geometry

removePunct <- function(v) gsub('\\)|\\(|\\/|\\.|\\?', '', v)
replaceSpace <- function(v, r="_") gsub(' ', r, v)
fixNames <- function(d) {
  names(d) <- replaceSpace(removePunct(toupper(names(d))))
  d
}
normalize <- function(x) {
  mu <- mean(x, na.rm=TRUE)
  sd <- sd(x, na.rm=TRUE)
  return( (x - mu) / sd )
}

setDataPaths <- function(dataname, verbose=TRUE) {
  # I know using <<- is dangerous but trust me... this is easier
  metapath <<- here::here(datapath, dataname, "meta")
  procpath <<- here::here(datapath, dataname, "processed")
  rawpath <<- here::here(datapath, dataname, "raw")
  if (verbose) {
    cat('resetting data paths to...\n')
    cat("  metapath: ", gsub(".*KAZA", "", metapath), '\n')
    cat("  rawpath:  ", gsub(".*KAZA", "", rawpath), '\n')
    cat("  procpath: ", gsub(".*KAZA", "", procpath), '\n')
  }
}

# get antelope data and temperature ranges
setDataPaths('antelope')
antfile <- file.path(procpath, 'antelope_collar_data.RData')
load(antfile)
temp.lims <- collar.df %>% nog() %>% 
  group_by(SEASON, TOD) %>% 
  summarize(min=quantile(TEMP_DEG_C, 0.05, na.rm=TRUE),
            max=quantile(TEMP_DEG_C, 0.95, na.rm=TRUE)) %>% 
  mutate(ID = paste0(SEASON, TOD)) %>% 
  column_to_rownames("ID")
temp.range <- seq(min(temp.lims$min), max(temp.lims$max), length.out=50)

message("   ...Basic functions loaded.")

# ******************************************************************************
#                             PLOTTING FUNCTIONS
# ******************************************************************************

# load protected areas maps
setDataPaths('khaudum_geography', verbose=FALSE)
load(here::here(procpath, "geographicData.RData"))
load(here::here(procpath, "fenceline.Rdata"))

khauPlot <- function(fill="transparent", lwd=1, ...) {
  ggplot(data=khau) + 
    geom_sf(fill=fill, linewidth=lwd, ...) + 
    plot.theme
}
kazaPlot <- function(fill="transparent", lwd=1, ...) {
  ggplot(data=kaza) + 
    geom_sf(fill=fill, linewidth=lwd, ...) + 
    plot.theme
}

message("   ...Plotting functions loaded.")


# ******************************************************************************
#                          OTHER CUSTOM FUNCTIONS
# ******************************************************************************

# Find distance in a row of collar.df
findDists <- function(row, df) {
  lon <- row$LON
  lat <- row$LAT
  latdist.m <- abs(df$latitude - lat) * 111.32 * 1000
  londist.m <- abs(df$longitude - lon) * 40075 * 1000 * cos(lat) / 360
  dist.m <- sqrt(londist.m^2 + latdist.m^2)
  return(min(dist.m))
}

# find landcover classes
findLC <- function(val) {
  inx <- ifelse(is.na(val), NA, 
                which.min(abs(lands.meta$category - val)))
  return( inx )
}
pullBufferCovariateData <- function(enc, collar.df, nhours=24, burnin=6) {
  
  getCovariates <- function(collar_i, range) {
    # get covariate data
    collar_i_sub <- collar_i[ collar_i$DATETIME %within% range, ]
    return(list(evi=collar_i_sub$EVI,
                water=collar_i_sub$DIST_NAT))
  }
  
  getCovariateStats <- function(covar_list) {
    # spit out means
    evi.mean = mean(covar_list$evi, na.rm=TRUE)
    water.mean = mean(covar_list$water, na.rm=TRUE)
    
    # spit out sds
    evi.sd = sd(covar_list$evi, na.rm=TRUE)
    water.sd = sd(covar_list$water, na.rm=TRUE)
    
    return(list(evi_mu=evi.mean, water_mu=water.mean,
                evi_sd=evi.sd, water_sd=water.sd))
  }
  
  # set up data needed for the FOR loop
  # message('Finding covariate data within ', nhours, ' hour window before and ',
  #         'after each crossing. Burnin (removed window): ', burnin, ' hours.')
  message('Window: ', nhours, 'hours; burnin: ', burnin, 'hours')
  th  <- hours(nhours)
  tb  <- hours(burnin)
  nrep = 2
  crosses <- enc %>% filter(EVENT %in% c("Cross", 'Jag'))
  time.names <- c('before', 'after')
  nrow = nrow(crosses)
  crosses$INX <- 1:nrow
  df <- data.frame(
    id = rep(crosses$AnimalID, each=nrep),
    burst = paste(sp, fixrate, rep(crosses$INX, each=nrep)),
    season="", 
    evi_mu=NA, evi_sd=NA, evi_sig=NA,
    water_mu=NA, water_sd=NA, water_sig=NA,
    time = factor(time.names, levels=time.names),
    side = ""
  )
  
  for (i in 1:nrow) {
    
    # choose encounter and pull id and times
    enc_i <- enc[i,]
    id <- enc_i$AnimalID
    t0 <- enc_i$start_time
    t1 <- enc_i$end_time
    
    # makes more sense to exclude points where they actually are 
    # crossing, since those will be similar EVI and water etc. 
    r0 <- interval(t0-th, t0-tb)
    r1 <- interval(t1+tb, t1+th)
    
    # get activity budgets
    collar_i <- collar.df %>% filter(COLLAR_ID == id) 
    
    # get before and after data
    before <- getCovariates(collar_i, r0)
    before_stats <- getCovariateStats(before)
    after  <- getCovariates(collar_i, r1)
    after_stats  <- getCovariateStats(after)
    
    # save to main df
    j <- (nrep*(i-1))+1
    inx <- j:(j+nrep-1)
    df$season[inx] <- enc_i$SEASON
    df$evi_mu[inx] <- c(before_stats$evi_mu, after_stats$evi_mu)
    df$evi_sd[inx] <- c(before_stats$evi_sd, after_stats$evi_sd)
    df$water_mu[inx] <- c(before_stats$water_mu, after_stats$water_mu)
    df$water_sd[inx] <- c(before_stats$water_sd, after_stats$water_sd)
    df$side[inx] <- c(enc_i$start_side, enc_i$end_side)
    
    # Welch's t tests
    if (sum(before$water) > 0) {
      wt.water <- t.test(before$water, after$water, alternative="two.sided", var.equal=FALSE)
      df$water_sig[inx] <- wt.water$p.value
    }
    wt.evi <- t.test(before$evi, after$evi, alternative="two.sided", var.equal=FALSE)
    df$evi_sig[inx] <- wt.evi$p.value
  }
  df <- df %>%
    mutate(season = factor(season),
           time=factor(time),
           evi_min = evi_mu-evi_sd,
           evi_max = evi_mu+evi_sd,
           water_min = water_mu-water_sd,
           water_max = water_mu+water_sd)
  return(df)
}



message("   ...Custom functions loaded.")


# ******************************************************************************
#                          ADVANCED HOUSEKEEPING
# ******************************************************************************

# Set of functions to list all custom functions and variables
listFuns <- function() {
  df <- my.funs
  types <- levels(df$ftype)
  for (type in types) {
    cat(type, "\n", df$fname[df$ftype == type], "\n\n")
  }
}
listVars <- function() {
  df <- my.vars
  types <- levels(df$vtype)
  for (type in types) {
    cat(type, "\n", df$vname[df$vtype == type], "\n\n")
  }
}
.listVars <- function() {
  lst <- ls(envir = .GlobalEnv)
  inx <- sapply(lst,function(var) any(class(get(var))!='function'))
  vs <- lst[inx]
  
  getType <- function(val) {
    val <- tolower(val)
    type <- "OTHER"
    if (grepl("path", val)) { type = "PATH"
    } else if (grepl("theme|scale|crs", val)) { type = "THEME"
    } else if (grepl("khau|kaza|waters", val)) { type = "MAPDATA"
    } 
    type = factor(type, levels=c("PATH", "THEME", "MAPDATA", "OTHER"))
  }
  vtypes <- unlist(lapply(vs, getType))
  df <- data.frame(vtype=vtypes, vname=vs) %>% arrange(vtype)
  df
}
.listFuns <- function() {
  lst <- ls(envir = .GlobalEnv)
  inx <- sapply(lst,function(var) any(class(get(var))=='function'))
  fs <- lst[inx]
  
  getType <- function(val) {
    val <- tolower(val)
    type <- "OTHER"
    if (grepl("as.", val)) { type = "LAZY"
    } else if (grepl("plot", val)) { type = "PLOTTING"
    } else if (grepl("list", val)) { type = "LISTING"
    } else if (grepl("fix|^re", val)) { type = "FIXING"
    } else if (grepl("set", val)) { type = "HOUSEKEEPING"
    } 
    type = factor(type, levels=c("HOUSEKEEPING", "FIXING", "LISTING", 
                                 "PLOTTING", "LAZY", "OTHER"))
  }
  ftypes <- unlist(lapply(fs, getType))
  df <- data.frame(ftype=ftypes, fname=fs) %>% arrange(ftype)
  df
}

# only create list of custom vars and funs once, via housekeeping switch
if ( !exists('housekeeping_switch') ) {
  my.vars <- .listVars()
  my.funs <- .listFuns()
  housekeeping_switch = TRUE
} 

message("   ...Advanced housekeeping functions loaded.\n")

# ******************************************************************************
message("All utility functions and parameters loaded.")
message("Use listFuns() to list all available utility functions, 
        and listVars() to list all base parameters.")
# ******************************************************************************
