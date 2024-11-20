# resampleCollarData.R
# Created 22 Apr 2023
# Margaret Swift <mes114@duke.edu>

# Resampling and shaving collar data to regularize

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(lubridate)
load(antfile)


# ******************************************************************************
#                             FIXING UP AND CLEANING DATA
# ******************************************************************************

# Special case: SAT1515
d <- collar.df %>% filter(COLLAR_ID == "SAT1515")
inx <- d$DATETIME < "2015-04-02"
a.inx <- match(d$INX[inx], collar.df$INX)
b.inx <- match(d$INX[!inx], collar.df$INX)
collar.df$COLLAR_ID[a.inx] <- "SAT1515a"
collar.df$COLLAR_ID[b.inx] <- "SAT1515b"

# giving collar types
collar.df <- collar.df %>% 
  mutate(COLLAR_TYPE = ifelse(COLLAR_TYPE == "SAT", "SAT", "SAT+"))


# ******************************************************************************
#                                     FUNCTIONS
# ******************************************************************************

# resamples data to a certain fix rate
resampleMe <- function(fixmax, diffmin, collartype=NULL) {
  if (is.null(collartype)) { df <- collar.df
  } else { df <- collar.df %>% filter(COLLAR_TYPE == collartype) }
  
  dates.freq <- df %>%
    nog() %>%
    group_by(DATE) %>% 
    group_by(COLLAR_ID, DATE) %>%
    summarize(n=n()) %>% 
    dcast(DATE ~ COLLAR_ID) %>% 
    column_to_rownames(var="DATE")
  ids <- colnames(dates.freq)
  rem.inx <- c()

  for ( i in 1:length(ids)) {
    id <- ids[i]
    message("resampling ", id)
    
    # pull data and find days with too many fixes
    f <- data.frame(fixes=dates.freq[,i])
    f$inx <- 1:nrow(f)
    row.names(f) <- row.names(dates.freq)
    errors <- which(f$fixes > fixmax)
    
    # loop over each burst of errors
    bursts <- split(errors, cumsum(c(1, diff(errors) != 1)))
    for (j in 1:length(bursts)) {
      burst <- bursts[[j]]
      message('burst: ', j)
      if (length(burst)) {
        # get dates and filter df to match
        dates <- as.Date(row.names(f)[burst])
        data <- df %>%
          filter(DATE %in% dates, COLLAR_ID == id) %>%
          nog()
        
        # Iteratively remove rows and check each next row's datetime
        # against the new previous until fixes are ~5hr apart.
        d <- data
        for (k in 2:nrow(data)) {
          d.inx <- which(d$INX == data$INX[k])[1]
          if (d.inx > 1) {
            cur <- d$DATETIME[d.inx]
            prev <- d$DATETIME[d.inx-1]
            diff <- difftime(cur, prev, units="mins")
            if ( diff < diffmin ) d <- d[-d.inx,]
          }
        }
        # save original INX of data points to remove
        rem.inx <- c(rem.inx, data$INX[ which(!(data$INX %in% d$INX)) ])
      }
    }
  }
  # in case inx doesn't match INX
  to.rem <- match(rem.inx, df$INX)
  df <- df[-to.rem,]
  df
}

# shaves data down on the edges and after too many errors
shaveMe <- function(df,
                    x, # appropriate number of fixes per day
                    t, # we want x fixes per day, but we will accept x+t or x-t
                    maxerr, #max errors allowed
                    shave) {#how many data points from start and end to shave
  ideal <- (x-t):(x+t)
  
  dates.freq <- df %>%
    st_drop_geometry() %>%
    group_by(DATE) %>% 
    group_by(COLLAR_ID, DATE) %>%
    summarize(n=n()) %>% 
    dcast(DATE ~ COLLAR_ID) %>% 
    column_to_rownames(var="DATE")
  
  cuts <- data.frame(ID=names(dates.freq), first=0, last=0)
  rownames(cuts) <- cuts$ID
  for (i in 1:length(dates.freq)) {
    # pull data
    data <- data.frame(fixes=dates.freq[,i])
    data$inx <- 1:nrow(data)
    row.names(data) <- row.names(dates.freq)
    id <- cuts$ID[i]
    message("shaving ", id)
    
    # shave first and last bits
    reals <- which(!is.na(data$fixes))
    first=reals[1]+shave
    last=reals[length(reals)]-shave
    data <- data %>% filter(inx %in% first:last)
    
    # find data that is outside of our bounds (or missing)
    x <- !(data$fixes %in% ideal)
    errors <- cumsum(replace_na(x, 0))
    
    # remove data after a certain number of 'errors' (cutoff)
    dates <- row.names(data)
    overs <- errors>=maxerr
    LD <- length(dates)
    cut.inx <- ifelse(!any(overs), LD, which(overs)[1])
    cuts$first[i] <- dates[1]
    cuts$last[i] <- dates[cut.inx]
  }
  cuts$interval <- interval(cuts$first, cuts$last)
  df <- df %>%
    filter(DATE %within% cuts[COLLAR_ID, 'interval']) %>% 
    mutate(FIRST = COLLAR_ID != lag(COLLAR_ID),
           FIRST = ifelse(is.na(FIRST), TRUE, FIRST))
  df
}

# histogram for data
histMe <- function(data, type=NULL, ...) {
  p <- data %>% 
    st_drop_geometry() %>%
    group_by(DATE) %>%
    ggplot() +
    geom_histogram(aes(x=DATE), binwidth=1, ...) +
    plot.theme
  if (!is.null(type)) p <- p + facet_wrap(~ANIMAL_ID)
  p
}


# ******************************************************************************
#                           RESAMPLING AND SHAVING
# ******************************************************************************

# five-hour fixes can resample the 1hour data down to an appropriate level
collar.df.5h <- resampleMe(fixmax=7, diffmin=275)
histMe(collar.df.5h, "faceted", fill="orange")
collar.df.5h <- shaveMe(df=collar.df.5h, x=5, t=2, maxerr=15, shave=5)
histMe(collar.df.5h, "faceted", fill="#008080")

# one-hour fixes need to just get rid of the 5h data
collar.df.1h <- resampleMe(fixmax=25, diffmin=50, collartype="SAT+")
histMe(collar.df.1h, "faceted", fill="orange")
collar.df.1h <- shaveMe(df=collar.df.1h, x=24, t=8, maxerr=30, shave=5)
histMe(collar.df.1h, "faceted", fill="#008080")

# ******************************************************************************
#                                   SAVING
# ******************************************************************************

# save
save(collar.df, 
     collar.meta, 
     collar.df.1h,
     collar.df.5h, 
     fake.df, file=antfile)
