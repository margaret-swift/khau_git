# analyzemicroclimates.R
# Created 02 July 2023
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am("02_scripts/analysis/microclimates/analyzeMicroclimates.R")
source(here::here("02_scripts", "utilities.R"))
pacman::p_load('ggspatial')
load(antfile)
setDataPaths('temperature')
load(here(procpath, 'era5estimates.rds'))

### FUNCTIONS ###
sumTemps <- function(v)  {
  v %>% summarize(oryxmin=quantile(Gemsbok, 0.05), 
                  oryxmax=quantile(Gemsbok, 0.95),
                  roanmin=quantile(Roan, 0.05), 
                  roanmax=quantile(Roan, 0.95))
}
ttestMe <- function(v) {
  print( sumTemps(v) )
  print( t.test(v$Roan, v$Gemsbok, paired=TRUE) )
  effsize::cohen.d(v$Roan, v$Gemsbok, paired=TRUE)
}


# ******************************************************************************
#                   MODELING MICROCLIMATE TEMPERATURE BY ERA5
# ******************************************************************************

# running simple LM on the two temps
lm_temp <- lm(TEMP_DEG_C ~ SEASON * ERA5_TEMP_C + 
                SPECIES * ERA5_TEMP_C, data=collar.df.1h)
summary(lm_temp)

################################################
## FUNCTIONS ##
################################################
modelMe <- function(sp, breaks) {
  sub <- collar.df.1h %>% filter(SPECIES == sp)
  x <- sub$ERA5_TEMP_C
  sub$categ <- ifelse(x<breaks[1], 'first', 
                 ifelse(x<breaks[2], 'second', 'third'))
  m <- lm(TEMP_DEG_C ~ ERA5_TEMP_C*categ + SEASON, data=sub)
  return(m)
}
plotModel <- function(sp, m, breaks) {
  res <- resid(m)
  plot(fitted(m), res)
  abline(0,0, col='blue')
  
  x=sub$ERA5_TEMP_C[sub$SPECIES==sp]
  y=sub$TEMP_DEG_C[sub$SPECIES==sp]
  plot(x,y, main=sp, xlab='ERA5', ylab='collar temperature')
  newx = data.frame(x=seq(0, 40, by=0.5))
  newx$categ=ifelse(newx$x<breaks[1], 'first', 
                    ifelse(newx$x<breaks[2], 'second', 'third'))
  lines(newx$x, predict.lm(m, newdata=newx), col='red', lwd=2)
}
################################################
################################################
################################################

# Running model for gemsbok
breaks.o <- c(13.5, 29.5)
m.o <- modelMe("Gemsbok", breaks.o)
plotModel("Gemsbok", m.o, breaks.o)

### best breaks for roan: 15.5 and 31
breaks.r <- c(15.5, 31)
m.r <- modelMe("Roan", breaks.r)
plotModel("Roan", m.r, breaks.r)

save(m.o, m.r, breaks.o, breaks.r, file=here(procpath, 'breaksAndModels.Rdata'))



# ******************************************************************************
#                   TEMPERATURE DATA -- COMPARING AVG TO HOT
# ******************************************************************************

diffs.df <- collar.df.1h %>% nog() %>% 
  filter(!is.na(TEMP_DEG_C),
         COLLAR_ID != "SAT1728") %>% #this one roan is out of bounds
  filter(ERA5_TEMP_C > 35) %>%
  dplyr::select(SPECIES, DATE, DAYTYPE, SEASON, TOD, 
                TEMP_DEG_C, ERA5_TEMP_C) %>% 
  mutate(diff = TEMP_DEG_C - ERA5_TEMP_C)

myMean <- function(x, by=2) round(mean(x, na.rm=TRUE), by)
mySD <- function(x, by=2) round(sd(x, na.rm=TRUE), by)
myCohen <- function(x) paste(round(x$estimate, 2), '(', as.c(x$magnitude), ')')
ttestMax <- function(sp, szn) {
  message(sp, ' ', szn)
  data = diffs.df %>% filter(SEASON==szn, SPECIES==sp)
  avg <- data$diff[data$DAYTYPE == "AVG"]
  hot <- data$diff[data$DAYTYPE == "HOT"]
  message('mean (avg): ', myMean(avg), ' +/- ', mySD(avg))
  message('mean (hot): ', myMean(hot), ' +/- ', mySD(hot))
  test <- t.test(avg, hot)
  d <- effsize::cohen.d(avg, hot)
  message('p value: ', test$p.value)
  message("cohen's d: ",myCohen(d))
}
ttestMax("Gemsbok", "DRY")
ttestMax("Gemsbok", "WET")
ttestMax("Roan", "DRY")
ttestMax("Roan", "WET")





# ******************************************************************************
#                   TEMPERATURE DATA -- SPECIFIC DIFFERENCES
# ******************************************************************************

# calculating species differences in experienced microclimate temperatures
diffs <- collar.df.1h %>% nog() %>% 
  group_by(SPECIES, SEASON, DATE, TOD, HOUR) %>% 
  summarize(TEMP = mean(TEMP_DEG_C)) %>% 
  dcast(SEASON + DATE + TOD + HOUR ~ SPECIES) %>% 
  mutate(diff = Roan - Oryx) 
wet <- diffs %>% filter(SEASON == "WET") %>% filter(!is.na(diff))
dry <- diffs %>% filter(SEASON == "DRY") %>% filter(!is.na(diff))


# comparing species differences in experienced microclimate temperatures
ttestMe(dry)
ttestMe(wet)

# JUST during the daytime
ttestMe(dry %>% filter(TOD=="DAY"))
ttestMe(wet%>% filter(TOD=="DAY"))

# Daily max/min
diffs.daily <- collar.df.1h %>% nog() %>% 
  group_by(SPECIES, SEASON, DATE, COLLAR_ID) %>% 
  filter(!is.na(TEMP_DEG_C)) %>% 
  summarize(min = min(TEMP_DEG_C, na.rm=TRUE),
            max = max(TEMP_DEG_C, na.rm=TRUE),
            n=n()) %>%
  filter(!is.na(min), !is.na(max),
         !is.infinite(min), !is.infinite(max),
         n > 15) %>% 
  mutate(range=max-min)
  group_by(SPECIES, SEASON, DATE) %>% 
  summarize(m_range = mean(range))

kwTest <- function(sp, seas="DRY") {
  dat <- diffs.daily %>% 
    ungroup() %>% 
    filter(SPECIES == sp, SEASON == seas) %>% 
    dplyr::select(DAYTYPE, m_range)
  
  # run a levene test
  LT <- leveneTest(m_range ~ DAYTYPE, data=dat)
  print(LT)
  
  # if levene test H is rejected, variance is not equal
  # we can't use ANOVA, so we use Kruskal-Wallis test
  KT <- kruskal.test(m_range ~ DAYTYPE, data=dat)
  print(KT)
  
  # Effect size
  ES <- kruskal_effsize(formula=m_range ~ DAYTYPE, data=dat)
  print(ES)
}

kwTest("Roan", "DRY")
kwTest("Gemsbok", "DRY")

kwTest("Roan", "WET")
kwTest("Gemsbok", "WET")





# ******************************************************************************
#                   TEMPERATURE DATA -- ERA5 CORRELATIONS
# ******************************************************************************

######### CHECKING ERA5 CORRELATION ETC
collar.df.1h$ERA5_TEMP_C = collar.df.1h$ERA5_TEMP_K - 273.15
sub <- collar.df.1h %>% filter(!is.na(ERA5_TEMP_K))
m <- lm( TEMP_DEG_C ~ ERA5_TEMP_K + SPECIES, data=sub)
GGally::ggpairs(sub %>% filter(SPECIES=="Roan") %>%  
                  dplyr::select(TEMP_DEG_C, ERA5_TEMP_C))
GGally::ggpairs(sub %>% filter(SPECIES=="Gemsbok") %>%  
                  dplyr::select(TEMP_DEG_C, ERA5_TEMP_C))


### Correlation by collar ID
ids <- unique(collar.df.1h$COLLAR_ID)
for (id in ids) {
  
  dat <- collar.df.1h %>% 
    filter(!is.na(ERA5_TEMP_K),
           !is.na(TEMP_DEG_C),
           COLLAR_ID == id) %>% 
    nog()
  corr <- round(cor(dat$TEMP_DEG_C, dat$ERA5_TEMP_K), 3)
  message(id, " : ", corr)
}


