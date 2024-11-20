# analyzeResourceUse.R
# Created 24 April 2023
# Margaret Swift <mes114@duke.edu>

#   - DO ROAN CROSS TO ACCESS BETTER RESOURCES?


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here("02_scripts", "utilities.R"))
load(antfile)
load(here("03_output", "barriers", "encounters.RData"))
enc <- encounters

time.names <- factor(c("before", "after"), levels=c('before', 'after'))
time.colors<- c('#C12D3A', 'black')
side.names <- c("NAM", "BOTS")
side.colors <- c("#C12D3A", "#7BA7CE")


# ******************************************************************************
#                                     FUNCTIONS
# ******************************************************************************
getCovariates <- function(collar_i, range) {
  # get covariate data
  collar_i_sub <- collar_i[ collar_i$DATETIME %within% range, ]
  
  # spit out means
  evi.mean = mean(collar_i_sub$EVI, na.rm=TRUE)
  water.mean = mean(collar_i_sub$DIST_NAT, na.rm=TRUE)
  
  # spit out sds
  evi.sd = sd(collar_i_sub$EVI, na.rm=TRUE)
  water.sd = sd(collar_i_sub$DIST_NAT, na.rm=TRUE)
  return(list(evi_mu=evi.mean, water_mu=water.mean,
              evi_sd=evi.sd, water_sd=water.sd))
}
pullBufferCovariateData <- function(hrs=24) {
  
  # set up the two datasets
  crosses <- encounters %>% filter(EVENT == "Cross")
  buffer <- hours(hrs)
  nrep = 2
  nrow = nrow(crosses)
  crosses$INX <- 1:nrow
  df <- data.frame(
    id = rep(crosses$AnimalID, each=nrep),
    burst = rep(crosses$INX, each=nrep),
    season="", evi_mu=NA, water_mu=NA, evi_sd=NA, water_sd=NA,
    time = time.names,
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
    # let's try a separate buffer of 6 hours as well.
    h <- hours(6)
    r0 <- interval(t0-buffer, t0-h)
    r1 <- interval(t1+h, t1+buffer)
    
    # get activity budgets
    collar_i <- collar.df %>% filter(COLLAR_ID == id) 
    before <- getCovariates(collar_i, r0)
    after  <- getCovariates(collar_i, r1)
    
    # save to main df
    j <- (nrep*(i-1))+1
    inx <- j:(j+nrep-1)
    df$season[inx] <- enc_i$SEASON
    df$evi_mu[inx] <- c(before$evi_mu, after$evi_mu)
    df$evi_sd[inx] <- c(before$evi_sd, after$evi_sd)
    df$water_mu[inx] <- c(before$water_mu, after$water_mu)
    df$water_sd[inx] <- c(before$water_sd, after$water_sd)
    df$side[inx] <- c(enc_i$start_side, enc_i$end_side)
  }
  df <- df %>%
    mutate(evi_min = evi_mu-evi_sd,
           evi_max = evi_mu+evi_sd,
           water_min = water_mu-water_sd,
           water_max = water_mu+water_sd,)
  return(df)
}

# ******************************************************************************
#                     DO ROAN CROSS TO ACCESS BETTER RESOURCES?
# ******************************************************************************

df.roan <- pullBufferCovariateData() 

# categorize whether botswana value is larger
after.better <- df.roan %>% 
  group_by(burst) %>% 
  reframe(evi_mu = evi_mu[time=="after"] > evi_mu[time=="before"],
          water_mu = water_mu[time=="after"] < water_mu[time=="before"],
  )

# attach larger botswana flag to main sub dataset
inx <- match(df.roan$burst, after.better$burst)
df.roan$post.evi.better <- after.better$evi_mu[inx]
df.roan$post.water.better <- after.better$water_mu[inx]

# testing to show there's no difference before and after
testSides <- function(seas) {
  message(seas, " SEASON PAIRED T-TESTS")
  before <- df.roan %>% filter(time=="before", season==seas)
  after <- df.roan %>% filter(time=="after", season==seas)
  message('EVI') 
  print( t.test(before$evi_mu,
                after$evi_mu,
         paired=TRUE) )
  
  message('WATER')
  print( t.test(before$water_mu,
                after$water_mu,
         paired=TRUE) )
}
testSides('WET')
testSides('DRY')


lm <- lme4::glmer(side ~ water_mu + evi_mu + season + (1 | id), 
                  data=df.roan, 
                  family="binomial")
print(summary(lm))

# ******************************************************************************
#                                     PLOTTING!
# ******************************************************************************
time.colors<- c('#C12D3A', 'black')

# EVI doesn't change on either side
df.roan$offset <- runif(nrow(df.roan), 0, 1)
df.roan <- df.roan %>% mutate(time2 = ifelse(time=="after", offset+5, offset))

p1 <- df.roan %>% 
  # ggplot(mapping=aes(x=time, y=evi_mu, group=burst, color=post.evi.better)) +
  # --------
  # uncomment for sd segments
  ggplot(mapping=aes(x=time2, y=evi_mu, group=burst, color=post.evi.better)) +
  geom_segment(mapping=aes(xend=time2,y=evi_min, yend=evi_max), alpha=0.3, linewidth=1) +
  # -------
  geom_point(size=2) + 
  geom_line(linewidth=1) + 
  theme(text=element_text(size=16)) +
  facet_wrap(~season) +
  guides(color="none") + scale_color_manual(values=time.colors) +
  ylab("mean EVI") + xlab('')

p2 <- df.roan %>% 
  # ggplot(aes(x=time, y=water_mu, group=burst, color=post.water.better)) + 
  # --------
  # uncomment for sd segments
  ggplot(aes(x=time2, y=water_mu, group=burst, color=post.water.better)) + 
  geom_segment(mapping=aes(xend=time2,y=water_min, yend=water_max), alpha=0.3, linewidth=1) +
  #   # --------
  geom_point(size=2) + geom_line(linewidth=1) + 
  facet_wrap(~season) +
  theme(text=element_text(size=16)) +
  guides(color="none") + scale_color_manual(values=time.colors) + 
  ylab('mean distance from water (m)') + xlab('')

p2 + p1 + 
  plot_annotation('24-hour mean resource access, before and after fence crossing',
                          theme=theme(plot.title=element_text(hjust=0.5, size=24)))
                          
# ggsave(file=here("03_output","barriers","resource_results_seasonal.png"),
       # width=12, height=7)
ggsave(file=here("03_output","barriers","resource_results_seasonal_jitter.png"),
       width=12, height=7)

# LC doesn't change on either side
(p4 <- df.roan.sub %>% 
    ggplot(aes(x=side, y=lc, group=burst)) + 
    geom_point(size=2) + geom_line(linewidth=1, alpha=0.1,position=position_jitter(w=0, h=0.2)) + 
    theme(text=element_text(size=18)) +
    guides(color="none") + scale_color_manual(values=time.colors) +
    ylab("landcover class"))


# ******************************************************************************
#                             SENSITIVITY FOR RESOURCES
# ******************************************************************************
getSums <- function(hrs) {
  sums <- combineEncounterActivities('Roan', roan, hrs=hrs) %>% 
    filter(btype == "Cross", season=="WET", state=="ruminating") %>% 
    dplyr::select(-n, -state) %>% 
    mutate(side=factor(side, levels=side.names)) %>% 
    group_by(side) %>% 
    dplyr::select(side, evi, dist_nat, dist_riv) %>% 
    summarize_all(.funs=list('mean'=mean, 'sd'=sd)) %>% 
    melt() %>% 
    mutate(measure=gsub("_mean|_sd", "", variable),
           stat=gsub('.*_', "", variable),
           timestep=hrs) %>% 
    dplyr::select(side, timestep, measure, stat, value)
  sums
}

# this takes a minute
t0 <- proc.time()[['elapsed']]
hrs <- c(8, 10, 16, 24, 48)
for (i in 1:length(hrs)) {
  hr <- hrs[i]
  if (i == 1) sums <- getSums(hr)
  else sums <- rbind(sums, getSums(hr))
}
proc.time()[['elapsed']]-t0
head(sums)

getData <- function(meas) {
  x <- sums %>% 
    filter(measure == meas) %>% 
    group_by(side, timestep) %>% 
    summarize(mean=first(value),
              sd=last(value)) %>% 
    mutate(min=mean-sd, max=mean+sd,
           scale=ifelse(side=="Botswana", 0, 0.5))
  x
}

plotData <- function(meas) {
  data <- getData(meas)
  p <- ggplot() +
    geom_segment(data = data,
                 aes(x=timestep+scale, xend=timestep+scale, 
                     y=min, yend=max, color=side), alpha=0.2, linewidth=2) +
    geom_point(data=data,
               aes(x=timestep+scale, y=mean, color=side)) +
    geom_line(data = data,
              aes(x=timestep+scale, y=mean, color=side)) +scale_color_manual(values=side.colors) + 
    plot.theme + 
    xlab('hour buffer')
  p
}

p1 <- plotData('dist_riv') + ylab('mean dist from river (m)')
p2 <- plotData('dist_nat') + ylab('mean dist from waterholes (m)')
p3 <- plotData('evi') + ylab('mean EVI')
p1 + p2 + p3 + plot_layout(guides="collect")

outfile <- here('03_output', 'barriers', 'sensitivity_testing', 'TEST_hourbuffer.png')
ggsave(outfile, width=16.5, height=7)
