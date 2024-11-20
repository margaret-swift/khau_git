# combineEncountersActivities.R
# Created 24 April 2023
# Margaret Swift <mes114@duke.edu>

# Script to combine encounters and activities

# ******************************************************************************
#                               DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here("02_scripts", "utilities.R"))
load(antfile)
load(here('03_output', 'hmm', 'hmm1h5hUTM.rdata'))
load(here("03_output", "barriers", "encounters.RData"))
state.names <- c('ruminating', 'exploring', 'foraging')
time.names <- c("before", "after")
time.names <- factor(time.names, levels=time.names)


# ******************************************************************************
#                                     FUNCTIONS
# ******************************************************************************
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
Mean <- function(x) {
  return(mean(x, na.rm=TRUE))
}
# getActivityBudget <- function(hmm_i, collar_i, range) {
#   # get covariate data
#   collar_i_sub <- collar_i[ collar_i$DATETIME %within% range, ]
#  
#   # summarize HMM states
#   hmm_i_sub <- hmm_i[ hmm_i$DATETIME %within% range, ]
#   states <- hmm_i_sub %>% 
#     summarize(n=n()) %>% 
#     complete(STATE) %>%
#     mutate(n=ifelse(is.na(n), 0, n)) %>% 
#     dplyr::select(n) %>% 
#     unlist() %>% as.n()
#   
#   # spit out means
#   ls <- list(states=states,
#              alt=Mean(collar_i_sub$ALT_M),
#              lc=Mode(collar_i_sub$LC_CATEG),
#              evi=Mean(collar_i_sub$EVI),
#              water=Mean(collar_i_sub$DIST_NAT),
#              step=Mean(hmm_i_sub$step),
#              angle=Mean(hmm_i_sub$angle)
#   ) 
#   return(ls)
# }

getCovariates <- function(collar_i, hmm_i,range) {
  # get covariate data
  collar_i_sub <- collar_i[ collar_i$DATETIME %within% range, ]
  
  # summarize HMM states
  hmm_i_sub <- hmm_i[ hmm_i$DATETIME %within% range, ]
  states <- table(hmm_i_sub$STATE)
  return(list(states=states, 
              evi=collar_i_sub$EVI,
              water=collar_i_sub$DIST_NAT))
}
getCovariateStats <- function(covar_list) {
  covar_list$water[covar_list$water == -99] = 7000
  
  # spit out means
  evi.mean = mean(covar_list$evi, na.rm=TRUE)
  water.mean = mean(covar_list$water, na.rm=TRUE)
  
  # spit out sds
  evi.sd = sd(covar_list$evi, na.rm=TRUE)
  water.sd = sd(covar_list$water, na.rm=TRUE)
  
  return(list(evi_mu=evi.mean, water_mu=water.mean,
              evi_sd=evi.sd, water_sd=water.sd))
}



combineEncounterActivities <- function(sp, fixrate, nhours=24, burnin=0) {
  # grab data name
  data.name = paste0(tolower(sp), '_', fixrate)
  message(data.name)
  
  # pull data from HMM and encounters for given species and fixrate
  hmm <- hmm.data[[data.name]] %>% group_by(STATE) %>% mutate(STATE=factor(STATE, levels=state.names))
  ids <- unique(hmm$ID)
  enc <- encounters %>% filter(AnimalID %in% ids)
  
  # set up parameters
  th  <- hours(nhours)
  tb  <- hours(burnin)
  NT = length(time.names)
  NS = length(state.names)
  nrep = NT
  nrow = nrow(enc)
  
  # set up empty DF to hold values
  df <- data.frame(
    id = rep(enc$AnimalID, each=nrep),
    sp = sp,
    fixrate=fixrate,
    period = nhours,
    burnin = burnin,
    burst = rep(1:nrow, each=nrep),
    season="", 
    btype = NA,
    time = time.names,
    side = "",
    evi_mu=NA, evi_sd=NA, evi_sig=NA, water_mu=NA, water_sd=NA, water_sig=NA,
    ruminating = NA, exploring=NA, foraging=NA, states_sig = NA
  )
  
  for (i in 1:nrow) {
    if (i %% 25 == 0) message(i, ' of ', nrow)
    
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
    hmm_i <- hmm %>% filter(ID == id) 
    collar_i <- collar.df %>% filter(COLLAR_ID == id) 
    # before <- getActivityBudget(hmm_i, collar_i, r0)
    # after <- getActivityBudget(hmm_i, collar_i, r1)
    
    # get before and after data
    before <- getCovariates(collar_i, hmm_i, r0)
    before_stats <- getCovariateStats(before)
    after  <- getCovariates(collar_i, hmm_i, r1)
    after_stats  <- getCovariateStats(after)
    
    # get indices to slot into main df
    j <- (nrep*(i-1))+1
    inx <- j:(j+nrep-1)
    
    # save to main df
    df$btype[inx] <- enc_i$EVENT
    df$season[inx] <- enc_i$SEASON
    df$side[inx] <- c(enc_i$start_side, enc_i$end_side)
    
    df[inx, names(before$states)] <- c(before$states, 
                                       after$states)
    
    df$evi_mu[inx] <- c(before_stats$evi_mu, after_stats$evi_mu)
    df$evi_sd[inx] <- c(before_stats$evi_sd, after_stats$evi_sd)
    df$water_mu[inx] <- c(before_stats$water_mu, after_stats$water_mu)
    df$water_sd[inx] <- c(before_stats$water_sd, after_stats$water_sd)
    
    # chisq test
    cs <- chisq.test(before$states, after$states)
    df$states_sig[inx] <- cs$p.value
    
    # Welch's t tests
    identifier = paste0(enc_i$AnimalID, ', fixrate ', fixrate, ', burst ', df$burst[j])
    if (length(before$water) > 2 && length(after$water) > 2) {
      if (before_stats$water_mu == after_stats$water_mu) {
        df$water_sig[inx] <- 1.00
      } else {
        wt.water <- t.test(before$water, after$water, alternative="two.sided", var.equal=FALSE)
        df$water_sig[inx] <- wt.water$p.value
      }
    } else {
      message('not enough data to run water comparison for ', identifier)
    }
    if (length(before$evi) > 2 && length(after$evi) > 2) {
      wt.evi <- t.test(before$evi, after$evi, alternative="two.sided", var.equal=FALSE)
      df$evi_sig[inx] <- wt.evi$p.value
    } else {
      message('not enough data to run evi comparison for ', identifier)
    }
  }
  return(df)
}


# ******************************************************************************
#                               COMBINE AND SAVE
# ******************************************************************************
#   CREATE COMBINED DATA
nhours <- 100; burnin <- 48;
roan.hmm.enc.1h <- combineEncounterActivities('Roan', '1h', nhours, burnin)
oryx.hmm.enc.1h <- combineEncounterActivities('Oryx', '1h', nhours, burnin)
roan.hmm.enc.5h <- combineEncounterActivities('Roan', '5h', nhours, burnin)
oryx.hmm.enc.5h <- combineEncounterActivities('Oryx', '5h', nhours, burnin)

hmm.enc.data.ls <- list(roan.hmm.enc.1h, oryx.hmm.enc.1h, 
                        roan.hmm.enc.5h, oryx.hmm.enc.5h)
names(hmm.enc.data.ls) <- names(hmm.data)
hmm.enc.data <- rbind(roan.hmm.enc.1h, oryx.hmm.enc.1h,
                      roan.hmm.enc.5h, oryx.hmm.enc.5h)

#   SAVE
save(hmm.enc.data.ls, hmm.enc.data,
     file=here("03_output", "barriers", "encounters_activities_combined.rdata"))
