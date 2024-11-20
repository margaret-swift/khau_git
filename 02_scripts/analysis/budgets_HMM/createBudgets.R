# hiddenMarkovModel_momentuHMM.R
# Created 22 March 2023
# Margaret Swift <mes114@duke.edu>

# HIDDEN MARKOV MODEL FITTING WITH momentuHMM
# PAPER: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12578
# VIGNETTE: https://cran.r-project.org/web/packages/momentuHMM/vignettes/momentuHMM.pdf

# Reasons to use momentuHMM instead:
#   Although I don't have irregular, high-error, or multivariate data,
#   one benefit of momentuHMM over moveHMM is the former's ability to
#   take factorial environmental data.


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

pacman::p_load(momentuHMM, rgdal, tictoc, here, doRNG)
source(here::here("02_scripts", "utilities.R"))
load(antfile)
collar.meta$ID = paste0('SAT', collar.meta$COLLAR_ID)


# ******************************************************************************
#                                SET UP FOR HMM
# ******************************************************************************

setUpData <- function(d, ts) {
  ctype = ifelse(ts=="1 hour", "SAT+", "SAT")
  df <- d %>%
    st_transform(crs = 3578) %>%
    mutate(ID = COLLAR_ID, 
           DATE = date(DATETIME),
           X = sf::st_coordinates(.)[,1],
           Y = sf::st_coordinates(.)[,2],
           diff = (DATETIME - lag(DATETIME)) / 60) %>% 
    nog() %>% 
    filter(SEX == "F", 
           COLLAR_TYPE == ctype, 
           diff > 30 # remove duplicate rows & ones that are too often!
    ) 
  
  df <- df[,c('INX', 'ID', 'DATETIME', 'X', 'Y')]
  tic()
  # singular: 1718, 1722, 1723
  ln.prior <- function(theta) dnorm(theta[2],-4,2,log=TRUE)
  crw <- crawlWrap(obsData=df, 
            timeStep=ts,
            coord = c('X', 'Y'),
            Time.name = "DATETIME",
            fillCols=TRUE,
            attempts=20,
            theta=c(6.855, -0.007), 
            fixPar=c(NA,NA),
            prior=ln.prior)
  toc()
  tic()
  pd <- prepData( crw, type="UTM",
              covNames=c('INX'))
  toc()
  pd <- pd %>% left_join(collar.meta[,c('ID', 'ANIMAL')], by='ID')
  pd
}

# fix up collar.df.1h
inx.rm.1h <- c(29, 11491, 11492, 49323, 53986, 102457, 152736)
collar.df.1h <- collar.df.1h[-inx.rm.1h]

# fix up collar.df.5h
inx.rm.5h <- c( 23974 )
collar.df.5h <- collar.df.5h[-inx.rm.5h]

data.1h <- setUpData(collar.df.1h, ts="1 hour")
data.5h <- setUpData(collar.df.5h, ts="5 hours")
save(data.1h, data.5h, file=here::here('03_output', 'hmm', 'pdata1h5hUTM.rdata'))

# ******************************************************************************
#                                HMM PREP DATA
# ******************************************************************************

# split roan and oryx
hmm.data.roan.1h <- data.1h %>% filter(ANIMAL == "Roan")
hmm.data.oryx.1h <- data.1h %>% filter(ANIMAL == "Oryx")
hmm.data.roan.5h <- data.5h %>% filter(ANIMAL == "Roan")
hmm.data.oryx.5h <- data.5h %>% filter(ANIMAL == "Oryx")

# ******************************************************************************
#                                PARAMETERS
# ******************************************************************************

# label states
nbStates <- 3
stateNames <- c("ruminating", "foraging", "exploring")

# distributions for observation processes
dist = list(step = "gamma", angle = "vm")

# step parameter priors
mu0.1h <- c(10, 100, 500)
sd0.1h <- c(10, 100, 500)
mu0.5h <- c(100, 500, 2000)
sd0.5h <- c(100, 500, 2000)
zm0 <- c(0.01, 0.05, 0.00005)

# angle parameter priors
angleMean0 <- c(0, pi, 0) # initial means (one for each state) 
angleCon0 <- c(0.001, 1.0, 1.0) # initial concentrations (one for each state)

#put it all together
Par0.1h = list(step = c(mu0.1h, sd0.1h, zm0), angle = c(angleMean0, angleCon0))
Par0.5h = list(step = c(mu0.5h, sd0.5h), angle = c(angleMean0, angleCon0))

# ******************************************************************************
#                             FITTING HMM - SIMPLE
# ******************************************************************************

# set up model structure
hmmModel <- function( data, pars) {
  momentuHMM::fitHMM(
    data = data, 
    nbStates = nbStates, 
    dist = dist, 
    Par0 = pars,
    estAngleMean = list(angle=TRUE), 
    stateNames = stateNames,
    formula=formula( ~ 1))
}

# fit model
tic()
m.roan.1h <- hmmModel(hmm.data.roan.1h, Par0.1h)
m.oryx.1h <- hmmModel(hmm.data.oryx.1h, Par0.1h)
m.roan.5h <- hmmModel(hmm.data.roan.5h, Par0.5h)
m.oryx.5h <- hmmModel(hmm.data.oryx.5h, Par0.5h)
toc()

# decoding states and adding to data
hmm.data.roan.1h$STATE <- stateNames[viterbi(m.roan.1h)]
hmm.data.oryx.1h$STATE <- stateNames[viterbi(m.oryx.1h)]
hmm.data.roan.5h$STATE <- stateNames[viterbi(m.roan.5h)]
hmm.data.oryx.5h$STATE <- stateNames[viterbi(m.oryx.5h)]

# save data and models
hmm.models <- list(m.roan.1h, m.oryx.1h, m.roan.5h, m.oryx.5h)
hmm.data <- list(hmm.data.roan.1h, hmm.data.oryx.1h,
                 hmm.data.roan.5h, hmm.data.oryx.5h)
names(hmm.data) <- names(hmm.models) <- c('roan_1h', 'oryx_1h', 'roan_5h', 'oryx_5h')

save(hmm.models, hmm.data, 
     file=here("03_output", "hmm", "hmm1h5hUTM.rdata"))

# ******************************************************************************
#                           FITTING HMM - WITH COVARS
# ******************************************************************************

# # # decoding states and adding to data
# roan$STATE <- stateNames[viterbi(m.roan.simp)]
# oryx$STATE <- stateNames[viterbi(m.oryx.simp)]
# 
# # adding "time since state change" to data
# calcTimeChange <- function(states) {
#   stateflip <- states == lag(states)
#   stateflip <- ifelse(is.na(stateflip), FALSE, stateflip)
#   flips <- ave(stateflip, cumsum(stateflip==0), FUN=cumsum)
#   flips
# }
# roan$TIMESINCECHANGE <- calcTimeChange(roan$STATE)
# oryx$TIMESINCECHANGE <- calcTimeChange(oryx$STATE)

# formula
form <- formula( ~ #TIMESINCECHANGE + 
                   TOD * DAYTYPE * SEASON )
                   #DIST_ART + DIST_NAT + DIST_RIV )

# fit model
m.roan <- momentuHMM::fitHMM (
  data = roan,
  nbStates = nbStates,
  dist = dist,
  Par0 = Par0,
  estAngleMean = list(angle=TRUE),
  stateNames = stateNames,
  formula=form
  )

m.oryx <- fitHMM (
  data = oryx,
  nbStates = nbStates,
  dist = dist,
  Par0 = Par0,
  estAngleMean = list(angle=TRUE),
  stateNames = stateNames,
  formula=form
  )


# ******************************************************************************
#                                 SAVING OUTPUTS
# ******************************************************************************

outfile = here('03_output', 'hmm', 'HMMOutputALL_TEMPC.RData')
save(pdata, 
     roan, oryx, 
     m.roan, m.oryx,
     m.roan.simp, m.oryx.simp, 
     file=outfile)
# load(outfile)


# ******************************************************************************
#                                PLOTTING TIME!
# ******************************************************************************

# # plot data
plot(roan, compact=T)
plot(oryx, compact=T)
# 
# # show model outputs
# m.roan
# m.oryx
# 
# # plot model outputs
plot(m.roan.simp, plotCI=TRUE)
plot(m.oryx.simp, plotCI=TRUE)
plotStates(m.roan.simp, animals="SAT1731")

# plotStates(m, animals="elk-115")
# plotStationary(m.oryx, plotCI=TRUE)

rbind(roan, oryx) %>%
  group_by(STATE) %>%
  summarize(step=mean(DIST, na.rm=TRUE),
            angle=mean(angle, na.rm=TRUE),
            n=n(),
            p=n/nrow(.))
m %>%
  ggplot() +
  geom_histogram(aes(x=angle)) +
  facet_wrap(~STATE)






