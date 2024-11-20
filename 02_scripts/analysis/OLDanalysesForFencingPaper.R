# analysesForFencingPaper.R
# Created 21 March 2024
# Margaret Swift <mes473@cornell.edu>

# One script to run all analyses for the fencing chapter publication 


# ******************************************************************************
#                           DATA & LIBRARY LOADING
# ******************************************************************************
here::i_am('02_scripts/analysis/analysesForFencingPaper.R')
source(here::here("02_scripts", "utilities.R"))
pacman::p_load(lme4, nlme, glmm, jtools, reshape2, sjPlot, units)
load(antfile)
no.y.theme <- theme(axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.title.y=element_blank())

# ******************************************************************************
#                         1 - SEASONAL MOVEMENT INCREASES
# ******************************************************************************
# We tested for interspecific differences in seasonal movement increases using a 
# GLMM with each individual-season hourly distance as the dependent variable. 
# Each year was taken as a replicate, with season and species as fixed effects 
# antelope collar ID as a random effect.
load(here('03_output', 'barriers', 'range_hulls.RData'))
test.hulls <- hulls %>% nog() %>% 
  mutate(LABEL=paste(COLLAR_ID, YEAR, SEASON, sep="_")) %>% 
  mutate(AREA = set_units(set_units(AREA, 'm2'), 'ha'))
test.hulls.k <- hulls.k %>% nog() %>% ungroup() %>% 
  mutate(LABEL=paste(COLLAR_ID, YEAR, SEASON, sep="_")) %>% 
  mutate(AREA.K = set_units(AREA, 'ha')) %>% 
  dplyr::select(LABEL, AREA.K)
# buffs
area.b <- st_area(buffs)
test.buffs <- buffs %>% nog() %>% 
  mutate(COLLAR_ID = as.n(gsub('SAT', '', COLLAR_ID))) %>% 
  mutate(LABEL=paste0("SAT", COLLAR_ID, '_', YEAR, '_', SEASON)) %>% 
  left_join(collar.meta %>% dplyr::select(SEX, ANIMAL, COLLAR_ID), by="COLLAR_ID")
test.buffs$AREA.B = area.b
test.buffs <- test.buffs %>% ungroup() %>% dplyr::select(LABEL, AREA.B)


test <- left_join(test.hulls, test.hulls.k, by='LABEL') %>% 
  left_join(test.buffs, by='LABEL') %>% 
  filter(SEX=="F", !ENC_FENCE) %>% ungroup() %>% 
  dplyr::select(COLLAR_ID, SPECIES, SEASON, YEAR, AREA, AREA.K, AREA.B)

# LMM - D.V.by species
m0k<-lmer(AREA.K~ SPECIES + (1|COLLAR_ID), data=test)
# m0a<-lmer(AREA~ SPECIES + (1|COLLAR_ID), data=test)
# m0b<-lmer(AREA.B~ SPECIES + (1|COLLAR_ID), data=test)
m0d<-lmer(DIST_WATER ~ SPECIES + (1|COLLAR_ID), data=collar.df)

# LMM - D.V. by season
m1k<-lmer(AREA.K~ SEASON + (1|COLLAR_ID), data=test)
# m1a<-lmer(AREA~ SEASON + (1|COLLAR_ID), data=test)
# m1b<-lmer(AREA.B~ SEASON + (1|COLLAR_ID), data=test)
m1d<-lmer(DIST_WATER ~ SEASON + (1|COLLAR_ID), data=collar.df)

# LMM - D.V. by season + species
m2k<-lmer(AREA.K~ SEASON + SPECIES + (1|COLLAR_ID), data=test)
# m2a<-lmer(AREA~ SEASON + SPECIES + (1|COLLAR_ID), data=test)
# m2b<-lmer(AREA.B~ SEASON + SPECIES + (1|COLLAR_ID), data=test)
m2d<-lmer(DIST_WATER ~ SEASON + SPECIES + (1|COLLAR_ID), data=collar.df)

# LMM - D.V. by season * species
m3k<-lmer(AREA.K~ SEASON*SPECIES + (1|COLLAR_ID), data=test)
# m3a<-lmer(AREA~ SEASON*SPECIES + (1|COLLAR_ID), data=test)
# m3b<-lmer(AREA.B~ SEASON*SPECIES + (1|COLLAR_ID), data=test)
m3d<-lmer(DIST_WATER ~ SEASON*SPECIES + (1|COLLAR_ID), data=collar.df)

# LMM - Species modeled separately
r <- test %>% filter(SPECIES == "Roan", !is.na(AREA.K)) %>% mutate(AREA.K = as.n(AREA.K))
mod.r.area<-lmer(AREA.K~ SEASON + (1|COLLAR_ID), data=r)
mod.r.d<-lmer(DIST_WATER~ SEASON + (1|COLLAR_ID), data=collar.df %>% filter(SPECIES=="Roan"))
g <- test %>% filter(SPECIES != "Roan", !is.na(AREA.K)) %>% mutate(AREA.K = as.n(AREA.K))
mod.o.d<-lmer(DIST_WATER~ SEASON + (1|COLLAR_ID), data=collar.df %>% filter(SPECIES!="Roan"))
mod.o.area<-lmer(AREA.K~ SEASON + (1|COLLAR_ID), data=g)


# testing things out
tab_model(m1k, m2k, m3k)
tab_model(m1d, m2d, m3d)
tab_model(mod.o.area, mod.r.area)
tab_model(mod.o.d, mod.r.d)

# plots
test %>% 
  mutate(SPECIES = ifelse(SPECIES == "Roan", "Roan antelope", "Gemsbok")) %>% 
  ggplot() + 
  geom_boxplot(aes(x=SPECIES, y=AREA.K, fill=SEASON), 
               alpha=0.8) + 
  scale_fill_brewer(palette="Dark2", direction=-1) + 
  theme(text=element_text(size=16))

collar.df %>% 
  mutate(SPECIES = ifelse(SPECIES == "Roan", "Roan antelope", "Gemsbok")) %>% 
  ggplot() + 
  geom_boxplot(aes(x=SPECIES, y=DIST_WATER, fill=SEASON), 
               alpha=0.8) + 
  scale_fill_brewer(palette="Dark2", direction=-1) + 
  theme(text=element_text(size=16)) + 
  ylab('distance from water (m)')





# ******************************************************************************
#              2 - CATEGORIZING ROAN AND GEMSBOK FENCE ENCOUNTERS
# ******************************************************************************
# For each individual, we calculated the number of fence encounters per total 
# number of fixes (encounter-per-fix), and the number of fence crossings per 
# total number of encounters (cross-per-encounter). We then used two GLMMs, 
# one with encounters-per-fix as a binomially distributed response, and one 
# with crosses-per-encounter as a binomially distributed response, to determine 
# whether differences in these two responses were driven significantly by either 
# season (wet and dry) or species (roan antelope versus gemsbok).
load(here("03_output", "barriers", "encounter_data.rdata"))
edata <- encounter.data %>% 
  filter(LON > 20.9, LON < 21.1)

# GLMM - encounters per fix
# THIS ONE TAKES A WHILE!
mod.epf <- glmer(IS_ENC ~ SPECIES*SEASON + (1|COLLAR_ID),
                 family=binomial, data=edata)
tab_model(mod.epf)

c.keep <- edata %>% 
  group_by(SPECIES, SEASON, COLLAR_ID) %>% 
  summarize(nfix = n()) %>% 
  mutate(keep = nfix > 500) %>% 
  group_by(SPECIES, COLLAR_ID) %>% 
  summarize(s=sum(keep)) %>% 
  filter(s == 2)

pdata <- edata %>% 
  filter(COLLAR_ID %in% c.keep$COLLAR_ID) %>% 
  group_by(SPECIES, SEASON, COLLAR_ID) %>% 
  summarize(nfix = n(),
            n=sum(IS_ENC==TRUE),
            epf=n/nfix) 
ggplot(pdata) + 
  geom_bar(aes(x=SPECIES, y=epf, fill=SEASON), 
           stat='identity', position='dodge') +
  scale_fill_brewer(palette="Dark2", direction=-1) + 
  theme(text=element_text(size=16)) + 
  ylab('number of encounters per GPS fixes') 
ggplot(pdata) + 
  geom_bar(aes(x=COLLAR_ID, y=epf, fill=SEASON), 
           stat='identity', position='dodge') +
  scale_fill_brewer(palette="Dark2", direction=-1) + 
  theme(text=element_text(size=16)) + 
  ylab('number of encounters per GPS fixes') + 
  facet_wrap(~SPECIES, nrow=2, scales='free_x')

# GLMM - crosses per encounter
mod.cpe <- glmer(IS_CROSS ~ SEASON + (1|COLLAR_ID), 
                 family=binomial, 
                 data=encounter.data %>% filter(IS_ENC))


encounter.data %>% 
  filter(IS_ENC) %>% 
  group_by(SEASON, COLLAR_ID) %>% 
  summarize(nenc = n(),
            n=sum(IS_CROSS==TRUE),
            cpe=n/nenc) %>% 
  ggplot() + 
    geom_bar(aes(x=COLLAR_ID, y=cpe, fill=SEASON), 
             stat='identity', position='dodge') +
    scale_fill_brewer(palette="Dark2", direction=-1) + 
    theme(text=element_text(size=16)) + 
    ylab('number of crosses per fence encounter') 
tab_model(mod.cpe)

# ******************************************************************************
#              3 - RESOURCE ACCESS BEFORE AND AFTER CROSSING
# ******************************************************************************
# To test whether there were statistically significant increases in EVI (more 
# forage) or water density (more water access) after crossing the border fence,
# we calculated the density of “water” pixels and median EVI for each GPS fix in 
# the 24 hours pre- and post-crossing. We used a GLMM with fence side (‘before’ 
# or ‘after’ crossing) as a binary response; EVI, season, and water density as 
# fixed effects; and antelope ID as a random effect.
load(here('03_output', 'barriers', 'encounters.RData'))

# Let's have a plot where we see the difference in covariate values for 
#  before/after, with various window and burnin sizes.
nhours <- seq(24, 24*7, by=8); NH = length(nhours);
nburns <- seq(0, 48, by=12); NB = length(nburns);
windows <- data.frame(hours  = rep(nhours, each=NB),
                      burnin = rep(nburns, times=NH)) %>% 
  filter(hours > burnin)
NW <- nrow(windows)
covars <- c('evi', 'water')
events <- c('Cross', 'Jag')
cols <- paste(covars, rep(events, each=2), sep="_")
covar.dat <- data.frame(matrix(0, nrow=NW, ncol=length(covars)*2))
names(covar.dat) <- cols
time.names <- factor(c('before', 'after'), levels=c('before', 'after'))

pullBufferCovariateData <- function(enc, collar.df, eventlist=NULL, 
                                    nhours=24, burnin=6) {
  
  getCovariates <- function(collar_i, range) {
    collar_i_sub <- collar_i[ collar_i$DATETIME %within% range, ]
    evi.mean = mean(collar_i_sub$EVI, na.rm=TRUE)
    water.mean = mean(collar_i_sub$DIST_NAT, na.rm=TRUE)
    return(list(evi_mu=evi.mean, water_mu=water.mean))
  }
  
  th=hours(nhours); tb=hours(burnin); nrep=2;
  NE <- nrow(enc)
  cols <- c('evi_mu', 'water_mu', 'time', 'side', 
            'season', 'id', 'sex', 'species')
  dat <- data.frame(matrix(data=NA, nrow=NE*nrep, ncol=length(cols)))
  names(dat) <- cols
  dat$time <- time.names
  if (!is.null(eventlist)) enc = enc %>% filter(EVENT %in% eventlist)
  for (i in 1:NE) {
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
    before <- getCovariates(collar_i, r0)
    after  <- getCovariates(collar_i, r1)
    
    # save to main data
    j <- (nrep*(i-1))+1
    inx <- j:(j+nrep-1)
    dat$evi_mu[inx] <- c(before$evi_mu, after$evi_mu)
    dat$water_mu[inx] <- c(before$water_mu, after$water_mu)
    dat$side[inx] <- c(enc_i$start_side, enc_i$end_side)
    
    # add encounter info
    dat$season[inx] = enc_i$SEASON
    dat$id[inx] = enc_i$AnimalID
    dat$sex[inx] = enc_i$SEX
    dat$species[inx] = enc_i$SPECIES
  }
  
  dat <- dat   %>% 
    mutate(water=ifelse(water_mu == -99, NA, water_mu), 
           evi=ifelse(water_mu == -99, NA, evi_mu)) %>% 
    mutate_at(covars, ~(scale(.) %>% as.vector))
  
  dat$INX <- rep(1:(nrow(dat)/2), each=2)
  INX.rm <- c(dat$INX[is.na(dat$water)], dat$INX[is.na(dat$evi)])
  dat <- dat[!(dat$INX %in% INX.rm),]
  
  dat <- dat %>% dplyr::select(id, sex, species, 
                               time, side, water, evi, season)
  return(dat)
}

for (i in 1:NW) {
  row <- windows[i,]
  message(paste(row, collapse='-'))
  for (event in events) {
    cat(event, '-')
    r <- pullBufferCovariateData(encounters, collar.df, 
                                 eventlist=event, 
                                 nhours=row$hours, 
                                 burnin=row$burnin)
    
    after <- r[r$time=="after",covars] 
    before <-r[r$time=="before",covars] 
    evi.stat <- cor(before$evi, after$evi)
    water.stat <- cor(before$water, after$water)
    my_col <- cols[grepl(event,cols)]
    covar.dat[i,my_col] <- c(evi.stat, water.stat)
  }
}

# fix it up
agg.data <- cbind(windows, covar.dat) %>% 
  filter(burnin %in% c(0, 12, 24, 36, 48)) %>% 
  reshape2::melt(id.vars=c('hours', 'burnin'), 
                 value.name='correlation', 
                 variable.name='covariate') %>% 
  mutate(covar = toupper(gsub('_.*', '', covariate)),
         event = toupper(gsub('.*_', '', covariate)),
         burnin= factor(burnin))

# plot it
ggplot(data=agg.data,
       mapping=aes(x=hours, color=burnin, group=burnin), 
       linewidth=1) + 
  geom_line(mapping=aes(y=correlation), linewidth=1) + 
  facet_wrap(~covar+event, scales="free_x") + big.theme + 
  scale_color_brewer(palette='Greens', direction=1, name="burn-in\n(hours)") + 
  theme_nice() + 
  xlab('covariate averaging window (hours)') + 
  ylab('correlation between values before and \nafter fence encounter')

# now with the chosen nhours and burnin buffer
resource.use <- pullBufferCovariateData(encounters, collar.df, 
                                        eventlist=c('Cross'),
                                        nhours=100, burnin=48)

# GLMM - fence side as a function of EVI, water density, and season
## boundary (singular) fit: see help('isSingular')
# glmm <- glmer(time ~ water + evi + season + (1 | id), 
#                       data=resource.use,
#                       family="binomial")
# tab_model(glmm)

# GLM - fence side as a function of EVI, water density, and season
glm <- glm( time ~ water + evi + season, 
            data=resource.use, 
            family="binomial")
tab_model(glm)

# stats
resource.use %>% 
  group_by(season, time) %>% 
  summarize(evi.mu = mean(evi, na.rm=TRUE),
            evi.sd = sd(evi, na.rm=TRUE),
            water.mu = mean(water, na.rm=TRUE),
            water.sd = sd(water, na.rm=TRUE))


# ******************************************************************************
#                   4 - ACTIVITY BUDGETS BEFORE CROSSING
# ******************************************************************************
# To determine whether activity budgets changed in the 24 hours before a cross, 
# when compared to other fence encounter types, we used a GLMM with a binomial 
# response. A “success” occurred when a GPS point was classed as an energy-
# spending activity (exploring) and a “failure” when classed as an energy-
# gaining activity (ruminating and foraging). The type of fence encounter 
# (bounce, cross, jag, average, and trace) was used as a categorical fixed 
# effect in the GLMM, interacting with species as a fixed effect, and animal 
# ID as a random effect. 

load(here("03_output", "hmm", "hmm1h5hUTM.rdata"))
load(here("03_output", "barriers", "encounters_activities_combined.rdata"))

# stats for step lengths and turning angles
rbind(hmm.data[['roan_1h']], hmm.data[['oryx_1h']]) %>% 
  group_by(STATE) %>% 
  summarize(step.mu=mean(step, na.rm=TRUE),
            step.sd=sd(step, na.rm=TRUE),
            an.mu=mean(angle, na.rm=TRUE),
            an.sd=sd(angle, na.rm=TRUE))
rbind(hmm.data[['roan_15']], hmm.data[['oryx_5h']]) %>% 
  group_by(STATE) %>% 
  summarize(step.mu=mean(step, na.rm=TRUE),
            step.sd=sd(step, na.rm=TRUE),
            an.mu=mean(angle, na.rm=TRUE),
            an.sd=sd(angle, na.rm=TRUE))

# GLMM - 24 hours BEFORE an encounter
# do.call(rbind, hmm.enc.data.ls) %>% filter(time=="before") %>% 
#   uncount(n) %>% 
#   mutate(is_bounce = btype == "Bounce",
#          enc_type = factor(btype, levels=c('Average', 'Trace', 'Jag', 'Cross', 'Bounce')),
#          energy_using = ifelse(state!="exploring", "not exploring", state),
#          energy_using = factor(energy_using, levels=c('not exploring', 'exploring'))
#   )
before <- do.call(rbind, hmm.enc.data.ls) %>% filter(time=="before") 
totals <- before %>% group_by(burst) %>% summarize(tot=sum(n))
before.data <- before %>% 
  filter(state == "exploring", sp=="Roan") %>% 
  left_join(totals, by='burst') %>% 
  mutate(prop.explore=n/tot,
         enc_type = factor(btype, levels=c('Average', 'Trace', 'Jag', 'Cross', 'Bounce')),
         is_cross = btype == "Cross")
before.mod <- glmer(is_cross ~ prop.explore + (1|id), 
                    family=binomial, 
                    data=before.data)
tab_model(before.mod)

# stats for fence encounter activity budget portions
giveStats <- function(dat) {
  tots <- dat %>% 
    group_by(enc_type, sp) %>% 
    summarize(ntot=n())
  stat <- dat %>% 
    filter(sp == "Roan") %>%
    group_by(sp, enc_type, energy_using) %>% 
    summarize(n=n()) %>%
    left_join(tots, by=c('enc_type', 'sp')) %>%
    mutate(p=n/ntot,
           label=paste0(round(p,2), 
                        ifelse(enc_type=="Cross", "*", "")))
  print(stat)
  p <- ggplot(data=stat, mapping=aes(x=enc_type, y=p)) + 
    geom_bar(mapping=aes(fill=energy_using), 
             stat="identity", position="stack") + 
    theme(text=element_text(size=18)) +
    ylab('activity budget proportion \n for 24 hours before encounter') +
    xlab('fence encounter type') +
    geom_text(data = stat %>% filter(energy_using == "exploring"),
              mapping=aes(label=label), 
              vjust = -1.0, 
              size=6,
              color = "black") +
    scale_fill_manual(values=c('darkgray', 'black'))
  p
}


# plain jane
p0 = do.call(rbind, hmm.data) %>% 
  # sample_n(10000) %>% 
  mutate(energy_using = ifelse(STATE!="exploring", "not exploring", STATE),
         energy_using = factor(energy_using, levels=c('not exploring', 'exploring')),
         enc_type="none", sp=ANIMAL) %>% 
  giveStats() + xlab('') + no.y.theme
  
# 24 hours before
pb = giveStats(before.data %>% 
                 mutate(energy_using = ifelse(state!="exploring", "not exploring", state),
                        energy_using = factor(energy_using, levels=c('not exploring', 
                                                                                   'exploring'))
))

# put them together
pb + p0 + patchwork::plot_layout(widths=c(4,1), guides="collect")
outfile = here::here(outpath, "barriers", 'results', 'prob_explore.png')
ggsave(filename=outfile, width=12, height=8)
