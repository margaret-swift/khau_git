# analysesForFencingPaper.R
# Created 21 March 2024
# Margaret Swift <mes473@cornell.edu>

# One script to run all analyses for the fencing chapter publication 


# ******************************************************************************
#                           DATA & LIBRARY LOADING
# ******************************************************************************
here::i_am('02_scripts/analysis/0_statsForFencingPaper.R')
source(here::here("02_scripts", "utilities.R"))
pacman::p_load(lme4, nlme, glmm, jtools, reshape2, sjPlot, units)
load(antfile)

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
  mutate(SPECIES = ifelse(SPECIES=="Oryx", "Gemsbok", "Roan")) %>% 
  filter(LON > 20.9, LON < 21.1)

# redefining seasons
edata$IS_LATE_DRY <- month(edata$DATETIME) %in% 9:11

# GLMM - Is the fix an encounter?
# THIS ONE TAKES A WHILE!
mod.epf <- glmer(IS_ENC ~ SPECIES*IS_LATE_DRY + (1|COLLAR_ID),
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
  group_by(SPECIES, IS_LATE_DRY, COLLAR_ID) %>% 
  summarize(nfix = n(),
            n=sum(IS_ENC==TRUE),
            epf=n/nfix) 
ggplot(pdata) + 
  geom_bar(aes(x=SPECIES, y=epf, fill=IS_LATE_DRY), 
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

# GLMM - is the encounter a cross?
mod.cpe <- glmer(IS_CROSS ~ IS_LATE_DRY + (1|COLLAR_ID), 
                 family=binomial, 
                 data=edata %>% filter(IS_ENC))
tab_model(mod.cpe)

edata %>% 
  filter(IS_ENC) %>% 
  group_by(IS_LATE_DRY, COLLAR_ID) %>%  #SEASON
  summarize(nenc = n(),
            n=sum(IS_CROSS==TRUE),
            cpe=n/nenc) %>% 
  ggplot() + 
  geom_bar(aes(x=COLLAR_ID, y=cpe, fill=IS_LATE_DRY), 
           stat='identity', position='dodge') +
  scale_fill_brewer(palette="Dark2", direction=-1) + 
  theme(text=element_text(size=16)) + 
  ylab('number of crosses per fence encounter') 
