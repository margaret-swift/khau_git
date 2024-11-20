# analyzeActivitiesResources.R
# Created 24 April 2023
# Margaret Swift <mes114@duke.edu>

#   - WHAT ARE ROAN DOING BEFORE AND AFTER ENCOUNTERS?


# ******************************************************************************
#                               DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here("02_scripts", "utilities.R"))
load(antfile)
load(here("03_output", "barriers", "encounters_activitiesFINAL.RData"))

# ******************************************************************************
#                                   FUNCTIONS
# ******************************************************************************

plotProportions <- function(df, title="") {
  df %>% 
    mutate(n=as.n(n)) %>% 
    group_by(state, time, btype) %>% 
    summarize(nsum=sum(n, na.rm=TRUE)) %>% 
    ggplot() + 
    geom_bar(aes(x=btype, y=nsum, fill=state), 
             position="fill",
             stat="identity",
             color='black') +
    scale_fill_brewer(palette="Accent") +
    facet_wrap(~time) + 
    theme(text=element_text(size=18)) + 
    ylab('proportion of activities') + 
    xlab("movement type") +
    ggtitle(title)
}

# chi squared etc function
chiSqMe <- function(df, t) {
  uc <- uncount(df, n) %>% filter(time == t)
  tab <- table( uc$btype, uc$state )
  print( tab )
  print( chisq.test( tab ) )
}

# ******************************************************************************
#                 WHAT ARE ROAN DOING BEFORE AND AFTER ENCOUNTERS?
# ******************************************************************************

# ROAN
df.roan <- hmm.enc.data %>% filter(sp=="Roan")
plotProportions(df.roan, "Roan")
# ggsave(here("03_output", "barriers", "figs", "roan_before_after_enc.png"),
#        width=13.1, height=6.19)
chiSqMe(df.roan, "before")
chiSqMe(df.roan, "after")

df.roan %>% 
  filter(time=="before") %>% 
  mutate(is.cross = btype%in% c("Cross")) %>% 
  dplyr::select(burst, state, is.cross, n) %>% 
  group_by(is.cross) %>% 
  mutate(total = sum(n)) %>% 
  ungroup() %>% 
  group_by(is.cross, state) %>% 
  summarize(n=sum(n),
            total=first(total)) %>% 
  mutate(p=n/total)

# ORYX
plotProportions(df.oryx, "Oryx")
chiSqMe(df.oryx, "before")
chiSqMe(df.oryx, "after")

# post-hoc series of chisquare tests to show which movement type has 
# a significantly different proportion of states than the others in "before"
uc <- uncount(df.roan, n) %>% filter(time == 'before')
tab <- table( uc$btype, uc$state )
# adjusted bonferroni p-value: 0.05/10 = 0.005
chisq.test(tab[c(1,2),])$p.value
chisq.test(tab[c(1,3),])$p.value #
chisq.test(tab[c(1,4),])$p.value
chisq.test(tab[c(1,5),])$p.value
chisq.test(tab[c(2,3),])$p.value #
chisq.test(tab[c(2,4),])$p.value
chisq.test(tab[c(2,5),])$p.value
chisq.test(tab[c(3,4),])$p.value #
chisq.test(tab[c(3,5),])$p.value #
chisq.test(tab[c(4,5),])$p.value 

#kruskal-wallis test
kruskal.test(state ~ btype, data = uc)




# let's just be sure Gemsbok aren't doing anything weird
plotProportions(df.oryx, "Oryx")
chiSqMe(df.oryx, "before")
chiSqMe(df.oryx, "after")





# ******************************************************************************
#                                   GLMM
# ******************************************************************************

pacman::p_load(mixcat, brms, R2jags, posterior, lme4)

totals <- df.roan %>% 
  group_by(burst, time) %>% 
  summarize(ntot=sum(n))
df <- left_join(df.roan, totals, by=c('burst', 'time')) %>% 
  mutate(id=factor(id),
         p = n / ntot,
         state=factor(state, levels=c('ruminating', 'foraging', 'exploring')),
         time=factor(time, levels=c('before', 'after'))) %>% 
  dplyr::select(id, season, state, btype, time, n, ntot, p)
df.exp <- df %>% filter(state=='exploring')
  
# model 1: no random effect, but quasibinomial
mod1 <- glm(n/ntot ~ btype*time, #(1|id),
             weights=ntot,
             data=df.exp,
           family="quasibinomial")
performance::check_overdispersion(mod1)
summary(mod1)

# model 2: random effect, but binomial (overdispersed)
mod2 <- glmer(n/ntot ~ btype*time + (1|id),
           weights=ntot,
           data=df.exp,
           family="binomial")
performance::check_overdispersion(mod2)
summary(mod2)

# model 3: brms
df.expanded <- df[rep(row.names(df), df$n),] %>% 
  dplyr::select(id, season, state, btype, time)
attach(df.expanded)
bayes.brms <- brm(state ~ btype*time + (1 | id),
                  family = categorical("logit"),
                  data = df.expanded,
                  chains = 2, # nb of chains
                  iter = 5000, # nb of iterations, including burnin
                  warmup = 1000, # burnin
                  thin = 1) # thinning
detach(df)
# plot(bayes.brms)
summary(bayes.brms)

library(reshape2)
library(tidyverse)
df <- data.frame(sex=c('male', 'female'),
                 river = c(14.5, 10.1),
                 road = c(25.8, 15.3),
                 fence = c(0, 3.5)) %>% 
  melt(id.var="sex")
ggplot(df, aes(x=variable, y=value*0.01, fill=sex)) + 
  geom_bar(stat='identity', position='dodge') + 
  theme(text=element_text(size=16)) + 
  xlab('linear feature') + ylab('chance of crossing') + 
  scale_fill_manual(values=c('darkgray', 'black'))
  
