# analysesForFencingPaper.R
# Created 21 March 2024
# Margaret Swift <mes473@cornell.edu>

# One script to run all analyses for the fencing chapter publication 


# ******************************************************************************
#                         PACKAGES, UTILITIES, DATA                             
# ******************************************************************************
here::i_am('02_scripts/analysis/0_analysesForFencingPaper.R')
source(here::here("02_scripts", "utilities.R"))
pacman::p_load(#lme4, nlme, glmm, jtools, 
               sjPlot, units)
paper.theme = theme_bw() +
  theme(text=element_text(size=22), 
        panel.grid=element_blank(), 
        legend.position="bottom")

# ******************************************************************************
# LOAD & SET UP DATA FOR ANALYSES
# antelope collar data
load(antfile)
### hmm results
load(here('03_output', 'hmm', 'hmm1h5h.rdata'))
# Encounter data
load(here('03_output', 'barriers', 'encounters.RData'))
# encounters + hmm activities
load(here("03_output", "barriers", "encounters_activities_combined.RData"))
# ele crossing sites
load(here("01_data", "ele", "processed", "geographicData.RData"))

# ******************************************************************************
### TESTING RESULTS FUNCTION ###
testGLMResults <- function(mod) {
  require(DHARMa)
  # model results
  tab = tab_model(mod)
  print(tab)
  
  # check for overdispersion
  disp = testDispersion(mod)
  print(disp)
  
  ## Checking residuals
  modelresiduals = residuals(mod)
  plot(modelresiduals, pch=19)
  abline(0, 0, lty=2, col='red')
  
  #Deviance residuals
  res = residuals(mod, type = 'deviance')
  hist(res, breaks=20)
  
  # Deviance goodness-of-fit
  p.value = pchisq(mod$deviance, df=mod$df.residual, lower.tail=FALSE)
  message('chi-squared p value: ', round(p.value, 4))
  # Fail to reject; looks like a good fit.
  
  ## DHARMa check model fit
  # https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
  check_gamma_model <- simulateResiduals(fittedModel=mod, n=500)
  plot(check_gamma_model)
}










# ******************************************************************************
#                             1 - SEASONAL CHANGES
#       Are antelope more likely to encounter/cross fences during times of 
#             greater stress (end of dry season; September-November?)
# ******************************************************************************

#   Q1: Do roan cross the fence more often in the late dry season? 
#   H0: Antelope do not cross the fence more often in the late dry season. 

#   Q2: Does seasonality drive fence crossing? 
#   H0: Effect size and p-value for late-dry-season are not large/significant 
#       when determining fence crossing potential.

# ******************************************************************************
# ******************************************************************************
## SET UP DATA
## 
enc.data <- encounters %>% 
  rename_all(tolower) %>% 
  mutate(date=as.Date(gsub('[0-9]{4}', '2010', start_time)), 
         month=month(date),
         stress_season = month %in% 9:11 ) %>% 
  dplyr::select(species, animalid, event, date, start_time, month, stress_season)
enc.data.lambda = enc.data %>% 
  arrange(animalid, start_time) %>% 
  mutate(lambda = as.n(start_time - lag(start_time)),
         lambda = ifelse(animalid != lag(animalid), NA, lambda))

# ******************************************************************************
# # PLOT FUNCTION
plotEncounterTimeline <- function(spname, title=spname) {
  df <- enc.data %>% filter(species == spname)
  start = as.Date("2010-01-01"); end = as.Date('2010-12-31');
  p1 = ggplot(df, aes(x=date, fill=stress_season)) + 
    geom_histogram(bins=50, color='white') + 
    scale_x_date(date_labels="%b",
                 date_breaks="1 month",
                 limits = c(start, end)) +
    xlab('month of year') + 
    ylab('# fence encounters') +
    ggtitle(title) +
    scale_fill_manual(values=c('darkgray', 'black')) +
    paper.theme
  p2 = ggplot(enc.data.lambda %>% filter(species==spname)) + 
    geom_histogram(aes(x=lambda, fill=stress_season), color='white') + 
    scale_fill_manual(values=c('darkgray', 'black')) +
    xlim(c(0, 20000)) + xlab('hours between fence encounters') +
    ggtitle(title) +
    paper.theme
  p1 + p2
}

# ******************************************************************************
## PLOTTING
p1 <- plotEncounterTimeline('Roan', "Roan antelope")
p2 <- plotEncounterTimeline('Oryx', 'Gemsbok')
( p <- p1 / p2 + plot_layout(guides='collect') & theme(legend.position="bottom"))

outfile = here::here(outpath, "fencing_paper", 'fence_encounter_timeline.png')
ggsave(file=outfile, width=18, height=10, units="in")

# ******************************************************************************
# tabulate
test.data %>% group_by(species) %>% summarize(n=n())
test.data %>% group_by(species) %>% dplyr::select(stress_season) %>% table()
test.data %>% group_by(species, stress_season) %>% 
  summarize(mean=mean(lambda, na.rm=TRUE),
            sd=sd(lambda, na.rm=TRUE))


# ******************************************************************************
## GLM
## run the model only using relevant encounter results.
test.data = enc.data.lambda %>% filter(lambda < 20000) 
mod1 <- glm( lambda ~ stress_season + species, 
             data = test.data, 
             family = Gamma )
mod2 <- glm( lambda ~ month*species, 
             data = test.data, 
             family = Gamma )
testGLMResults(mod1)
testGLMResults(mod2)
tab_model(mod1)
tab_model(mod2)

tab2 = test.data %>% 
  group_by(species, stress_season) %>% 
  summarize(L=mean(lambda, na.rm=TRUE), n=n())
print(tab2)

# ******************************************************************************
# chi-squared to show difference in encounter type ratios

chitab <- encounters %>% group_by(SPECIES) %>% dplyr::select(EVENT) %>% table() 
chisq.test(chitab)

# running many bootstrapping
nrep=100
oryx <- encounters %>% filter(SPECIES == "Oryx")
p.vals <- vector(mode="numeric", length=nrep)
for (i in 1:nrep) {
  roan <- encounters %>% filter(SPECIES == "Roan") %>% sample_n(75)
  df <- rbind(roan, oryx) 
  tab <- df %>% group_by(SPECIES) %>% dplyr::select(EVENT) %>% table() 
  p.vals[i] <- chisq.test(tab)$p.value
}




# ******************************************************************************
#                             2 - WATER / VEG ACCESS
#           Do roan antelope cross the fence to access more resources?
# ******************************************************************************

#   Q1: Is distance-from-water significantly smaller post-crossing than pre-crossing? 
#   H0: Distance is the same before and after a fence encounter. 
#   Q2: Do distance values significantly drive fence crossing? 
#   H0: Effect size and p-value for distance are not large/significant when determining 
#       fence crossing potential.
#   Q1: Is EVI significantly greater post-crossing than pre-crossing? 
#   H0: EVI is the same before and after a fence encounter. 
#   Q2: Do EVI values significantly drive fence crossing? 
#   H0: Effect size and p-value for EVI are not large/significant when determining fence crossing potential.

#   Q3: Do roan antelope increase their foraging percentage after a fence crossing compared to before? 
#   H0: Antelope do not change their foraging percentage after a fence crossing.
#   H0: Antelope do increase their foraging, but this increase is seen after other fence encounter types, 
#       or this increase is also seen after gemsbok encounter a fence.

# To test whether there were statistically significant increases in EVI and foraging 
# percentage (more forage) or water density (more water access) after crossing the fence,
# we calculated the density of “water” pixels and median EVI for each GPS fix in 
# w hours pre- and post-crossing, with burn-in b removed. We used a GLMM with fence side (‘before’ 
# or ‘after’ crossing) as a binary response; EVI, season, and water density as 
# fixed effects; and antelope ID as a random effect.


# ******************************************************************************
## SET UP DATA

# resource use data
totals <- hmm.enc.data %>% dplyr::select(ruminating, foraging, exploring) %>% rowSums()
resource.data <- cbind(hmm.enc.data, tot=totals) %>% mutate(p.forage=foraging/tot)
resource.data$water_norm <- resource.data$water_mu %>% normalize()
resource.data$evi_norm <- resource.data$evi_mu %>% normalize()
resource.data$forage_norm <- resource.data$p.forage %>% normalize()


# ******************************************************************************
# PLOTTING FUNCTION
plotResourceUse <- function(res, name=res, title=res) {
  if (res == "forage") {
    dat <- resource.data[,c('id', 'time', 'burst', 'p.forage', 'states_sig')] %>% 
      rename(mu = p.forage, sig=states_sig) %>% 
      mutate(sd=0)
  } else {
    cols <- names(resource.data)[grepl(res, names(resource.data))]
    dat <- resource.data[,c('id', 'time', 'burst', cols)]
    names(dat) <- gsub(".*_", "", names(dat))
  }
  
  before <- dat[dat$time == 'before',]
  after  <- dat[dat$time == 'after',]
  
  # see if resource access is better after crossing
  if (grepl('water', res)) { better <- after$mu < before$mu
  } else { better <- after$mu > before$mu }
  
  dat$alpha <- ifelse(dat$sig < 0.05, 1, 0.4)
  dat$color <- rep(better, each=2)
  dat$id <- gsub("SAT", "", dat$id)
  
  totals = dat %>% 
    filter(time=="before") %>% 
    group_by(id) %>% 
    summarize(n=n()) %>% 
    arrange(n)
  sortby = rev(totals$id)
  
  # fix it up
  dat <- dat %>% mutate(id = factor(id, levels=sortby))
  
  ggplot(data=dat,
         mapping=aes(x=time, y=mu, group=burst, color=color, alpha=alpha)) + 
    geom_line() +
    geom_segment(aes(y=mu-sd, xend=time, yend=mu+sd),
                 alpha=0.07, linewidth=1) +
    geom_point(size=2) + 
    facet_wrap(~id, nrow=1) + 
    scale_color_manual(values=c('#FFA622', '#AE00A8', 'lightgray')) + 
    labs(colour="resource use is greater after crossing") +
    ylab(name) + xlab('fence crossing side') + 
    paper.theme +
    theme(axis.ticks.x = element_blank(),
          axis.line.x = element_blank()) + 
    guides(alpha="none") + ggtitle(title)
}

# ******************************************************************************
## PLOTTING
p1 <- plotResourceUse('water', 'distance from water (m)', 'distance from surface water') + 
  xlab('') + theme(legend.position="none")
p2 <- plotResourceUse('evi', 'EVI', "density of vegetation") + theme(legend.position="none")
p3 <- plotResourceUse('forage', '% of activity budget', 
                      "time spent foraging") + theme(legend.position="bottom")
( p <- p1 / p2 / p3 )

# save
outfile = here::here(outpath, "fencing_paper", 'resource_use_comparison.png')
ggsave(p, file=outfile, width=15, height=15, units="in")


# ******************************************************************************
# Tabulate

# hmm output
names(hmm.models) <- names(hmm.data) <- names(hmm.enc.data.ls)
type = "roan_1h"
my.hmm <- hmm.models[[type]]$mle
par.step <- my.hmm$step
par.angle <- my.hmm$angle

# resourc uses
resource.data %>% 
  group_by(time) %>% 
  summarize(m_evi = mean(evi_mu, na.rm=TRUE),
            sd_evi = sd(evi_mu, na.rm=TRUE),
            m_water = mean(water_mu, na.rm=TRUE),
            sd_water = sd(water_mu, na.rm=TRUE),
            m_forage = mean(p.forage, na.rm=TRUE),
            sd_forage = sd(p.forage, na.rm=TRUE)
            )


# ******************************************************************************
# ANALYSIS
## ========== RANDOM EFFECTS ARE TOO SMALL =========
# GLMM - using animal ID as a random effect.
# glmm <- glmer(time ~ water_norm + evi_norm + (1|id),
#               data=resource.use, family="binomial")
# ===================================================

# GLM - fence side as a function of EVI, water density, and season, 
#   with NO random effect
glm.resource <- glm( time ~ water_norm + evi_norm*p.forage, #normalized covariates
            data = resource.data, 
            family = binomial(link="logit"))
tab_model(glm.resource)
testGLMResults(glm.resource)

glm.forage <- glm( time ~ forage_norm, #normalized covariates
                     data = resource.data, 
                     family = binomial(link="logit"))
tab_model(glm.forage)
testGLMResults(glm.forage)


glm.no.forage <- glm( time ~ water_norm + evi_norm, #normalized covariates
                   data = resource.data, 
                   family = binomial(link="logit"))
tab_model(glm.no.forage)
testGLMResults(glm.no.forage)







# ******************************************************************************
#                           2 - FORAGING INCREASES
#       Do roan antelope cross the fence to access more vegetation?
# ******************************************************************************

#   Q3: Do roan antelope increase their foraging percentage after a fence crossing compared to before? 
#   H0: Antelope do not change their foraging percentage after a fence crossing.
#   H0: Antelope do increase their foraging, but this increase is seen after other fence encounter types, 
#       or this increase is also seen after gemsbok encounter a fence.
# 
# # ******************************************************************************
# # PLOTTING FUNCTION
# plotForagingComparison <- function(spname, title=spname) {
#   hmm.enc.data %>% 
#     mutate(n=as.n(n),
#            btype=factor(btype),
#            state=factor(state, levels=c('exploring', 'ruminating', 'foraging'))) %>% 
#     group_by(sp, state, time, btype, .drop=FALSE) %>% 
#     summarize(nsum=sum(n, na.rm=TRUE)) %>% 
#     filter(sp==spname) %>% 
#   ggplot() + 
#     geom_bar(aes(x=time, y=nsum, fill=state), 
#              stat="identity",
#              color=NA) +
#     scale_fill_manual(values=c('lightgray', 'darkgray', 'darkgreen')) +
#     facet_wrap(~btype, nrow=1) + 
#     paper.theme +
#     ylab('# time steps in state') + 
#     xlab("fence encounter time") + ggtitle(title)
# }
# # ******************************************************************************
# # PLOTTING
# p1 <- plotForagingComparison('Roan') + xlab('') + theme(legend.position="none")
# p2 <- plotForagingComparison('Oryx', 'Gemsbok')
# ( p <- p1 / p2 )  
# 
# # save it
# outfile = here::here(outpath, "fencing_paper", 'foraging_comparison.png')
# ggsave(p, file=outfile, width=18, height=10, units="in")
# 
# 
# # ******************************************************************************
# # ANALYSIS
# dat <- do.call(rbind, hmm.enc.data.ls) %>%
#   uncount(n) %>%
#   mutate(enc_type = factor(btype, levels=c('Average', 'Trace', 'Jag', 'Cross', 'Bounce')),
#          is_cross = btype == "Cross",
#          is_foraging = state=="foraging")
# 
# 
# # GLMM - AFTER a fence encounter, do they forage more?
# forage.mod1 <- glm(is_foraging ~ enc_type*sp + time,
#                     family='binomial', data=dat)
# tab_model(forage.mod1)
# testGLMResults(forage.mod1)
# 
# # GLMM - AFTER a fence crossing, do they forage more?
# forage.mod2 <- glm(is_foraging ~ is_cross*time,
#                   family='binomial', data=dat %>% filter(sp=="Roan"))
# tab_model(forage.mod2)
# testGLMResults(forage.mod2)
# 
# ## GLM - foraging as a count
# dat <- do.call(rbind, hmm.enc.data.ls)











# ******************************************************************************
#                             4 - ELEPHANT FACILITATION
#                 Do elephant crossings facilitate roan crossings?
# ******************************************************************************

#   Q1: What’s the average distance from an elephant crossing point for these roan?
#   Q2: How many roan crossings are made within 100m of an elephant crossing?
#       (100m allows for some noise from GPS jitter and hour-long fix rates)

# ******************************************************************************
## SET UP DATA
crs <- st_crs(fence)
ele.sf <- st_transform(ele, crs=crs) %>% 
  mutate(date=as.Date(DATE))

x <- st_as_sf(encounters, coords=c('easting', 'northing'), crs=crs) %>% 
  mutate(date=as.Date(start_time), 
         iscross = EVENT == "Cross",
         iscrossjag = EVENT %in% c("Cross", "Jag"))
x$distm_fence <- as.numeric(st_distance(x, fence))
hist( x$distm_fence, breaks=20 )
fence.buff <- st_buffer(fence, dist=max(x$distm))

# restrict dates
date.range <- c(max(min(x$date), min(ele.sf$date, na.rm=TRUE)), 
                min(max(x$date), max(ele.sf$date, na.rm=TRUE)))
windows <- seq(date.range[1], date.range[2], length.out=20)

# distance from ele crossings -- make sure crossing exists
y <- x %>% filter(date %within% interval(date.range[1], date.range[2]))
y$distm_ele <- NA
for (i in 2:length(windows)) {
  w <- windows[i]
  w0 <- windows[i-1]
  print(i)
  print(w)
  ele.sub <- ele.sf %>% filter(date <= w)
  y.inx <- (y$date <= w) + (y$date >= w0) == 2
  y.sub <- y[y.inx,]
    
  # distance from existing ele crossings
  if (nrow(y.sub)) {
    dists <- st_distance(y.sub, ele.sf)
    Min <- function(e) min(e, na.rm=TRUE)
    y$distm_ele[y.inx] <- apply(dists, 1, Min)
  } else { message('none')}
}

# ******************************************************************************
# plotting
fplot<- ggplot() + 
  geom_sf(data=fence.buff, fill='lightgreen', color=NA, alpha=0.4) + 
  geom_sf(data=fence) + 
  geom_sf(data=ele.sf, color='red', size=3, 
          shape=3, alpha=0.2) + 
  geom_sf(data=x, mapping=aes(color=iscross),
          size=3, #color="purple", 
          alpha=0.5) + 
  # geom_sf(data=point[2,], color='green', size=6, 
          # shape=0) +
  coord_sf(ylim=c(-19.01, -18.3),
           xlim=c(20.99, 21.1)) + 
  scale_color_manual(values=c('gray', 'purple')) +
  paper.theme + theme(axis.text.x=element_blank())+ 
  theme(legend.position="none")
fplot

# subplot location
ylim <- c(-18.61, -18.69)
xlim <- c(20.992, 21.005)
point <- st_as_sf(data.frame(x=xlim, 
                             y=ylim), 
                  coords=c('x', 'y'),
                  crs=st_crs(fence))
fplot_sub = ggplot() + 
  geom_sf(data=fence.buff, fill='lightgreen', alpha=0.4) + 
  # geom_sf(data=fence, linetype='dashed', linewidth=1) +
  geom_sf(data=ele.sf, color='red', size=4, 
          shape=3) + 
  geom_sf(data=y, #color='purple',
          mapping=aes(color=iscross), 
          size=4, alpha=0.5) + 
  scale_color_manual(values=c('gray', 'purple')) +
  scale_x_continuous(breaks=seq(20.992, 21.005, length.out=3)) +
  coord_sf(ylim=ylim, xlim=xlim) + 
  paper.theme + theme(axis.text.x=element_blank()) + 
  theme(legend.position="none")
fplot_sub

## HISTOGRAMS ##
plotHistogramEleCrossings <- function(events, fill='gray', title) {
  hlims <- xlim(c(0, 800))
  xlab = "distance (m) from nearest elephant crossing site"
  hdata <- nog(y) %>% filter(EVENT %in% events) %>% mutate()
  mean <- mean(hdata$distm_ele)
  sd <- sd(hdata$distm_ele)
  my_hist <- ggplot(hdata) +
    geom_histogram(aes(x=distm_ele, fill=iscross), bins=40, 
                   color='black', fill=fill) + 
    geom_vline(xintercept=mean, color="red", linetype='dotted', linewidth=1) +
    geom_vline(xintercept=mean+sd, color="gray", linetype='dashed') +
    geom_vline(xintercept=mean-sd, color="gray", linetype='dashed') +
    paper.theme + hlims + ylab('frequency') +
    ggtitle(title) + 
    xlab(xlab) 
  my_hist
}

events <- c("Jag", "Trace", "Bounce", "Average")
h1 <- plotHistogramEleCrossings("Cross", fill="purple", title="antelope crossings")
h2 <- plotHistogramEleCrossings(events, fill="gray", title="antelope non-crossings")
h1 / h2

# save it
fplot + fplot_sub + (h1 / h2) + plot_layout(nrow=1)
outfile = here::here(outpath, "fencing_paper", 'ele_crossing_map.png')
ggsave(file=outfile, width=18, height=10, units="in")

# ******************************************************************************
# TABULATE
y %>% nog() %>% group_by(SPECIES) %>% 
  summarize(m=mean(distm_ele, na.rm=TRUE),
            sd=sd(distm_ele, na.rm=TRUE))
y %>% nog() %>% filter(SPECIES == "Roan") %>% 
  group_by(iscross) %>% 
  summarize(m=mean(distm_ele, na.rm=TRUE),
            sd=sd(distm_ele, na.rm=TRUE))
y %>% nog() %>% group_by(SPECIES, iscrossjag) %>% 
  summarize(m=mean(distm_ele, na.rm=TRUE),
            sd=sd(distm_ele, na.rm=TRUE))

# ******************************************************************************
## ANALYSIS
glm.crossjag <- glm( iscrossjag ~ distm_ele + SPECIES, 
            data = nog(y), 
            family = "binomial" )
tab_model(glm.crossjag)
testGLMResults(glm.crossjag)

glm.cross <- glm( iscross ~ distm_ele, 
                     data = nog(y) %>% filter(SPECIES == "Roan"), 
                     family = "binomial" )
tab_model(glm.cross)
testGLMResults(glm.cross)