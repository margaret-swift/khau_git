# testingWindowSize.R
# Created 30 April 2024
# Margaret Swift <mes473@cornell.edu>

# To test whether there were statistically significant increases in resource
# access after crossing the border fence, we need to decide what window size
# is the most informative. This script compares different combinations of before
# and after crossing window sizes and burnin (throwaway) periods immediately
# before and after the crossing. This will allow us to remove some of the 
# spatial autocorrelation inherent to comparing habitats before and after an event.


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
here::i_am('02_scripts/analysis/testingWindowSize.R')
source(here::here("02_scripts", "utilities.R"))
pacman::p_load(lme4, nlme, glmm, jtools, tictoc, reshape2, sjPlot, units)
load(here('03_output', 'barriers', 'encounters.RData'))


# ******************************************************************************
#                               COMPARING WINDOWS
# ******************************************************************************

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

for (i in 1:NW) {
  tic()
  row <- windows[i,]
  message(paste(row, collapse='-'))
  for (event in events) {
    cat(event, '-')
    r <- pullBufferCovariateData(encounters, collar.df, 
                                 eventlist=event, 
                                 nhours=row$hours, 
                                 burnin=row$burnin)
    
    # get before/after data
    after <- r[r$time=="after",covars] 
    before <-r[r$time=="before",covars] 
    
    # calculate correlations
    evi.stat <- cor(before$evi, after$evi)
    water.stat <- cor(before$water, after$water)

    # save
    my_col <- cols[grepl(event,cols)]
    covar.dat[i,my_col] <- c(evi.stat, water.stat)
  }
  toc()
}


# ******************************************************************************
#                              PLOTTING RESULTS
# ******************************************************************************

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
p <- ggplot(data=agg.data,
       mapping=aes(x=hours, color=burnin, group=burnin), 
       linewidth=1) + 
  geom_line(mapping=aes(y=correlation), linewidth=1) + 
  facet_wrap(~covar+event, scales="free_x") + big.theme + 
  scale_color_brewer(palette='Greens', direction=1, name="burn-in\n(hours)") + 
  theme_nice() + 
  xlab('covariate averaging window (hours)') + 
  ylab('correlation between values before and \nafter fence encounter')

# save it
ggsave(p, filename=here('03_output', "windowTesting.png"))


# ******************************************************************************
#                                     EOF
# ******************************************************************************