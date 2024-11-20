# plotWaterDists.R
# Created 22 April 2023
# Margaret Swift <mes114@duke.edu>


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
load(here::here("03_output", "waterdists", "allModelsDENSITY.RData"))
load(antfile)


# ******************************************************************************
#                             PLOT RELIANCE BY MONTH
# ******************************************************************************

LV <- data.frame(LV=c('seasonal','permanent'))
LVS <- c("WATER_DENS_300m", "DIST_ART")
row.names(LV) <- LVS

# set up the data
data <- collar.df.1h %>% 
  filter( IN_BORDER ) %>% 
  nog() %>% 
  mutate(YEAR = ifelse(as.n(MONTH)<5, 2017, 2016),
         DATEX = as_date(paste(YEAR, MONTH, "01", sep="-")),
         DATEX = as.Date(paste(DATEX, "00"), format="%Y-%m-%d %H")
  ) %>% 
  group_by(SPECIES, MONTH, DATEX)
data.m <- data %>% 
  summarize_at(.vars=c('DIST_ART', 'WATER_DENS_300m'), .funs='mean') %>% 
  melt(measure.vars=row.names(LV), value.name="D", variable.name="WATER_TYPE") %>% 
  mutate(WATER_TYPE = factor(LV[WATER_TYPE,"LV"], levels=LV$LV))

# the plot
ggplot() + 
  geom_line(data=data.m, aes(x=DATEX, y=D, 
                             group=SPECIES, linetype=SPECIES), 
            linewidth=1) +
  facet_wrap(~ WATER_TYPE, scales="free", nrow=3) +
  xlab('month') + ylab('average monthly\ndistance from water (m)') +
  scale_x_date(date_breaks="1 month", date_labels="%m") +
  plot.theme + 
  guides(linetype="none")

outfile <- here::here("03_output", "waterdists", "waterDistLinesMonthly.png")
ggsave(file=outfile, width=7, height=7)


# ******************************************************************************
#                                 PLOT MODEL RESULTS
# ******************************************************************************

createLinesTemp <- function(model, data) {
  nrep=100
  # get EVI
  evi <- data %>% group_by(SEASON) %>% 
    summarize(evi=mean(EVI_300m, na.rm=TRUE)) %>% 
    column_to_rownames("SEASON")
  # get temps
  temps <- temp.lims %>% group_by(SEASON) %>% 
    summarize(min=min(min), max=max(max)) %>% column_to_rownames("SEASON")
  temps.dry <- seq(temps['DRY', 'min'], temps['DRY', 'max'], length.out=nrep)
  temps.wet <- seq(temps['WET', 'min'], temps['WET', 'max'], length.out=nrep)
  # put them into one df
  newdat <- data.frame(SEASON = rep(c("DRY", "WET"), each=nrep),
             TEMP_DEG_C = c(temps.dry, temps.wet),
             EVI_300m = rep(unlist(evi), each=nrep))
  
  p <- predict(model, newdat, 
                   se.fit = TRUE, interval = "confidence", level = 0.95)
  lines <- data.frame(SEASON = newdat$SEASON,
                      TEMP_DEG_C = newdat$TEMP_DEG_C,
                      p$fit)
  return(lines)
}

plotFit <- function(response) {
  if (response == "art") { mod.o <- mod.art.o; mod.r <- mod.art.r
  } else { mod.o <- mod.nat.o; mod.r <- mod.nat.r }
  set.seed(1000)
  lines.o <- createLinesTemp(mod.o, oryx)
  lines.r <- createLinesTemp(mod.r, roan)
  lines <- rbind(lines.o %>% mutate(SPECIES = "gemsbok"),
                 lines.r %>% mutate(SPECIES = "roan"))
  xlims = c(0, 50); xbreaks = seq(xlims[1], xlims[2], by=10);
  ggplot(data=lines, 
         mapping=aes(x=TEMP_DEG_C, color=SEASON)) +
    geom_line(mapping=aes(y=fit), linewidth=1) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.1) +
    xlab('') + ylab('') +
    guides(color="none") +
    scale_color_manual(values=c("#E28B07", "#06A081")) +
    scale_x_continuous(breaks=xbreaks, limits=xlims) +
    plot.theme + guides(color="none") + 
    facet_wrap(~SPECIES, nrow=1)
}

# Plotting seasonal water
p1 <- plotFit('nat')
p2 <- plotFit('art')
p1 / p2
ggsave(file=here::here("03_output", "waterdists", "allDistsTempLM.png"),
       width=9, height=9)

