# analyzemicroclimates.R
# Created 02 July 2023
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am("02_scripts/analysis/microclimates/analyzeMicroclimates.R")
source(here::here("02_scripts", "utilities.R"))
pacman::p_load(ggspatial, car, rstatix)
load(antfile)
setDataPaths('temperature')
load(here(procpath, 'era5estimates.rds'))

mydata <- collar.df.1h %>% nog() %>% 
  filter(!is.na(TEMP_DEG_C),
         COLLAR_ID != "SAT1728") #this one roan is out of bounds


# ******************************************************************************
#             TEMPERATURE DATA PLOTTING BY SPECIES AND BY MONTH
# ******************************************************************************

# by hour and species
plotByMonth <- function(sp) {
  p <- mydata %>% 
    filter(SPECIES == sp) %>% 
    ggplot(aes(x=hour(DATETIME), y=TEMP_DEG_C)) +
    geom_point(aes(y=TEMP_DEG_C, color=TEMPGROUP)) +
    geom_smooth(color='black') +
    facet_wrap(~MONTH, nrow=1) +
    plot.theme + 
    theme(strip.text=element_blank()) +
    scale_color_brewer(palette="RdBu", direction=-1, name="Temp (C)") +
    ylab('collar temperature (C)') + 
    ylim(0, 50) +
    ggtitle(sp) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
  p
}
p1 <- plotByMonth("Roan")
p2 <- plotByMonth("Gemsbok")


temp.diffs <- mydata %>% 
  group_by(SPECIES, MONTH, HOUR) %>% 
  summarize(TEMP=mean(TEMP_DEG_C)) %>% 
  dcast(MONTH + HOUR ~ SPECIES) %>% 
  filter(!is.na(Gemsbok), !is.na(Roan)) %>% 
  mutate(DIFFERENCE = Roan - Gemsbok,
         TEMPGROUP = factor(plyr::round_any(DIFFERENCE, 0.75)))
p3 <- ggplot( temp.diffs, aes(x=HOUR, y=DIFFERENCE, fill=TEMPGROUP)) + 
  geom_bar(stat="identity") +
  facet_wrap(~MONTH, nrow=1) +
  plot.theme + 
  ggtitle("Roan minus Gemsbok average temperature difference") +
  theme(strip.text=element_blank()) +
  scale_fill_brewer(palette="RdBu", direction=-1, name="Temp (C)") +
  scale_x_discrete(breaks=seq(0, 23, by=5)) +
  xlab('hour') + ylab('difference in collar\n temperature (C)')

p1 / p2 / p3 + plot_annotation(tag_level='a')
outfile=here::here('03_output', 'temperature', 'temperatureBySpecies.png')
ggsave(file=outfile)




################################################
### PLOTTING ERA5
################################################
set.seed(10001)
sub <- mydata %>% sample_n(20000)
sub$HMS <- format(as.POSIXct(sub$DATETIME), format = '%H:%M:%S')
sub$HMS <- round(period_to_seconds(lubridate::hms(sub$HMS))/60)

# Plotting collar temp vs era5
ggplot(data=sub,
       mapping=aes(x=ERA5_TEMP_C, y=TEMP_DEG_C, color=SPECIES)) +
                   # color=gsub("SAT", "", COLLAR_ID))) + 
  geom_point(alpha=0.2) + 
  facet_wrap(~SPECIES) + 
  guides(color=guide_legend("Collar ID")) +
  xlab("ERA5 reanalysis temperature (C)") + 
  ylab("Collar temperature (C)") +
  geom_abline(slope=1, intercept=0, color='black', 
              linetype="dashed", linewidth=1) +
  plot.theme + 
  scale_color_manual(values=c('#f72a4c', 'darkgray')) +
  guides(color="none") +
  theme(strip.text=element_text(size=20))
ggsave(here("03_output", "temperature", 'microvsERA5.png'),
       width=9, height=9)


# plotting era5 and collar temp vs hour of day
levels=c("ERA5 (reanalysis)\n temperature (C)", 
         "Collar temperature (C)")
hourly <- sub %>% 
  melt(measure.vars=c('ERA5_TEMP_C', 'TEMP_DEG_C'),
       value.name="temp") %>% 
  mutate(variable=ifelse(grepl("ERA5", variable), 
                         levels[1], levels[2]),
         variable = factor(variable, levels=levels))

ggplot(data=hourly, 
       mapping=aes(x=HMS, y=temp, color=SPECIES)) + 
  geom_point(alpha=0.01) +
  geom_smooth(method=loess) +
  facet_wrap(~variable, nrow=1) +
  scale_color_manual(values=c("#f72a4c", 'darkgray')) + 
  ylab('temperature (C)') +
  xlab('minutes since day start') +
  plot.theme + 
  theme(strip.text=element_text(size=20))
ggsave(here("03_output", "temperature", 'dailyTempsComparison.png'), 
         width=8, height=6)
  



################################################
### DAILY MAX COMPARISON
################################################
max.daily <- mydata %>% 
  group_by(SPECIES, SEASON, DATE, MONTH, DAYTYPE) %>% 
  summarize(max_era5 = max(ERA5_TEMP_C, na.rm=TRUE),
            max_temp = max(TEMP_DEG_C, na.rm=TRUE),
            n=n()) %>%
  filter(!is.na(max_temp), !is.infinite(max_temp), n > 15) %>% 
  ungroup() %>% 
  group_by(SPECIES, SEASON, DAYTYPE, MONTH, DATE) %>% 
  summarize(m_max_temp = mean(max_temp),
            m_max_era5 = mean(max_era5))

# get expected temperatures given species and season
lm_temp <- lm(m_max_temp ~ SEASON + m_max_era5, data=max.daily)
p <- predict(lm_temp, newdata = max.daily)
max.daily$m_max_exp = p

### BOXPLOT ###
ggplot() + 
  geom_boxplot(data=max.daily, 
               mapping=aes(x=DAYTYPE, y=m_max_exp),
               fill='#bdd6ff', color='#a4c6fc',
               outlier.shape=NA, coef=0) +
  geom_boxplot(data=max.daily, 
               mapping=aes(x=DAYTYPE, y=m_max_temp, fill=SPECIES),
               position=position_dodge(0.7), width=0.6) +
  facet_wrap(~SEASON) +
  ylab('maximum daily \nexperienced temp (C)') +
  xlab('daily heat category') +
  scale_fill_manual(values=c('#f72a4c', 'darkgray')) +
  theme(strip.text=element_text(size=20),
        text=element_text(size=20))
ggsave(here("03_output", "temperature", 'maxTempDAYTYPE.png'),
       width=9, height=9)


# collar vs era5 differences
diffs.df <- mydata %>%
  filter(ERA5_TEMP_C > 35) %>%
  dplyr::select(SPECIES, DATE, DAYTYPE, SEASON, TOD, 
                TEMP_DEG_C, ERA5_TEMP_C) %>% 
  mutate(diff=TEMP_DEG_C - ERA5_TEMP_C)

ggplot(diffs.df, 
       aes(x=DAYTYPE, y=diff, fill=SPECIES)) + 
  geom_boxplot() + 
  facet_wrap(~SPECIES+SEASON, nrow=1) + 
  theme(strip.text = element_blank()) +
  scale_fill_manual(values=c("#f72a4c", 'darkgray')) +
  # scale_fill_manual(values=c("#E28B07", "#06A081")) +
  theme(text=element_text(size = 22)) + 
  xlab('daily heat category') + 
  ylab('collar minus ambient temperature')
ggsave(here("03_output", "temperature", 'tempDiffs.png'), 
       width=9, height=6)
