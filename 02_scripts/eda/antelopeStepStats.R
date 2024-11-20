# antelopeWaterUsage.R
# Created 25 March 2023
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
load(antfile)
pacman::p_load(rgeos)

blank.grid.theme <- 
  theme( panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank())
mf <- collar.df.resamp %>% st_drop_geometry() %>% 
  group_by(COLLAR_ID, SEX) %>% 
  summarize(SPECIES = first(SPECIES), SEX=first(SEX), NAME = first(ANIMAL_ID)) %>% 
  column_to_rownames('COLLAR_ID')
mf$ID <- paste0(mf$SEX, 1:nrow(mf))



# ******************************************************************************
#                             DATA TIMELINE BY INDIV
# ******************************************************************************
# melt dataset for plotting
lvls <- with(collar.meta, unique(COLLAR_ID[order(DATE_START)]))
collar.m <- collar.meta %>%
  filter(COLLAR_ID %in% gsub("SAT", "", unique(collar.df.1h$COLLAR_ID))) %>% 
  mutate(COLLAR_ID = factor(COLLAR_ID, levels=lvls)) %>% 
  melt(measure.vars=c('DATE_START', 'DATE_END'), value.name="DATE") %>% 
  mutate(DATE = as.POSIXct.Date(DATE),
         ANIMAL = ifelse(ANIMAL == "Oryx", "Gemsbok", ANIMAL))

# Plotting lengths of time for each animal
adj.theme <- theme( text = element_text(size=18), 
                    title = element_text(size=20),
                    axis.text.x = element_text(angle=45, hjust=1),
                    strip.text = element_text(size=20))
collar.m %>%
  ggplot(aes(x=DATE, y=COLLAR_ID), color="black") +
  geom_line(linewidth=2) +
  facet_wrap( ~ ANIMAL, nrow=2, scales = "free_y" ) +
  scale_x_datetime(date_breaks="1 month") +
  xlab("date") + ylab("Collar ID") + 
  adj.theme 
ggsave(file=here("03_output", "eda", "CollarTimeline.png"), width=12, height=7)

# plotting herd size for each animal
collar.m %>%
  mutate(HERD_SIZE = as.numeric(HERD_SIZE)) %>%
  filter(variable=="DATE_START") %>%
  ggplot() +
  geom_col(mapping=aes(x=HERD_SIZE, y=COLLAR_ID), fill="black") +
  adj.theme +
  facet_wrap( ~ ANIMAL, nrow=2, scales = "free_y" ) +
  xlab("Herd Size") 
ggsave(file=here("03_output", "eda", "HerdSizes.png"), width=5, height=7)



# ******************************************************************************
#                             NUMBER OF FIXES BY INDIVIDUAL
# ******************************************************************************

collar.df %>% 
  nog() %>% 
  group_by(COLLAR_ID) %>% 
  summarize(n=n()) %>% 
  print(n=30)




# ******************************************************************************
#                                INDIVIDUAL MAPS
# ******************************************************************************

data <- collar.df %>% 
  group_by(SPECIES, SEX, ANIMAL_ID) %>% 
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")

my.af <- world %>% 
  filter(continent == "Africa", !is.na(iso_a2),
         name_long %in% c('Namibia', 'Botswana')) %>% 
  st_transform(crs=st_crs(khau))

khauPlot(color='transparent', fill='gray') + 
  geom_sf(data=fence, color="red", linewidth=1.5) +
  geom_sf(data=data, aes(color=ANIMAL_ID), alpha=0.4) +
  facet_wrap(~SPECIES + SEX, nrow=1) +
  guides(color="none") +
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())
ggsave(file=here::here("03_output", "eda", "pathsBySpeciesSex.png"))

# ******************************************************************************
#                       STEP LENGTH / VELOCITY SUMMARIES
# ******************************************************************************

sum.df <- collar.df %>%
  st_drop_geometry() %>%
  group_by(ANIMAL_ID, DATE, SEASON) %>%
  summarize(KPD = sum(DIST)) %>%
  mutate(KPD = as.numeric(KPD / 1000)) %>% 
  filter(!is.na(KPD)) %>% 
  ungroup() %>% 
  mutate(SPECIES=mf$SPECIES[match(ANIMAL_ID, mf$NAME)])
my.sum <- sum.df %>% 
  filter(KPD<20) %>% 
  mutate(SEASON = factor(SEASON, levels=c("WET", "DRY")))
qlist <- my.sum %>% 
  group_by(SPECIES, SEASON) %>% 
  filter(!is.na(KPD)) %>% 
  summarize(avg = mean(KPD),
            q25=quantile(KPD, 0.25),
            q50=quantile(KPD, 0.50),
            q75=quantile(KPD, 0.75),
            q95=quantile(KPD, 0.95))
ggplot() + 
  geom_histogram(data=my.sum, aes(x=KPD, fill=SEASON),
                 alpha=0.6, bins=50, position="identity") +
  # geom_density(data=my.sum, aes(x=KPD, color=SEASON), 
               # position="identity", linewidth=1) +
  geom_vline(data=qlist, aes(xintercept=q50, linetype=SEASON), linewidth=0.7) +
  geom_vline(data=qlist, aes(xintercept=q95, linetype=SEASON), linewidth=0.7) +
  facet_wrap(~SPECIES, nrow=2) +
  scale_x_continuous(breaks=seq(0,20,by=1)) +
  scale_fill_brewer(palette="Dark2", direction=1) +
  scale_color_brewer(palette="Dark2", direction=1) +
  xlab('kilometers per day') +
  plot.theme
ggsave(here::here('03_output', 'eda', 'km_per_day_hist.png'))


# ******************************************************************************
#                            CHECKING FIX FREQUENCY
# ******************************************************************************

# set up data
data <- collar.df %>%
  st_drop_geometry() %>%
  mutate(
    TYPE = collar.meta$TYPE[match(gsub("SAT", "", ID), collar.meta$COLLAR_ID)],
    TYPE = ifelse(TYPE == "SAT", "SAT", "SAT+")) %>% 
  filter(SEX == "F")

my_bins = 24*60*60 # number of seconds in a day
ggplot(data=data %>% filter(TYPE == "SAT")) + 
  geom_histogram( aes(x=DATETIME), fill='orange', binwidth=my_bins ) + 
  facet_wrap(~ANIMAL_ID, nrow=3) + plot.theme
ggplot(data=data %>% filter(TYPE == "SAT+")) + 
  geom_histogram( aes(x=DATETIME), fill='darkgreen', binwidth=my_bins ) + 
  facet_wrap(~ANIMAL_ID, nrow=3) + plot.theme

# document fix numbers per day
fix.freq <- collar.df %>%
  group_by(ANIMAL_ID, DATE) %>%
  st_drop_geometry() %>%
  summarize(n=n())

fix.freq %>% filter(n>6) 
fix.freq %>% filter(n<4) 
dates.freq <- dcast(fix.freq, DATE ~ n)



# ******************************************************************************
#                       WATER DISTANCE BY SEASON
# ******************************************************************************

levels <- c('ephemeral','river','artificial')

# plot distance from waterholes by date
data.m <- collar.df.1h %>% 
  filter( IN_BORDER ) %>% 
  st_drop_geometry() %>% 
  mutate(HOUR = as.n(gsub(":.*", "", TIME))) %>% 
  group_by(SPECIES, MONTH) %>% 
  summarize(DIST_ART = mean(DIST_ART),
            DIST_NAT = mean(DIST_NAT),
            DIST_RIV = mean(DIST_RIV)) %>% 
  melt(measure.vars=c("DIST_ART", "DIST_NAT", "DIST_RIV"), 
       value.name="DISTANCE",
       variable.name="WATER_TYPE") %>% 
  mutate(WATER_TYPE = ifelse(WATER_TYPE == "DIST_NAT", levels[1], 
                             ifelse(WATER_TYPE == "DIST_RIV", levels[2], 
                                    levels[3])),
         WATER_TYPE = factor(WATER_TYPE, levels=levels))
         # SEASON = tolower(SEASON))
m <- fake.df %>% 
  st_drop_geometry() %>% 
  dplyr::select(DIST_NAT, DIST_RIV, DIST_ART) %>% 
  colMeans() %>% 
  melt(value.name="DISTANCE")
apathy.m <- data.frame(
  SPECIES = rep(c("Roan", "Oryx"), each=nrow(m)),
  WATER_TYPE = rep(row.names(m), times=2),
  DISTANCE = rep(m$DISTANCE, times=2)
)
apathy.m <- apathy.m %>% 
  mutate(WATER_TYPE = ifelse(WATER_TYPE == "DIST_NAT", levels[1], 
                             ifelse(WATER_TYPE == "DIST_RIV", levels[2], 
                                    levels[3])),
         WATER_TYPE = factor(WATER_TYPE, levels=levels))
colors <- c('orange', 'blue')
base <- ggplot(data=data.m) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, 
                ymin=-Inf, ymax=Inf),
            fill=colors[2], alpha=0.25, color='transparent') +
  geom_rect(aes(xmin=5, xmax=10, 
                ymin=-Inf, ymax=Inf),
            fill=colors[1], alpha=0.25, color='transparent')
ybreaks <- seq(0, 8500, by=1000)
base + geom_hline(yintercept=ybreaks, color='lightgray') +
  geom_line(data=data.m, aes(x=as.n(MONTH), y=DISTANCE, 
                group=SPECIES, linetype=SPECIES), linewidth=1) +
  geom_hline(data=apathy.m,
             aes(yintercept=DISTANCE),
             color="black", linetype="dotted") +
  facet_wrap(~ WATER_TYPE) +
  ylim(0, 8500) +
  xlab('month') + 
  ylab('distance from water (m)') +
  scale_x_continuous(breaks=seq(1, 12, by=1)) +
  plot.theme

outfile <- here::here("03_output", "eda", "waterDistLinesMonthly.png")
ggsave(file=outfile)

# to make insets, add free_y and remove ylim()
outfile <- here::here("03_output", "eda", "waterDistLinesInsets.png")
ggsave(file=outfile)

# ******************************************************************************
#                   WATERHOLE DISTANCE BY SEASON -- BOXPLOT
# ******************************************************************************

data.m.box <- collar.df %>%
  st_drop_geometry() %>% 
  dplyr::select(SPECIES, SEASON, DIST_ART, DIST_NAT, DIST_RIV) %>% 
  melt(id.vars=c("SPECIES", "SEASON"), 
       value.name="DISTANCE",
       variable.name="WATER_TYPE") %>% 
  mutate(
    WATER_TYPE = ifelse(WATER_TYPE == "DIST_NAT", levels[1], 
                        ifelse(WATER_TYPE == "DIST_RIV", levels[2], 
                               levels[3])),
    WATER_TYPE = factor(WATER_TYPE, levels=levels))


### BOXPLOT ###
ggplot(data=data.m.box, 
       aes(x=SEASON,
           y=DISTANCE, 
           group=interaction(WATER_TYPE, SEASON),
           fill=SEASON)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~SPECIES+WATER_TYPE, nrow=1) +
  ylab('distance from water (m)') +
  scale_y_continuous(breaks=seq(0, 13000, by=1000),
                     limits=c(0, 13000)) +
  scale_fill_brewer(palette="Dark2", name="Season", direction=-1) +
  plot.theme

outfile <- here::here("03_output", "eda", "waterDistBoxplot.png")
ggsave(file=outfile)

### VIOLIN ###
ggplot(data=data.m.box, 
       aes(x=SEASON,
           y=DISTANCE, 
           group=interaction(WATER_TYPE, SEASON),
           fill=SEASON)) +
  geom_violin(draw_quantiles=0.5) +
  facet_wrap(~SPECIES+WATER_TYPE, nrow=1) +
  ylab('distance from water (m)') +
  scale_y_continuous(breaks=seq(0, 13000, by=1000),
                     limits=c(0, 13000)) +
  scale_fill_brewer(palette="Dark2", name="Season", direction=-1) +
  plot.theme

outfile <- here::here("03_output", "eda", "waterDistViolin.png")
ggsave(file=outfile)

### DENSITY ###
ggplot(data=data.m.box %>% filter(DISTANCE < 20000), 
       aes(x=DISTANCE, 
           fill=SEASON)) +
  geom_density(alpha=0.7) +
  facet_wrap(~SPECIES + WATER_TYPE, nrow=2, scales='free_x') +
  xlab('distance from water (m)') +
  ylab('density') +
  scale_fill_brewer(palette="Dark2", name="Season", direction=-1) +
  plot.theme +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  coord_flip()

outfile <- here::here("03_output", "eda", "waterDistDensityFlip.png")
ggsave(file=outfile)



# ******************************************************************************
#                          EVI BY SEASON -- BOXPLOT
# ******************************************************************************

evi.vals <- data.frame(month=1:12, evi=c(0.42, 0.40, 0.41, 0.37, 0.29, 0.25, 0.23, 0.22, 
                       0.22, 0.23, 0.30, 0.39))

# plot mean EVI by season
data.m <- collar.df %>% 
  st_drop_geometry() %>% 
  mutate(MONTH = as.numeric(MONTH)) %>% 
  group_by(SPECIES, MONTH) %>% 
  summarize(EVI = mean(EVI, na.rm=TRUE),
            DIST_WATER = mean(DIST_WATER, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(EVI_DIFF = EVI - evi.vals$evi[MONTH])
lasts <- data.m[data.m$MONTH == 1,] %>% mutate(MONTH = 13)
data.m <- rbind(data.m, lasts)

colors <- RColorBrewer::brewer.pal("BrBG", n=8)[c(4,5)]

base <- ggplot(data=data.m) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, 
                ymin=-Inf, ymax=Inf),
            fill=colors[2], color='transparent') +
  geom_rect(aes(xmin=5, xmax=10, 
                ymin=-Inf, ymax=Inf),
            fill=colors[1], color='transparent')

ybreaks <- seq(round(min(data.m$EVI_DIFF), 2),
               max(data.m$EVI_DIFF), by=0.01)
base + 
  geom_hline(yintercept=ybreaks, color='lightgray') +
  geom_hline(yintercept=0, color='darkgray') +
  geom_line(aes(x=MONTH, 
                y=EVI_DIFF, 
                group=SPECIES,
                linetype=SPECIES), 
            linewidth=1) +
  scale_x_continuous(breaks=1:13) +
  scale_y_continuous(breaks=ybreaks) +
  plot.theme + 
  blank.grid.theme +
  ylab('EVI difference from monthly mean')
outfile <- here::here("03_output", "eda", "EVIUsage.png")
ggsave(file=outfile)

ybreaks <- seq(250, 1050, by=100)
base + geom_hline(yintercept=ybreaks, color='lightgray') +
  geom_line(aes(x=MONTH, 
                y=DIST_WATER, 
                group=SPECIES,
                linetype=SPECIES), 
            linewidth=1) +
  scale_x_continuous(breaks=1:12) +
  scale_y_continuous(breaks=ybreaks) +
  scale_fill_manual(values=c('white', 'lightgray')) +
  plot.theme + 
  blank.grid.theme
outfile <- here::here("03_output", "eda", "EVIUsage.png")
ggsave(file=outfile)



# ******************************************************************************
#                          More water distance
# ******************************************************************************
data <- collar.df.1h
base <- ggplot(data=data %>% sample_n(10000),
               aes(x=TEMP_DEG_C, 
                   color=SEASON)) +
  facet_wrap(~SPECIES) +
  plot.theme + 
  xlab('temperature (C)') +
  scale_color_brewer(palette="Dark2", direction=-1)

p1 <- base +
  geom_point(aes(y=DIST_ART), alpha=0.1) +
  geom_smooth(aes(y=DIST_ART), method="lm", linewidth=2) +
  ylab('dist. from artificial water') + guides(color="none")
p2 <- base +
  geom_point(aes(y=DIST_NAT), alpha=0.1) +
  geom_smooth(aes(y=DIST_NAT), method="lm", linewidth=2) +
  ylab('dist. from ephemeral natural water') + guides(color="none")
p3 <- base +
  geom_point(aes(y=DIST_RIV), alpha=0.1) +
  geom_smooth(aes(y=DIST_RIV), method="lm", linewidth=2) +
  ylab('dist. from rivers') + guides(color="none")

ggsave(p1, file=here::here("03_output", "eda", "LMArtDist.png"))
ggsave(p2, file=here::here("03_output", "eda", "LMNatDist.png"))
ggsave(p3, file=here::here("03_output", "eda", "LMRivDist.png"))
