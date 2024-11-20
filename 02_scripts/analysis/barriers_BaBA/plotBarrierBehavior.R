# plotBarrierBehavior.R
# Created 14 April 2023
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
here::i_am("02_scripts/analysis/barriers_BaBA/plotBarrierBehavior.R")
source(here::here("02_scripts", "utilities.R"))
load(antfile)
base <- here('03_output', 'barriers')
load(here(base, "encounters.RData"))

# ******************************************************************************
#                                   PARAMETERS
# ******************************************************************************

# Plotting parameters
my.theme <- 
  theme(text = element_text(size = 18),
        axis.title.x=element_blank())
colors = c(
 '#AEABAB', #avg mvmt
 '#F5C242', #bounce
 '#BF3288', #cross
 '#E8837A', #jag
 "#A0CE63", #trace or back and forth
 '#C0AEDC' #trapped
)
base.plot <- ggplot(encounters) +
  my.theme + 
  scale_fill_manual(values=colors) + 
  # guides(fill="none") +
  ylab('# encounters')


# massaging data for %encounter type
n.enc <- encounters %>% 
  group_by(SPECIES) %>% 
  summarize(n=n()) %>% 
  as.data.frame()
row.names(n.enc) <- paste(n.enc$SPECIES, n.enc$SEASON, sep="_")
tab <- encounters %>% 
  group_by(SPECIES, EVENT) %>% 
  summarize(n=n()) %>% 
  mutate(perc=n/n.enc[SPECIES, 'n'])
tab <- encounters %>% 
  group_by(SPECIES, SEASON) %>% 
  summarize(n=n())
###

# ******************************************************************************
#                             PLOTTING COVERAGE
# ******************************************************************************

lines2 <- collar.df %>% 
  group_by(SPECIES, SEX, ANIMAL_ID) %>% 
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")
xmin=20.4; xmax=21.15; 
ymin=-19.2; ymax=-18;
ggplot() + 
  geom_sf(data=khau, fill='#EBFBA6', linewidth=1.5) +
  geom_sf(data=lines2, 
          mapping=aes(color=ANIMAL_ID), alpha=0.5) + 
  facet_wrap(~SPECIES) + plot.theme + 
  geom_sf(data=fence, color='black', linewidth=2) +
  geom_sf(data=cbpp, color='darkgray', linewidth=1.5) +
  coord_sf(xlim=c(xmin, xmax), 
           ylim=c(ymin, ymax))

ggsave(here('03_output', 'eda', 'all_paths.png'))
#









# ******************************************************************************
#                        ENCOUNTER PROPORTIONS BY SPP
# ******************************************************************************

## Plot encounters by sex
base.plot + 
  geom_bar(aes(x=SEX, fill=EVENT), 
           position=position_dodge(preserve = "single")) +
  facet_wrap(~SPECIES)
ggsave(file.path(base, 'encounters_by_type.png'), width=7, height=10)

# Plot encounters by individual
base.plot + 
  geom_bar(aes(x=gsub('SAT', '', AnimalID), fill=EVENT), 
           position=position_dodge(preserve = "single")) +
  facet_wrap(~SPECIES, scales = "free_x") +
  guides(fill="none")
ggsave(file.path(base, 'encounters_by_indiv.png'), width=9, height=7.5)

# ******************************************************************************
#                        ENCOUNTER HISTOGRAM OVER TIME
# ******************************************************************************

# Plot encounters over time
rects <- data.frame(xstart = as.POSIXct(paste0(2013:2016, "-05-01")), 
                    xend = as.POSIXct(paste0(2013:2016, "-11-01")))
col <- collar.df %>% 
  mutate(COLLAR_TYPE = ifelse(COLLAR_TYPE == "SAT", "SAT", "SAT+"),
         SCALE = ifelse(SPECIES == "Oryx", 20, 10))
histbase <- ggplot(encounters) +
  my.theme + 
  ylab('# encounters') + 
  geom_rect(xmin=as.POSIXct(min(col$DATETIME)), xmax=Inf, 
            ymin=-Inf, ymax=Inf, fill="#bce0d3") +
  geom_rect(data = rects, mapping=aes(xmin = xstart, xmax = xend, 
                              ymin = -Inf, ymax = Inf),
                              fill = '#fadea5') +
  xlim(min(col$DATETIME), max(col$DATETIME))

p1 <- histbase + 
  geom_histogram(data=encounters %>% filter(SPECIES == "Oryx"), aes(x=start_time, fill=SPECIES), bins=80, fill="black") +
  geom_density(data=col %>% filter(SPECIES == "Oryx"), mapping=aes(x=DATETIME, y=..scaled..*20), color="red")
p2 <- histbase + 
  geom_histogram(data=encounters %>% filter(SPECIES == "Roan"), aes(x=start_time, fill=SPECIES), bins=80, fill="black") +
  geom_density(data=col %>% filter(SPECIES == "Roan"), mapping=aes(x=DATETIME, y=..scaled..*15), color="red")
p1 / p2
ggsave(file.path(base, 'encounters_over_time.png'), width=9, height=6, units="in")

# Plot encounters by season
base.plot + 
  geom_bar(aes(x=SEASON, fill=EVENT), 
           position=position_dodge(preserve = "single")) +
  facet_wrap(~SPECIES, scales = "free_x")
ggsave(file.path(base, 'encounters_by_season.png'), width=6, height=6, units="in")

# Plot % encounters by species
base.plot + 
  geom_bar(aes(x=SPECIES, fill=EVENT), position="fill") +
  facet_wrap(~SPECIES, nrow=1, scales="free_x") +
  ylab('% encounters') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(file.path(base, 'encounters_by_perc.png'), width=4, height=8, units="in")


# ******************************************************************************
#                                   MAPPING
# ******************************************************************************

# Plot the fence borderline and all 'encountering' individuals
ids <- unique(encounters$AnimalID)
data <- collar.df %>% 
  group_by(COLLAR_ID, SPECIES, SEX) %>% 
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING") %>% 
  filter(COLLAR_ID %in% ids) %>% 
  mutate(ID=paste0(gsub("SAT", "", COLLAR_ID), " (", SEX, ")"))
xmin=20.9; xmax=21.15; ymin=-19.0; ymax=-18.3
no.axes <- theme(axis.text.y = element_blank(),
                 axis.ticks.y=element_blank())
base.plot.2 <- ggplot() + 
  geom_sf(data=khau, fill='#abed9f', color='#71de5d', linewidth=2) + 
  geom_sf(data=rivers %>% filter(RIV_ORD <= 6), color='lightblue',
          linewidth=1, alpha=0.6) +
  geom_sf(data=fence, color='black', linewidth=2) +
  geom_sf(data=cbpp, color='black', linewidth=2) +
  theme(panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=1),
        text=element_text(size=18)) + guides(color="none") +
  scale_x_continuous(breaks=c(20.9, 21, 21.1))

p1=base.plot.2 +
  geom_sf(data=data%>%filter(SPECIES == "Oryx"), linewidth=1,
          mapping=aes(color=ID, group=ID), alpha=0.7) +
  coord_sf(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) + 
  scale_color_manual(values=c('#2891ed', '#1448d9', '#590ec9', '#a40ec9', '#de04cf'))
p2=base.plot.2 + 
  geom_sf(data=data%>%filter(SPECIES == "Roan"), linewidth=1,
          mapping=aes(color=ID, group=ID), alpha=0.7) + 
  no.axes +
  coord_sf(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) + 
  scale_color_manual(values=c('#fc3112', '#e65e10', '#e67205', '#ffbb1c', 
                              '#f4f716', '#e6052b'))
p1 + p2
ggsave(file.path(base, 'encounters_mapped.png'), width=8.5, height=12.5)







# ******************************************************************************
#                             MAPPING SEASONAL RANGES
# ******************************************************************************


load(here("03_output", "barriers", "range_hulls.RData"))

szn.lvls <- c("WET", "DRY")

id <- "SAT1728"
h <- hulls %>% filter(COLLAR_ID == id) %>% mutate(SEASON=factor(SEASON, levels=szn.lvls))
centers <- data.frame(SEASON=h$SEASON, YEAR=h$YEAR, LON=h$LON, LAT=h$LAT) %>% 
  st_as_sf(coords=c("LON", "LAT"), crs=4326)
pts <- collar.df.1h %>% filter(COLLAR_ID==id) %>% 
  mutate(YEAR=year(DATETIME), SEASON=factor(SEASON, levels=szn.lvls))
p1 <- ggplot() + 
  geom_sf(data=h, mapping=aes(color=SEASON), linewidth=1.5, fill=NA) + 
  geom_sf(data=pts, mapping=aes(color=SEASON), size=0.4) + 
  geom_sf(data=centers, color='red', size=6) +
  geom_sf(data=centers, mapping=aes(color=SEASON), size=4) +
  scale_color_manual(values=c('black', 'darkgray'), name="season") + 
  facet_wrap(~YEAR, nrow=1) +
  plot.theme + ggtitle(paste(gsub("SAT", "", id), "(Roan)"))


id <- "SAT1718"
h <- hulls %>% filter(COLLAR_ID == id) %>% 
  mutate(SEASON=factor(SEASON, levels=szn.lvls)) %>% 
  filter(YEAR != 2017)
centers <- data.frame(SEASON=h$SEASON, YEAR=h$YEAR, LON=h$LON, LAT=h$LAT) %>% 
  st_as_sf(coords=c("LON", "LAT"), crs=4326)
pts <- collar.df.1h %>% filter(COLLAR_ID==id) %>% 
  mutate(YEAR=year(DATETIME), SEASON=factor(SEASON, levels=szn.lvls)) %>% 
  filter(YEAR != 2017)
p2 <- ggplot() + 
  geom_sf(data=h, mapping=aes(color=SEASON), linewidth=1.5, fill=NA) + 
  geom_sf(data=pts, mapping=aes(color=SEASON), size=0.4) + 
  geom_sf(data=centers, color='red', size=6) +
  geom_sf(data=centers, mapping=aes(color=SEASON), size=4) +
  scale_color_manual(values=c('black', 'darkgray'), name="season") + 
  facet_wrap(~YEAR, nrow=1) +
  plot.theme + ggtitle(paste(gsub("SAT", "", id), "(Gemsbok)"))


p2 + p1
ggsave(file=here('03_output', 'barriers', 'rangeexamples.png'),
       width=22, height=6.8)



# EOF




# General stats
# encounters %>% group_by(SPECIES, SEASON) %>% summarize(n=n())
# m <- encounters %>% filter(SPECIES == "Oryx")
# chisq.test(table(m$SEASON, m$EVENT))
# 
# m <- encounters %>% 
#   filter(SPECIES == "Roan") %>% 
#   group_by(SPECIES, SEX, EVENT) %>%
#   summarize(n=n())
# chisq.test(table(m$SEX, m$EVENT))
# 
# 
# # looking at Trapped animals
# enc.trap <- encounters %>% filter(EVENT == "Trapped")
# int <- interval(enc.trap$start_time - hours(24),
#                enc.trap$end_time + hours(24))
# df <- data.frame(id=enc.trap$AnimalID,
#                  int=int)
# for (i in 1:nrow(df)) {
#   row <- df[i,]
#   d <- collar.df %>% 
#     filter(COLLAR_ID == row$id, 
#            DATETIME %within% row$int) %>% 
#     st_union() %>% st_cast("LINESTRING") 
#   bb <- st_bbox(d)
#   p <- ggplot(d) +
#     geom_sf() + 
#     ggtitle(row$id) +
#     geom_sf(data=fence, color='red') + 
#     coord_sf(xlim=c(bb['xmin'],bb['xmax']),
#              ylim=c(bb['ymin'],bb['ymax']))
#   print(p)
# }
