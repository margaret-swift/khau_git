# spaceUseSeasonal.R
# Created 19 April 2023
# Margaret Swift <mes114@duke.edu>

# For now, running a simple convex hull to see the space used by each 
# individual in each season for the study period.


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(lme4, nlme, glmm)
load(antfile)
load(here("03_output", "barriers", "encounters.RData"))


# ******************************************************************************
#                             INDIVIDUAL-YEAR STATS
# ******************************************************************************

tots <- collar.df %>% 
  st_drop_geometry() %>% 
  filter(SEX=="F") %>% 
  group_by(ANIMAL_ID, SPECIES, DATE) %>% 
  summarize() %>% 
  group_by(SPECIES) %>% 
  summarize(nday=n()) %>%  
  mutate(perc=nday/sum(nday)) %>% 
  column_to_rownames("SPECIES")

# ******************************************************************************
#                             SET UP SPACE USE RANGES
# ******************************************************************************

ids <- unique(collar.df$COLLAR_ID)
hulls <- collar.df %>% 
  mutate(YEAR = year(DATETIME)) %>% 
  group_by(SPECIES, SEX, COLLAR_ID, YEAR, SEASON) %>% 
  summarize() %>% 
  st_cast("MULTIPOINT") %>% 
  st_convex_hull() %>% 
  mutate(AREA=as.n(st_area(geometry)),
         ENC_FENCE=COLLAR_ID %in% encounters$AnimalID)

centers <- st_centroid(hulls) %>% st_coordinates()
hulls$LAT <- centers[,2]
hulls$LON <- centers[,1]
LON.means <- hulls %>% 
  group_by(COLLAR_ID) %>% 
  st_drop_geometry() %>% 
  summarize(ML=mean(LON)) %>% 
  column_to_rownames("COLLAR_ID")
hulls <- hulls %>% 
  mutate(LON.chg = LON - LON.means[COLLAR_ID, "ML"])




# ******************************************************************************
#                                 FUNCTIONS
# ******************************************************************************

plotHulls <- function(data, type, split=FALSE) {
  colors <- RColorBrewer::brewer.pal(name="BrBG", n=4)[c(1,2,4,3)]
  data = data %>% 
    filter(SEX=="F") %>% 
    mutate(ENC = ifelse(ENC_FENCE, "FENCE", "NO.FENCE"),
           ENC_SZN = paste(SEASON, ENC, sep='.') )
  if (type == "area") {
    data$YVAR <- data$AREA / 10000
    ylab <- "area (ha)"
  } else {
    data$YVAR <- data$LON.chg
    ylab <- "change in longitude"
  }
  if (split) data$XVAR <- data$ENC_SZN
  else data$XVAR <- data$SEASON
  
  p <- ggplot(data) +
    geom_boxplot(aes(x=XVAR, y=YVAR, fill=XVAR)) +
    facet_wrap(~SPECIES) +
    ylab(ylab) +
    plot.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank())
  if (split) {
    p  <- p + scale_fill_manual(values=colors, name="season x encounter")
  } else {
    p <- p + 
      scale_fill_manual(values=c('darkgray', 'black'), name="season")
  }
  p
}




# ******************************************************************************
#                                 EXPANSION HYPOTHESIS
# ******************************************************************************

# checking seasonal area shift
hulls.sub <- hulls %>% 
  filter(SEX=="F", !ENC_FENCE) %>% 
  st_drop_geometry()
  
mod.r.area <- lme(AREA ~ SEASON,
             random = ~1|COLLAR_ID,
             data=hulls.sub %>% filter(SPECIES == "Roan"))
anova(mod.r.area)
mod.o.area <- lme(AREA ~ SEASON, 
             random = ~1|COLLAR_ID,
             data=hulls.sub %>% filter(SPECIES == "Oryx"))
anova(mod.o.area)




# ******************************************************************************
#                               MIGRATION HYPOTHESIS
# ******************************************************************************

# checking seasonal centroid eastward shift
mod.r.long <- lme(LON.chg ~ SEASON, 
             random = ~1|COLLAR_ID,
             data=hulls.sub %>% filter(SPECIES == 'Roan'))
anova(mod.r.long)
mod.o.long <- lme(LON.chg ~ SEASON, 
             random = ~1|COLLAR_ID,
             data=hulls.sub %>% filter(SPECIES == 'Oryx'))
anova(mod.o.long)

# plot it
plotHulls(hulls.sub, "long") + guides(fill="none")
outfile <- here::here("03_output", "barrier_analysis", "migration.png")
ggsave(file=outfile, width=7, height=5)




# ******************************************************************************
#                                   PLOTTING
# ******************************************************************************
# Plotting 
p0 <- khauPlot() + geom_sf(data=hulls %>% filter(SEX=="F"), 
                           aes(color=SEASON, linetype = ENC_FENCE),
                           alpha=0.4, linewidth=1) +
  scale_color_manual(values=c('darkgray', 'black')) + 
  guides(color="none", linetype="none") +
  facet_wrap(~SPECIES+SEASON, ncol=2) + 
  blank.theme
p1 <- plotHulls(hulls.sub, type="longs") + guides(fill="none")
p2 <- plotHulls(hulls.sub, type="area") + guides(fill="none")
p0 + p1 + p2
outfile <- here::here("03_output", "barriers", "result_figs", 
                      "migration_expansion.png")
ggsave(file=outfile, width=20, height=10.5)




# ******************************************************************************
#                     ARE ORYX RANGES SMALLER NEAR FENCES?
# ******************************************************************************
mod.r <- lm(AREA ~ ENC_FENCE + SEASON, data=hulls %>% filter(SPECIES=="Roan"))
summary(mod.r)
mod.o <- lm(AREA ~ ENC_FENCE + SEASON, data=hulls %>% filter(SPECIES=="Oryx"))
summary(mod.o)


# checking seasonal area shift
area.m <- hulls %>% 
  filter(SEX=="F") %>% 
  st_drop_geometry() %>% 
  group_by(SPECIES, ENC_FENCE, SEASON) %>% 
  summarize(MAREA = mean(AREA) / 10000, n=n())

hulls %>% 
  mutate(AREA=AREA/10000) %>% 
  filter(SEX=="F") %>% 
  st_drop_geometry() %>% 
  ggplot() + 
  geom_boxplot(aes(x=interaction(ENC_FENCE, SEASON), y=AREA, fill=SEASON)) +
  facet_wrap(~SPECIES) + 
  plot.theme + 
  scale_fill_manual(values=c('darkgray', 'black'), name="season") + 
  # scale_alpha_manual(values=c(0.6, 1.0), name="fence encounter") +
  ylab("mean area (ha)") + 
  xlab("encountered fence")
outfile <- here::here("03_output", "barriers", "result_figs", 
                      "hulls_comparison.png")
ggsave(file=outfile)




# ******************************************************************************
#                     IS EVI/WATER IN WET SEASON RANGE BETTER?
# ******************************************************************************

hulls.r <- hulls %>% 
  filter(SPECIES == "Roan") %>% 
  arrange(COLLAR_ID, YEAR, SEASON) %>% 
  mutate(ID = paste0(COLLAR_ID, YEAR)) %>% 
  relocate(ID) %>% 
  st_as_sf()

dry <- hulls.r %>% filter(SEASON=="DRY")
wet <- hulls.r %>% filter(SEASON=="WET")

for (i in 1:nrow(wet)) {
  i=i+1
  pysanky <- dry$geometry[ which(dry$COLLAR_ID == wet$COLLAR_ID[i]) ]
  pysanka <- st_union(pysanky)
  donut <- st_difference(wet$geometry[i], pysanka)
  ggplot() + 
    geom_sf(data=pysanka, fill="orange", alpha=0.4) +
    geom_sf(data=wet[i,] %>% st_as_sf(), fill='blue', alpha=0.4) +
    geom_sf(data=donut %>% st_as_sf(), color='green', fill="transparent") +
    ggtitle(wet$ID[i])
  # st_geometry(wet[i,]) <- donut
}

ggplot() + 
  geom_sf(data=pysanka, fill="orange", alpha=0.4) +
  geom_sf(data=wet[i,] %>% st_as_sf(), fill='blue', alpha=0.4) +
  geom_sf(data=donut %>% st_as_sf(), color='green', fill="transparent")



# ******************************************************************************
#                     DOES HERD SIZE AFFECT CROSSINGS?
# ******************************************************************************

# create encounters dataset
en <- encounters %>% 
  mutate(COLLAR_ID = gsub('SAT', '', AnimalID),
         HERD_SIZE = collar.meta$HERD_SIZE[match(COLLAR_ID, collar.meta$COLLAR_ID)],
         HERD = factor(ifelse(HERD_SIZE == 1, "solitary", "herd"), levels=c('solitary', 'herd')),
         HERD_TYPE = ifelse(HERD_SIZE == 1, "solitary", 
                            ifelse(HERD_SIZE < 10, "small (<10)", 
                                   ifelse(HERD_SIZE < 25, "medium (<25)", "large"))),
         ISCROSS = ifelse(EVENT %in% c('Cross', "Jag"), "Cross", "Not Cross")) %>% 
  filter(EVENT != "Trapped", SPECIES=='Roan')


# Herd size
en %>% 
  filter(HERD_SIZE != "HERD") %>% # filter size we don't know
  dplyr::select(HERD_TYPE, ISCROSS) %>% 
  table() %>% 
  summary()
(p1 <- ggplot() + 
  geom_bar(data=en, 
           aes(x=HERD_SIZE, fill=ISCROSS), 
           position="fill", color='black') + 
  theme(text=element_text(size=18)) + 
  xlab('herd size') + ylab('frequency') +
  scale_fill_manual(values = c('#6e6e6e','#ffffff'), name="event type") + 
  facet_wrap(~SPECIES))

# Solitary status
en %>% 
  dplyr::select(HERD, ISCROSS) %>% 
  table() %>% 
  summary()
(p2 <- ggplot() + 
    geom_bar(data=en, aes(x=HERD, fill=ISCROSS), 
             position="fill", color='black') + 
    theme(text=element_text(size=18)) + 
    xlab('solitary status') + ylab('frequency') +
    scale_fill_manual(values = c('#6e6e6e','#ffffff'), name="event type") + 
    facet_wrap(~SPECIES))

# note though if you repeat the "HERD" filter you get a 
# p-value more like 0.022 since 1510 accounts for a lot of data.
p1 + p2 + plot_layout(guides = 'collect', widths = c(2, 1))
outfile <- here::here("03_output", "barriers", "result_figs", 
                      "herd_size_fig.png")
ggsave(file=outfile, width=9.5, height=8)


collar.meta %>% 
  filter(HERD_SIZE != "HERD") %>% 
  mutate(HERD_SIZE = as.n(HERD_SIZE)) %>% 
  group_by(ANIMAL) %>% 
  summarize(mean=mean(HERD_SIZE),
            med=median(HERD_SIZE),
            max=max(HERD_SIZE))


# ******************************************************************************
#         DO ROAN CROSSINGS CORRELATE WITH DAYS SINCE WET SEASON START?
# ******************************************************************************

# create encounters dataset
nov1 <- yday("2015-11-01")
may1 <- yday("2015-05-01")
en <- encounters %>% 
  mutate(ISCROSS = ifelse(EVENT %in% c('Cross', "Jag"), "Cross", "Not Cross"),
         DOY = yday(start_time),
         DSS = ifelse(DOY>=nov1, DOY-nov1, DOY-nov1+365)) %>% 
  filter(SPECIES=='Roan', EVENT != "Trapped")
col <- collar.df %>% 
  mutate(AnimalID=COLLAR_ID,
         DOY = yday(DATETIME),
         DSS = ifelse(DOY>=nov1, DOY-nov1, DOY-nov1+365)) %>% 
  filter(COLLAR_ID %in% unique(en$AnimalID))


colors=c('red', 'black')
base <- ggplot(data=en, mapping=aes(x=DSS)) +
  facet_wrap(~AnimalID, ncol=1, scales='free_y') + 
  xlab('days since wet season start') + 
  geom_vline(xintercept=(365-(nov1-may1)), linetype='dashed') + 
  scale_fill_manual(values=colors)

base + geom_histogram(data=col, mapping=aes(x=DSS), bins=50)
base + geom_density(mapping=aes(fill=ISCROSS), alpha=0.5)
base + geom_histogram(mapping=aes(fill=ISCROSS), bins=50, position='dodge')

ggplot() +
  geom_histogram(data=col, mapping=aes(x=DSS)) +
  facet_wrap(~COLLAR_ID, ncol=1, scales='free_y') + 
  xlab('days since wet season start') + 
  geom_vline(xintercept=(365-(nov1-may1)), linetype='dashed') + 
  scale_fill_manual(values=colors)
glm()



# ******************************************************************************
#                   DO CROSSINGS CORRELATE WITH ELEPHANT GAPS?
# ******************************************************************************
crosses <- encounters %>% 
  filter(EVENT == "Cross") %>% #, year(start_time) > 2015) %>% 
  st_as_sf(coords=c('easting', 'northing'), crs=st_crs(khau)) %>% 
  mutate(mindist=NA)
ele <- ele %>% 
  filter(year(DATE) <= 2017)

fakes <- st_sample(
    fence,
    1000,
    type = "random"
  ) %>% 
  st_cast("POINT") %>% 
  as.data.frame()
coords <- st_coordinates(fakes$geometry)
fakes$LON <- coords[,1]
fakes$LAT <- coords[,2]
fakes <- fakes %>% filter(LON < 21.2)
for (i in 1:nrow(fakes)) {
  fakes$mindist[i] <- min( st_distance(ele, fakes$geometry[i]) )
  fakes$mindist2[i] <- min( st_distance(crosses, fakes$geometry[i]) )
}
for (i in 1:nrow(crosses)) {
  cross <- crosses[i,]
  gaps <- ele#[ which(ele$DATE < cross$start_time), ]
  crosses$mindist[i] <- min( st_distance(gaps, cross) )
}

summary(fakes$mindist2)

# distance from other crossings
dists <- st_distance(crosses, crosses) %>% as.data.frame()
minNoZero <- function(x) {
  inx <- as.numeric(x) == 0
  x <- x[!inx]
  return(min(x))
}
sort(unlist( apply(dists, 1, minNoZero) ))

ggplot() +
  geom_sf(data=fence) +
  geom_sf(data=crosses, color='red', size=5) +
  geom_sf(data=gaps %>% filter(year(DATE) < 2017), size=0.5, color='blue') + 
  coord_sf(xlim=c(20.95, 21.05), ylim=c(-19, -18.6))
