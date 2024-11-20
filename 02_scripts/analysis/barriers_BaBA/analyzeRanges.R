# analyzeRanges.R
# Created 19 April 2023
# Margaret Swift <mes114@duke.edu>

# For now, running a simple convex hull to see the space used by each 
# individual in each season for the study period. Answering some questions on 
# species ranges

#   - MIGRATION HYPOTHESIS
#   - EXPANSION HYPOTHESIS
#   - ARE ORYX RANGES SMALLER NEAR FENCES?
#   - IS EVI/WATER IN WET SEASON RANGE BETTER?

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(lme4, nlme, glmm, adehabitatHR)
load(antfile)
load(here("03_output", "barriers", "encounters.RData"))
sf_use_s2(FALSE)

### SAVE THAT SHITIIITITITT! ###
# save(hulls, hulls.k, buffs, lines,
     # file=here("03_output", "barriers", "range_hulls.RData"))



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

collar.df %>% nog() %>% 
  mutate(YEAR=year(DATE)) %>% 
  group_by(COLLAR_ID, SPECIES, SEX, YEAR, SEASON) %>% 
  summarize(n=n()) %>% 
  print(n=50)



# ******************************************************************************
#                             Simple Convex Hulls
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
#                           k-Nearest Neighbor : k check
# ******************************************************************************

# figuring out best value of k
hullMe <- function(pts, k) {
  hull <- LoCoH.k(pts %>% as("Spatial"),
                  k=k, 
                  unin = "m",
                  unout = "ha",
                  duplicates="random", amount = NULL) %>% 
    st_as_sfc()
  
  ggplot() +
    geom_sf(data=hull, fill="black", alpha=0.1) +
    geom_sf(data=pts, color="red") + 
    ggtitle(paste(k, "nearest neighbors")) + 
    plot.theme + blank.theme
}

# good examples: 
# - 1020 dry season 2013
# - 1510 wet season 2016
season="WET"; year=2016; id="SAT1510"
pts <- collar.df %>% 
  mutate(YEAR=year(DATE)) %>% 
  filter(COLLAR_ID == id, YEAR==year, SEASON==season) 
title <- paste(id, season, "season", year)
p00 <- pts %>% st_union() %>% st_convex_hull() %>% 
  ggplot() + geom_sf(fill='pink', alpha=0.5) +
  geom_sf(data=pts, color="red") +
  ggtitle(title) + plot.theme 
p1 <- hullMe(pts, 5)
p2 <- hullMe(pts, 15)
p3 <- hullMe(pts, 25)
(p00+p1) / (p2 + p3)
f <- paste0("k_testing_", id, season, year, ".png")
ggsave(filename=here('03_output', 'barriers', 'sensitivity_testing', f),
       width=14.2, height=11.7)


# ******************************************************************************
#                           k-Nearest Neighbor Convex Hulls
# ******************************************************************************

createHull <- function(pts, k) {
  hull <- LoCoH.k(pts %>% as("Spatial"),
                  k=k, unin = "m", unout = "ha", duplicates="random", amount = NULL) %>% 
    st_as_sfc() %>% 
    st_union()
  return(hull)
}

hulls.k <- collar.df %>% nog() %>% 
  mutate(YEAR=year(DATE)) %>% 
  group_by(COLLAR_ID, SPECIES, SEX, YEAR, SEASON) %>% 
  summarize() %>% 
  mutate(AREA=NA, LON=NA, LAT=NA)

### WARNING THIS TAKES A LONG TIME! ###
for (i in 1:nrow(hulls.k)) {
  row <- hulls.k[i,]
  message(row$COLLAR_ID, row$YEAR, row$SEASON)
  pts <- collar.df %>% mutate(YEAR=year(DATE)) %>% 
    filter(COLLAR_ID==row$COLLAR_ID, 
           YEAR==row$YEAR, SEASON==row$SEASON) 
  
  #ERROR HANDLING
  hull <- tryCatch(
    createHull(pts, k=25),
    error=function(e) e
  )
  
  if(inherits(hull, "error")) {
    print("ERROR!")
    next
  }
  
  #REAL WORK
  hulls.k$AREA[i] <- st_area(hull) %>% units::set_units("ha")
  centers <- st_centroid(hull) %>% st_coordinates()
  hulls.k$LAT[i] <- centers[,2]
  hulls.k$LON[i] <- centers[,1]
}

LON.means <- hulls.k %>% 
  group_by(COLLAR_ID) %>% 
  st_drop_geometry() %>% 
  summarize(ML=mean(LON)) %>% 
  column_to_rownames("COLLAR_ID")
hulls.k <- hulls.k %>% 
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
#                           Caterpillars + MC
# ******************************************************************************
library(units)
lines <- collar.df %>% 
  mutate(ENC_FENCE=COLLAR_ID %in% encounters$AnimalID,
         YEAR = year(DATETIME)) %>% 
  group_by(SPECIES, SEX, COLLAR_ID, ENC_FENCE, YEAR, SEASON) %>% 
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")

# this took about 20 minutes
buffs <- lines %>% st_buffer(dist=.01)
areas <- buffs %>% st_area() %>% units::set_units('ha')
lines$AREA <- areas

ggplot() + 
  geom_sf(data=buffs[2,]) + 
  geom_sf(data=lines[2,])


# Seasonal shift for CH
buffs.sub <- lines %>% 
  filter(SEX=="F", !ENC_FENCE) %>% 
  st_drop_geometry() %>% 
  mutate(AREA = as.n(AREA))

# models
mod.r.area <- lme(AREA ~ SEASON,
                  random = ~1|COLLAR_ID,
                  data=buffs.sub %>% filter(SPECIES == "Roan"))
anova(mod.r.area)

mod.o.area <- lme(AREA ~ SEASON, 
                  random = ~1|COLLAR_ID,
                  data=buffs.sub %>% filter(SPECIES == "Oryx"))
anova(mod.o.area)



save()










# ******************************************************************************
#                               EXPANSION HYPOTHESIS
# ******************************************************************************

# Seasonal shift for CH
hulls.sub <- hulls %>% 
  filter(SEX=="F", !ENC_FENCE) %>% 
  st_drop_geometry()
  
# area change tab
hulls.sub %>% group_by(SPECIES, SEASON) %>% summarize(mean(AREA)/10000)

# models
mod.r.area <- lme(AREA ~ SEASON,
             random = ~1|COLLAR_ID,
             data=hulls.sub %>% filter(SPECIES == "Roan"))
anova(mod.r.area)
mod.o.area <- lme(AREA ~ SEASON, 
             random = ~1|COLLAR_ID,
             data=hulls.sub %>% filter(SPECIES == "Oryx"))
anova(mod.o.area)

p1 <- ggplot(hulls.sub) +
  geom_boxplot(aes(x=SEASON, y=AREA/10000000, fill=SEASON)) +
  facet_wrap(~SPECIES) +
  plot.theme +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  guides(fill="none") +
  ylab("area (hectare x 1000)")+
  scale_fill_manual(values=c('darkgray', 'black'), name="season")

# Seasonal shift for k-NNCH
hulls.k.sub <- hulls.k %>% 
  mutate(ENC_FENCE=COLLAR_ID %in% encounters$AnimalID) %>% 
  filter(SEX=="F", !ENC_FENCE, !is.na(AREA) )

# area change tab
hulls.k.sub %>% group_by(SPECIES, SEASON) %>% summarize(mean(AREA))

# models
mod.r.area <- lme(AREA ~ SEASON,
                  random = ~1|COLLAR_ID,
                  data=hulls.k.sub %>% filter(SPECIES == "Roan"))
anova(mod.r.area)
mod.o.area <- lme(AREA ~ SEASON, 
                  random = ~1|COLLAR_ID,
                  data=hulls.k.sub %>% filter(SPECIES == "Oryx"))
anova(mod.o.area)

p2 <- ggplot(hulls.k.sub) +
  geom_boxplot(aes(x=SEASON, y=AREA/1000, fill=SEASON)) +
  facet_wrap(~SPECIES) +
  plot.theme +
  ylab("area (hectare x 1000)")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  guides(fill="none") +
  scale_fill_manual(values=c('darkgray', 'black'), name="season")


# ******************************************************************************
#                               MIGRATION HYPOTHESIS
# ******************************************************************************


# lon change tab
hulls.sub %>% group_by(SPECIES, SEASON) %>% summarize(mean(LON))

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
p3 <- plotHulls(hulls.sub, "long") + guides(fill="none") +
  geom_hline(yintercept=0, color='red', alpha=0.5) + 
  ylab('change in longitude from the mean')

p1 + p2 + p3
outfile <- here::here("03_output", "barriers", "migration_expansion.png")
ggsave(file=outfile, width=16, height=7)



# ******************************************************************************
#                                 PLOTTING HULLS
# ******************************************************************************
# Plotting simple CH
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


