# analyzeCrossings.R
# Created 01 June 2023
# Margaret Swift <mes114@duke.edu>

# Analyzing what affects crossing rates and locations
#   - DOES HERD SIZE AFFECT CROSSINGS?
#   - DO ROAN CROSSINGS CORRELATE WITH DAYS SINCE WET SEASON START?
#   - DO CROSSINGS CORRELATE WITH ANTELOPE/ELEPHANT GAPS?


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here("02_scripts", "utilities.R"))
load(antfile)
load(here("03_output", "barriers", "encounters.RData"))



# ******************************************************************************
#                         ENCOUNTER:FIX AND CROSS:ENCOUNTER
# ******************************************************************************

ids <- unique(encounters$AnimalID)
fixes <- collar.df %>% 
  nog() %>% 
  filter(COLLAR_ID %in% ids) %>% 
  group_by(COLLAR_ID, SPECIES, SEX) %>% 
  summarize(n.fix=n())
data <- encounters %>% 
  nog() %>% 
  rename(COLLAR_ID=AnimalID) %>% 
  mutate(iscross=EVENT=="Cross") %>% 
  group_by(COLLAR_ID, SPECIES, SEX) %>% 
  summarize(n.enc=n(),
            n.cross=sum(iscross)) %>% 
  left_join(fixes) %>% 
  mutate(r1 = n.enc / n.fix,
         r2 = n.cross / n.enc) %>% 
  dplyr::relocate(n.fix, .after=SEX)

tab=data %>% 
  group_by(SPECIES) %>% 
  summarize(enc=sum(n.enc),
            not.enc=sum(n.fix)-sum(n.enc)) %>% 
  column_to_rownames("SPECIES")
chisq.test(tab)

tab=data %>% 
  group_by(SPECIES) %>% 
  summarize(cross=sum(n.cross),
            not.cross=sum(n.enc)-sum(n.cross)) %>% 
  column_to_rownames("SPECIES")
chisq.test(tab)



# ******************************************************************************
#                         DOES SEASON AFFECT CROSSINGS?
# ******************************************************************************

ids <- unique(encounters$AnimalID)
fixes <- collar.df %>% 
  nog() %>% 
  filter(COLLAR_ID %in% ids) %>% 
  group_by(COLLAR_ID, SPECIES, SEX, SEASON) %>% 
  summarize(n.fix=n())
data <- encounters %>% 
  nog() %>% 
  rename(COLLAR_ID=AnimalID) %>% 
  mutate(iscross=EVENT=="Cross") %>% 
  group_by(COLLAR_ID, SPECIES, SEX, SEASON) %>% 
  summarize(n.enc=n(),
            n.cross=sum(iscross)) %>% 
  left_join(fixes) %>% 
  mutate(r1 = n.enc / n.fix,
         r2 = n.cross / n.enc) %>% 
  dplyr::relocate(n.fix, .after=SEASON)



tab=data %>% 
  filter(SPECIES == "Oryx") %>% 
  group_by(SEASON) %>% 
  summarize(enc=sum(n.enc),
            not.enc=sum(n.fix)-sum(n.enc)) %>% 
  column_to_rownames("SEASON")
chisq.test(tab)

tab=data %>% 
  filter(SPECIES == "Oryx") %>% 
  group_by(SEASON) %>% 
  summarize(cross=sum(n.cross),
            not.cross=sum(n.enc)-sum(n.cross)) %>% 
  column_to_rownames("SEASON")
chisq.test(tab)

# ******************************************************************************
#                         DOES HERD SIZE AFFECT CROSSINGS?
# ******************************************************************************

# create encounters crossing dataset
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

# collar.meta %>% 
#   filter(HERD_SIZE != "HERD") %>% 
#   mutate(HERD_SIZE = as.n(HERD_SIZE)) %>% 
#   group_by(ANIMAL) %>% 
#   summarize(mean=mean(HERD_SIZE),
#             med=median(HERD_SIZE),
#             max=max(HERD_SIZE))


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
#               DO CROSSINGS CORRELATE WITH ANTELOPE/ELEPHANT GAPS?
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
