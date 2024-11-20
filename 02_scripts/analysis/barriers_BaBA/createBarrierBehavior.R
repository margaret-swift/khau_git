# barrierBehavior.R
# Created 14 April 2023
# Margaret Swift <mes114@duke.edu>

# Xu, W., Dejid, N., Herrmann, V., Sawyer, H. & Middleton, A. D. 
#   Barrier Behaviour Analysis (BaBA) reveals extensive effects of fencing on 
#   wide‐ranging ungulates. J Appl Ecol 58, 690–698 (2021). 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2664.13806
#   https://github.com/wx-ecology/BaBA/

# Path straightness is the ratio between the displacement distance and the 
# accumulated step length of a trajectory, ranging from 0 (sinuous) to 1 (straight) 

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
source(here::here("02_scripts", "analysis", "barriers_BaBA", "_BaBA_fix.R"))
load(antfile)
pacman::p_load(sp, BaBA)

seasons <- collar.df %>% nog() %>% 
  group_by(DATE) %>% 
  summarize(SEASON = first(SEASON)) %>% 
  column_to_rownames("DATE")

# ******************************************************************************
#                             SET UP DATA TO FIT BABA
# ******************************************************************************

# Data points mutation
points_1h <- collar.df %>% 
  mutate(date=DATETIME, Animal.ID=COLLAR_ID,
         diff=ifelse(FIRST, NA, difftime(DATETIME, lag(DATETIME)))/60) %>% 
  filter(diff <= 70, diff >= 50)

points_5h <- collar.df.resamp %>% 
  filter(!(INX %in% points_1h$INX)) %>% 
  mutate(date=DATETIME, Animal.ID=COLLAR_ID)

# getting m/f data and #indivs
mf <- collar.df.resamp %>% st_drop_geometry() %>% 
  group_by(COLLAR_ID, SEX) %>% 
  summarize(SPECIES = first(SPECIES), SEX=first(SEX), NAME = first(ANIMAL_ID)) %>% 
  column_to_rownames('COLLAR_ID')
mf$ID <- paste0(mf$SEX, 1:nrow(mf))
nindiv <- collar.df %>% st_drop_geometry() %>% 
  group_by(SPECIES, SEX) %>% 
  summarize(n_distinct(ANIMAL_ID))

# ******************************************************************************
#                                 BABA PARAMETERS
# ******************************************************************************

b_time <- 10 #max #hours for behavior = short barrier interaction (bounce or quick)
p_time <- 100 #min #hours for behavior = prolonged barrier interaction (i.e. trapped)
w <- 10 #moving window size, in days.
max_cross <- 3 #max crosses allowed in trace and back-and-forth behavior
d <- 80 #distance considered "fence encounter"

# ******************************************************************************
#                             BARRIER ANALYSIS RUN
# ******************************************************************************
# source(here::here("02_scripts", "analysis", "BaBA_fix.R"))
baba_1h <- .BaBA(
  animal = points_1h, barrier = fence, 
  d = d, b_time = b_time, p_time = p_time, w = w, max_cross = max_cross,
  interval = 1, units = "hours", tolerance = 10
)
baba_5h <- .BaBA(
  animal = points_5h, barrier = fence, 
  d = d, b_time = b_time, p_time = p_time, w = w, max_cross = max_cross,
  interval = 5, units = "hours", tolerance = 25
)

# combine the two data frames
ambiguous = c("Trace", "Back_n_forth", "Trace_OR_Back_n_Forth")
encounters <- rbind(baba_1h$classification, baba_5h$classification) %>% 
  mutate(SPECIES = mf[AnimalID, "SPECIES"],
         SEX = mf[AnimalID, "SEX"],
         ID = mf[AnimalID, "ID"],
         NAME = mf[AnimalID, "NAME"],
         eventTYPE = ifelse(eventTYPE %in% ambiguous, "Trace", eventTYPE),
         eventTYPE = gsub("_", " ", eventTYPE),
         SEASON = seasons[as.c(date(start_time)),'SEASON'],
         start_side = "",
         end_side = "") %>% 
  relocate(start_side, end_side, .after=end_time) 

# decide what side the beginning and end points are on
# this takes ~2min
nb <- geodata::gadm("Namibia", level=0, path='tmp') %>% st_as_sf()
pb <- txtProgressBar(style = 3)
NE <- nrow(encounters)
for (i in 1:NE) {
  
  # set progress bar and pull matching data
  setTxtProgressBar(pb, i/NE)
  id <- encounters$AnimalID[i]
  data_i <- collar.df %>% filter(COLLAR_ID == id)
  
  # find start and finish points +-1 buffer
  t0 <- encounters$start_time[i]
  t1 <- encounters$end_time[i]
  inx.start <- which(data_i$DATETIME == t0) - 1
  inx.end <- which(data_i$DATETIME == t1) + 1
  
  # then intersection with Namibia to find sides
  ints <- data_i[c(inx.start,inx.end),] %>% 
    st_intersects(nb) %>% as.n()
  sides = ifelse(is.na(ints), "Botswana", "Namibia")
  encounters$start_side[i] <- sides[1]
  encounters$end_side[i] <- sides[2]
}

# assign "jag" and "cross" and rename eventTYPE to EVENT
encounters <- encounters %>% 
  mutate(EVENT = ifelse(eventTYPE == "Cross", 
                        ifelse(start_side == end_side, "Jag", "Cross"), 
                        eventTYPE))

outfile <- here::here("03_output", "barriers", "encounters.RData")
save(encounters.cbpp, file=outfile)
# load(outfile)



# ******************************************************************************
#                        SETTING UP DATA FOR PAPER ANALYSIS
# ******************************************************************************
encounter.data <- collar.df %>% nog() %>% 
  dplyr::select(COLLAR_ID, SPECIES, SEX, SEASON, DATETIME, LAT, LON)
ints <- encounters %>% mutate(int = interval(start_time, end_time))
encounter.data$EVENT = NA
NR = nrow(collar.df)
doubles <- c()
# THIS TAKES A WHILE
for (i in 1:NR) {
  # i = 17613
  if ((i %% 10000) == 0) message(i, ' / ', NR)
  row <- encounter.data[i,]
  c1 <- encounters$AnimalID == row$COLLAR_ID
  c2 <- encounters$SEASON == row$SEASON
  inx1 <- c1 + c2 == 2
  if (sum(inx1)) {
    inx2 <- row$DATETIME %within% ints$int
    inx3 <- (inx1 + inx2) == 2
    enc <- encounters[inx3,]
    if (nrow(enc) == 1) { encounter.data$EVENT[i] <- enc$EVENT
    } else if (nrow(enc)>1) {
      message("TWO MATCHES FOR ROW ", i)
      doubles <- c(doubles, i)
      encounter.data$EVENT[i] <- enc$EVENT[2]
    }
  }
}
encounter.data$IS_ENC = !is.na(encounter.data$EVENT)
encounter.data$IS_CROSS = ifelse(is.na(encounter.data$EVENT), NA,
                                 encounter.data$EVENT == "Cross")
save(encounter.data, file=here("03_output", "barriers", "encounter_data.rdata"))


# ******************************************************************************
#                             CBPP ANALYSIS RUN
# ******************************************************************************
iko <- 
samo <- cbpp %>% filter(grepl("Samochima", Name)) %>% dplyr::select(geometry)

findEncounters <- function(time, fencename) {
  fencename = toupper(fencename)
  barrier = cbpp %>% filter(grepl(fencename, toupper(Name))) %>% dplyr::select(geometry)
  if (time == "1h") {interval=1; tolerance=10; animal=points_1h
  } else {interval=5; tolerance=25; animal=points_5h;}
  baba <- .BaBA(
        animal = animal, barrier = barrier, 
        d = d, b_time = b_time, p_time = p_time, w = w, max_cross = max_cross,
        interval = interval, units = "hours", tolerance = tolerance
      ) 
  df <- baba$classification %>% as.data.frame() %>% mutate(FENCE=fencename)
  return(df)
}
baba_1h_sam <- findEncounters("1h", "Samochima")
baba_1h_iko <- findEncounters("1h", "Ikoga")
baba_5h_sam <- findEncounters("5h", "Samochima")
baba_5h_iko <- findEncounters("5h", "Ikoga")

# combine the two data frames
ambiguous = c("Trace", "Back_n_forth", "Trace_OR_Back_n_Forth")
encounters.cbpp <- rbind(baba_1h_iko, 
                         baba_5h_sam
                         ) %>% 
  mutate(SPECIES = mf[AnimalID, "SPECIES"],
         SEX = mf[AnimalID, "SEX"],
         ID = mf[AnimalID, "ID"],
         NAME = mf[AnimalID, "NAME"],
         eventTYPE = ifelse(eventTYPE %in% ambiguous, "Trace", eventTYPE),
         eventTYPE = gsub("_", " ", eventTYPE),
         SEASON = seasons[as.c(date(start_time)),'SEASON'])


for (i in 1:nrow(encounters.cbpp)) {
  start <- encounters.cbpp$start_time[i]
  end <- encounters.cbpp$end_time[i]
  df <- collar.df %>%
    filter(COLLAR_ID == encounters.cbpp$AnimalID[i],
           DATETIME %within% interval(start-hours(24), end+hours(24))
    ) %>%
    st_union() %>%
    st_cast("LINESTRING")
  p = ggplot() +
    geom_sf(data=cbpp, color='red') +
    geom_sf(data=df) + 
    ggtitle(i)
  print(p)
}

encounters.cbpp$eventTYPE[c(7,8)] <- "Jag"
encounters.cbpp$eventTYPE[c(13)] <- "Cross"
encounters.cbpp <- encounters.cbpp[-c(1,2,16)]

outfile <- here::here("03_output", "barriers", "encounters_cbpp.RData")
save(encounters.cbpp, file=outfile)
# load(outfile)

# ******************************************************************************
#                             OPTIMIZE BUFFER d
# ******************************************************************************

# By applying fence buffer distances every 10 meters from 50m - 150m, we
# define the optimal buffer as the distance at which the number of quick cross
# events begins to level off (< 1 % increase).
prepareQC <- function(points, interval, tolerance) {
  d.list <- seq(30, 110, by=10)
  qc.df <- data.frame(d=d.list, nqc=0, ntot=0, perc=0)
  
  for (i in 1:nrow(qc.df)) {
    d <- d.list[i]
    message("BaBA for ", d, "m")
    baba.out <- .BaBA(
      animal = points,
      barrier = fence,
      d = d,
      b_time = b_time,
      p_time = p_time,
      w = w,
      interval = interval,
      units = "hours",
      tolerance = tolerance,
      max_cross = max_cross
    )
    classes <- baba.out$classification 
    if (nrow(classes)) {
      qc.df$nqc[i] <- sum(classes$eventTYPE == "Cross")
      qc.df$ntot[i] <- nrow(classes)
    }
  }
  qc.df$perc <- qc.df$nqc / qc.df$ntot
  return(qc.df)
}

qc.df.1h <- prepareQC(points_1h, 1, 10)
qc.df.5h <- prepareQC(points_5h, 5, 50)

pt.1h <- data.frame(x=80, y=qc.df.1h$perc[qc.df.1h$d==80] )
pt.5h <- data.frame(x=80, y=qc.df.5h$perc[qc.df.5h$d==80] )
ggplot() +
  geom_point(data=qc.df.1h, aes(x=d, y=perc), color='black') +
  geom_line(data=qc.df.1h, aes(x=d, y=perc), color='black') +
  geom_point(data=qc.df.5h, aes(x=d, y=perc), color='darkgray') +
  geom_line(data=qc.df.5h, aes(x=d, y=perc), color='darkgray') +
  geom_point(data=pt.1h, mapping=aes(x=x, y=y), color="red", size=3) +
  geom_point(data=pt.5h, mapping=aes(x=x, y=y), color="red", size=3) +
  xlab("distance (m)") +
  ylab("proportion 'cross'") + plot.theme + 
  scale_x_continuous(breaks=d.list)
ggsave(here('03_output', 'barriers', 'sensitivity_testing', 'cross_dists.png'), 
       width=10, height=3.5)
