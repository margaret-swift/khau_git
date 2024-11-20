# CTMM.R
# Created 07 Jul 2023
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(ctmm, tidyverse, RSpectra)
load(antfile)

test <- collar.df.1h %>% 
  filter(COLLAR_ID == 'SAT1717') %>% 
  rename(individual.local.identifier=COLLAR_ID, 
         timestamp=DATETIME,
         location.long=LON,
         location.lat=LAT) %>% 
  nog() %>% 
  as.telemetry()
test
plot(test)
vg <- variogram(test)
variogram.fit(vg, interactive=TRUE)
plot(vg)
fitted.mods <- ctmm.select(test, CTMM=GUESS, verbose=TRUE)
