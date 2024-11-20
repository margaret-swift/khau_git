# analyzeWaterDists.R
# Created 29 March 2023
# Margaret Swift <mes114@duke.edu>


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am("./02_scripts/analysis/waterdists/analyzeWaterDists.R")
source(here::here("02_scripts", "utilities.R"))
load(antfile)

# ******************************************************************************
#                                   FUNCTIONS
# ******************************************************************************

getFullEffect <- function(model) {
  cfs <- coef(model)
  wet <- cfs[["TEMP_DEG_C"]] + cfs[["TEMP_DEG_C:SEASONWET"]]
  dry <- cfs[["TEMP_DEG_C"]]
  response <- all.vars(model$call)[1]
  sigs <- summary(model)$coefficients[,4]
  
  message('coefficients:')
  print(cfs)
  
  message('significances:')
  print(sigs)
}


# ******************************************************************************
#                     ROAN AND DIST FROM ART WATER BY TEMP
# ******************************************************************************

roan <- collar.df.1h %>% nog() %>% filter(SPECIES=="Roan")
oryx <- collar.df.1h %>% nog() %>% filter(SPECIES=="Oryx")


# *************************************************************
#                             ORYX 
# *************************************************************
mod.nat.o <- lm(DIST_NAT ~ TEMP_DEG_C*SEASON + EVI, data=oryx)
mod.art.o <- lm(DIST_ART ~ TEMP_DEG_C*SEASON + EVI, data=oryx)

getFullEffect(mod.nat.o)
getFullEffect(mod.art.o)

# *************************************************************
#                             ROAN 
# *************************************************************

mod.nat.r <- lm(DIST_NAT ~ TEMP_DEG_C*SEASON + EVI, data=roan)
mod.art.r <- lm(DIST_ART ~ TEMP_DEG_C*SEASON + EVI, data=roan)

getFullEffect(mod.nat.r)
getFullEffect(mod.art.r)


# ******************************************************************************
#                                   SAVE
# ******************************************************************************
outfile <- here::here("03_output", "waterdists", "allModelsUPDATED.rdata")
save(mod.nat.r, mod.art.r,
     mod.nat.o, mod.art.o,
     roan, oryx,
     file=outfile
)

ttestDist <- function(sp, var) {
  message(sp, ' ', var, ' seasonal comparison')
  vals <- collar.df.1h %>% nog() %>% filter(SPECIES == sp)
  dry <- vals[vals$SEASON=="DRY", var]
  wet <- vals[vals$SEASON=="WET", var]
  message("wet season mean: ", 
          round(mean(wet, na.rm=TRUE), 3),
          " +/- ",
          round(sd(wet, na.rm=TRUE)))
  message("dry season mean: ", 
          round(mean(dry, na.rm=TRUE), 3),
          " +/- ",
          round(sd(dry, na.rm=TRUE)))
  print(t.test(wet, dry))
  print(effsize::cohen.d(wet, dry))
}
ttestAllDists <- function(sp) {
  ttestDist(sp, "DIST_NAT")
  ttestDist(sp, "DIST_ART")
}
ttestAllDists("Oryx")
ttestAllDists("Roan")



