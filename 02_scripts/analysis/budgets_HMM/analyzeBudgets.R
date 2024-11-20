# analyzeBudgets.R
# Created 18 April 2023
# Margaret Swift <mes114@duke.edu>

# Using results of the simple HMM model (momentuHMM), what can we reveal about 
# each species' behavior patterns?

# Useful UCLA: https://stats.oarc.ucla.edu/r/dae/multinomial-logistic-regression/
# book chapter: https://bookdown.org/chua/ber642_advanced_regression/multinomial-logistic-regression.html
# book chapter: https://bookdown.org/egarpor/PM-UC3M/app-ext-multinomialreg.html

# cramer's V: 
# for 4 degrees of freedom, small=0.05	med=0.15	large=0.25
# webpage: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426219/#B2
# book: https://www.google.com/books/edition/Statistical_Power_Analysis_for_the_Behav/rEe0BQAAQBAJ


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(nnet, foreign, DescTools, lmtest, effects)
load(here('03_output', 'hmm', 'HMMOutputSimpleV2.Rdata'))
state.names <- c('ruminating', 'exploring', 'foraging')

# ******************************************************************************
#                                   FUNCTIONS
# ******************************************************************************

chiSqMe <- function(sp, by="SEASON", tod=NULL, season=NULL, hour=NULL) {
  if (sp == "Roan") { data <- roan
  } else { data <- oryx }
  if (!is.null(tod)) data <- data %>% filter(TOD==tod)
  if (!is.null(season)) data <- data %>% filter(SEASON == season)
  if (!is.null(hour)) data <- data %>% filter(hour(DATETIME) == hour)
  if (by == "SEASON") { 
    tab <- table( data$STATE, data$SEASON )
  } else { 
    tab <- table( data$STATE, data$DAYTYPE ) 
  }
  # ptab <- tab
  # ptab[,"DRY"] <- round( ptab[,"DRY"] / sum(tab[,'DRY']), 2 )
  # ptab[,"WET"] <- round( ptab[,"WET"] / sum(tab[,'WET']), 2 )
  mytest <- chisq.test(tab)
 
  print(tab)
  # print(ptab)
  print( mytest )
  rcompanion::cramerV(tab)
}





# ******************************************************************************
#               TESTING PROPORTION DIFFS FOR RUMINATION PEAKS
# ******************************************************************************


chiSqMe(sp="Oryx", by="SEASON", hour=3)
chiSqMe(sp="Oryx", by="SEASON", hour=12)
chiSqMe(sp="Oryx", by="SEASON", hour=20)


chiSqMe(sp="Roan", by="SEASON", hour=3)
chiSqMe(sp="Roan", by="SEASON", hour=12)
chiSqMe(sp="Roan", by="SEASON", hour=20)



# ******************************************************************************
#           TESTING PROPORTION DIFFS FOR SEASONS AND DAYTYPES
# ******************************************************************************

# fig 7A: dry-wet
chiSqMe("Oryx", by="SEASON")
chiSqMe("Roan", by="SEASON")

# fig 7B: mild-avg-hot
chiSqMe("Oryx", by="DAYTYPE")
chiSqMe("Roan", by="DAYTYPE")


# ******************************************************************************
#                 TESTING PROPORTION DIFFS BY DAY/NIGHT
# ******************************************************************************

chiSqMe("Oryx", by="DAYTYPE", tod="DAY", season="DRY")
chiSqMe("Oryx", by="DAYTYPE", tod="NIGHT", season="DRY")
chiSqMe("Oryx", by="DAYTYPE", tod="DAY", season="WET")
chiSqMe("Oryx", by="DAYTYPE", tod="NIGHT", season="WET")

chiSqMe("Roan", by="DAYTYPE", tod="DAY", season="DRY")
chiSqMe("Roan", by="DAYTYPE", tod="NIGHT", season="DRY")
chiSqMe("Roan", by="DAYTYPE", tod="DAY", season="WET")
chiSqMe("Roan", by="DAYTYPE", tod="NIGHT", season="WET")

my.df <- rbind(roan, oryx) %>%
  mutate(COMBO = paste(SPECIES, SEASON, DAYTYPE, TOD, sep="_")) %>% 
  group_by(SPECIES, SEASON, DAYTYPE, TOD, COMBO, STATE, .drop=FALSE) %>%
  summarize(n=n())

df.n <- my.df %>% 
  group_by(COMBO) %>% 
  summarize(n=sum(n))

my.df <- my.df %>% 
  mutate(total=df.n$n[match(COMBO, df.n$COMBO)],
         perc=n/total) %>% 
  ungroup() %>% 
  dplyr::select(-COMBO, -total) %>% 
  filter(TOD %in% c("DAY", "NIGHT"))

diffs.df <- my.df %>% 
  group_by(SPECIES, SEASON, TOD) %>% 
  filter(DAYTYPE != "AVG") %>% 
  arrange(SPECIES, SEASON, TOD, STATE) %>% 
  mutate(diff = perc - lag(perc)) %>% 
  filter(DAYTYPE == "HOT", STATE == "ruminating")










