# analyzeTransitions.R
# Created 14 June 2023
# Margaret Swift <mes114@duke.edu>

# HIDDEN MARKOV MODEL FITTING WITH momentuHMM
# PAPER: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12578
# VIGNETTE: https://cran.r-project.org/web/packages/momentuHMM/vignettes/momentuHMM.pdf


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(momentuHMM, rgdal)
load(antfile)
load(here('03_output', 'hmm', 'HMMOutputALL.Rdata'))

# ******************************************************************************
#                           TRANSITION REGRESSIONS
# ******************************************************************************

getProbabilitiesDAYTYPE <- function(t, seas, sp) {
  
  message("getting probabilities for ", sp, " ", seas, ' ', t)
  
  # pull data and set params
  row <- paste0(seas, t)
  tod <- paste0("TOD", t)
  ref.tod <- "DAWN"
  ref.day <- "MILD"
  
  # pull betas from model
  if (sp == "Oryx") { m <- m.oryx
  } else { m <- m.roan }
  coefs.df <- m$mle$beta %>% as.data.frame()
  
  # extract estimates for each transition type
  extractEstimate <- function(type, heat) {
    message('  extracting for transition ', type, ' : ', heat)
    coefs <- coefs.df %>% dplyr::select_at(type)
    wet.flag <- seas=="WET"
    day.heat <- paste0("DAYTYPE", heat)
    
    # set up estimate with baseline (intercept, temp, time since change)
    int <- coefs['(Intercept)',]
    wet.part <- coefs['SEASONWET',] * wet.flag
    
    # add DAYTYPE levels if it isn't MILD (reference level)
    if (heat != ref.day) {
      day.part <- 
        coefs[day.heat,] + 
        coefs[paste0(day.heat, ":SEASONWET"),] * wet.flag
    }
    
    # add TOD-associated values if it isn't DAWN (reference level)
    if (t != ref.tod) {
      tod.part <- 
        coefs[tod,] + 
        coefs[paste0(tod, ":SEASONWET"),] * wet.flag
    }
    
    # add both if both aren't ref level
    if ( t != ref.tod && heat != ref.day ) {
      day.tod.part <- 
        coefs[paste0(tod, ":", day.heat),] +
        coefs[paste(tod, day.heat, "SEASONWET", sep=":"),] * wet.flag
    }
    
    # calculate estimate baseline (without)
    est <- int + 
      wet.part + 
      day.part + 
      tod.part + 
      day.tod.part
      
    return(exp(est))
  }
  extractEstimates <- function(type) {
    ests <- c(extractEstimate(type, "MILD"),
              extractEstimate(type, "MED"),
              extractEstimate(type, "HOT"))
    ests
  }
  
  # get estimates for all transitions
  vals <- 1:3
  probs.df <- data.frame(matrix(nrow=3, ncol=9))
  names(probs.df) <- paste(rep(vals, each=3), vals, sep=" -> ")
  for (i in vals) {
    types <- paste(i, vals[vals != i], sep=" -> ")
    
    d2 <- extractEstimates(types[1])
    d3 <- extractEstimates(types[2])
    denom <- d2 + d3 + 1
    
    p2 <- d2 / denom 
    p3 <- d3 / denom
    p1 <- 1 - ( p2 + p3 )
    
    probs <- data.frame( p1, p2, p3 )
    names(probs) <- c(paste(i, i, sep=" -> "), types)
    probs.df[,names(probs)] <- probs
  }
  
  # put it together
  df <- cbind(data.frame(SEASON = seas, TOD = t, 
                         DAYTYPE=c("MILD", "MED", "HOT")), probs.df)
  df
}
getSeasonProbabilitiesDAYTYPE <- function(seas, sp) {
  df <- rbind(
    getProbabilitiesDAYTYPE("DAWN", seas, sp),
    getProbabilitiesDAYTYPE("DAY",  seas, sp),
    getProbabilitiesDAYTYPE("DUSK", seas, sp),
    getProbabilitiesDAYTYPE("NIGHT",seas, sp)
  )
  df <- df %>% 
    melt(variable.name="TRANSITION", 
         value.name="PROB", 
         id.vars=c('SEASON', 'TOD', 'DAYTYPE')) %>% 
    mutate(FROM = gsub(" ->.*", "", TRANSITION),
           TO = gsub(".*-> ", "", TRANSITION))
  nrep <- nrow(df)/3
  df$DAYTYPE <- rep(c("MILD", "MED", "HOT"), times=nrep)
  df
}

plotProbsDAYTYPE <- function(seas, sp) {
  colors <- c('#59aebd', '#76e2f5', "#f5bc62", "#85632d" )
  df <- getSeasonProbabilitiesDAYTYPE(seas, sp)
  ggplot() +
    geom_bar(data=df, 
             mapping=aes(x=DAYTYPE, y=PROB, fill=TOD, group=TOD),
             stat="identity", position="dodge") +
    facet_wrap(~TRANSITION) + 
    scale_fill_manual(values=colors) + 
    plot.theme + xlab("") + ylab("") + 
    theme(plot.title=element_blank()) +
    guides(color="none") 
}
plotProbs2 <- function(seas, sp) {
  df <- getSeasonProbabilitiesDAYTYPE(seas, sp) %>% 
    mutate(TOD = factor(TOD, levels=c("DAY", 'DUSK', 'NIGHT', 'DAWN')),
           DAYTYPE = factor(DAYTYPE, levels=c("MILD", "MED", "HOT")),
           TO = ifelse(TO == 1, "ruminating", ifelse(TO == 2, "foraging", "exploring")),
           FROM = ifelse(FROM == 1, "ruminating", ifelse(FROM == 2, "foraging", "exploring")),
           TRANTYPE = factor(ifelse(FROM == TO, "REMAIN", "SWITCH"), levels=c("SWITCH", "REMAIN")))
  
  # line version
  ggplot() +
    geom_line(data=df, 
              mapping=aes(x=DAYTYPE, y=PROB, 
                          color=FROM, group=FROM),
              linewidth=1) +
    facet_wrap(~TO + TOD, ncol=4) + 
    scale_color_brewer(palette="Accent") +
    plot.theme + xlab("") + ylab("") + 
    theme(plot.title=element_blank())
  ggplot() +
    
    # bar version
    geom_bar(data=df, 
              mapping=aes(x=DAYTYPE, y=PROB, 
                          fill=FROM, group=FROM),
              stat="identity", position="dodge") +
    facet_wrap(~TO + TOD, ncol=4) + 
    scale_fill_brewer(palette="Accent") +
    plot.theme + xlab("") + ylab("") + 
    theme(plot.title=element_blank())
}

p1 <- plotProbs("DRY", "Oryx")
p2 <- plotProbs("DRY", "Roan")
p3 <- plotProbs("WET", "Oryx")
p4 <- plotProbs("WET", "Roan")

p1 <- plotProbs2("DRY", "Oryx")
p2 <- plotProbs2("DRY", "Roan")
p3 <- plotProbs2("WET", "Oryx")
p4 <- plotProbs2("WET", "Roan")

