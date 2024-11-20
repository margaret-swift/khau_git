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
load(here('03_output', 'hmm', 'HMMOutputALL_TEMPC.Rdata'))

# ******************************************************************************
#                           TRANSITION REGRESSIONS
# ******************************************************************************

getProbabilities <- function(t, seas, sp) {

  message("getting probabilities for ", sp, " ", seas, ' ', t)

  # pull data and set params
  row <- paste0(seas, t)
  tod <- paste0("TOD", t)
  NT  <- length(temp.range)
  ref.lev <- "DAWN"

  # temperature should be standardized
  tmp <- (temp.range - mean(temp.range))/sd(temp.range)

  # pull betas from model
  if (sp == "Oryx") { m <- m.oryx
  } else { m <- m.roan }
  coefs.df <- m$mle$beta %>% as.data.frame()

  # extract estimates for each transition type
  extractEstimate <- function(type) {
    message('  extracting for transition ', type)
    coefs <- coefs.df %>% dplyr::select_at(type)
    wet.flag <- seas=="WET"

    # set up estimate with baseline (intercept, temp, time since change)
    int <- coefs['(Intercept)',]
    tmp.part <- coefs["TEMP_DEG_C",]
    wet.part <- coefs['SEASONWET',]
    wet.tmp.part <- coefs["TEMP_DEG_C:SEASONWET",]

    # calculate estimate baseline (without)
    est <- int +
      wet.part * wet.flag +
      tmp.part * tmp +
      wet.tmp.part * wet.flag * tmp

    # add TOD-associated values if it isn't DAWN (reference level)
    if (t != ref.lev) {
      est <- est +
        coefs[tod,] +
        coefs[paste0(tod, ":SEASONWET"),] * wet.flag +
        coefs[paste0(tod, ":TEMP_DEG_C"),] * tmp +
        coefs[paste0(tod, ":TEMP_DEG_C:SEASONWET"),] * wet.flag * tmp
    }

    return(exp(est))
  }

  # get estimates for all transitions
  vals <- 1:3
  probs.df <- data.frame(matrix(nrow=NT, ncol=9))
  names(probs.df) <- paste(rep(vals, each=3), vals, sep=" -> ")
  for (i in vals) {
    types <- paste(i, vals[vals != i], sep=" -> ")

    d2 <- extractEstimate(types[1])
    d3 <- extractEstimate(types[2])
    sum <- d2 + d3

    p2 <- d2 / ( 1 + sum )
    p3 <- d3 / ( 1 + sum )
    p1 <- 1 - ( p2 + p3 )

    probs <- data.frame( p1, p2, p3 )
    names(probs) <- c(paste(i, i, sep=" -> "), types)
    probs.df[,names(probs)] <- probs
  }

  # put it together
  df <- cbind(data.frame(SEASON = seas, TOD = t, TEMP=tmp), probs.df)
  df
}
getSeasonProbabilities <- function(seas, sp) {
  df <- rbind(
    getProbabilities("DAWN", seas, sp),
    getProbabilities("DAY",  seas, sp),
    getProbabilities("DUSK", seas, sp),
    getProbabilities("NIGHT",seas, sp)
  )
  df <- df %>%
    melt(variable.name="TRANSITION",
         value.name="PROB",
         id.vars=c('SEASON', 'TOD', 'TEMP')) %>%
    mutate(FROM = gsub(" ->.*", "", TRANSITION),
           TO = gsub(".*-> ", "", TRANSITION))
  nrep <- nrow(df)/length(temp.range)
  df$TEMP <- rep(temp.range, times=nrep)
  df <- df %>% 
    mutate(TOD = factor(TOD, levels=c("DUSK", "DAY", 'NIGHT', 'DAWN')),
           TO = ifelse(TO == 1, "ruminating", ifelse(TO == 2, "foraging", "exploring")),
           FROM = ifelse(FROM == 1, "ruminating", ifelse(FROM == 2, "foraging", "exploring")),
           TRANTYPE = factor(ifelse(FROM == TO, "REMAIN", "SWITCH"), levels=c("SWITCH", "REMAIN")))
  df
}
limitTemps <- function(df) {
  df <- df %>%
    filter(TEMP >= temp.lims[paste0(SEASON, TOD), 'min'],
           TEMP <= temp.lims[paste0(SEASON, TOD), 'max'])
  df
}
plotProbs <- function(seas, sp, plotType) {
  df <- getSeasonProbabilities(seas, sp)
  dfbold <- limitTemps(df)
  
  if (plotType == "TOD") {
    colors <- c('#59aebd', '#76e2f5', "#f5bc62", "#85632d" )
    p <- ggplot() +
      geom_line(data=df, mapping=aes(x=TEMP, y=PROB, group=TOD),
                alpha=0.25, color='gray', linewidth=0.75) +
      geom_line(data=dfbold, mapping=aes(x=TEMP, y=PROB, color=TOD, group=TOD), linewidth=1) +
      facet_wrap(~TRANSITION) +
      scale_color_manual(values=colors)
  } else if (plotType == "FROM") {
    p <- ggplot() +
      geom_line(data=df, 
                mapping=aes(x=TEMP, y=PROB, 
                            color=FROM, group=FROM),
                alpha=0.25, color='gray', linewidth=0.75) +
      geom_line(data=dfbold, mapping=aes(x=TEMP, y=PROB, 
                                         color=FROM, group=FROM), 
                linewidth=1) +
      facet_wrap(~TO + TOD, ncol=4) + 
      scale_color_brewer(palette="Accent")
  }
  p <- p + 
    xlim(c(min(temp.range)-5, max(temp.range)+5)) +
    plot.theme + xlab("") + ylab("") + 
    theme(plot.title=element_blank())
  return(p)
}

p1 <- plotProbs("DRY", "Oryx", "TOD")
p2 <- plotProbs("DRY", "Roan", "TOD")
p3 <- plotProbs("WET", "Oryx", "TOD")
p4 <- plotProbs("WET", "Roan", "TOD")
p <- (p1 + p2) / (p3 + p4)
ggsave(filename=here('03_output', 'hmm', "Transitions.png"),
       plot=p,
       width=15, height=15)

plotProbs2("WET", "Roan", "FROM")
plotProbs2("WET", "Oryx", "FROM")
plotProbs2("DRY", "Roan", "FROM")
plotProbs2("DRY", "Oryx", "FROM")

