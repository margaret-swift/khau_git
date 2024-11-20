# plotBudgets.R
# Created 22 April 2023
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
pacman::p_load(nnet, foreign, DescTools, lmtest, momentuHMM)
load(here('03_output', 'hmm', 'HMMOutputSimpleV2.Rdata'))
# load(here('03_output', 'hmm', 'HMMbudgetModels.Rdata'))

makePretty <- function(p) {
  p + 
    scale_fill_brewer(palette="Accent") +
    scale_color_brewer(palette="Accent") +
    guides(fill="none", color="none") + 
    xlab("") + ylab("") +
    theme(text = element_text(size = 18),
          strip.background = element_blank(),
          strip.text = element_blank() ) 
}



# ******************************************************************************
#                      ( FIG X ) 1h COLLAR DATA PLOTTED
# ******************************************************************************


lines <- collar.df.1h %>% 
  mutate(YEAR = year(DATETIME), ID=gsub("SAT", "", COLLAR_ID)) %>% 
  group_by(SPECIES, SEX, ID, YEAR, SEASON) %>% 
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")
no.bg <- theme(rect = element_rect(fill = "transparent"))

khauPlot() + 
  geom_sf(data=lines, aes(color=ID)) + 
  geom_sf(data=water.artificial, color='black', size=4) +
  plot.theme + 
  guides(color="none") + 
  no.bg +
  blank.theme

ggsave(file=here('03_output', 'eda', 'collarLocations.png'),
       width=13, height=12)



# ******************************************************************************
#                        ( FIG 6 ) DAILY ACTIVITY BUDGETS PLOT
# ******************************************************************************

# time budget percents
percents <- rbind(roan, oryx) %>% 
  mutate(HOUR=hour(DATETIME),
         DATEX=as_date(ifelse(HOUR<=5, Sys.Date()+1, Sys.Date())),
         DATETIMEX = as.POSIXct(paste(DATEX, HOUR), format="%Y-%m-%d %H")) %>% 
  group_by(SPECIES, SEASON, HOUR, DATEX, DATETIMEX, STATE) %>% 
  summarize(COUNT = n()) %>% 
  mutate(PER = COUNT/sum(COUNT))

# FIG 6A
makePretty(
  percents %>%   
    group_by(SPECIES, SEASON, STATE) %>% 
    summarize(COUNT = sum(COUNT)) %>% 
    ggplot() +
    geom_bar(aes(x=SEASON, fill=STATE, y=COUNT), position="fill", stat="identity") +
    facet_wrap(~SPECIES)
)
outfile = here('03_output', 'hmm', "TimeBudgetBasic.png")
ggsave(file=outfile, width=3, height=5)

# FIG 6B
makePretty(
  percents %>% 
  ggplot() + 
  geom_bar(aes(x=DATETIMEX, fill=STATE, y=COUNT), position="fill", stat="identity") +
  facet_wrap(~SEASON+SPECIES) +
  scale_x_datetime(date_breaks="2 hours", date_labels="%H")
)
outfile = here('03_output', 'hmm', 'TimeBudget.png')
ggsave(file=outfile, width=12, height=7)


#FIG 6C
percents %>% 
  ungroup() %>% 
  filter(HOUR %in% c(3,12,21)) %>%
  mutate(PER=round(PER*100)) %>% 
  arrange(SPECIES, SEASON, STATE) %>% 
  dplyr::select(-DATEX, -DATETIMEX, -COUNT) %>% 
  print(n=36)
percents %>% 
  ungroup() %>% 
  filter(HOUR %in% 14:20, STATE != "exploring",
         SEASON=="DRY", SPECIES == "Roan") %>%
  arrange(STATE) %>% 
  mutate(PER=round(PER*100)) %>% 
  dplyr::select(-DATEX, -DATETIMEX, -COUNT) %>% 
  print(n=24)

makePretty(
  ggplot(percents, aes(x=DATETIMEX, y=PER, group=STATE, color=STATE)) + 
  facet_wrap(~SEASON + SPECIES) + 
  geom_line(linewidth=1) + 
  plot.theme + xlab('hour') + 
  scale_x_datetime(date_breaks="2 hour", date_labels="%H")
  )
outfile = here('03_output', 'hmm', 'LineBudget.png')
ggsave(file=outfile, width=12, height=7)


# ******************************************************************************
#                      ( FIG 7 ) DAILY BUDGET BY DAY TEMPERATURE
# ******************************************************************************

# Time budgets plots
plot.temp <- makePretty(
  rbind(roan, oryx) %>% 
  mutate(SPECIES = ifelse(SPECIES == "Oryx", "Gemsbok", "Roan"),
         TOD = factor(TOD, levels=c("DAY", "DUSK", "NIGHT", "DAWN"))) %>% 
  ggplot() + 
  geom_bar(aes(x=DAYTYPE, fill=STATE), position="fill")
)

# FIG 7A (SUPPLEMENT)
plot.temp + facet_wrap(~SPECIES)
ggsave(file=here('03_output', 'hmm', "TimeBudgetTEMPSBasic.png"),
       width=4, height=5)

# FIG 7B
makePretty(
  rbind(roan, oryx) %>% 
    mutate(SPECIES = ifelse(SPECIES == "Oryx", "Gemsbok", "Roan"),
           TOD = factor(TOD, levels=c("DAY", "DUSK", "NIGHT", "DAWN"))) %>% 
    filter(TOD %in% c("DAY", "NIGHT")) %>% 
    ggplot() + 
      geom_bar(aes(x=DAYTYPE, fill=STATE), position="fill") +
      facet_wrap(~SEASON + TOD + SPECIES, nrow=2)
)
ggsave(file=here('03_output', 'hmm', "TimeBudgetTEMPS.png"),
       width=8, height=9.5)

# FIG 7C
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

ggplot(diffs.df) +
  geom_bar(aes(x=interaction(TOD, SEASON, DAYTYPE), 
               y=diff, fill=SPECIES),
           stat="identity", position="dodge") + 
  geom_hline(yintercept=0, color='black', linewidth=1) +
  facet_wrap(~SEASON, nrow=1, scales="free_x") +
  plot.theme + 
  xlab("") + ylab("change in proportion ruminating") + 
  theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c('#F4A958', '#F4C99B'))
ggsave(file=here('03_output', 'hmm', "RuminDiffs.png"),
       width=6.2, height=8.2)



ggplot(diffs.df) +
  geom_bar(aes(x=interaction(TOD, DAYTYPE), y=diff, fill=SPECIES),
           stat="identity", position="dodge") + 
  geom_hline(yintercept=0, color='black', linewidth=1) +
  facet_wrap(~SEASON, nrow=2) +
  plot.theme + guides(fill="none") +
  xlab("") + ylab("") + theme(axis.text.x=element_blank()) + 
  scale_fill_manual(values=c('black', 'darkgray'))





# ******************************************************************************
#                      ( FIG SUPP ) HMM STEP AND ANGLE PARAMS
# ******************************************************************************

pacman::p_load(ConnMatTools, reshape2, geostats)

createData <- function(mod, stat, states=NULL, melt=FALSE) {
  dVM <- function( col ) unlist(vonMises(x, mu=col[1], kappa=col[2])) 
  dGamma <- function( col ) {
    params <- unlist(gammaParamsConvert(mean=col[1], sd=col[2]))
    return( unlist(dgamma(x, shape=params[3], scale=params[4])) )
  }
  
  # determine which function and params to use
  if (stat=="step") {
    data <- data.frame(mod$mle$step);  xmin=0; xmax=2; fun=dGamma;
  } else if (stat == "angle") {
    data <- data.frame(mod$mle$angle); xmin=-pi; xmax=pi; fun=dVM;
  } else { message('stat not found'); return(NA);}
  if (!is.null(states)) colnames(data) <- states
  
  # create distribution and melt if necessary
  x <- seq(xmin, xmax, length.out=1000)
  y <- apply(data, 2, fun)
  df <- data.frame(cbind(x, y))
  if (melt) df <- melt(df, id.vars='x', variable.name="state")
  return(df)
}

# get fitted parameters
stateNames <- c('ruminating', 'foraging', 'exploring')
dgam <- createData(m.roan.simp, 'step', stateNames, melt=TRUE) 
dvm  <- createData(m.roan.simp, 'angle', stateNames, melt=TRUE)

# plotting
plotStat <- function(data, lab) {
  # if (grepl('step', lab)) scales="free_y"
  scales = NULL
  ggplot(data=data, aes(x=x, y=value, color=state, group=state)) + 
    geom_line(linewidth=1) + 
    facet_wrap(~state, scales=scales, ncol=1) + 
    scale_color_brewer(palette="Dark2") + 
    xlab(lab) + ylab('density') + plot.theme 
}
plotStat(dgam, 'step length (km)') + plotStat(dvm, 'turning angle (rad)') + plot_layout(guides='collect')
ggsave(here("03_output", 'hmm', 'StateParameters.png'),
       width=13, height=9)



ggplot(oryx, aes(x=step, fill=STATE, group=STATE)) + 
  geom_histogram() + 
  facet_wrap(~STATE)
