# calculateAverageDistances.R
# Created 21 Feb 2023
# Margaret Swift <mes114@duke.edu>


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here("02_scripts", "utilities.R"))
load(antfile)
load()


# ******************************************************************************
#                        CALCULATING CHARACTERISTICS
# ******************************************************************************
# for each individual, what's the balance of time spent near water and not near
# water, on a seasonal and time-of-day basis? What about favored terrain, and
# average distance traveled from the previous point?

data <- collar.df.1h %>% 
  nog() %>% 
  mutate(HOUR=hour(DATETIME),
         MONTH=month(DATETIME),
         YEAR = year(DATETIME),
         MSS = interval("2013-08-01", DATE) %/% months(1))

makeTable <- function(df, groups=NULL) {
  tab <- df %>% 
    nog() %>% 
    group_by_at(vars(one_of(groups))) %>% 
    summarize(
            nat_distm = mean(WATER_DENS_300m),
            art_distm = mean(DIST_ART))
  return(tab)
}

budget_tbl_1 <- makeTable(data, c("SPECIES", "SEASON"))
budget_tbl_2 <- makeTable(data, c("SPECIES", "SEASON", "HOUR"))
budget_tbl_3 <- makeTable(data, c("SPECIES", "SEASON", "MONTH", "YEAR", "MSS"))



### AVERAGE DISTANCE FROM ARTIFICIAL WATERHOLES
p1 <- ggplot() + 
  geom_line(data=budget_tbl_2, 
            aes(x=HOUR, y=nat_distm, color=SEASON, group=SEASON), 
            linewidth=1) + 
  facet_wrap(~SPECIES) + plot.theme + 
  ylab('% time within 100m\n of natural water') + 
  ggtitle("natural waterhole visitation budget") +
  xlab('hour of day') +
  # xlab('') + guides(x="none") +
  scale_color_brewer(palette="Dark2", direction=-1)


### AVERAGE DISTANCE FROM ARTIFICIAL WATERHOLES
p3 <- budget_tbl_2 %>%
  ggplot(aes(x=HOUR, y=art_distm, color=SEASON, group=SEASON)) + 
  geom_line(linewidth=1) + 
  facet_wrap(~SPECIES) + plot.theme + 
  ylab('distance (m) from artificial water') + 
  xlab("hour of day") +
  ggtitle("artificial waterhole distances") +
  # scale_y_reverse() +
  scale_color_brewer(palette="Dark2", direction=-1)
p3