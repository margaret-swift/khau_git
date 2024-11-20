# getting metadata from file and summarizing

pacman::p_load(lubridate, tidyverse)
setDataPaths('antelope')
file <- file.path(metapath, "Roan_Oryx data.csv")
data <- read.csv(file)
size <- as.numeric( unlist( lapply(data$Remarks, function(e) 
           str_extract_all(e, "[0-9]+")[[1]][1]) ) )
size <- ifelse(grepl('single', tolower(data$Remarks)), 1, size)

collar.meta <- data %>% 
  mutate(Date.Start = dmy(Date.Deployed),
         Date.End = dmy(Date.stopped),
         N.Months = interval(Date.Start, Date.End) %/% months(1),
         Herd.Size = ifelse(is.na(size), "HERD", size),
         Sp = ifelse(Animal == "Oryx", "Og", "He"),
         Animal.ID = paste(Sp, Collar.ID, sep="_"),
         Is.Solitary = Herd.Size == 1 ) %>%
  dplyr::select(Collar.ID, Animal.ID, Type, Frequency, Animal, Sex, Length.Age,
                Date.Start, Date.End, N.Months, Herd.Size, Is.Solitary,
                Replaced.with., Area.Deployed)
names(collar.meta) <- toupper(gsub('\\.', '_', names(collar.meta)))

#fix the two animals who had collars replaced
olds.inx <- which(!is.na(collar.meta$REPLACED_WITH_))
news.inx <- match(collar.meta$REPLACED_WITH_[olds.inx], collar.meta$COLLAR_ID)
collar.meta$ANIMAL_ID[olds.inx] <- collar.meta$ANIMAL_ID[news.inx]

outfile <- file.path(procpath, 'Collar_Metadata.RData')
save(collar.meta, file=outfile)


### Trait data ###
file <- "./TraitDataPantheria.csv"
traits <- read.csv(file) %>%
  mutate(adult_mass_kg = round(adult_mass_g / 1000, 0),
         female_maturity_m = round(female_maturity_d / 30, 2),
         male_maturity_m = round(male_maturity_d / 30, 2),
         gestation_length_m = round(gestation_length_d / 30, 2),
         interbirth_interval_m = round(interbirth_interval_d / 30, 2),
         weaning_age_m = round(weaning_age_d / 30, 2),
         ) %>%
  dplyr::select(order, family, genus, species, 
                adult_mass_kg, brain_mass_g,
                female_maturity_m, male_maturity_m,
                gestation_length_m, interbirth_interval_m,
                litter_size_n, litters_per_year_n, weaning_age_m,
                dispersal_km, social_group_n, 
                upper_elevation_m, biogeographical_realm,
                )
write.csv(traits, file='./TraitDataPantheriaClean.csv')
