# Script to compile 2022 and 2023 data and clean it
# Sabrina McNew
# Updated May 2025
library(tidyr)
library(vegan)
library(dplyr)
library(lubridate) # modify dates

# Load data ---------------------------------------------------------------

dat_2022 <- read.csv("input_data/base_de_datos_aves2022.csv") %>% rename_with(tolower)
dat_2023 <- read.csv("input_data/Campo2023_base_datos.csv") %>% rename_with(tolower)
data <- bind_rows(dat_2022, dat_2023)


head(data)
dim(data) #5552 total rows

# Process data  -----------------------------------------------------------

data <- data %>%
  mutate(taxon = case_when(species == "CRA" ~ "Platyspiza crassirostris", # add latin name
                           species == "FOR" ~ "Geospiza fortis",
                           species == "FUL" ~ "Geospiza fuliginosa",
                           species == "GAMO" ~ "Mimus parvulus",
                           species == "MAG" ~ "Geospiza magnirostris",
                           species == "MELA" ~ "Coccyzus melacoryphus",
                           species == "MYI" ~ "Myiarchus magnirostris",
                           species == "OLI" ~ "Certhidea olivacea",
                           species == "PAL" ~ "Camarhynchus pallidus",
                           species == "PET" ~ "Setophaga petechia",
                           species == "PSI" ~ "Camarhynchus psittacula",
                           species == "SCA" ~ "Geospiza scandens",
                           species == "ZEN" ~ "Zenaida galapagoensis",
                           species == "PAR" ~ "Camarhynchus parvulus",
                           species == "GALLARETA" ~ "Neocrex erythrops",
                           species == "ANI" ~ "Crotophaga ani")) %>%
  mutate(type_of_site = tolower(type_of_site)) %>%
  mutate(type_of_site = case_when(type_of_site == "humeda" ~ "bosque", #synonymize habitat types
                                  type_of_site == "" ~ "agricola",
                                  site == "Maria Elena Guerra" ~ "agricola", # per Ilke 2025
                                  site == "Reina Oleas" ~ "agricola", # per Ilke 2025
                                  TRUE ~ type_of_site)) %>%
  mutate(site = case_when(site == "El Barranco" ~ "Barranco", # synonymize micro-sites within CDRS
                          site == "FCD" ~ "Barranco",
                          site == "Corrales" ~ "Barranco",
                          TRUE ~ site)) %>%
  filter(!is.na(band)) %>%
  filter(!is.na(pox_iur)) %>%
  filter(!is.na(species)) %>%
  group_by(band, species) %>%
  mutate(ncaps = n()) %>% # number of total captures for this individual
  mutate(capturen = 1:n()) %>% # which capture was this for the individual
  mutate(n_estados = n_distinct(pox_iur)) %>% # number of pox states, i.e. if captured I and then U = 2
  as.data.frame()

# Add additional site information --------------------------------------------

# Shannon's diversity per site
# Create a species x column data frame

shannons <- table(data$site, data$species) %>%
  diversity(., index = "shannon")
shannon_df <- data.frame(site = names(shannons), shannons = as.numeric(shannons))

# Add additional information about type of site
sites <- read.csv("input_data/sites.csv") %>% left_join(., shannon_df)
data <- left_join(data, sites)

boxplot(shannons ~ habitat, data = sites)
# Process data to create net hours
# Use lubridate package to turn hours from character into time format

data <- data %>%
  mutate(rain_time = ifelse(is.na(rain_time), "00:00", rain_time)) %>%
  mutate(rain_time = ifelse(rain_time == "", "00:00", rain_time)) %>%
  mutate(across(c(net_open, net_close, rain_time), hm)) %>%
  mutate(operating_hours = time_length(net_close - net_open - rain_time, unit = "hour")) %>%
  mutate(net_hours = operating_hours * net_no)

table(data$net_hours)

# calculate total captures per day
bird_dens <- data %>%
  group_by(date) %>%
  summarize(total_daily_caps = n())


# add daily captures to data
data <- left_join(data, bird_dens)


# calculate density as daily captures / net hours
data <- mutate(data, bird_dens = total_daily_caps / net_hours)


# create a column called "banded" that shows whether the bird was banded or not

data <- data %>%
  mutate (
    banded = case_when(
      grepl("FUL", band) ~ 0,
      grepl("FUL", band) ~ 0,
      grepl("OLI", band) ~ 0,
      grepl("PET", band) ~ 0,
      grepl("MYI", band) ~ 0,
      grepl("ANI", band) ~ 0,
      grepl("PAR", band) ~ 0,
      grepl("GALLARETA", band) ~ 0,
      grepl("NA-GALL", band) ~ 0,
      grepl("NO ANILLADO", band) ~ 0,
      grepl("UNBANDED", band) ~ 0,
      band == "" ~ 0,
      TRUE ~ 1)
          )

# Fix date formatting,
# Make sure data are sorted by date
data <-
  data %>%
  mutate(date = as.Date(paste(day, month, year, sep = "/"), format = "%d/%m/%Y")) %>%
  arrange(date) %>%
  mutate(year_month = floor_date(date, unit = "month")) %>% # add id column for year x month
  mutate(year_week = floor_date(date, unit = "week")) %>%  # add id column for week of year
  mutate(day_of_year = yday(date))  # correct small errors in day of year


# Add a column that indicates inf_stat; 1 = infected, 0 is uninfected
data <-
  data %>% mutate (inf_stat = case_when(pox_iur == "U" ~ 0,
                                        pox_iur == "R" ~ 0,
                                        pox_iur == "I" ~ 1))




# Calculate Body Condition ---------------------------------------------------------------
# First remove Anis from study
data <- data %>% filter(species != "ANI")

smi_slope <- function(thedata) {
  # Calculate Beta coefficient from OLS regression
  Bols <- lm(log(mass) ~ log(tarsus), data = thedata)$coefficients[2]
  # Calculate R coefficient
  pearsons_r <- cor(
    x = log(thedata$tarsus),
    y = log(thedata$mass),
    use = "complete.obs",
    method = "pearson"
  )
  smi_slope = Bols / pearsons_r # calculate smi_slope
  smi_slope = as.vector(smi_slope)
  return(smi_slope)
}

# Calculate a slope for each species
bsmi <- data %>% group_by(species) %>%
  do(smis = smi_slope(.)) %>%
  mutate(smis = unlist(smis)) %>%
  as.data.frame

# Calculate mean tarsus for each species and add that to species-slopes
meant <- data %>% group_by(species) %>%
  summarize(meant = mean(tarsus, na.rm = T))
bsmi <- left_join(bsmi, meant)

# Add those slopes on to the big data frame
data <- left_join(data, bsmi)

# Calculate individual body condition
# condition = mass * (mean_tarsus/tarsus)^smi

data <- data %>%
  mutate(condition = mass *((meant/tarsus)^smis))




# Export data -------------------------------------------------------------


# exclude strange fuliginosa with 60g mass
data <- data %>% filter(band != "SM827")

write.csv(data, "formatted_data/combined_2022_2023_finch.csv", row.names = FALSE)


# Compile locality info ---------------------------------------------------

locals <- data %>% select(site, latitude, longitude, habitat) %>%
  distinct(site, .keep_all = TRUE) %>%
  mutate(habitat = case_when(site == "Garrapatero" ~ "arid", TRUE ~ habitat)) %>%
  na.omit()
write.csv(locals, "formatted_data/locals.csv", row.names = F)


# Check tarsus and mass  --------------------------------------------------
# Subsequent analysis of tarsus and mass suggested that some entries should be
# checked in the data. Code below looks for outliers in these two metrics and flags
# entries for checking. March 2024: exported file for checking and subsequently
# corrected a few (~25 entries) based on paper data sheets.


# Procedure:
# find a way to check the expected mass and tarsus for each species and
# flag outliers
# calculate species-specific quantiles for mass and tarsus
# find "upper" and "lower" bounds for each species equal to 1.5 x the
# interquartile range, subtracted or added to the first or third quartile

# mass_bounds <- aggregate(mass ~ species, data = data, fivenum) %>%
#   do.call("data.frame", .) %>%
#   mutate(m.lowerbound = mass.2 - (1.5 * (mass.4 - mass.2))) %>%
#   mutate(m.upperbound = mass.4 + (1.5 * (mass.4 - mass.2)))
#
# tarsus_bounds <- aggregate(tarsus ~ species, data = data, fivenum) %>%
#   do.call("data.frame", .) %>%
#   mutate(t.lowerbound = tarsus.2 - (1.5 * (tarsus.4 - tarsus.2))) %>%
#   mutate(t.upperbound = tarsus.4 + (1.5 * (tarsus.4 - tarsus.2)))
#
# both_bounds <- left_join(mass_bounds[,c("species", "m.lowerbound", "m.upperbound")],
#                          tarsus_bounds[,c("species", "t.lowerbound", "t.upperbound")])
#
# # add bounds onto the main data frame, create a "flag mass" and "flag tarsus"
# # column to id problem entries
# data <- left_join(data, both_bounds) %>%
#         mutate(flag_mass = ifelse(mass < m.lowerbound | mass > m.upperbound, 1, 0)) %>%
#         mutate(flag_tarsus = ifelse(tarsus < t.lowerbound | tarsus > t.upperbound, 1,0))
#
# data %>% filter(flag_mass == 1 ) %>% dim
# data %>% filter(flag_tarsus == 1 ) %>% dim



