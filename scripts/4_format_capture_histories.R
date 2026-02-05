# Script to format capture histories for multi-state mark recap
# analysis in R and Rmark

library(tidyr)
library(dplyr)

# Load data
captures <- read.csv("formatted_data/combined_2022_2023_finch.csv") %>%
  filter(sex != "N") %>%
  filter(banded == 1) %>%
  mutate(date = as.Date(date))

# Split data set in to highlands (2022) and lowlands (2023-2024)
d22 <- captures %>% filter(year == 2022)
d23 <- captures %>%
  filter(year == 2023 | year == 2024) %>%
  filter(zone == "Puerto Ayora" )

# Check and make sure that date column is in date order
# alldates <- captures$date %>% unique
# alldates == sort(alldates) # all true

# How many days apart were our sampling points?
# Take off the last NA from the difference in time between the last day and the following non-existent day.
# time_int_all <- (lead(alldates) - alldates) %>%
#  head(., -1) %>%
#  as.vector

# Capture histories infected, recovered, uninfected ----------------------------

# 2023, including species, and including pox coded as I/U/R
capt_hist23 <- d23 %>%
  select(date, band, species, pox_iur) %>%
  distinct(date, band, date, .keep_all = TRUE) %>%
  pivot_wider(names_from = date, values_from = pox_iur, values_fill = "0") %>%
  group_by(band) %>%
  unite("ch", 3:ncol(.), sep = "") %>%
  mutate(species = as.factor(species)) %>%
  ungroup %>%
  dplyr::select(-band) %>%
  as.data.frame()

#2022 including species, and including pox coded as I/U/R
capt_hist22 <- d22 %>%
  select(date, band, species, pox_iur) %>%
  distinct(date, band, .keep_all = TRUE) %>%
  pivot_wider(names_from = date, values_from = pox_iur, values_fill = "0") %>%
  group_by(band) %>%
  unite("ch", 3:ncol(.), sep = "") %>%
  mutate(species = as.factor(species)) %>%
  ungroup %>%
  dplyr::select(-band) %>%
  as.data.frame()
head(capt_hist22)


# Function to format capture histories 2 levels (inf vs. uninf) ------------------------------------

format_cap_hist <- function(dataz) {
    dataz %>%
    mutate(inf_stat = ifelse(inf_stat == 0, "A", "B")) %>%
    select(date, band, species, inf_stat) %>%
    distinct(date, band, .keep_all = TRUE) %>%
    pivot_wider(names_from = date, values_from = inf_stat, values_fill = "0") %>%
    group_by(band) %>%
    unite("ch", 3:ncol(.), sep = "") %>%
    mutate(species = as.factor(species)) %>%
    ungroup %>%
    dplyr::select(-band) %>%
    as.data.frame()

}


format_time_int <- function(data) {
  dates <- data$date %>% unique
  ints <- (lead(dates) - dates) %>%
    head(., -1) %>%
    as.vector
  return(ints)
}

#
## Identify high recapture/pox times  -------------------

# 2022
pox_sum22 <- table(d22$year_month, d22$pox_iur) %>%
  as.data.frame.matrix %>%
  mutate(prev = I / (I + U + R)) %>%
  tibble::rownames_to_column("year_month") %>%
  mutate(year_month = as.numeric(year_month))

#2023
pox_sum23 <- table(d23$year_month, d23$pox_iur) %>%
  as.data.frame.matrix %>%
  mutate(prev = I / (I + U + R)) %>%
  tibble::rownames_to_column("year_month") %>%
  mutate(year_month = as.numeric(year_month))

# When were recapture rates highest?
pox_sum_22 <- aggregate(recaptured ~ year_month, data = d22, sum) %>% left_join(pox_sum22)
pox_sum_23 <- aggregate(recaptured ~ year_month, data = d23, sum) %>% left_join(pox_sum23)
30/(96+76+175)
# Create infected / uninfected capture histories ---------------------------
capt_hist22_simp <- d22 %>% format_cap_hist
capt_hist23_simp <- d23 %>% format_cap_hist
time_int_2022 <- d22 %>% format_time_int()
time_int_2023 <- d23 %>% format_time_int()

# Best of 2022
capt_hist_22_limited <- d22 %>%
  filter(as.character(year_month) %in% c("2022.58333333333", "2022.66666666667")) %>%
  format_cap_hist()

capt_hist_22_ints <- d22 %>%
  filter(as.character(year_month) %in% c("2022.58333333333", "2022.66666666667")) %>%
  format_time_int()

d22 %>%
  select(year_month) %>% distinct %>% mutate(year_month = as.character(year_month))

# Best of 2023/2024
capt_hist24_simp <- d23 %>% filter(year == 2024 ) %>% format_cap_hist()
time_int_2024_simp <- d23 %>% filter(year == 2024) %>% format_time_int()

 d23 %>% filter(year == 2024 ) %>%
  filter(species == "FOR") %>% count(ncaps, inf_stat)
 d23 %>%
   filter(species == "GAMO") %>%
   mutate(inf_stat = ifelse(inf_stat == 0, "A","B")) %>%
   count(ncaps, inf_stat)

 d22 %>%
   filter(species == "PAR") %>%
   distinct(band, .keep_all = TRUE) %>%
aggregate(ncaps ~ inf_stat, ., mean)


 # Export data -------------------------------------------------------------
write.csv(capt_hist22, "capture_histories/capt_hist_22.csv", row.names = F)
write.csv(capt_hist23, "capture_histories/capt_hist_23.csv", row.names = F)
write.csv(capt_hist22_simp, "capture_histories/capt_hist_22_simp.csv", row.names = F)
write.csv(capt_hist23_simp, "capture_histories/capt_hist_23_simp.csv", row.names = F)
write.csv(time_int_2022, "capture_histories/time_int_2022", row.names = F)
write.csv(time_int_2023, "capture_histories/time_int_2023", row.names = F)
write.csv(capt_hist24_simp, "capture_histories/capt_hist24_simp.csv", row.names = F)
write.csv(time_int_2024, "capture_histories/time_int_2024", row.names = F)


