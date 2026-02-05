### mark
### Trying multi-state models in RMark....
# Notes about package installation:
# Followed instructions here for installing MARK
# http://www.phidot.org/software/mark/downloads/
#
library(RMark)
library(tidyr)
library(dplyr)


# load data  --------------------------------------------------------------

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
alldates <- captures$date %>% unique
alldates == sort(alldates) # all true

# Create vectors of time intervals: how many days apart were our sampling points?
# Take off the last NA from the difference in time between the last day and the following non-existent day.
time_int_all <- (lead(alldates) - alldates) %>%
  head(., -1) %>%
  as.vector

# 2022
dates_2022 <- d22$date %>% unique
time_int_2022 <- (lead(dates_2022) - dates_2022) %>%
  head(., -1) %>%
  as.vector
# 2023
dates_2023 <- d23$date %>% unique
time_int_2023 <- (lead(dates_2023) - dates_2023) %>%
  head(., -1) %>%
  as.vector

# Format capture histories ------------------------------------------------
d23 %>% head
capt_hist23 <- d23 %>%
  select(date, band, species, pox_iur) %>%
  distinct(band, .keep_all = TRUE) %>%
  pivot_wider(names_from = date, values_from = pox_iur, values_fill = "0") %>%
  group_by(band) %>%
  unite("ch", 3:ncol(.), sep = "") %>%
  mutate(species = as.factor(species)) %>%
  ungroup %>%
  dplyr::select(-band) %>%
  as.data.frame()

mstrata.processed=process.data(capt_hist23, model="Multistrata",
                               strata.labels = c("I", "R", "U"),# labels for strata
                               #groups=c("InitAge"),# set initial ages
                               time.intervals=time_int_2023,# set intervals (prop. of mean)
                               nocc=length(time_int_2023)) # number of capture occasions

mstrata.ddl = make.design.data(mstrata.processed)
# examine data
head(mstrata.ddl$S)
head(mstrata.ddl$p)
head(mstrata.ddl$Psi)

# survival probability
S.site.time = list(formula = ~ time * stratum)

# detection probability
p.site.time = list(formula = ~ time * stratum)

# transition probs
Psi.site.time = list(formula = ~ time * stratum)

S.sitetime.p.sitetime.psi.sitetime = mark(mstrata.processed, mstrata.ddl,
                                          model.parameters = list(S = S.site.time, p = p.site.time, Psi = Psi.site.time),
                                          output=TRUE, delete=TRUE)





