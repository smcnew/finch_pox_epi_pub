#install.packages("marked") # first time using the package
#
library(marked)
library(tidyr)
library(dplyr)
library(magrittr)
#library(RMark)


# Load data ---------------------------------------------------------------
# Keep only banded adults
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
# For just 2022
capt_hist22<- d22 %>%
  select(date, band, species, inf_stat) %>%
  distinct() %>% #remove duplicates for the day
  mutate(detect = 1) %>% # add a column that shows each bird was detected on that day
    spread(date, detect, fill = 0 ) %>% # fill in 0s to show days that birds were not detected
    group_by(band) %>%
  unite("ch", 4:ncol(.), sep = "") %>%# unite columns 4:end,
  mutate(inf_stat = as.factor(inf_stat)) %>%
  ungroup %>%
  dplyr::select(-band) %>%
  as.data.frame

capt_hist23 <- d23 %>%
  select(date, band, species, inf_stat) %>%
  distinct() %>% #remove duplicates for the day
  mutate(detect = 1) %>% # add a column that shows each bird was detected on that day
  spread(date, detect, fill = 0 ) %>% # fill in 0s to show days that birds were not detected
  group_by(band) %>%
  unite("ch", 4:ncol(.), sep = "") %>% # unite columns 4:end, creates the ch vector of 01 detection history
  mutate(inf_stat = as.factor(inf_stat)) %>%
  ungroup %>%
  dplyr::select(-band) %>%
  as.data.frame()



# CJS wrapper functions ---------------------------------------------------
# Create a function that compares all models
fit.bird.cjs.models <- function (data, ddl) {
  # Phi = apparent survival
  Phi.pox = list(formula=~inf_stat) # survival differs between pox status
  Phi.dot = list(formula=~1) # constant survival
  # p = detection probability
  p.pox = list(formula=~inf_stat) #detection differs between infection status
  p.dot = list(formula=~1) # constant detection rates

  cml <- create.model.list(c("Phi", "p"))
  results <- crm.wrapper(cml,
                         data = data,
                         ddl = ddl,
                         external = FALSE,
                         accumulate = FALSE,
                         hessian = TRUE)
  return(results)
}

# Create another function, throw the data processing and make design steps in there
fit.bird.cjs.modelsb <- function (capture_history, time_int) {
  data <- process.data(capture_history, group = "inf_stat", time.intervals = time_int) # group by pox
  ddl <- make.design.data(data)

  # Phi = apparent survival
  Phi.pox = list(formula=~inf_stat) # survival differs between pox status
  Phi.dot = list(formula=~1) # constant survival
  # p = detection probability
  p.pox = list(formula=~inf_stat) #detection differs between infection status
  p.dot = list(formula=~1) # constant detection rates

  cml <- create.model.list(c("Phi", "p"))
  results <- crm.wrapper(cml,
                         data = data,
                         ddl = ddl,
                         external = FALSE,
                         accumulate = FALSE,
                         hessian = TRUE)
  return(results)
}


# Data Analysis -------------------------------------------------------

# process data for 2022
d22.proc <- process.data(capt_hist22, group = "inf_stat", time.intervals = time_int_2022) # group by pox
d22.ddl <- make.design.data(d22.proc)

# Create model set function
year_22mods <- fit.bird.cjs.models(data = d22.proc, ddl = d22.ddl)
year_22mods[[3]]
predict(year_22mods[[3]])

# try creating bayesian model by hand
Phi.pox = list(formula=~inf_stat)
p.pox = list(formula=~inf_stat)

y22_mcmc <- crm(
  model = "probitCJS",
  d22.proc,
  d22.ddl,
  model.parameters = list(Phi = Phi.pox,
                          p = p.pox),
  accumulate = FALSE, burnin=1000,iter=5000)

predict(y22_mcmc)

# process data for 2023
d23.proc <- process.data(capt_hist23, group = "inf_stat", time.intervals = time_int_2023) # group by pox
d23.ddl <- make.design.data(d23.proc)
year_23mods <- fit.bird.cjs.models(data = d23.proc, ddl = d23.ddl)
year_23mods[[4]]
predict(year_23mods[[4]])

# comparing JS vs CJS vs mcmcs
y23_JS <- crm(
  model = "JS",
  d23.proc,
  d23.ddl,
  model.parameters = list(Phi = Phi.pox,
                          p = p.pox),
  accumulate = FALSE)

y23_CJS <- crm(
  model = "CJS",
  d23.proc,
  d23.ddl,
  model.parameters = list(Phi = Phi.pox,
                          p = p.pox),
  accumulate = FALSE)


y23_mcmc <- crm(
  model = "probitCJS",
  d23.proc,
  d23.ddl,
  model.parameters = list(Phi = Phi.pox,
                          p = p.pox),
  burnin=1000,iter=5000)



data(dipper)
Flood=matrix(rep(c(0,1,1,0,0,0),each=nrow(dipper)),ncol=6)
colnames(Flood)=paste("Flood",1:6,sep="")
dipper=cbind(dipper,Flood)
design.parameters=list(Phi=list(time.varying="Flood"))
model.parameters=list(Phi=list(formula=~Flood), p=list(formula=~time+sex))
MCMCfit=crm(dipper,model="probitCJS",
          model.parameters=model.parameters,
          design.parameters=design.parameters,
          burnin=1000,iter=5000)

Phi.pox = list(formula=~inf_stat) # survival differs between pox status
Phi.dot = list(formula=~1) # constant survival
# p = detection probability
p.pox = list(formula=~inf_stat) #detection differs between infection status
p.dot = list(formula=~1) # constant detection rates

design.parameters=list(Phi=list(~1))
model.parameters=list(Phi=list(formula=~inf_stat), p=list(formula=~inf_stat))

# this takes about 15 min on the laptop
MCMC23=crm(capt_hist23,model="probitCJS",
            model.parameters=model.parameters,
            design.parameters=design.parameters,
            burnin=100,iter=500)

#"JS" and "CJS" seem to be the same;

y23_JS
predict(y23_JS)
y23_CJS
predict(y23_CJS)
MCMC23
predict(MCMC23) #predicted survival for Phi is lower. Still don't know what Occ means

# subsets of animals  -----------------------------------------------------

#Finches
finches <- filter(capt_hist23, species %in% c("FOR", "FUL", "CRA", "SCA","PAR"))
finch_mods <- fit.bird.cjs.modelsb (finches)
finch_mods[[3]]

finch.proc <- process.data(finches)
finch.ddl <- make.design.data(finch.proc)
finch_mods2 <- fit.bird.cjs.models(data = finch.proc, ddl = finch.ddl)

finches <- finches %>% mutate(species = as.factor(species))
finch.proc <- process.data(finches, group = "species")
finch.ddl <- make.design.data(finch.proc)
finch_mods3 <- fit.bird.cjs.models(data = finch.proc, ddl = finch.ddl)
predict(finch_mods3[[3]])


# Models for individual species -------------------------------------------
# GAMO
gamomods <- filter(capt_hist23, species == "GAMO") %>%
            fit.bird.cjs.modelsb(time_int = time_int_2023)
gamomods[[3]] %>% predict

# FOR
formods <- filter(capt_hist23, species == "FOR") %>%
          fit.bird.cjs.modelsb(time_int = time_int_2023)
formods[[2]] %>% predict

# FUL
fulmods <- filter(capt_hist23, species == "FUL") %>%
  fit.bird.cjs.modelsb(time_int = time_int_2023)
fulmods[[3]] %>% predict

# SCA: probably not enough data?
scamods <- filter(capt_hist23, species == "SCA") %>%
  fit.bird.cjs.modelsb(time_int_2023)
scamods[[4]] %>% predict

# PAR: also not a ton of data
parmods <- filter(capt_hist23, species == "PAR") %>%
  fit.bird.cjs.modelsb(time_int_2023)
parmods[[3]] %>% predict

#CRA
cramods <- filter(capt_hist23, species == "CRA") %>%
  fit.bird.cjs.modelsb(time_int_2023)
cramods
# Models with species as covariate/group ----------------------------------

fit.bird.cjs.modelsc <- function (capture_history) {
  data <- process.data(capture_history) # group by pox
  ddl <- make.design.data(data)

  # Phi = apparent survival
  Phi.pox = list(formula=~inf_stat*species) # survival differs between pox status
  Phi.dot = list(formula=~1) # constant survival
  # p = detection probability
  p.pox = list(formula=~inf_stat*species) #detection differs between infection status
  p.dot = list(formula=~1) # constant detection rates

  cml <- create.model.list(c("Phi", "p"))
  results <- crm.wrapper(cml,
                         data = data,
                         ddl = ddl,
                         external = FALSE,
                         accumulate = FALSE,
                         hessian = TRUE)
  return(results)
}

# test with a small dataset

test_finches <-  filter(capt_hist23, species %in% c("GAMO", "CRA", "FOR")) %>%
                  mutate(species = as.factor(species))

test_sppcomp <- fit.bird.cjs.modelsc(test_finches)

test_sppcomp[[4]] %>% predict
# Think about infection as time varying  ----------------------------------
d23 %>%
  select(date, band, species, inf_stat) %>%
  distinct() %>% #remove duplicates for the day
  mutate(detect = 1)  %>% head


# add a column that shows each bird was detected on that day
  spread(date, detect, fill = 0 ) %>% # fill in 0s to show days that birds were not detected
  group_by(band) %>%
  unite("ch", 4:ncol(.), sep = "") %>% # unite columns 4:end, creates the ch vector of 01 detection history
  mutate(inf_stat = as.factor(inf_stat)) %>%
  ungroup %>%
  dplyr::select(-band) %>%
  as.data.frame()
