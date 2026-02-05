# Script to run multi-state models in RMark
# Notes about package installation:
# Followed instructions here for installing MARK
# http://www.phidot.org/software/mark/downloads/
library(RMark)
library(tidyr)
library(dplyr)


# Load capture histories --------------------------------------------------
# States = A(uninfected); B(infected)
capt_hist23_simp <- read.csv("capture_histories/capt_hist_23_simp.csv")

# Intervals between captures.
time_int_2023 <- read.csv("capture_histories/time_int_2023") %>% pull


# Process data  -----------------------------------------------------------
# Work with one species as an example: FOR

for23 <- filter(capt_hist23_simp, species == "FOR") %>% select(-species)
for.23.processed <- process.data(for23, model = "Multistrata",
                                 strata.labels = c("A","B"),
                                 time.intervals=time_int_2023)
for.23.ddl <- make.design.data(for.23.processed)

all_highcaps.processed <- process.data(capt_hist23_simp, model = "Multistrata",
                                       strata.labels = c("A","B"))

all_highcaps.ddl <- make.design.data(all_highcaps.processed)

fit_model_for <- mark(all_highcaps.processed, all_highcaps.ddl,
                      model.parameters = list(S = S.stratum,
                                              p = p.stratum,
                                              Psi = Psi.stratum), delete = TRUE)


#### ---- Fit the models ---- ####
#### Create the formula(s) for each model component
# Detection (formulas must use a lowercase "p")
# Decided to just model intercept and stratum, not time for initial model
p.intercept = list(formula = ~1)
p.stratum = list(formula = ~stratum)

# State transitions (formulas must use a captial "P")
Psi.intercept = list(formula = ~ 1)
Psi.stratum = list(formula = ~ stratum)

# Survival (formulas must use a captial "S" )
S.intercept = list(formula= ~1)
S.stratum = list(formula= ~stratum)

#### Run the MARK function to execute the model
# You can run a single model using the "mark" function and specifying the model formulas above, like so:

#Intercept only model
fit_model_1 = mark(for.23.processed, for.23.ddl,
                   model.parameters = list(S = S.intercept,
                                           p = p.intercept,
                                           Psi = Psi.intercept), delete = TRUE)
# Stratum model
fit_model_for <- mark(for.23.processed, for.23.ddl,
                   model.parameters = list(S = S.stratum,
                                           p = p.stratum,
                                           Psi = Psi.stratum), delete = TRUE)





