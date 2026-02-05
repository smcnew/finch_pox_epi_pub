# Script for multi-state mark recapture analysis

# Data --------------------------------------------------------------------
library(marked)
library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)

# Data --------------------------------------------------------------------
# Two states: uninfected (a) and infected (b)
capt_hist22_simp <- read.csv("capture_histories/capt_hist_22_simp.csv")
capt_hist23_simp <- read.csv("capture_histories/capt_hist_23_simp.csv")
capt_hist24_simp <- read.csv("capture_histories/capt_hist_24.csv")

# Other data frames that have 3 states (IUR)
# capt_hist23 <- read.csv("capture_histories/capt_hist_23.csv")


# Intervals between captures
time_int_2023 <- read.csv("capture_histories/time_int_2023") %>% pull
time_int_2022 <- read.csv("capture_histories/time_int_2022") %>% pull

#
# Model parameters --------------------------------------------------------
# Probability of detection
p.pox = list(formula=~time) # Pox status
p.int = list(formula=~1) # intercept

# Probability of transitions between states
Psi.pox = list(formula=~stratum) # Pox status
Psi.int = list(formula=~1) # intercept

# Probability of survival
S.pox = list(formula=~stratum)
S.int = list(formula=~1)

# create model list with all combos of parameters
cml <- create.model.list(c("Psi","p", "S"))

#
# Step-by-step functions-----------------------------------


d24.proc <- process.data(capt_hist24,model="Mscjs",
                         strata.labels = c("A","B"), time.intervals = capt_hist24_ints) # group by pox
d24.ddl <- make.design.data(d24.proc)

# Basic model
mod2 <- crm(d24.proc, d24.ddl, method="nlminb",
            model.parameters = list(S = S.pox, p = p.pox, Psi = Psi.pox),
            hessian = FALSE, use.tmb =  TRUE)
predict(mod2, newdata = data.frame(pox = c("A", "B")), se = T)

# Run all combinations of parameters
results1 <- crm.wrapper(cml,
                        data = d24.proc, ddl = d24.ddl,
                        hessian = FALSE, use.tmb =  TRUE, external = FALSE)

results1
# Model selection function ---------------------------------------------------
# Create a function to do model selection for different combos of data

model_select <- function(data, time_ints) {
  d.proc <- process.data(data,
                         model="Mscjs",
                         strata.labels = c("A","B"),
                         time.intervals = time_ints) # group by pox
  d.ddl <- make.design.data(d.proc)
  results<- crm.wrapper(cml,
             data = d.proc, ddl = d.ddl,
             hessian = FALSE, use.tmb =  TRUE, external = FALSE)
  return(results)
}

model <- function(capture_history, time_ints, S = S.int, Psi = Psi.int, p = p.int, hessian = FALSE) {
  d.proc <- process.data(capture_history,model="Mscjs",
                         strata.labels = c("A","B"), time.intervals = time_ints) # group by pox
  d.ddl <- make.design.data(d.proc)
  model <- crm(d.proc, d.ddl, method="nlminb",
               model.parameters = list(S = S, p = p, Psi = Psi),
               hessian = hessian, use.tmb =  TRUE)
  return(model)
}

pull_coefs <- function(model){
  predictions <- predict(model, newdata = data.frame(pox = c("A", "B")), se = T)
  return(predictions$S)
}

# 2023 variations ---------------------------------------------------------

FOR23 <- capt_hist23_simp %>% filter(species == "FOR")
FUL23 <- capt_hist23_simp %>% filter(species == "FUL")
GAMO23 <- capt_hist23_simp %>% filter(species == "GAMO")
PAR23 <- capt_hist23_simp %>% filter(species == "PAR")
CRA23 <- capt_hist23_simp %>% filter(species == "CRA")

formods23 <- model_select(FOR23, time_int_2023)
fulmods23 <- model_select(FUL23, time_int_2023)
gamomods23 <- model_select(GAMO23, time_int_2023)
parmods23 <- model_select(PAR23, time_int_2023)
cramods23 <- model_select(CRA23, time_int_2023)

formod <- model(FOR23, time_int_2023, S=S.pox, p = p.pox, hessian = TRUE)


d.proc <- process.data(GAMO23, # data for mockingbirds
                       model="Mscjs",
                       strata.labels = c("A","B"), # infected = A; uninf = B
                       time.intervals = time_int_2023) # time interval
d.ddl <- make.design.data(d.proc)

# All model params = list(formula=~stratum)
gamomod <- crm(d.proc, d.ddl,
             model.parameters = list(S = S.int, p = p.pox, Psi = Psi.int),
             hessian = TRUE, use.tmb =  TRUE, useHess = TRUE, debug = TRUE)

gamomod$results





#
# 2022 variations ---------------------------------------------------------

FOR22 <- capt_hist22_simp %>% filter(species == "FOR")
FUL22 <- capt_hist22_simp %>% filter(species == "FUL")
GAMO22 <- capt_hist22_simp %>% filter(species == "GAMO")
PAR22 <- capt_hist22_simp %>% filter(species == "PAR")
CRA22 <- capt_hist22_simp %>% filter(species == "CRA")
OLI22 <- capt_hist22_simp %>% filter(species == "OLI")

formods22 <- model_select(FOR22, time_int_2022)
fulmods22 <- model_select(FUL22, time_int_2022)
gamomods22 <- model_select(GAMO22, time_int_2022)
parmods22 <- model_select(PAR22, time_int_2022)
cramods22 <- model_select(CRA22, time_int_2022)
olimods22 <- model_select(OLI22, time_int_2022)

gamomods22[[4]] %>% predict(newdata = data.frame(pox = c("A", "B"), se = T))

#
#

# 2024 variations ---------------------------------------------------------
FOR24 <- capt_hist24 %>% filter(species == "FOR")
FUL24 <- capt_hist24 %>% filter(species == "FUL")
GAMO24 <- capt_hist24 %>% filter(species == "GAMO")
PAR24 <- capt_hist24 %>% filter(species == "PAR")
CRA24 <- capt_hist24 %>% filter(species == "CRA")
#SCA24 <- capt_hist24 %>% filter(species == "SCA")
all24 <- capt_hist24

formods24 <- model_select(FOR24, capt_hist24_ints)
fulmods24 <- model_select(FUL24, capt_hist24_ints)
gamomods24 <- model_select(GAMO24, capt_hist24_ints)
parmods24 <- model_select(PAR24, capt_hist24_ints)
cramods24 <- model_select(CRA24, capt_hist24_ints)
#scamods24 <- model_select(SCA24, capt_hist24_ints)
allmods24 <- model_select(all24, capt_hist24_ints)

m1 <- model(FOR24, capt_hist24_ints, S=S.int, p = p.pox)
m2 <- model(FUL24, capt_hist24_ints, S=S.pox)
m3 <- model(CRA24, capt_hist24_ints, S=S.pox)
m4 <- model(PAR24, capt_hist24_ints, S=S.pox)
m5 <- model(capt_hist24, capt_hist24_ints, S=S.pox)
pull_coefs(m5)
betas <- pull_coefs(m1) %>% mutate(species = "FOR")
betas <- pull_coefs(m2) %>% mutate(species = "FUL") %>% bind_rows(., betas)
betas <- pull_coefs(m3) %>% mutate(species = "CRA") %>% bind_rows(., betas)
betas <- pull_coefs(m5) %>% mutate(species = "ALL") %>% bind_rows(., betas)

# make some sort of figure
betas %>% filter(stratum == "B") %>%
  ggplot(aes(x = species, y = estimate, color = stratum)) +
  geom_point() +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.2)



# Run other variations with limited data from 2022

formods22lim <- capt_hist_22_limited %>% filter(species =="FOR" ) %>% model_select(., capt_hist_22_ints)
fulmods22lim <- capt_hist_22_limited %>% filter(species == "FUL") %>% model_select(., capt_hist_22_ints)
olimods22lim <- capt_hist_22_limited %>% filter(species == "OLI") %>% model_select(., capt_hist_22_ints)
parmods22lim <- capt_hist_22_limited %>% filter(species == "PAR" ) %>% model_select(., capt_hist_22_ints)
palmods22lim <- capt_hist_22_limited %>% filter(species == "PAL") %>% model_select(., capt_hist_22_ints)
table(capt_hist_22_limited$species)

m6 <- capt_hist_22_limited %>% model(., capt_hist_22_ints, hessian = TRUE, S = S.pox)
pull_coefs(m6)


# run a test with deer  ---------------------------------------------------

deer <- read.csv("input_data/deer.csv", header = T)
# This data describes observations of individuals in two possible states over time: breeder and non-breeder
# column headers must be "ch" and "freq"

#### Create processed data frame
deer.processed <- process.data(deer,
                               model="Mscjs",
                               strata.labels=c("B","N"))

#### Make design data frame from processed data set
deer.ddl <- make.design.data(deer.processed)

p.intercept = list(formula = ~1)
p.time = list(formula = ~time) # RMark has a series of columns already built in to the design matrix that you can take advantage of
p.stratum = list(formula = ~stratum)
p.stratum.time = list(formula = ~stratum + time)

# State transitions (formulas must use a captial "P")
Psi.intercept = list(formula = ~ 1)
Psi.time = list(formula = ~ time)
Psi.stratum = list(formula = ~ stratum)
Psi.stratum.time = list(formula = ~ stratum + time)

# Survival (formulas must use a captial "S" )
S.intercept = list(formula= ~1)
S.time = list(formula= ~time)
S.stratum = list(formula= ~stratum)
S.stratum.time = list(formula= ~stratum + time)

#### Run the MARK function to execute the model
# You can run a single model using the "mark" function and specifying the model formulas above, like so:
deer = mark(deer.processed, deer.ddl,
            model.parameters = list(S = S.stratum,
                                    p = p.stratum,
                                    Psi = Psi.stratum), delete = TRUE)
deer_mscjs <- crm(
  model = "multistrata",
  deer.processed, deer.ddl,
 model.parameters = list(S = S.stratum,
       p = p.stratum,
       Psi = Psi.stratum), use.tmb =TRUE, hessian = TRUE,
  accumulate = FALSE)
predict(deer_mscjs, newdata = data.frame(strat = c("A", "B")), se = T)


# this example requires admb
# The same example is in the RMark package and it is included here to
# illustrate the differences in the handling of mlogit parameters between RMark
# and marked.  The MARK software handles parameters like Psi which must sum to 1
# by excluding one of the cells that is used as a reference cell and is computed by
# subtracting the other cell values from 1 so the total sums to 1.  This is often
# handled with an mlogit parameter in which the cell values are exp(beta) and the
# reference cell is set to 1 and the values are divided by the sum across the cells
# so the resulting values are probabilities that sum to 1. In marked, instead of removing
# one of the cells, all are included and the user must select which should be the
# reference cell by setting the value fix=1 for that cell and others are NA so they are
# estimated. For transition parameters like Psi, the default design data is setup so
# that the probability of remaining in the cell (stratum=tostratum) is the reference cell
# and fix set to 1.  Thus, this means 2 changes are needed to the script in RMark.
# The first is to remove the statement skagit.ddl$Psi$fix=NA because that over-rides
# the default fix values.  The other is to add
# skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="B"&
#  skagit.ddl$Psi$time==5]=0
# to change the value from 1 to 0 which forces movement from B to A in the interval 5 to 6. If
# this is not done then Psi B to B=Psi B to A=0.5 because each is 1 and when they are normalized
# they are divided by the sum which is 2 (1/2).
#if(!is(try(setup_admb("mscjs")),"try-error"))
#{
  data(skagit)
  skagit.processed=process.data(skagit,model="Mscjs",groups=c("tag"),strata.labels=c("A","B"))
  skagit.ddl=make.design.data(skagit.processed)
  head(skagit)
  #
  # p
  #
  # Can't be seen at 5A or 2B,6B (the latter 2 don't exist)
  skagit.ddl$p$fix=ifelse((skagit.ddl$p$stratum=="A"&skagit.ddl$p$time==5) |
                            (skagit.ddl$p$stratum=="B"&skagit.ddl$p$time%in%c(2,6)),0,NA)
  # Estimated externally from current data to allow estimation of survival at last interval
  skagit.ddl$p$fix[skagit.ddl$p$tag=="v7"&skagit.ddl$p$time==6&skagit.ddl$p$stratum=="A"]=0.687
  skagit.ddl$p$fix[skagit.ddl$p$tag=="v9"&skagit.ddl$p$time==6&skagit.ddl$p$stratum=="A"]=0.975
  #
  # Psi
  #
  # only 3 possible transitions are A to B at time interval 2 to 3 and
  # for time interval 3 to 4 from A to B and from B to A
  # rest are fixed values
  ############ change for RMark to marked; remove next line
  #skagit.ddl$Psi$fix=NA
  # stay in A for intervals 1-2, 4-5 and 5-6
  skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="A"&
                       skagit.ddl$Psi$tostratum=="B"&skagit.ddl$Psi$time%in%c(1,4,5)]=0
  # stay in B for interval 4-5
  skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"
                     &skagit.ddl$Psi$time==4]=0

  # leave B to go to A for interval 5-6
  skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"&
                       skagit.ddl$Psi$time==5]=1
  ############ change for RMark to marked; add next line to set B to B to 0 otherwise it has
  ############ been set to 1 by default which would make psi B to B = psi B to A = 0.5
  skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="B"&
                       skagit.ddl$Psi$time==5]=0
  # "stay" in B for interval 1-2 and 2-3 because none will be in B
  skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"&
                       skagit.ddl$Psi$time%in%1:2]=0
  #
  # S
  #
  # None in B, so fixing S to 1
  skagit.ddl$S$fix=ifelse(skagit.ddl$S$stratum=="B"&skagit.ddl$S$time%in%c(1,2),1,NA)
  skagit.ddl$S$fix[skagit.ddl$S$stratum=="A"&skagit.ddl$S$time==4]=1
  # fit model
  p.timexstratum.tag=list(formula=~time:stratum+tag,remove.intercept=TRUE)
  Psi.sxtime=list(formula=~-1+stratum:time)
  S.stratumxtime=list(formula=~-1+stratum:time)
  #
  mod1=crm(skagit.processed,skagit.ddl,
           model.parameters=list(S=S.stratumxtime,p= p.timexstratum.tag,Psi=Psi.sxtime),hessian=TRUE,
           use.tmb = TRUE)
  if(!is(mod1,"try-error")) mod1
}
