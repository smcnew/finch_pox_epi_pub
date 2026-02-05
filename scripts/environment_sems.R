# Script to analyze pox prevalence using partial structural equation models
# Author: Sabrina McNew
# Last update: June 2026

# For stats
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(vegan)
library(mgcv) # for gam

# For plotting
library(visreg)
library(ggplot2)
library(ggthemr)

# For wrangling
library(gt)
library(tidyr)
library(dplyr)

# Plot params
ggthemr("flat", layout = "clear", text_size = 20) # clean is no grid lines

# Load data ---------------------------------------------------------------
# mean_temp is for the month; mean_air_temp is daily


captures <- read.csv("formatted_data/combined_2022_2023_finch_climate.csv") %>%
  mutate(year_month = as.Date(year_month)) %>%
  mutate(date = as.Date(date)) %>%
  filter(sex != "N") %>% # get rid of nestlings
  #filter(banded == 1) %>% # keep only banded birds
  #filter(species != "ZEN") %>% # remove non-passerines
  #filter(species != "MELA") %>%
  #mutate(pox_iur = factor(pox_iur, levels = c("U", "I", "R"),
  #                        labels = c("Uninfected", "Infected", "Recovered"))) %>% # relevel pox scale
  #mutate(pox_scale = case_when(is.na(pox_scale) ~ "U", TRUE ~ pox_scale)) %>%
  #mutate(pox_scale = factor(pox_scale, levels = c("U","A","B","C","D"))) %>%
  #mutate(pox_size_lesion = case_when(is.na(pox_size_lesion) ~ 0, TRUE ~ pox_size_lesion)) %>% # replace NA with 0
  #mutate(type_of_site = relevel(factor(type_of_site), ref = "bosque")) %>% # relevel type of site
  mutate(habitat = relevel(factor(habitat), ref = "natural")) %>% # relevel habitat
  mutate(taxon = factor(taxon, levels = c("Zenaida galapagoensis" ,"Coccyzus melacoryphus",
                                          "Myiarchus magnirostris", "Mimus parvulus",
                                          "Certhidea olivacea", "Platyspiza crassirostris",
                                          "Camarhynchus parvulus", "Camarhynchus psittacula",
                                          "Camarhynchus pallidus", "Geospiza fortis",
                                          "Geospiza fuliginosa", "Geospiza scandens",
                                          "Setophaga petechia"))) %>% # taxonomic order
  droplevels()


# Calculate daily Shannon's Diversity
dpd <- table(captures$date, captures$species) %>% diversity(., index = "shannon")
dpd_df <- data.frame(date = as.Date(names(dpd)), shannon_dpd = as.numeric(dpd))
captures <- left_join(captures, dpd_df) # Add to main data frame.

glimpse(captures)
# Select columns to model, and remove NAs. Scale the predictor columns

fdat <- captures %>%
  select(templ2, templ4, templ6, templ8, templ10,
         precipl2, precipl4, precipl6, precipl8, precipl10,
         humidityl2, humidityl4, humidityl6, humidityl8, humidityl10, humidity,
         habitat, site, mean_air_temp, week_precip, precipitation, year,
         inf_stat, bird_dens, species, date, shannon_dpd, elevation) %>%
  mutate(across(c(templ2, templ4, templ6, templ8, templ10,
                  precipl2, precipl4, precipl6, precipl8, precipl10, humidity,
                  humidityl2, humidityl4, humidityl6, humidityl8, humidityl10,
                  mean_air_temp, week_precip, precipitation,
                  bird_dens, shannon_dpd, elevation),
                  ~ as.numeric(scale(.x)),
                .names = "scale_{.col}")) %>%
  rename(scale_elev = scale_elevation,
         scale_dens = scale_bird_dens,
         scale_shannons = scale_shannon_dpd,
         scale_temp = scale_mean_air_temp,
         scale_precip = scale_week_precip) %>%
  na.omit()


# Mutate data so that each day is a row
# Select rows that we'll use to model,
# Pivot wider to create two columns- infected and uninfected totals per day
ddat <-  fdat %>%
  pivot_wider(names_from = inf_stat, values_from = inf_stat, values_fn = length,
              values_fill = 0, names_prefix = "inf") %>%
  as.data.frame %>%
  mutate (total = inf0 + inf1) %>%
  mutate (prop_inf = inf1/total) %>%
  mutate(habitat_num = as.numeric(habitat))

# variation in shannon's among habitat  -----------------------------------

head(ddat)
ddat <- ddat %>% mutate(numeric_date = as.numeric(date))
ddat22 <- ddat[ddat$year==2022,]
ddat23 <- ddat[ddat$year!=2022,]

# Models of habitat and date
diversity_habitat <- lmer(shannon_dpd ~ habitat + (1|site), data = ddat)
summary(diversity_habitat)
# Date 2022
diversity_date22 <- gam(shannon_dpd ~ s(numeric_date), data = ddat22)
summary(diversity_date22)

# Date 2023
diversity_date23 <- gam(shannon_dpd ~ s(numeric_date), data = ddat23)
summary(diversity_date23)

d1 <-visreg(diversity_habitat, gg = T) +
  labs(y = "Shannon's diversity index")

d2 <- visreg(diversity_date22, gg = T ) +
  labs(y = "Shannon's diversity index", x = "Date (2022)")
d3 <- visreg(diversity_date23, gg = T ) +
  labs(y = "Shannon's diversity index",  x = "Date (2023-2024)")
pdf("output_plots/diversity_habitats_time.pdf", height = 10, width = 10)
d1/(d2 + d3)
dev.off()

# check VIF elev habitat --------------------------------------------------

#r2 <- summary(lm(elevation ~ habitat, data = ddat))$r.squared
#vif_x1 <- 1 / (1 - r2)  # highly correlated

# PSEMs for pox prevalence ---------------------------------------------------

# Identify best time lag for precipitation
lag0 <- psem(
  glmer(prop_inf ~ scale_temp + scale_precip + scale_humidity + scale_dens + scale_shannons + habitat + (1|site),
        weights = total, data = ddat, family = "binomial"))

lag2 <- psem(
  glmer(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidity + scale_dens + scale_shannons + habitat + (1|site),
        weights = total, data = ddat, family = "binomial"))

lag4 <- psem(
  glmer(prop_inf ~ scale_templ4 + scale_precipl4 + scale_humidity + scale_dens + scale_shannons + habitat + (1|site),
        weights = total, data = ddat, family = "binomial"))

lag6 <- psem(
  glmer(prop_inf ~ scale_templ6 + scale_precipl6 + scale_humidity + scale_dens + scale_shannons + habitat + (1|site),
        weights = total, data = ddat, family = "binomial"))

lag8 <- psem(
  glmer(prop_inf ~ scale_templ8 + scale_precipl8 + scale_humidity + scale_dens + scale_shannons + habitat + (1|site),
        weights = total, data = ddat, family = "binomial"))

lag10 <- psem(
  glmer(prop_inf ~ scale_templ10 + scale_precipl10 + scale_humidity + scale_dens + scale_shannons + habitat + (1|site),
        weights = total, data = ddat, family = "binomial"))


AIC(lag0, lag2, lag4, lag6, lag8, lag10)# 2 week lag is the best
AIC(lag0, lag2, lag4, lag6, lag8, lag10) %>%
  as.data.frame %>%
  mutate(Lag_weeks = c(0,2,4,6,8,10)) %>%
  select(Lag_weeks, AIC,K,n) %>%
  gt %>%
  gtsave("output_tables/AIC_selection.rtf")




# # Inspect individual models
# summary(lag0)
# summary(lag2)
# summary(lag4)
# summary(lag6)
# summary(lag8)
# summary(lag10)
# summary(lagmonth)

#
# Full PSEMS --------------------------------------------------------------


# Habitat model
modd <- psem(
  glmer(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + habitat + scale_dens + scale_shannons + (1|site),
        weights = total, data = ddat, family = "binomial"),
  lmer(scale_dens ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + habitat +  (1|site) , data = ddat),
  lmer(scale_shannons ~ scale_templ2 + scale_precipl2 + scale_humidityl2+ scale_dens + habitat + (1|site) , data = ddat), # take out habitat to improve model convergence.
  data = ddat)
smodd <- summary(modd)
smodd

# Elevation model
modde <- psem(
  glmer(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_shannons + scale_elev + (1|site),
        weights = total, data = ddat, family = "binomial"),
  lmer(scale_dens ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_elev +  (1|site) , data = ddat),
  lmer(scale_shannons ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_elev + (1|site) , data = ddat), # take out habitat to improve model convergence.
  data = ddat)
smodde <- summary(modde)



# Format output tables for habitat model
smodd$coefficients %>%
  select(-9) %>%
  gt() %>%
  gtsave("output_tables/sem_environment_habitat_factor.rtf")


# Format output tables for elevation model
smodde$coefficients %>%
select(-9) %>%
  gt() %>%
  gtsave("output_tables/sem_environment_elevation.rtf")

plot(modde)

# # Create a model with habitat coded as a numerical variable to get standardized
# # estimates for beta coefficients of the other variables.
#
# modd2 <- psem(
#   glmer(prop_inf ~ habitat_num +scale_templ2 + scale_precipl2 + scale_dens + scale_shannons + (1|site),
#         weights = total, data = ddat, family = "binomial"),
#   lmer(scale_dens ~ scale_templ2 + scale_precipl2 + habitat_num + (1|site), data = ddat), # take out random effect bc it causes singularity error
#   lmer(scale_shannons ~ scale_templ2 + scale_precipl2 + scale_dens +habitat_num+ (1|site), data = ddat),
#   data = ddat)
# smodd2 <- summary(modd2)
# smodd2
#
#
# # Save coefficients
# smodd2$coefficients %>%
#   select(-9) %>%
#   gt %>%
#   gtsave("output_tables/sem_environment_habitat_num.rtf")
#
#
# # Note, psem plots won't write to a pdf IDK why.
# plot(modd2)



# Visualizations -----------------------------------------------------------------
# Use visreg to visualize differences among habitats

prev_model <- glmer(cbind(inf1, inf0) ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + habitat + (1|site), data = ddat, family = "binomial")


summary(prev_model)$coefficients

str(summary(prev_model))
visreg(prev_model, "habitat", scale = "response")
predicted_vals <- visreg(prev_model,"habitat", scale = "response")$fit
predicted_resids <- visreg(prev_model,"habitat", scale = "response")$res

# base visreg plot
#visreg(prev_model,"habitat", scale = "response", partial = TRUE, gg= T) +
#  labs(y = "Prevalence")


# Visreg plot just means and confidence bands
predicted_vals %>%
  ggplot(aes (x= habitat, y = visregFit)) +
  geom_tile(width = .9) +
  geom_tile(aes(height = visregUpr - visregLwr, width = .9),
  alpha = 0.5, show.legend = FALSE, fill = "#34495e") +
  labs(x = "Habitat", y = "Pox prevalence")


# Visreg plot with partial residuals on top
pdf("output_plots/visreg_habitatxprev_part_resids.pdf", width = 10)
ggplot() +
  geom_tile(data = predicted_vals, aes(
    x= habitat, y = visregFit), width = .9, height = .005) +
  geom_tile(data = predicted_vals, aes(
    x= habitat, y = visregFit, height = visregUpr - visregLwr, width = .9), alpha = 0.5, show.legend = FALSE, fill = "#34495e") +
  geom_jitter(data = predicted_resids, aes(x = habitat, y = visregRes), width = .1) +
  labs(x = "Habitat", y = "Pox Prevalence") +
  scale_x_discrete(labels = c("natural" = "Forest", "arid" = "Coast/Urban",
                              "coffee" = "Coffee", "fruit_vegetable" = "Fruit/Veg",
                              "pasture" = "Pasture"))
dev.off()

  # Visreg plot with raw data on top (less good)
# pdf("output_plots/visreg_habitatxprev_rawdat.pdf", width = 10)
#   ggplot() +
#     geom_tile(data = predicted_vals, aes(
#       x= habitat, y = visregFit), width = .9, height = .005) +
#     geom_tile(data = predicted_vals, aes(
#       x= habitat, y = visregFit, height = visregUpr - visregLwr, width = .9), alpha = 0.5, show.legend = FALSE, fill = "#34495e") +
#      geom_jitter(data = ddat, aes(x = habitat, y = prop_inf), width = .1) +
#     labs(x = "Habitat", y = "Pox Prevalence") +
#     scale_x_discrete(labels = c("natural" = "Forest", "arid" = "Coast/Urban",
#                                 "coffee" = "Coffee", "fruit_vegetable" = "Fruit/Veg",
#                                 "pasture" = "Pasture"))
# dev.off()
#

# PSEMs of divided datasets -----------------------------------------------
#2022
ddat22 <- filter(ddat, year == 2022, habitat !="arid") %>% droplevels()
modd22 <- psem(
  glmer(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_shannons + (1|site),
        weights = total, data = ddat22, family = "binomial"),
  lmer(scale_dens ~ scale_templ2 + scale_humidityl2 + scale_precipl2 +  (1|site) , data = ddat22),
  lmer(scale_shannons ~ scale_templ2 + scale_humidityl2 + scale_precipl2 + scale_dens + (1|site) , data = ddat22),
  data = ddat22)
smodd22 <- summary(modd22)
smodd22

smodd22$coefficients %>%
  select(-9) %>%
  gt() %>%
  gtsave("output_tables/sem_environment_habitat_2022.rtf")


# 2023
ddat23 <- filter(ddat, year != 2022, habitat == "arid")
modd23 <- psem(
  glm(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_shannons ,
        weights = total, data = ddat23, family = "binomial"),
  lm(scale_dens ~ scale_templ2 + scale_humidityl2 + scale_precipl2 , data = ddat23),
  lm(scale_shannons ~ scale_templ2 + scale_humidityl2 + scale_precipl2 + scale_dens, data = ddat23),
  data = ddat23)
smodd23 <- summary(modd23)


smodd23$coefficients %>%
  select(-9) %>%
  gt() %>%
  gtsave("output_tables/sem_environment_habitat_2023.rtf")



# understand coefficients better
test <- filter(captures, species %in% c("SCA", "FOR", "CRA")) %>%
  mutate(species = factor(species, levels = c("SCA", "FOR", "CRA")))
lm(mass ~ species, data = test) %>% emmeans (~ species) %>% cld()
psem(
  lm(mass ~ species, data = test),
  data = test
) %>% summary
lm(mass ~ species, data = test) %>% summary

mod <- glmer(prop_inf ~ scale_templ2 + scale_precipl2 + habitat + scale_dens + scale_shannons + (1|site),
      weights = total, data = ddat, family = "binomial")
emmeans(mod, ~ habitat)






