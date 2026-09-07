# Script to analyze pox prevalence using partial structural equation models
# Author: Sabrina McNew
# Last update: Feb 2026

# For stats
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(vegan) # for calculating Shannon's diversity
library(mgcv) # for gam
#library(multcompView)
library(boot) # for bootstrapping psems

# For plotting
library(visreg)
library(ggplot2)
library(patchwork)
#devtools::install_github('Mikata-Project/ggthemr')
library(ggthemr)

# For wrangling
library(tidyr) # data manipulation
library(dplyr) # data manipulation

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
  mutate(site = factor(site)) %>%
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
  mutate(habitat_num = as.numeric(habitat)) %>%
  mutate(numeric_date = as.numeric(date))

ddat22 <- filter(ddat, year == 2022, habitat !="arid") %>% droplevels()
ddat23 <- filter(ddat, year != 2022, habitat == "arid") %>% droplevels()
#
# Variation in Shannon's Diversity among habitats  -----------------------------------


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

#
# Identify best fit time lag for environmental variables ---------------------------------------------------

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

# Format output tables for habitat model
smodd$coefficients %>%
  select(-9) %>%
  mutate(across(.cols = c(Response, Predictor), .fns =
                  ~recode(.x,
                          prop_inf = "Pox prevalence",
                          scale_dens = "Host abundance",
                          scale_shannons = "Host diversity",
                          scale_templ2 = "Temperature",
                          scale_elev = "Elevation",
                          scale_precipl2 = "Precipitation",
                          scale_humidityl2 = "Humidity"))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  write.csv("output_tables/sem_environment_habitat_factor.csv", row.names = FALSE)

# Elevation model
modde <- psem(
  glmer(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_shannons + scale_elev + (1|site), weights = total,
        data = ddat, family = "binomial"),
  lmer(scale_dens ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_elev +  (1|site) , data = ddat),
  lmer(scale_shannons ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_elev + (1|site) , data = ddat), # take out habitat to improve model convergence.
  data = ddat)
smodde <- summary(modde)
smodde


# Bootstrap models for 95CIs for elevation --------------------------------------


boot_elevation_sem <- function(data, indices) {
  d <- data[indices, , drop = FALSE]
  modde <- psem(
    glmer(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_shannons + scale_elev + (1|site), weights = total,
          data = d, family = "binomial",
          control = glmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 200000))),
    lmer(scale_dens ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_elev + (1|site), data = d),
    lmer(scale_shannons ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_elev + (1|site), data = d), # take out habitat to improve model convergence.
    data = d)
coefficients <- modde %>% summary %>% .$coefficients %>% .$Std.Estimate
return(coefficients)
}


full_elev_model_boot <- boot(data = ddat, statistic = boot_elevation_sem, R = 2000, strata = ddat$site)


# format boot results
boot_results <- data.frame(
  boot_std_estimate = full_elev_model_boot$t0,
  boot_bias = summary(full_elev_model_boot)$bootBias,
  boot_mean = apply(full_elev_model_boot$t, 2, function(x) mean(x)),
  Lower_CI = apply(full_elev_model_boot$t, 2, function(x) quantile(x, 0.025)),
  Upper_CI = apply(full_elev_model_boot$t, 2, function(x) quantile(x, 0.975))
)


#Some sanity checks
# test <- smodde$coefficients %>%
#   select(-9) %>%
#   cbind(boot_results)
#
# test$Std.Estimate == test$boot_std_estimate
# test$Upper_CI > test$Std.Estimate
# test$Lower_CI < test$Std.Estimate

# Format output tables for elevation model
smodde$coefficients %>%
select(-9) %>%
  cbind(boot_results) %>%
  select(-boot_std_estimate, -boot_mean, -boot_bias) %>%
  mutate(across(.cols = c(Response, Predictor), .fns =
                  ~recode(.x,
                          prop_inf = "Pox prevalence",
                          scale_dens = "Host abundance",
                          scale_shannons = "Host diversity",
                          scale_templ2 = "Temperature",
                          scale_elev = "Elevation",
                          scale_precipl2 = "Precipitation",
                          scale_humidityl2 = "Humidity"))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  write.csv("output_tables/sem_environment_elevation.csv", row.names = F)




#
# Visualizations -----------------------------------------------------------------
# Use visreg to visualize differences among habitats

prev_model <- glmer(cbind(inf1, inf0) ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + habitat + (1|site), data = ddat, family = "binomial")


summary(prev_model)$coefficients


visreg(prev_model, "habitat", scale = "response")
predicted_vals <- visreg(prev_model,"habitat", scale = "response", plot = FALSE)$fit
predicted_resids <- visreg(prev_model,"habitat", scale = "response", plot = FALSE)$res


# Visreg plot with partial residuals on top
pdf("output_plots/visreg_habitatxprev_part_resids.pdf", width = 10)
ggplot() +
  geom_tile(data = predicted_vals, aes(
    x= habitat, y = visreg_fit), width = .9, height = .005) +
  geom_tile(data = predicted_vals, aes(
    x= habitat, y = visreg_fit, height = visreg_upr - visreg_lwr, width = .9), alpha = 0.5, show.legend = FALSE, fill = "#34495e") +
  geom_jitter(data = predicted_resids, aes(x = habitat, y = visreg_res), width = .1) +
  labs(x = "Habitat", y = "Pox Prevalence") +
  scale_x_discrete(labels = c("natural" = "Forest", "arid" = "Coast/Urban",
                              "coffee" = "Coffee", "fruit_vegetable" = "Fruit/Veg",
                              "pasture" = "Pasture"))
dev.off()



# PSEMs of divided datasets -----------------------------------------------
#2022

modd22 <- psem(
  glmer(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_shannons + scale_elev + (1|site), weights = total,
        data = ddat22, family = "binomial",
        control = glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 200000))),
  lmer(scale_dens ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_elev + (1|site), data = ddat22),
  lmer(scale_shannons ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_elev + (1|site), data = ddat22), # take out habitat to improve model convergence.
  data = ddat22)
smodd22 <- summary(modd22)
smodd22

# Boot results

full_elev_model_boot22 <- boot(data = ddat22, statistic = boot_elevation_sem, R = 2000, strata = ddat22$site)

boot_results22 <- data.frame(
  boot_std_estimate = full_elev_model_boot22$t0,
  boot_bias = summary(full_elev_model_boot22)$bootBias,
  boot_mean = apply(full_elev_model_boot22$t, 2, function(x) mean(x)),
  Lower_CI = apply(full_elev_model_boot22$t, 2, function(x) quantile(x, 0.025)),
  Upper_CI = apply(full_elev_model_boot22$t, 2, function(x) quantile(x, 0.975))
)

# Do some sanity checks

# temp <- smodd22$coefficients %>%
#   select(-9) %>%
#   cbind(boot_results22)
# temp$Std.Estimate == temp$boot_std_estimate
# temp$Std.Estimate > temp$Lower_CI
# temp$Std.Estimate < temp$Upper_CI

# format results for export
smodd22$coefficients %>%
  select(-9) %>%
  cbind(boot_results22) %>%
  select(-boot_std_estimate, -boot_mean, -boot_bias) %>%
  mutate(across(.cols = c(Response, Predictor), .fns =
                  ~recode(.x,
                          prop_inf = "Pox prevalence",
                          scale_dens = "Host abundance",
                          scale_shannons = "Host diversity",
                          scale_templ2 = "Temperature",
                          scale_elev = "Elevation",
                          scale_precipl2 = "Precipitation",
                          scale_humidityl2 = "Humidity"))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  write.csv("output_tables/sem_environment_elevation_2022.csv", row.names = F)


# 2023

modd23 <- psem(
   glm(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_shannons ,
        weights = total, data = ddat23, family = "binomial"),
  lm(scale_dens ~ scale_templ2 + scale_humidityl2 + scale_precipl2 , data = ddat23),
  lm(scale_shannons ~ scale_templ2 + scale_humidityl2 + scale_precipl2 + scale_dens, data = ddat23),
  data = ddat23)
smodd23 <- summary(modd23)

# Boot results, make a new function without random effects

boot_elevation_sem23 <- function(data, indices) {
d <- data[indices, , drop = FALSE]
modd23 <- psem(
  glm(prop_inf ~ scale_templ2 + scale_precipl2 + scale_humidityl2 + scale_dens + scale_shannons ,
      weights = total, data = d, family = "binomial"),
  lm(scale_dens ~ scale_templ2 + scale_humidityl2 + scale_precipl2 , data = d),
  lm(scale_shannons ~ scale_templ2 + scale_humidityl2 + scale_precipl2 + scale_dens, data = d),
  data = d)
coefficients <- modd23 %>% summary %>% .$coefficients %>% .$Std.Estimate
return(coefficients)
}

full_elev_model_boot23 <- boot(data = ddat23, statistic = boot_elevation_sem23, R = 2000)


boot_results23 <- data.frame(
  boot_std_estimate = full_elev_model_boot23$t0,
  boot_bias = summary(full_elev_model_boot23)$bootBias,
  boot_mean = apply(full_elev_model_boot23$t, 2, function(x) mean(x)),
  Lower_CI = apply(full_elev_model_boot23$t, 2, function(x) quantile(x, 0.025)),
  Upper_CI = apply(full_elev_model_boot23$t, 2, function(x) quantile(x, 0.975))
)

# # Some sanity checks
# temp <- smodd23$coefficients %>%
#   select(-9) %>%
#   cbind(boot_results23)
# temp$Std.Estimate == temp$boot_std_estimate
# temp$Std.Estimate > temp$Lower_CI
# temp$Std.Estimate < temp$Upper_CI

# Format results for export
smodd23$coefficients %>%
  select(-9) %>%
  cbind(boot_results23) %>%
  select(-boot_std_estimate, -boot_mean, -boot_bias) %>%
  mutate(across(.cols = c(Response, Predictor), .fns =
                  ~recode(.x,
                          prop_inf = "Pox prevalence",
                          scale_dens = "Host abundance",
                          scale_shannons = "Host diversity",
                          scale_templ2 = "Temperature",
                          scale_precipl2 = "Precipitation",
                          scale_humidityl2 = "Humidity"))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  write.csv("output_tables/sem_environment_elevation2023.csv", row.names = F)






