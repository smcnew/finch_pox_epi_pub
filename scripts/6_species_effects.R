# For revision:
# Investigating the role of each species in contributing to overall pox prevalence
#

library(piecewiseSEM) # for SEMs
library(lme4) # for GLMs

library(ggthemr) # plot params
library(ggplot2) # plotting
library(patchwork) # plotting
library(sjPlot) # for saving tables

library(vegan) # Calculate shannon's diversity per day
library(lubridate) # wrangling dates

library(boot) # for calculating 95% CIs

library(tidyr)  # data wrangling
library(stringr) # data wrangling
library(dplyr) # data wrangling


# Plot parameters
ggthemr("flat", layout = "clear", text_size = 20) # clean is no grid lines


# Load data ---------------------------------------------------------------

captures <- read.csv("formatted_data/combined_2022_2023_finch_climate.csv") %>%
  mutate(year_month = as.Date(year_month)) %>%
  mutate(date = as.Date(date)) %>%
  mutate(year_week = as.Date(year_week)) %>%
  filter(sex != "N") %>% # get rid of nestlings
  filter(banded == 1) %>% # keep only banded birds
  filter(species != "ZEN") %>% # remove non-passerines
  filter(species != "MELA") %>%
  filter(species != "GALL") %>%
  filter(site != "Garrapatero") %>% # remove Garrapatero, only 13 captures
  mutate(pox_iur = factor(pox_iur, levels = c("U", "I", "R"),
                          labels = c("Uninfected", "Infected", "Recovered"))) %>% # relevel pox scale
  mutate(pox_scale = case_when(is.na(pox_scale) ~ "U", TRUE ~ pox_scale)) %>%
  mutate(pox_scale = factor(pox_scale, levels = c("U","A","B","C","D"))) %>%
  mutate(pox_size_lesion = case_when(is.na(pox_size_lesion) ~ 0, TRUE ~ pox_size_lesion)) %>% # replace NA with 0
  mutate(type_of_site = relevel(factor(type_of_site), ref = "bosque")) %>% # relevel type of site
  mutate(habitat = relevel(factor(habitat), ref = "natural")) %>% # relevel habitat
  mutate(taxon = factor(taxon, levels = c("Zenaida galapagoensis", "Coccyzus melacoryphus",
                                          "Myiarchus magnirostris", "Mimus parvulus",
                                          "Certhidea olivacea", "Platyspiza crassirostris",
                                          "Camarhynchus parvulus", "Camarhynchus psittacula",
                                          "Camarhynchus pallidus", "Geospiza magnirostris", "Geospiza fortis",
                                          "Geospiza fuliginosa", "Geospiza scandens",
                                          "Setophaga petechia"))) %>% # taxonomic order
  droplevels()

# Calculate shannon's diversity per day
dpd <- table(captures$date, captures$species) %>% diversity(., index = "shannon")
dpd_df <- data.frame(date = as.Date(names(dpd)), shannon_dpd = as.numeric(dpd))
captures <- left_join(captures, dpd_df) # Add to main data frame.



# Format data  ------------------------------------------------------------


# select columns for modeling, scale predictors, remove NAs
fdat <- captures %>%
  select(species,
         templ2,
         precipl2,
         humidityl2,
         site,
         inf_stat,
         bird_dens,
         species,
         shannon_dpd,
         elevation,
         date,
         habitat) %>%
         na.omit()


# Daily proportion of each species ----------------------------------------

prop <- captures %>%
  count(date, species) %>%
  group_by(date) %>%
  mutate(
    total_daily_caps = sum(n),
    prop = n / total_daily_caps
  ) %>%
  select(date, species, prop) %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = prop,
    names_prefix = "prop_",
    values_fill = 0
  )


# check that the proportion sums to 1.
apply(prop[,2:13], 1, sum)

# Add proportion information to main captures df
fdat <- left_join(fdat, prop) %>%
            mutate(across(where(is.numeric),
                ~ as.numeric(scale(.x)),
                .names = "scale_{.col}"))
head(fdat)
#
# Function to run psem for species ----------------------------------------
# list of all our species
spp <- captures$species %>% unique

# A function to run a psem for each species.
# Step 1. Filter dataset to exclude records for that species.
# Step 2. Run a psem to test whether the proportion of that species significantly
# affects the proportion of infection in other species.
#
# Component functions:
# Main function: prevalence ~ target species abundance + environment.
# Second function: target species abundance ~ environment.
psem_spp <- function (our_species, data) {
  ddat <-  data %>%
    filter(species!= our_species) %>%
    pivot_wider(names_from = inf_stat, values_from = inf_stat, values_fn = length,
                values_fill = 0, names_prefix = "inf") %>%
    as.data.frame %>%
    mutate (total = inf0 + inf1) %>%
    mutate (prop_inf = inf1/total)
  # compile column name for target species
  prop_var <- paste0("scale_prop_", our_species)

  f <- as.formula(
    paste(
      "prop_inf ~",
      prop_var,
      "+ scale_templ2 + scale_elevation + scale_precipl2 + scale_humidityl2 + ",
      " + (1|site)"
    )
  )

  mod1 <- glmer(
    f,
    data = ddat,
    weights = total,
    family = "binomial",
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)))

  f2 <- as.formula(
    paste0(
      prop_var,
      "~ scale_templ2 + scale_elevation + scale_precipl2 + scale_humidityl2  + ",
      " + (1|site)"
    )
  )

  mod2 <- lmer(f2, data = ddat)

  fitted_psem <- psem(
    mod1,
    mod2,
    data = ddat)
  return(fitted_psem)
}

fits_psem <- lapply(spp, function(x) psem_spp(our_species = x, data = fdat))


# summarize coefficients
psem_species_coefs <- lapply(fits_psem, function(x) {
  summary(x)$coefficients[, -9, drop = FALSE] %>%
    as.data.frame()
})


#
# Bootstrap to get 95% CIs on Std. Estimate -------------------------------




# make a function to process the data for each species
process_dat <- function(data, our_species) {
  data <-  data %>%
    filter(species!= our_species) %>%
    pivot_wider(names_from = inf_stat, values_from = inf_stat, values_fn = length,
                values_fill = 0, names_prefix = "inf") %>%
    as.data.frame %>%
    mutate (total = inf0 + inf1) %>%
    mutate (prop_inf = inf1/total) %>%
    mutate(site = as.factor(site)) %>%
  return(data)
}


# make a function to bootstrap the processed data
my_stat <- function(data, indices, our_species) {
  d <- data[indices, , drop = FALSE]
  # compile column name for target species
  prop_var <- paste0("scale_prop_", our_species)

  f <- as.formula(
    paste(
      "prop_inf ~",
      prop_var,
      "+ scale_templ2 + scale_elevation + scale_precipl2 + scale_humidityl2 + ",
      " + (1|site)"
    )
  )

    mod1 <- glmer(
    f,
    data = d,
    weights = total,
    family = "binomial",
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)))

  f2 <- as.formula(
    paste0(
      prop_var,
      "~ scale_templ2 + scale_elevation + scale_precipl2 + scale_humidityl2 + ",
      " + (1|site)"
    )
  )

  mod2 <- lmer(f2, data = d)

  fitted_psem <- psem(
    mod1,
    mod2,
    data = d)

  coefficients <- fitted_psem %>% summary %>% .$coefficients %>% .$Std.Estimate
  return(coefficients)
}


# Test functions on one species
process_dat(data = fdat, our_species = "FOR") %>%
  boot(data = ., statistic = my_stat, our_species = "FOR", R = 10, strata = .$site)

# Apply function across all 12 species

boot_all <- lapply(spp, FUN = function(x) process_dat(data = fdat, our_species = x) %>%
         boot(data = ., statistic = my_stat, our_species = x, R = 2000, strata = .$site))




# Calculate 95% CIs   -----------------------------------------


results_boot <- lapply(seq_along(boot_all), function(i) {

  x <- boot_all[[i]]

  data.frame(
    species = spp[i],
    boot_estimate = x$t0,
    std.error = summary(x)$bootSE,
    Bias = summary(x)$bootBias,
    Lower_CI = apply(x$t, 2, function(x) quantile(x, 0.025)),
    Upper_CI = apply(x$t, 2, function(x) quantile(x, 0.975))
  )
})




# Add bootstrap info to species coefs

psem_species_coefs_95 <- mapply(cbind, psem_species_coefs, results_boot, SIMPLIFY = FALSE) %>%
  bind_rows

# Sanity check that the model std. estimate and boot estimate are the same

psem_species_coefs_95$Std.Estimate == psem_species_coefs_95$boot_estimate # all T

# Write out model results for path diagram making


psem_species_coefs_95 %>%
  write.csv("output_tables/species_psem.csv", row.names = F)



# Plot effect of species on prevalence ------------------------------------


# Just pull the betas for species on infection
psem_species_coefs_plots <- psem_species_coefs_95 %>%
  filter(str_detect(Predictor, "scale_prop") ) %>%
  mutate(species = spp)



# Add right taxon information
taxon <- captures %>%
          select(taxon, species) %>%
          unique()

psem_species_coefs_plots <- left_join(psem_species_coefs_plots, taxon)


# Adjust those P values
psem_species_coefs_plots <- psem_species_coefs_plots %>%
  mutate(Adjusted.P = p.adjust(P.Value, method = "bonferroni"))

# Write out results for a table
psem_species_coefs_plots %>%
  select(taxon, Estimate, Std.Error, DF, Crit.Value, P.Value, Adjusted.P,
         Std.Estimate, Lower_CI, Upper_CI) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  rename("Taxon" = "taxon") %>%
  write.csv("output_tables/sem_species_effects.csv", row.names = F)


# Plot results  -----------------------------------------------------------

pdf("output_plots/effect_of_bird_on_pox.pdf", width = 10)
psem_species_coefs_plots %>%
  ggplot(aes(x = Std.Estimate, y = taxon)) +
  geom_linerange(
    aes(
      xmin = Lower_CI,
      xmax = Upper_CI
    ),
    linewidth = 0.75,
    orientation = "y",
    color = "#34495e"
  ) +

  # ± 1 SE
  geom_linerange(
    aes(
      xmin = Std.Estimate - std.error,
      xmax = Std.Estimate + std.error
    ),
    linewidth = 1.5,
    orientation = "y",
    color = "#34495e",

  ) +

  geom_point(size = 3, color = "#34495e") +

  labs(
    x = "Standardized Effect (β ± 95% CI)",
    y = "Species"
  ) +
  theme(axis.text.y = element_text(face = "italic")) +
scale_y_discrete(
    labels = c(
      "Geospiza fuliginosa" = "Geospiza fuliginosa*",
      "Geospiza scandens" = "Geospiza scandens*",
      "Geospiza magnirostris" = "Geospiza magnirostris*"
    ))

dev.off()



# Tiny bit of reformatting
coefs <- read.csv("output_tables/sem_species_effects.csv")

coefs <- coefs %>% mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(coefs, "output_tables/sem_species_effects.csv", row.names = F)
