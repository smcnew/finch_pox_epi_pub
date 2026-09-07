
#install.packages("piecewiseSEM")
#install.packages("effects")

library(piecewiseSEM)
library(effects)
library(semEff) # not helping.
library(boot)
library(lme4)
library(ggplot2) # Plotting
library(ggthemr)
library(gt) # to save output tables
library(dplyr) # wrangling

# Plot parameters
ggthemr("flat", layout = "clear", text_size = 20) # clean is no grid lines

#
# load data  --------------------------------------------------------------
# mean_temp is for the month; mean_air_temp is daily
captures <- read.csv("formatted_data/combined_2022_2023_finch_climate.csv") %>%
  mutate(year_month = as.Date(year_month)) %>%
  mutate(date = as.Date(date)) %>%
  filter(sex != "N") %>% # get rid of nestlings
  filter(banded == 1) %>% # keep only banded birds
  #filter(species != "ZEN") %>% # remove non-passerines
  #filter(species != "MELA") %>%
  mutate(pox_iur = factor(pox_iur, levels = c("U", "I", "R"),
                          labels = c("Uninfected", "Infected", "Recovered"))) %>% # relevel pox scale
  mutate(pox_scale = case_when(is.na(pox_scale) ~ "U", TRUE ~ pox_scale)) %>%
  mutate(pox_scale = factor(pox_scale, levels = c("U","A","B","C","D"))) %>%
  mutate(pox_size_lesion = case_when(is.na(pox_size_lesion) ~ 0, TRUE ~ pox_size_lesion)) %>% # replace NA with 0
  mutate(type_of_site = relevel(factor(type_of_site), ref = "bosque")) %>% # relevel type of site
  mutate(habitat = relevel(factor(habitat), ref = "natural")) %>% # relevel habitat
  mutate(taxon = factor(taxon, levels = c("Zenaida galapagoensis" ,"Coccyzus melacoryphus", #put birds in taxonomic order
                                          "Myiarchus magnirostris", "Mimus parvulus",
                                          "Certhidea olivacea", "Platyspiza crassirostris",
                                          "Camarhynchus parvulus", "Camarhynchus psittacula",
                                          "Camarhynchus pallidus", "Geospiza fortis",
                                          "Geospiza fuliginosa", "Geospiza scandens",
                                          "Setophaga petechia"))) %>%
  droplevels()


# create a z scaled variable of condition and mass
# per species so that effect sizes are comparable across species

captures <- captures %>%
  group_by(species) %>%
  mutate(condition_z = scale(condition)) %>%
  mutate(mass_z = scale(mass)) %>%
  ungroup()

# # Split data set in to highlands (2022) and lowlands (2023-2024)
# d22 <- captures %>%
#   filter(year == 2022) %>%
#   filter(site != "El Barranco") %>% # filter out 8 test birds from lowlands
#   mutate(type_of_site = relevel(factor(type_of_site), ref = "bosque")) %>%
#   group_by(species) %>% # exclude species with fewer than 10 records, (G. scandens)
#   filter(n() >= 10) %>%
#   ungroup() %>%
#   droplevels()


# focus on birds from just lowlands during second year of study
d23 <- captures %>%
  filter(year == 2023 | year == 2024) %>%
  filter(zone == "Puerto Ayora" ) %>%
  droplevels()

#
# SEM bird condition by species  ------------------------------------------
# Hypothesis: Bird condition may be affected by pox, but effects may vary among
# species. Temperature should be included in the model because it is potentially
# a common cause of both condition and pox lesions.
# Use just the six species for which we have good pox numbers.

# Format data for model, remove NAs
mdat <- d23 %>%
  filter(species %in% c("FOR","FUL","PAR","SCA", "GAMO", "CRA")) %>%
  select(condition, pox_size_lesion, species, mass, condition_z, mass_z, inf_stat,
         mean_air_temp, bander_id, taxon) %>%
  na.omit %>%
  mutate(scaled_temp = scale(mean_air_temp))

# PSM model
condition_spp_psm <- psem(
  lm(condition ~ pox_size_lesion + mean_air_temp, data = mdat),
  lm(pox_size_lesion ~ mean_air_temp, data = mdat),
  data = mdat)
summary(condition_spp_psm)

# Test for differences among species
# Species is a significant interactor between condition and lesion size
# Individual results significant just for SCA
multi_analysis <- multigroup(condition_spp_psm, group = "taxon")
multi_analysis
#summary(condition_spp_psm)


# Plot the DAG
plot(condition_spp_psm)


# Get the model outputs into a useful format
multi_analysis$group.coefs %>%
  bind_rows(, .id = "species") %>%
  rename("signif" = "...9", "constrained" = "...10") %>%
  gt(groupname_col = "species") %>%
  gtsave("output_tables/sem_multigroup_condition.rtf", )



# Plot with model estimates ---------------------------------

mod_coefs <- multi_analysis$group.coefs %>%
  bind_rows(, .id = "species") %>%
  filter(Predictor == "pox_size_lesion") %>%
  mutate(Estimate = as.numeric(Estimate)) %>%
  mutate(Std.Error = as.numeric(Std.Error)) %>%
  mutate(lwrSE = Estimate - (1.96*Std.Error)) %>%
  mutate(uprSE = Estimate + (1.96*Std.Error)) %>%
  mutate(species = factor(species, levels = c("Mimus parvulus",
                                            "Platyspiza crassirostris",
                                            "Camarhynchus parvulus",
                                            "Geospiza fortis",
                                            "Geospiza fuliginosa",
                                            "Geospiza scandens")))

pdf("output_plots/effect_size_condition_psem_2SE.pdf", width = 10)
mod_coefs %>%
  ggplot(aes(x = Estimate, y = species)) +
  geom_tile(aes(width = uprSE - lwrSE, height = 0.15),
            alpha = 0.5, fill = "#1abc9c", show.legend = F) +
  geom_point(size = 4) +
  labs(x = "Effect of lesion size on condition (β ± 95% CI)", y = "Species")  +
  theme(axis.text.y = element_text(face = "italic", size = 18),
        axis.title.x = element_text(size = 18))
  dev.off()



#
# Bootstrap to calculate 95% CIs THIS DOESN'T WORK ON A MULTIGROUP------------------------


# First, make a function that will resample data, run a psem, pull coefficients
my_stat <- function(data, indices) {
  # sample data, indices will be named from boot
  d <- data[indices, ]
  psem_mod <- psem(
    lm(condition ~ pox_size_lesion + mean_air_temp, data = d),
    lm(pox_size_lesion ~ mean_air_temp, data = d),
    data = d
  )
  # Compare results among species.
  multi_test <- multigroup(psem_mod, group = "species")
  #results <- sapply(1:6, function(x) multi_test$group.coefs[[x]][1,8], simplify = TRUE) %>% as.numeric
  if (length(multi_test$group.coefs) < 2) {
    results <- rep(NA,6)
  } else {
    results <- sapply(1:6, function(x) multi_test$group.coefs[[x]][1,8], simplify = TRUE) %>% as.numeric
  }
  return(results)

  # return(mean(d$condition))  # or extract model coefficient, etc.
}

boot_out <- boot(data = mdat, statistic = my_stat, R = 10)

my_stat <- function(data, indices) {
  d <- data[indices, ]
  mod <- multigroup(psem(
    lm(condition ~ pox_size_lesion + mean_air_temp, data = d),
    lm(pox_size_lesion ~ mean_air_temp, data = d),
    data = d
  ), group = "species") %>%
    .$group.coefs %>%
    .$FOR %>% .[1,8] %>% as.numeric
  stat <- mod
  return(stat)

}

boot_out <- boot(data = mdat, statistic = my_stat, R = 100)
boot_out
mean(mdat$condition)

psem(
  lm(condition ~ pox_size_lesion + mean_air_temp, data = mdat),
  lm(pox_size_lesion ~ mean_air_temp, data = mdat),
  data = mdat) %>% summary %>% .$coefficients %>% .[1,8]



# Manually bootstrap with a for loop and calculate CIs ------------------------------------

condition_spp_psm <- psem(
  lm(condition ~ pox_size_lesion + mean_air_temp, data = mdat),
  lm(pox_size_lesion ~ mean_air_temp, data = mdat),
  data = mdat)

multi_analysis <- multigroup(condition_spp_psm2, group = "species")


output_coefs <- multi_analysis$group.coefs %>%
  bind_rows(, .id = "species") %>%
  filter(Predictor == "pox_size_lesion") %>%
  arrange(species)

names_results[i,] <- output_coefs$species
results <- output$Std.Estimate

results <- sapply(1:6, function(x) multi_analysis$group.coefs[[x]][1,8])

# For loop to manually resample data
#boot_results <- data.frame(PAR = NA, CRA = NA, FUL = NA, SCA = NA, FOR = NA, GAMO = NA)
boot_results <- matrix(NA, nrow = 10000, ncol = 6) # 100 reps, 6 species
colnames(boot_results) <- c("CRA","FOR","FUL","GAMO", "PAR", "SCA")

names_results <- matrix(NA, nrow = 10000, ncol = 6)
# Loop to run bootstraps.
for (i in 1:10000) {
  # Sample the dataframe with replacement
  indices <- sample(1:nrow(mdat), replace = TRUE)
  d <- mdat[indices, ]

  # Create the psem model with the resampled dataset
  psem_mod <- psem(
    lm(condition ~ pox_size_lesion + mean_air_temp, data = d),
    lm(pox_size_lesion ~ mean_air_temp, data = d),
    data = d
  )
  # Compare results among species.
  multi_test <- multigroup(psem_mod, group = "species")

  # Some iterations will not find significant effect of species, so the output
  # will only include 1 data frame, rather than 6 (i.e., one per each spp.)
  # Check to see how long output is; if it's 1, return NA, if its 6, return coefs.

  if (length(multi_test$group.coefs) < 2) {
    results <- rep(NA,6)
    names <- rep(NA,6)
  } else {
    output_coefs <- multi_test$group.coefs %>%
      bind_rows(, .id = "species") %>%
      filter(Predictor == "pox_size_lesion") %>%
      arrange(species)

    names <- output_coefs$species
    results <- output_coefs$Std.Estimate

  }
  boot_results[i, ] <- as.numeric(results)
  names_results[i, ] <- names

}

# check that each column exclusively has the correct name
sapply(1:6, function(x) table(names_results[,x]))

# check that names match
names_results[1,] == colnames(boot_results)

boot_results <- as.data.frame(boot_results)
hist(boot_results$FOR)
abline(v = mean(boot_results$FOR, na.rm = T))

# calculate 95% CIs and save as data frame
CIs <- apply(boot_results, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)) %>%
  t() %>%
  as.data.frame %>%
  rename("lwrCI" = "2.5%", "uprCI" = "97.5%")
# add mean and standard error
CIs$SE <- apply(boot_results, 2, sd, na.rm = TRUE) %>% as.vector
CIs$mean <- apply(boot_results, 2, mean, na.rm = TRUE) %>% as.vector

# Add species
CIs$species <- rownames(CIs)
CIs <- mdat %>% select(species, taxon) %>%
  distinct %>%
  droplevels %>%
  left_join(CIs) %>%
  mutate(numeric_axis = as.numeric(taxon))


# Save results because bootstrapping takes too long
write.csv(CIs, "formatted_data/bootstrapped_species_psem_estimates.csv")
bootresults<- read.csv("formatted_data/bootstrapped_species_psem_estimates.csv")


# plot up the results

pdf("output_plots/effect_size_condition_bootstrappedpsem.pdf", width = 10)
CIs %>%
  ggplot(aes(x = mean, y = numeric_axis)) +
  geom_rect(aes(xmin = lwrCI, xmax = uprCI,
                ymin = numeric_axis-0.1, ymax = numeric_axis + .1),
                alpha = 0.5, fill = "#1abc9c", show.legend = F) +
  geom_point(size = 3, show.legend = FALSE) +
labs(x = "Effect of lesion size on condition (β ± 95% CI)", y = "Species")  +
  theme(axis.text.y = element_text(face = "italic", size = 18),
        axis.title.x = element_text(size = 18)) +
  scale_y_continuous(
    breaks = CIs$numeric_axis,
    labels = CIs$taxon)

dev.off()



#
# recaptures  -------------------------------------------------------------
# Does the number of captures vary among species? not really
head(captures)
firstcaps <- captures %>%
            filter(capturen == 1) %>%
            select(pox_iur, ncaps, site, taxon) %>%
  na.omit()

# how many captures per bird
mean(firstcaps$ncaps)

recap_psm <- psem(
  glmer(ncaps ~ pox_iur + (1|site), data = firstcaps, family = "poisson"),
  data = firstcaps)
glmer(ncaps ~ pox_iur + (1|site), data = firstcaps, family = "poisson") %>% summary


firstcaps2 <-
  captures %>%
  filter(capturen == 1) %>%
  select(pox_iur, ncaps, site, taxon, species, inf_stat) %>%
  filter(species %in% c("FOR", "GAMO", "PAR", "CRA", "CRA", "FUL")) %>%
  na.omit()


recap_psm <- psem(
  glm(ncaps ~ pox_iur , data = firstcaps2, family = "poisson"),
  data = firstcaps)

summary(recap_psm)
multigroup(recap_psm, group = "species", test.type = "LR")



captures %>%
  filter(species %in% c("PET","MYI")) %>%
  filter(pox_iur != "U") %>%
  group_by(species, habitat) %>% tally
