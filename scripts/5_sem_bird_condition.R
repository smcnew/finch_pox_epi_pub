
#install.packages("piecewiseSEM")
#install.packages("effects")

library(piecewiseSEM)
library(effects)
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

# investigating connection between temperature and condition:

cor(d23$condition, d23$templ2, use = "complete.obs")^2
lm(condition ~ templ2, data = d23) %>% summary

# Format data for model, remove NAs
mdat <- d23 %>%
  filter(species %in% c("FOR","FUL","PAR","SCA", "GAMO", "CRA")) %>%
  select(condition, pox_size_lesion, species, mass, condition_z, mass_z, inf_stat,
         mean_air_temp, bander_id, taxon, pox_iur, templ2) %>%
  na.omit %>%
  mutate(scaled_temp = scale(templ2)) %>%
  mutate(scaled_lesion = scale(pox_size_lesion))

table(mdat$taxon, mdat$pox_iur)

# PSM model
condition_spp_psm <- psem(
  lm(condition ~ scaled_lesion + scaled_temp, data = mdat),
  lm(scaled_lesion ~ scaled_temp, data = mdat),
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
  filter(Predictor == "scaled_lesion") %>%
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




# Plot relationship for each species  -------------------------------------

head(mdat)
mdat %>%
  filter(pox_iur == "Infected") %>%
  ggplot(aes(x = pox_size_lesion, y = condition)) +
  geom_point() +
  facet_wrap(~ taxon) +
  geom_smooth(method = "lm", se = TRUE)

head(multi_analysis)
condition_spp_psm %>% str
condition_spp_psm[2]

multi_analysis %>% str

#
