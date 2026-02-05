# Code for analysis question #3: Do species vary in susceptibility to pox?


library(lubridate)
library(ggthemr) # plot params
library(ggplot2) # visualization
library(patchwork) # plotting
library(sjPlot) # for saving tables

library(nnet) # multinomial model
library(emmeans) # posthoc comparisons
library(tidyr)
library(dplyr)


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


# Multinomial model -------------------------------------------------------
mult_data <- captures %>% select(pox_iur, taxon)

# Make a multinomial model of passerine species (i.e, those that get pox )
mod3 <- multinom(pox_iur ~ taxon, data = mult_data)
emm3 <- emmeans(mod3, ~ taxon | pox_iur, mode = "prob")
summary(mod3)

# Evaluate the species model compared to the null
null_model <- multinom(pox_iur ~ 1, data = mult_data)

# Compare with full model
anova(null_model, mod3, test = "Chisq")


# Figure out sample sizes for each species
table(captures$species)

# Plot of species prevalence
# plot(emm3, type = "response") # simple plot
# CI's are symmetric around the estimate so we can use geom_tile instead of geom_rect

pdf("output_plots/I_R_plot.pdf", width = 10)
emm3 %>%
  as.data.frame() %>%
  filter(pox_iur != "Uninfected") %>%
  ggplot(aes(x = prob, y = taxon, color = pox_iur)) +
  geom_tile(aes(width = upper.CL - lower.CL, height = 0.3, fill = pox_iur),
            position = position_dodge(width = 0.4), alpha = 0.5, show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.4), size = 2)+
  labs(x = "Estimated proportion", y = "Species", color = "Pox Status") +
  theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = c("Recovered" = "#f39c12", "Infected" = "#e74c3c" )) +
  scale_color_manual(values = c("Recovered" = "#f39c12", "Infected" = "#e74c3c")) +
  scale_y_discrete(labels = c("Camarhynchus psittacula" = "Camarhynchus psittacula*", "Geospiza magnirostris" = "Geospiza magnirostris*"))
dev.off()



