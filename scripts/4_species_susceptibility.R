# Code for analysis question #3: Do species vary in susceptibility to pox?


library(ggthemr) # plot params
library(ggplot2) # plotting
library(patchwork) # plotting
library(sjPlot) # for saving tables

#library(nnet) # multinomial model; deprecated in revision
#library(emmeans) # posthoc comparisons; deprecated in revision

library(brms) # multinomial model for revision
library(tidybayes) # plotting results of brms
library("modelr") # plotting results of brms

library(sjPlot) # Exporting model results
library(lubridate) # wrangling dates
library(tidyr) # data wrangling
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

mult_data <- captures %>% select(pox_iur, taxon, elevation, templ2, site)

# read in table of observed prevalence for supplement
summary_species <- read.csv("output_tables/prevalence_by_species.csv")


#
# NNet Multinomial model; DEPRECATED FOR REVISION------------------------------------------


# Make a multinomial model of passerine species (i.e, those that get pox )
mod3 <- multinom(pox_iur ~ taxon + elevation + templ2 , data = mult_data)
emm3 <- emmeans(mod3, ~ taxon | pox_iur, mode = "prob")
summary(mod3)

# Evaluate the species model compared to the null
null_model <- multinom(pox_iur ~ elevation + templ2, data = mult_data)

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




#
# revision, mult model in brms -------------------------------------------------
# for revision: more complex multinomial model:
fit <- brm(
  formula = pox_iur ~ taxon + elevation + templ2 + (1 | site),
  data = mult_data,
  family = categorical(link = "logit"),
  chains = 4,
  cores = 3,
  iter = 4000,
  warmup = 1000,
  seed = 123,
  save_pars = save_pars(all = TRUE)
)
summary(fit)
tab_model(
  fit,
  transform = "exp",
  show.ci = 0.95,
  show.r2 = FALSE,
  prefix.labels = "none",
  dv.labels = c(
    "Symptomatic",
    "Recovered"
  )
)

plot(fit)

fit_nospp <- brm(
  formula = pox_iur ~ elevation + templ2 + (1 | site),
  data = mult_data,
  family = categorical(link = "logit"),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 1000,
  seed = 123,
  save_pars = save_pars(all = TRUE)
)

# compare the fits

loo_compare(
  loo(fit),
  loo(fit_nospp)
)

# check observations
loo_fit <- loo(fit)

# Two records: an infected C. olivacea and a recovered M. magnirostris
# are kind of anomolous. But we have no reason to exclude them so there they are.
# which(loo_fit$diagnostics$pareto_k > 0.7)
# mult_data[c(1251, 2134), ]



# Figure of Bayesian results  ---------------------------------------------


# Make a figure
newdat <- expand.grid(
  taxon = levels(mult_data$taxon),
  elevation = mean(mult_data$elevation, na.rm = TRUE),
  templ2 = mean(mult_data$templ2, na.rm = TRUE)
)

pred <- newdat %>%
  add_epred_draws(
    fit,
    re_formula = NA) %>%
  filter(.category != "Uninfected")

pdf("output_plots/I_R_plot_revision_brms.pdf", width = 10)
pred %>%
  ggplot(aes(x = .epred, y = taxon, color = .category)) +
  stat_pointinterval(.width = c(.66, .95), position = position_dodge(width = 0.4)) +
  labs(x = "Estimated proportion", y = "Species", color = "Pox Status") +
  theme(axis.text.y = element_text(face = "italic")) +
  scale_color_manual(
    values = c(
      "Recovered" = "#f39c12",
      "Infected" = "#e74c3c"
    ),
    labels = c("Recovered" = "Recovered", "Infected" = "Symptomatic")
  ) +
  scale_y_discrete(
    labels = c(
      "Camarhynchus psittacula" = "Camarhynchus psittacula*",
      "Geospiza magnirostris" = "Geospiza magnirostris*"
    )
  ) +
scale_x_continuous(
    limits = c(0, 0.6),
    breaks = seq(0, 0.6, by = 0.1)
  )
dev.off()




# Format data for supplemental table --------------------------------------------------

pred_summary <- pred %>%
  median_qi(.epred, .width = 0.95) %>%
  select(-elevation, -templ2, -.row, -.width, -.point, -.interval) %>%
  pivot_wider(
    names_from = .category,
    values_from = c(.epred, .lower, .upper)
  ) %>%
  mutate(across(where(is.numeric),
                ~ round(.x, 3))) %>%
  mutate(Posterior_symptomatic = paste0(.epred_Infected, " [", .lower_Infected, "—", .upper_Infected, "]"),
         Posterior_recovered = paste0(.epred_Recovered, " [", .lower_Recovered, "—", .upper_Recovered, "]"))


summary_species

summary_species <- select(pred_summary, taxon, Posterior_symptomatic, Posterior_recovered) %>%
  left_join(summary_species, .) %>%
  mutate(prevalence = round(prevalence, 3)) %>%
  rename("Latin_name" = taxon,
         "Observed_prevalence" = prevalence,
         "Captures" = total,
         "Symptomatic" = I,
         "Recovered" = R,
         "Unaffected" = U)

# write output
write.csv(summary_species, "output_tables/species_prevalence_with_models.csv", row.names = FALSE)
