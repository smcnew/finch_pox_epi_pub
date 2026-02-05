# Script to investigate effects of pox on condition of birds
# Started April 2024
# Sabrina McNew

# Packages ----------------------------------------------------------------
# devtools::install_github('Mikata-Project/ggthemr')

library(lme4) # mixed models
library(lmerTest)
library(mgcv) # gams
library(lubridate)

library(ggthemr) # plot params
library(ggplot2) # visualization
library(visreg) # model visualization
library(patchwork) # plotting
library(sjPlot) # for saving tables

# pacakges for infection comparison between species
library(nnet) # multinomial model
library(emmeans) # posthoc comparisons
# library(multcomp)
# library(multcompView)
library(tidyr)
library(dplyr)



# Set plotting params

ggthemr("flat", layout = "clear", text_size = 20) # clean is no grid lines
# plot(swatch()) #check em out.

# Load data  --------------------------------------------------------------

captures <- read.csv("formatted_data/combined_2022_2023_finch_climate.csv") %>%
  mutate(date = as.Date(date)) %>%
  mutate(year_week = as.Date(year_week)) %>%
  mutate(year_month = as.Date(year_month)) %>%
  filter(sex != "N") %>% # get rid of nestlings
  #filter(banded == 1) %>% # keep only banded birds
  #filter(species != "ZEN") %>% # remove non-passerines
  #filter(species != "MELA") %>%
  filter(site != "Garrapatero") %>% # remove Garrapatero, only 13 captures
  mutate(pox_iur = factor(pox_iur, levels = c("U", "I", "R"),
                          labels = c("Uninfected", "Infected", "Recovered"))) %>% # relevel pox scale
  mutate(pox_scale = case_when(is.na(pox_scale) ~ "U", TRUE ~ pox_scale)) %>%
  mutate(pox_scale = factor(pox_scale, levels = c("U","A","B","C","D"))) %>% # order lesion scale in order of severity
  mutate(pox_size_lesion = case_when(is.na(pox_size_lesion) ~ 0, TRUE ~ pox_size_lesion)) %>% # replace NA with 0 for lesion size
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


# mean_temp is for the month; mean_air_temp is daily
# create a z scaled variable of condition and mass
# per species so that effect sizes are comparable across species

# captures <- captures %>%
#   group_by(species) %>%
#   mutate(condition_z = scale(condition)) %>%
#   mutate(mass_z = scale(mass)) %>%
#   ungroup()

# Split data set in to highlands (2022) and lowlands (2023-2024)
# d22 <- captures %>%
#   filter(year == 2022) %>%
#   filter(site != "El Barranco") %>% # filter out 8 test birds from lowlands
#   mutate(type_of_site = relevel(factor(type_of_site), ref = "bosque")) %>%
#   group_by(species) %>% # exclude species with fewer than 10 records, (G. scandens)
#   filter(n() >= 10) %>%
#   ungroup() %>%
#   droplevels()
#
# d23 <- captures %>%
#   filter(year == 2023 | year == 2024) %>%
#   filter(zone == "Puerto Ayora" ) %>%
#   droplevels()


# Summary stats -----------------------------------------------------------
# How many sites in 2022?
captures %>% filter (year == 2022) %>% dplyr::select(site, habitat) %>% unique %>% group_by(habitat) %>% tally

# How many captures in each year
captures %>% group_by(year) %>% tally
dim(captures)

# How many species each year
captures %>% group_by(year) %>% summarize(n_spp = length(unique(species)))
length(unique(captures$species))


# How many visits to the sites each year
df <- captures %>% select(date, year, site) %>% unique
table(df$year, df$site)

# Format site info for supplement
# Number of visits to each site per year
visits_per_year <-
  captures %>%
  select(site, pox_iur, year, date) %>%
  group_by(site, year) %>%
  summarize(visits = length(unique(date))) %>%
  ungroup() %>%
  pivot_wider(names_from = year, values_from = visits, values_fill = 0, names_prefix = "sampling_days_")

# Overall pox prevalence at each site and the number of captures
prevalence <-
  captures %>%
  group_by(site) %>%
  summarize(Prevalence = round(mean(pox_iur == "Infected"), 3),
            Total_captures = n())


#Lat/long and elevation per site
site_info <-
  captures %>%
  select(site, habitat, latitude, longitude, elevation) %>% distinct(site,, .keep_all = TRUE)

site_info <- left_join(site_info, visits_per_year) %>%
              left_join(., prevalence)

# Write out table for supplement
write.csv(site_info, "output_tables/Table_S1.csv", row.names = F)

# Prevalence by year
captures %>% group_by(year) %>% summarize(prevalene = mean(pox_iur == "Infected"))

# Prevalence by species

# Calculate prevalence for all species across years
prev_all <- captures %>%
  group_by(taxon) %>%
  summarize(prevalence = mean(pox_iur == "Infected"),
            total = n(),
            I = sum(pox_iur == "Infected"),
            U = sum (pox_iur == "Uninfected"),
            R = sum (pox_iur == "Recovered"))

# Write out prevalence for the supplement
write.csv(prev_all, "output_tables/prevalence_by_species.csv", row.names = F)


# correlation between pox scale and lesion size
aov(pox_size_lesion ~ pox_scale, data = captures) %>% summary


# Environmental patterns over the course of the year (PLOTS)  -----------------------------------

# make a data frame with just one record per day
single_day_d22 <- captures %>% filter(year == 2022) %>% distinct(date, .keep_all = TRUE)

# Plots of variation in covariates during 2022
dens22 <- single_day_d22 %>% ggplot(aes(x = date, y = bird_dens)) +
  geom_point(color = "#2ecc71") +
  geom_smooth() +
  labs(x = "Date (2022)", y = "Host density (captures/net hour)")

temp22 <- single_day_d22 %>% ggplot(aes(x = date, y = mean_air_temp)) +
  geom_point(color = "#2ecc71") +
  geom_smooth() +
  labs(x = "Date (2022)", y = "Mean daily air temp (C)")

precip22 <- single_day_d22 %>% ggplot(aes(x = year_week, y = week_precip)) +
  geom_point(color = "#2ecc71") +
  geom_smooth() +
  labs(x = "Month (2022)", y = "Weekly total precipitation (mm)")


humid22 <- single_day_d22 %>% ggplot(aes(x = date, y = humidity)) +
  geom_point(color = "#2ecc71") +
  geom_smooth() +
  labs(x = "Month (2022)", y = "Daily humidity (%)")



# Plots of covariates during 2023-2024
single_day_d23 <- captures %>% filter(year != 2022) %>% distinct(date, .keep_all = TRUE)

dens23 <- single_day_d23 %>% ggplot(aes(x = date, y = bird_dens)) +
  geom_point() +
  geom_smooth(fill = "#3498db", color = "#3498db") +
  labs(x = "Date (2023-2024)", y = "Host density (captures/net hour)")

temp23 <- single_day_d23 %>% ggplot(aes(x = date, y = mean_air_temp)) +
  geom_point() +
  geom_smooth(fill = "#3498db", color = "#3498db") +
  labs(x = "Date (2023-2024)", y = "Mean daily air temp (C)")

precip23 <-
  single_day_d23 %>% ggplot(aes(x = year_week, y = week_precip)) +
  geom_point() +
  geom_smooth(color = "#3498db", fill = "#3498db") +
  labs(x = "Month (2022)", y = "Weekly total precipitation (mm)")

humid23 <- single_day_d23 %>% ggplot(aes(x = date, y = humidity)) +
  geom_point(color = "#3498db") +
  geom_smooth(color = "#3498db", fill = "#3498db") +
  labs(x = "Month (2023)", y = "Daily humidity (%)")


pdf("output_plots/covariates_annual.pdf", width = 14, height = 20)
(temp22 + temp23)/(humid22 + humid23)/(precip22 + precip23)/(dens22 + dens23)+ plot_annotation(tag_levels = 'A')
dev.off()



 #
# Patterns over time raw data  DEPRECATED -----------------------------------------

# Make a figure showing I and R prevalence over months
# prev_by_date <- captures %>%
#   group_by(year_month) %>%
#   summarize(Infected = mean(pox_iur == "Infected"),
#             Recovered = mean(pox_iur == "Recovered")) %>%
#   pivot_longer(cols = c(Infected, Recovered), values_to = "Proportion",
#                names_to = "Category")
#
# prev_22 <- prev_by_date %>%
#   filter(year(year_month) == 2022) %>%
#   ggplot(aes(x = year_month, y = Proportion, color = Category, fill = Category)) +
#   geom_point() + geom_smooth() +
#
#   labs(x = "2023")
#
#
# prev_23 <- prev_by_date %>%
#   filter(year(year_month) != 2022) %>%
#   ggplot(aes(x = year_month, y = Proportion, color = Category, fill = Category)) +
#   geom_point() + geom_smooth() +
#   labs(x = "2023 - 2024")
#
# pdf("output_plots/prevalence_by_month_raw_data.pdf", width = 10)
# prev_22/prev_23
# dev.off()

# Make a figure showing variation in prevalence by week
prev_by_date <- captures %>%
  group_by(year_week) %>%
  summarize(
    Infected = mean(pox_iur == "Infected"),
    Recovered = mean(pox_iur == "Recovered")
  )  %>%
  pivot_longer(
    cols = c(Infected, Recovered),
    values_to = "Proportion",
    names_to = "Category"
  )

# Add columns by hand because date functions are difficult
prev_by_date$year <- year(prev_by_date$year_week)
prev_by_date$week <- week(prev_by_date$year_week) # week of year

# Add separate column to make 2024 years continue in sequence from 2023
prev_by_date <- prev_by_date %>%
  mutate(plot_week = case_when(year == 2024 ~ week + 52,
                               year !=2024 ~ week))

# Make a plot for 2022

prev_22 <- prev_by_date %>%
            filter(year == 2022)

axis_breaksa <- prev_22 %>% distinct(year_week, .keep_all = T) %>% pull(plot_week)
axis_labelsa <- prev_22 %>% distinct(year_week, .keep_all = T) %>% pull(year_week) %>% format(., "%d-%b")

prev_22p <- prev_22 %>%
  ggplot(aes(
    x = week,
    y = Proportion,
    color = Category,
    fill = Category
  )) +
  geom_point() + geom_smooth(method = "gam") +
  labs(x = "2022") +
  scale_x_continuous(breaks = axis_breaksa[seq(1, length(axis_breaksa), by = 5)], labels = axis_labelsa[seq(1, length(axis_labelsa), by = 5)]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Make a plot for 2023 - 2024
# First filter the data frame to the right years
prev_23 <- prev_by_date %>%
  filter(year != 2022)

# Format axis label vectors
axis_breaks <- prev_23 %>% distinct(year_week, .keep_all = T) %>% pull(plot_week)
axis_labels <- prev_23 %>% distinct(year_week, .keep_all = T) %>% pull(year_week) %>% format(., "%d-%b")

# make a plot
prev_23p <- prev_23 %>%
  ggplot(aes(
    x = plot_week,
    y = Proportion,
    color = Category,
    fill = Category
  )) +
  geom_point() + geom_smooth(method = "gam") +
  scale_x_continuous(breaks = axis_breaks[seq(1, length(axis_breaks), by = 5)], labels = axis_labels[seq(1, length(axis_labels), by = 5)]) +
  labs(x = "2023 - 2024") +
  theme(  axis.text.x = element_text(angle = 45, hjust = 1)
  )

# export the plots
pdf("output_plots/prevalence_by_week_raw_data.pdf", width = 10)
prev_22p / prev_23p
dev.off()


#
# Analysis 1: Prevalence over weeks (gam) ------------------------------------------------
# Weekly GAM

# format data by week
prev_by_week <- captures %>%
  select(year_week, year, pox_iur) %>%
  group_by(year_week, year) %>%
  summarize(inf = sum(pox_iur == "Infected"),
            rec = sum(pox_iur == "Recovered"),
            uninf = (sum(pox_iur == "Recovered") + sum(pox_iur == "Uninfected")),
            unrec = (sum(pox_iur == "Infected") + sum(pox_iur == "Uninfected")),
            n = n()) %>%
  mutate(num_week = as.numeric(year_week)) %>%
  mutate(prev_i = inf/n) %>%
  mutate(prev_r = rec/n)

# Make objects for 2022 and 2023 separately so visreg works
pbw22 <- prev_by_week %>% filter(year == 2022)
pbw23 <- prev_by_week %>% filter(year != 2022)

# Make some gams
# 2022 model infected
gam22i <-   gam(
  cbind(inf, uninf) ~ s(num_week),
  data = pbw22,
  family = binomial(link = "logit")
)

# 2023-2024 model infected
gam23i <-   gam(
  cbind(inf, uninf) ~ s(num_week),
  data = pbw23,
  family = binomial(link = "logit")
)

# 2022 model recovered
gam22r <-   gam(
  cbind(rec, unrec) ~ s(num_week),
  data = pbw22,
  family = binomial(link = "logit")
)


# 2023-2024 model recovered
gam23r <-   gam(
  cbind(rec, unrec) ~ s(num_week),
  data = pbw23,
  family = binomial(link = "logit")
)


summary(gam22i)
summary(gam23i)
summary(gam22r)
summary(gam23r)
# Make plots
# Get predicted values from visreg and add x axis date column
pred_obj22i <- visreg(gam22i, scale = "response") %>% .$fit %>%
  mutate(date = as.Date(num_week, origin = "1970-01-01")) %>%
  mutate(Status = "Infected")
pred_obj22r <- visreg(gam22r, scale = "response") %>% .$fit %>%
  mutate(date = as.Date(num_week, origin = "1970-01-01")) %>%
  mutate(Status = "Recovered")
pred_obj23i <- visreg(gam23i, scale = "response") %>% .$fit %>%
  mutate(date = as.Date(num_week, origin = "1970-01-01")) %>%
  mutate(Status = "Infected")
pred_obj23r <- visreg(gam23r, scale = "response") %>% .$fit %>%
  mutate(date = as.Date(num_week, origin = "1970-01-01"))%>%
  mutate(Status = "Recovered")



# Plot for 2022
prev_22p <- ggplot() +
  geom_ribbon(data = pred_obj22i,
              aes(x = date, ymin = visregLwr, ymax = visregUpr, fill = Status),
              alpha = 0.4) +
  geom_line(data = pred_obj22i, aes(x = date, y = visregFit, color = Status), linewidth = 1.2) +
  geom_ribbon(data = pred_obj22r,
              aes(x = date, ymin = visregLwr, ymax = visregUpr, fill = Status),
              alpha = 0.4) +
  geom_line(data = pred_obj22r, aes(x = date, y = visregFit, color = Status), linewidth = 1.2) +
  geom_point(data = pbw22, aes(x = year_week, y = prev_i), color = "#e74c3c") +
  geom_point(data = pbw22, aes(x = year_week, y = prev_r), color = "#f39c12") +
  scale_x_date(
    date_labels = "%d-%b",
    date_breaks = "4 week"
  ) +
  labs(x = "2022", y = "Proportion of captures") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Infected" = "#e74c3c", "Recovered" = "#f39c12")) +
  scale_color_manual(values = c("Infected" = "#e74c3c", "Recovered" = "#f39c12"))

prev_22p
# Plot for 2023 - 2024
prev_23p <- ggplot() +
  geom_ribbon(data = pred_obj23i,
              aes(x = date, ymin = visregLwr, ymax = visregUpr, fill = Status),
              alpha = 0.4) +
  geom_line(data = pred_obj23i, aes(x = date, y = visregFit, color = Status), linewidth = 1.2) +
  geom_ribbon(data = pred_obj23r,
              aes(x = date, ymin = visregLwr, ymax = visregUpr, fill = Status),
              alpha = 0.4) +
  geom_line(data = pred_obj23r, aes(x = date, y = visregFit, color = Status), linewidth = 1.2) +
  geom_point(data = pbw23, aes(x = year_week, y = prev_i), color = "#e74c3c") +
  geom_point(data = pbw23, aes(x = year_week, y = prev_r), color = "#f39c12") +
  scale_x_date(
    date_labels = "%d-%b",
    date_breaks = "4 week"
  ) +
  labs(x = "2023 – 2024", y = "Proportion of captures") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Infected" = "#e74c3c", "Recovered" = "#f39c12")) +
  scale_color_manual(values = c("Infected" = "#e74c3c", "Recovered" = "#f39c12"))


pdf("output_plots/prevalence_by_week_modeled_data.pdf", width = 12, height = 12)
prev_22p / prev_23p
dev.off()


#
# Daily variation in prevalence (gam) -------------------------------------

# Daily GAM
# Add a column that creates date as numeric
d22$num_date <- as.numeric(d22$date)
d23$num_date <- as.numeric(d23$date)

# Run model for 2022
gam22d <- gam(inf_stat ~ s(num_date),
             data = d22,
             family = binomial(link = "logit"))

# Save model summary for supplement
tab_model(gam22d, file = "output_tables/daily_prev_gam_22.html")


# Get predicted values using visreg
pred_obj22d <- visreg(gam22w, scale="response")
pred_obj22d$fit  <- pred_obj22d$fit %>% mutate(date = as.Date(num_date, origin = "1970-01-01"))


# Plot the model predictions
p_daily_gam_22d <-
pred_obj22d$fit %>%
ggplot(aes(x = date, y = visregFit)) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.4) +
  geom_line(linewidth = 1.2) +
  labs(x = "Date (2022)", y = "Infection prevalence") +
  scale_x_date(
    date_labels = "%d-%b",
    date_breaks = "4 week"
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Run model for 2023
gam23d <- gam(inf_stat ~ s(num_date),
             data = d23,
             family = binomial(link = "logit"))
tab_model(gam23d, file = "output_tables/daily_prev_gam_23.html")

# Get model predicted values
pred_obj23d <- visreg(gam23d, scale="response")
pred_obj23d$fit  <- pred_obj23d$fit %>% mutate(date = as.Date(num_date, origin = "1970-01-01"))

# Plot model
p_daily_gam_23d <- pred_obj23d$fit %>%
  ggplot(aes(x = date, y = visregFit)) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.4) +
  geom_line(linewidth = 1.2) +
  labs(x = "Date (2023-2024)", y = "Infection prevalence") +
  scale_x_date(
    date_labels = "%d-%b",
    date_breaks = "4 week"
  ) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plots
pdf("output_plots/daily_inf_prev_gam.pdf", height = 10, width = 8)
p_daily_gam_22d/p_daily_gam_23d
dev.off()


# Environmental predictors of prevalence DEPRECATED----------------------------------------------
# New approach- do this analysis with a PSEM (see environment_sems.R)

### good model; total precip and mean temp for the month; no random effect of site_month
### because causes convergence error, but estimates are the same either way.
mod1 <- glmer(inf_stat ~  scale(total_precip) + scale(mean_temp) + habitat +
        (1|species) + (1|site), data = d22, family = "binomial")

summary(mod1)
visreg(mod1,  "mean_temp", by = "habitat")

glmer(inf_stat ~  scale(total_precip) + scale(mean_temp)  +
        (1|species) , data = d23, family = "binomial") %>% summary()

glmer(inf_stat ~  scale(total_precip) + scale(mean_temp)  + habitat + bird_dens +
        (1|species)  + (1|year), data = captures, family = "binomial") %>% summary()

mod2 <- glmer(inf_stat ~  scale(total_precip) + scale(mean_temp)  + habitat +
        (1|species)  + (1|year), data = captures, family = "binomial")
visreg(mod2,  "habitat", scale="response")

d22 %>%
  group_by(habitat) %>%
  summarize(prevalence = mean(pox_iur == "Infected"),
            total = n())

boxplot(bird_dens ~ habitat, data = captures)
dens_model <- lm(bird_dens ~ habitat + day_of_year, data = captures)
visreg(dens_model, "habitat")
captures %>%
  ggplot(aes(x = date, y = bird_dens, color = as.factor(habitat))) +
  geom_point() +
  geom_smooth()



captures %>% filter(is.na(day_of_year))
plot(bird_dens ~ day_of_year, data = captures)


captures %>% filter(is.na(habitat))


#
# Analysis 3: Species differences -----------------------------------------------------



# Make a multinomial model
# Remove non passerine species (MELA and ZEN) that had no pox
mod3 <-
multinom(pox_iur ~ taxon, data = captures)
emm3 <- emmeans(mod3, ~ taxon | pox_iur, mode = "prob")
summary(mod3)

# Evaluate the species model compared to the null
null_model <- multinom(pox_iur ~ 1, data = captures)

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

# If one wants to compare species for significant differences in infection status.
# Update: species x species comparisons are kind of dumb because we have 12 different
# species.
# two x standard error
# twose <- emm3 %>% as.data.frame %>% pull(SE) %>% mean*2
# cld(emm3, delta = twose, adjust = "none")
#
# # Look for significant differences
#
# cld(emm3)

#
# Correlates of bird condition DEPRECATED IN FAVOR OF SEM-------------------------------------
prev_all %>% arrange(prevalence)

# Mass varies slightly among habitat type, but condition does not
lmer(mass ~ habitat + scale(mean_temp) +  (1|site) + (1|species), data = captures) %>% summary
lmer(condition ~ habitat + (1|site) + (1|species), data = captures) %>% summary()


# Focus just on variation in one site with high prevalence ()
# Both mass and condition associated with air temperature
lm(mass ~ scale(mean_air_temp), data = d23) %>% summary()
lm(condition ~ scale(mean_air_temp), data = d23) %>% summary()


# Write a function to pull beta coefficients and 95% CI from all models
con_mod <- function(spp, df) {
  model <- df %>% filter(species == spp) %>%
            lm(condition_z ~ scale(pox_size_lesion) + scale(mean_air_temp), data = .)
  beta <- model$coefficients[2]
  ci_lower <- confint(model)[2,1]
  ci_upper <- confint(model)[2,2]
  format_results <- data.frame(species = spp, beta = beta, ci_lower = ci_lower, ci_upper = ci_upper)
  return(format_results[1,])
}

results_condition <- data.frame()
for (i in c("GAMO", "SCA", "CRA", "FOR", "FUL", "PAR")) {
  tmp <- con_mod(i, df = "captures")
  results_condition <- rbind(results_condition, tmp)
}

# Add Latin name to the data frame
results_condition <- d23 %>%
  select(species, taxon) %>%
  distinct() %>%
  right_join(., results_condition)

pdf("output_plots/effect_size_condition2.pdf", width = 10)
results_condition %>%
  ggplot(aes(x = beta, y = taxon)) +
  geom_tile(aes(width = ci_upper - ci_lower, height = 0.15),
            alpha = 0.5, show.legend = FALSE, fill ="#1abc9c" ) +
  geom_point(size = 3, show.legend = FALSE) +
  labs(x = "Effect of lesion size on condition (β ± 95% CI)", y = "Species")  +
  theme(axis.text.y = element_text(face = "italic", size = 18),
        axis.title.x = element_text(size = 18))
dev.off()


# make a more complex model
captures %>% filter(species == "PAR") %>%
lmer(mass_z ~ scale(pox_size_lesion) + scale(mean_air_temp) + (1|site) + (1|year), data = .)

# %>%

con_mod2 <- function(spp, df) {
  model <- df %>% filter(species == spp) %>%
    lmer(condition_z ~ scale(pox_size_lesion) + scale(mean_air_temp) + (1|site) + (1|bander_id), data = .)
  beta <- model@beta[2]
  ci_lower <- confint(model)[5,1]
  ci_upper <- confint(model)[5,2]
  format_results <- data.frame(species = spp, beta = beta, ci_lower = ci_lower, ci_upper = ci_upper)
  return(format_results[1,])
}

results_condition <- data.frame()
for (i in c("FOR", "FUL", "PAR", "CRA", "SCA","GAMO", "PAL", "PSI")) {
  tmp <- con_mod2(i, captures)
  results_condition <- rbind(results_condition, tmp)
}


results_condition <- captures %>%
  select(species, taxon) %>%
  distinct() %>%
  right_join(., results_condition)

pdf("output_plots/effect_size_z_condition_all_spp.pdf", width = 10)
results_condition %>%
  ggplot(aes(x = beta, y = taxon)) +
  geom_tile(aes(width = ci_upper - ci_lower, height = 0.15),
            alpha = 0.5, show.legend = FALSE, fill ="#1abc9c" ) +
  geom_point(size = 3, show.legend = FALSE) +
  labs(x = "Effect of lesion size on mass (β ± 95% CI)", y = "Species")  +
  theme(axis.text.y = element_text(face = "italic", size = 18),
        axis.title.x = element_text(size = 18))
dev.off()


# Mass ~ pox lesion size (DEPRECATED) --------------------------------------------
aggregate(elevation ~ habitat, captures, mean)

# with condition
d23 %>% filter(species == "GAMO") %>%
  lmer(mass ~ scale(pox_size_lesion) + scale(mean_air_temp) + (1|bander_id) , data = .) %>% summary()

d23 %>% filter(species == "CRA") %>%
  lmer(mass ~ scale(pox_size_lesion) + scale(mean_air_temp) + (1|bander_id) , data = .) %>% summary()

d23 %>% filter(species == "SCA") %>%
  lmer(mass ~ scale(pox_size_lesion) + scale(mean_air_temp) + (1|bander_id) , data = .) %>% summary()

d23 %>% filter(species == "PAR") %>%
  lmer(mass ~ scale(pox_size_lesion) + scale(mean_air_temp) + (1|bander_id) , data = .) %>% summary()

d23 %>% filter(species == "FUL") %>%
  lmer(mass ~ scale(pox_size_lesion) + scale(mean_air_temp) + (1|bander_id) , data = .) %>% summary()

d23 %>% filter(species == "FOR") %>%
  lmer(mass ~ scale(pox_size_lesion) + scale(mean_air_temp) + (1|bander_id) , data = .) %>% summary()



#
