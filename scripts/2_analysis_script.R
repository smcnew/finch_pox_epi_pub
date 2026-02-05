# Analyze some data from finches
# This script is analyzing trends over the study period for the birds
# e.g. how many birds, species and abundance over time, mass over time
# packages and themes  ----------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggthemr)
library(khroma)

# set themes

ggthemr("flat", layout = "clean", text_size = 20)

# add more colors from Paul tol colors
muted <- color("muted")
muted(9)
set_swatch(c(swatch(), muted(9)))
plot(swatch()) #check em out.

# load data ---------------------------------------------------------------
# Just keep adults
captures <- read.csv("formatted_data/combined_2022_2023_finch.csv") %>%
            filter(sex != "N") #take out nestlings

d22 <- captures %>% filter(year == 2022)
d23 <- captures %>% filter(year == 2023 | year == 2024)
# a vector of all the months we worked
month_year_vector <- captures %>% select(month, year) %>%
  unique %>% mutate (monthyear = paste (month, year, sep = "-")) %>% pull (monthyear)

# 1. summary statistics ---------------------------------------------------

dim(captures) # 5596 individuals
dim(d22) # 3025 in 2022
dim(d23) # 2471 in 2023-2024

# how many individuals?
captures %>% select(band) %>% unique %>% dim # 4,741 individual birds


# 2. population trends ----------------------------------------------------

# First look at 2022:
# 14 species, ranging from 3 captures to 1494
filter(captures, year == 2022) %>% group_by(species) %>% tally
filter(captures, year == 2022) %>% group_by(type_of_site) %>% tally # all in highlands

# Make a plot by month for 2022
# Create a data frame with captures per month for each species
caps_per_m_22 <- captures %>%
  filter(year == 2022) %>%
  filter(species %in% c("CRA", "FOR", "FUL", "MYI", "OLI", "GAMO",
                        "PAL", "PAR", "PET", "PSI")) %>%
  group_by(species, year_month) %>%
  tally() %>%
  mutate(year_month = as.numeric(as.character(year_month))) %>%
  as.data.frame

# Create a total number of captures
all_caps_22 <- table(captures$year_month) %>%
  as.data.frame %>%
  rename(year_month = Var1, total = Freq) %>%
  mutate(year_month = as.numeric(as.character(year_month))) %>%
  left_join(caps_per_m_22, .) %>%
  mutate(prop_caps = n/total)

# Create a plot where y is "proportion of captures"
# x is "month" and each line is a species
pdf("output_plots/capturas_species_month_22.pdf", width = 10)
all_caps_22 %>%
  ggplot(aes(x = year_month, y = prop_caps, color = species)) +
  geom_line(linewidth = 1.3) +
  labs(x = "Month", y = "Proportion of captures", color = "") +
  scale_x_continuous(breaks = sort(unique(all_caps_22$year_month)),
                     labels = c(3:10))
dev.off()


# 2023
# how many species were captured each month?
filter(captures, year != 2022) %>% group_by(species) %>% tally
filter(captures, year != 2022) %>% group_by(type_of_site) %>% tally # most in lowlands
# Make a plot:
# First, aggregate captures by species for each month
caps_per_m_23 <- captures %>%
  filter(year != 2022) %>%
  filter(species %in% c("CRA", "FOR", "FUL", "GAMO", "MAG", "PAR", "SCA")) %>%
  group_by(species, year_month) %>%
  tally() %>%
  mutate(year_month = as.numeric(as.character(year_month))) %>%
  as.data.frame

# Add number of captures per month
all_caps_23 <- table(captures$year_month) %>%
  as.data.frame %>%
  rename(year_month = Var1, total = Freq) %>%
  mutate(year_month = as.numeric(as.character(year_month))) %>%
  left_join(caps_per_m_23, .) %>%
  mutate(prop_caps = n/total)


pdf("output_plots/capturas_species_month_23.pdf", width = 10)
all_caps_23 %>%
  ggplot(aes(x = year_month, y = prop_caps, color = species)) +
  geom_line(linewidth = 1.3) +
  labs(x = "Month", y = "Proportion of captures", color = "") +
  scale_x_continuous(breaks = sort(unique(all_caps_23$year_month)),
                     labels = c(2:12,1:3))
dev.off()


# How many birds per month?
capture_summary <-
captures %>%
  group_by(year_month, type_of_site) %>%
  summarize(total_caps = n(),
            net_hours = sum(net_hours),
            month = month) %>%
  mutate(caps_per_h = total_caps / net_hours)

capture_summary %>%
  ggplot(aes(x = year_month, y = caps_per_h, color= type_of_site)) +
          geom_line(lwd = 2) +
          scale_x_continuous(breaks = sort(unique(capture_summary$year_month)),
          labels = month_year_vector) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.3))


# Mass effects over the year  ---------------------------------------------

head(captures)
# do a test case with just ful
captures %>%
  filter(species == "FOR") %>%
  filter(year != "2022") %>%
  ggplot(aes(x = month, y = mass, group = month)) + geom_boxplot()

captures %>%
  filter(species == "FUL") %>%
  filter(year == "2022") %>%
  ggplot(aes(x = month, y = mass, group = month)) + geom_boxplot()

captures %>%
  filter(species == "FUL") %>%
  filter(year == "2022") %>%
  ggplot(aes(x = day_of_year, y = mass)) + geom_point() + geom_smooth()


lm(mass ~ day_of_year + as.factor(year), data = captures) %>% summary()


# gamo males --------------------------------------------------------------

head(captures)

captures %>% filter(species == "GAMO") %>%
  filter(sex != "J") %>%
  filter(sex !="U") %>%
  ggplot(aes(x = wing, y = beak_length, color = sex)) + geom_point()
