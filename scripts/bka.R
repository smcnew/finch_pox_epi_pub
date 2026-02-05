# BKA analysis
library(lme4)
library(lmerTest)
library(ggthemr) # plot params
library(ggplot2)
library(dplyr)

# Set plotting params

ggthemr("flat", layout = "clean", text_size = 14)


# Data --------------------------------------------------------------------


bka <- read.csv("input_data/2024GalapagosBKAsFinal.csv") %>%
  rename_all(tolower) %>%
  filter(species != "MAG") %>%
  filter(species != "GAMO")

bka %>%
  filter(age == "A") %>%
#  filter(pox_iur != "I") %>%
  ggplot(aes(x = species, y = killedperplas)) +
  geom_boxplot() +
  geom_point(color = "#34495e") +
  scale_x_discrete(labels = c("CRA" = "Vegetarian", "FOR" = "Med. ground", "FUL" = "Sm. ground",
                              "PAR" = "Sm. tree", "SCA" = "Cactus")) +
  labs(y = "Bacterial killing capacity", x = "")


table(bka$species)
               swatch() %>% plot
