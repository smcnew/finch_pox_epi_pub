# Script compare immune varation among species
#
library(lme4)
library(lmerTest)
library(ggthemr) # plot params
library(ggplot2)
library(patchwork)
library(dplyr)

# Set plotting params

ggthemr("flat", layout = "clean", text_size = 12)

# load data  --------------------------------------------------------------

leukocytes <- read.csv("input_data/leukocytes_diana.csv", fill = 0) %>%
          mutate(date = as.Date(date, "%d-%b-%y")) %>%
          filter(band != "SMA1361") %>%
          filter(band != "JP5867") %>%
          filter(exclude != 1)

head(captures)

captures <- read.csv("formatted_data/combined_2022_2023_finch.csv") %>%
  mutate(date = as.Date(date)) %>%
  mutate(pox_iur = factor(pox_iur, levels = c("U", "I", "R"))) %>%
  select(species, band, date, pox_iur, condition, mass, inf_stat)

leuks <- left_join(leukocytes, captures) %>%
  filter(!is.na(mass)) %>% # take out problem samples with ambiguous data
  filter(species %in% c("CRA","FOR","FUL","PAR", "SCA")) #core finch species

# explore data ------------------------------------------------------------

head(leuks)
# Calculate leukocytes per 10,000 erythrocytes
leuks <- mutate(leuks, tenKerth = total_erythrocytes/1000) %>%
          mutate(across(heterophil:total_leukocyte, ~./tenKerth, .names = "{col}_per10k")) %>%
          mutate(hl_ratio = heterophil_per10k/lymphocyte_per10k) %>%
          mutate(eos_bas = eosinophil_per10k + basophil_per10k)


# Make some plots
table(leuks$species)
leuks %>%
  ggplot(aes(x = species, y = total_leukocyte_per10k)) +
  geom_boxplot()


leuks %>%
  ggplot(aes(x = species, y = lymphocyte_per10k, fill = species)) +
  geom_boxplot()

leuks %>%
  ggplot(aes(x = species, y = log(monocyte_per10k))) +
  geom_boxplot()

pdf("output_plots/hlratio.pdf")
lplot <- leuks %>% filter(pox_iur == "U") %>%
  filter(species %in% c("CRA","FOR","OLI","PAL","PAR","SCA")) %>%
  mutate(species = factor(species, levels = c("OLI","CRA", "PAR", "PAL", "FUL", "FOR","SCA"))) %>%
  ggplot(aes(x = species, y = lymphocyte_per10k, fill = species)) +
  geom_boxplot() +
  geom_point(color = "#34495e") +
  scale_x_discrete(labels = c("CRA" = "Vegetarian", "FOR" = "Med. ground", "FUL" = "Sm. ground",
                              "PAR" = "Sm. tree", "SCA" = "Cactus")) +
  labs(y = "Lymphocytes per 1000 RBCs", x = "")
dev.off()


bka <- read.csv("input_data/2024GalapagosBKAsFinal.csv") %>%
  rename_all(tolower) %>%
  filter(species != "MAG") %>%
  filter(species != "GAMO")

bplot <- bka %>%
  filter(species %in% c("CRA","FOR","OLI","PAL","PAR","SCA")) %>%
  mutate(species = factor(species, levels = c("OLI","CRA", "PAR", "PAL", "FUL", "FOR","SCA"))) %>%

  filter(age == "A") %>%
  #  filter(pox_iur != "I") %>%
  ggplot(aes(x = species, y = killedperplas,fill = species)) +
  geom_boxplot() +
  geom_point(color = "#34495e") +
  scale_x_discrete(labels = c("CRA" = "Vegetarian", "FOR" = "Med. ground", "FUL" = "Sm. ground",
                              "PAR" = "Sm. tree", "SCA" = "Cactus")) +
  labs(y = "Bacterial killing capacity", x = "")

pdf("output_plots/hlplot_bkaplot.pdf")
lplot /bplot
dev.off()




