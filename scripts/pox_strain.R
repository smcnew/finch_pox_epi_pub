# Analyzing pox sequencing results
#

# Load your packages. Do every time you start R/Rstudio
library(ggthemr) # plot params
library(ggplot2) # for making plots
library(dplyr) # for manipulating data

# My easy way to make good looking graphs
# Read more here: https://github.com/Mikata-Project/ggthemr
ggthemr("flat", layout = "clean", text_size = 20)
plot(swatch())

# Load data. Note: I put my data in a folder called "input data" within my project
# folder. Adjust the path to wherever the "pox_sequence_info" file is.
pox <- read.csv("input_data/pox_sequence_info.csv") %>%
  rename_all(tolower)

# plot 'em up.
# Are any of the strains species-specific? Doesn't look like it
pox %>%
  ggplot(aes(x = strain, fill = species)) + geom_bar()

# Compare size of lesions in strain 1 vs. strain 2.
pox %>%
  filter(!is.na(pox_scale)) %>%
  ggplot(aes(x = strain, y = lesion_size)) +
  geom_boxplot()

pdf("output_plots/pox_strain2.pdf")
pox %>%
  filter(!is.na(pox_scale)) %>%
  ggplot(aes(x = pox_scale, fill = strain)) +
  geom_bar(position = "dodge")+
  labs(x = "Infection severity", y = "Individuals", fill = "Pox strain") +
  scale_fill_manual(values = c("#d35400", "#9b59b6"),
    labels = c("Gal 1", "Gal 2")) +
  scale_x_discrete(labels = c("Mild", "Moderate", "Severe", "Very severe"))
dev.off()
# Make the figure look better
pdf("output_plots/pox_strain.pdf")
pox %>%
  filter(!is.na(pox_scale)) %>%
  ggplot(aes(x = strain, fill = pox_scale)) +
  geom_bar() +
  labs(x = "Pox Strain", y = "Individuals", fill = "Infection scale") +
  scale_fill_discrete(labels = c("Mild","Moderate","Severe","Very severe")) +
  scale_x_discrete(labels = c("Strain 1", "Strain 2"))
  dev.off()
plot(swatch())
# Statistics --------------------------------------------------------------


pox_table <- table(pox$pox_scale, pox$strain) %>%
              as.data.frame.matrix %>%
              select(-"mixed infection")
chisq.test(pox_table)  # keeping the four categories gives us p = 0.06

# Lets group things into AB vs CD (this is more complicated code)
pox <- pox %>%
  mutate(group = case_when(pox_scale == "A" ~ "mild",
                           pox_scale == "B" ~ "mild",
                           pox_scale == "C" ~ "severe",
                           pox_scale == "D" ~ "severe"))
pox_table2 <- table(pox$group, pox$strain) %>%
  as.data.frame.matrix %>%
  select(-"mixed infection")

# check it out
pox_table2

# run a chisq test for small samples in a 2x2 (aka fisher's exact test)

fisher.test(pox_table2) # P = 0.019. Sweet!
chisq.test(pox_table2)

head(pox)
pox %>%
  filter(!is.na(pox_scale)) %>%
  t.test(lesion_size ~ strain, data = .)

t.test(lesion_size ~ strain, data = pox)

aggregate(lesion_size ~ strain, data = pox, mean)
