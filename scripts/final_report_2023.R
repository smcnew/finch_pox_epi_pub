# Analysis end of 2023
#devtools::install_github('Mikata-Project/ggthemr')
library(dplyr)
library(ggplot2)
library(sjPlot)
library(ggthemr)
library(lmerTest)
library(lme4)
library(patchwork)

#install.packages("patchwork")
# data --------------------------------------------------------------------
ggthemr("flat", text_size = 22, layout = "clean")
plot(swatch())
rain <- read.csv("raw_data/precip_2023.csv") %>%
  na.omit %>%
  mutate(year_month = (month-1)/12 + year)
captures <- read.csv("raw_data/Campo2023_base_datos_feb.csv") %>%
  rename_all(tolower) %>%
  mutate(pox_iur = case_when(pox_iur == "u" ~ "U", TRUE ~ pox_iur)) %>%
  mutate(species = case_when(species == "FUL " ~"FUL", TRUE ~ species)) %>%
  mutate(taxon = case_when(species == "CRA" ~ "Platyspiza crassirostris",
                           species == "FOR" ~ "Geospiza fortis",
                           species == "FUL" ~ "Geospiza fuliginosa",
                           species == "GAMO" ~ "Mimus parvulus",
                           species == "MAG" ~ "Geospiza magnirostris",
                           species == "MELA" ~ "Coccyzus melacoryphus",
                           species == "MYI" ~ "Myiarchus magnirostris",
                           species == "OLI" ~ "Certhidea olivacea",
                           species == "PAL" ~ "Camarhynchus pallidus",
                           species == "PET" ~ "Setophaga petechia",
                           species == "PSI" ~ "Camarhynchus psittacula",
                           species == "SCA" ~ "Geospiza scandens",
                           species == "ZEN" ~ "Zenaida galapagoensis",
                           species == "PAR" ~ "Camarhynchus parvulus")) %>%
  filter(!is.na(band)) %>%
  filter(!is.na(pox_iur)) %>%
  filter(!is.na(species)) %>%
  filter(!species %in% c("ZEN", "PSI", "MELA","MYI", "PET", "PAL","OLI")) %>%
  filter(!is.na(month)) %>%
  mutate(year_month = (month-1)/12 + year) %>%
  mutate(year_month = as.numeric(as.character(year_month))) %>%
  group_by(band, species) %>%
  mutate(ncaps = n()) %>%
  mutate(capturen = 1:n()) %>%
  mutate(n_estados = n_distinct(pox_iur)) %>% # numero de estados de pox, ej. si se capturan I y despues U, seria 2
  as.data.frame

# check out
dim(captures)
table(captures$taxon) %>% as.data.frame %>% pull(Freq) %>% sum
table(captures$month)

# analysis ----------------------------------------------------------------

tab_xtab(var.col = captures$pox_iur,
         var.row = captures$taxon,
         show.summary = FALSE,
         wrap.labels = 30)

tab_xtab(var.col = captures$pox_iur,
         var.row = captures$site,
         show.summary = FALSE,
         wrap.labels = 30)
head(captures)

table(captures$pox_iur, captures$species)

pdf("output_plots/pox_incidence_23_focalspp.pdf", width = 11)
captures %>% ggplot(aes(x = taxon, fill = pox_iur)) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18)) +
  labs(x = "", y = "Proporción", fill = "Estado") +
  scale_fill_manual(values =  swatch()[c(5,8,6)],
                    labels = c("Infectado", "Recuperado", "No infectado"))
dev.off()
table(captures$type_of_site)
# captures over time  -----------------------------------------------------

caps_per_m <- captures %>%
  #filter(species %in% c("PAR", "FOR", "FUL", "MAG", "SCA", "GAMO",
  #                     "CRA", "MYI")) %>%
  filter(type_of_site == "arida") %>%
  group_by(taxon, year_month) %>%
  tally() %>%
  mutate(year_month = as.numeric(as.character(year_month))) %>%
  as.data.frame


all_caps <- table(captures$year_month) %>%
  as.data.frame %>%
  rename(year_month = Var1, total = Freq) %>%
  mutate(year_month = as.numeric(as.character(year_month))) %>%
  left_join(caps_per_m, .) %>%
  mutate(prop_caps = n/total)



pdf("output_plots/capturas_species_month.pdf", width = 10)
all_caps %>%
  ggplot(aes(x = year_month, y = prop_caps, color = taxon)) +
  geom_line(linewidth = 1.7) +
  labs(x = "Mes", y = "Proporción de capturas", color = "") +
  scale_x_continuous(breaks = sort(unique(all_caps$year_month)),
                     labels = c(2:12,1:2))
dev.off()

pdf("output_plots/capturas_month_prevalence2.pdf", width = 8)
pprev <- captures %>%
  group_by(year_month, pox_iur) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(pox_iur == "I") %>%
  ggplot(aes(x = year_month, y = freq)) +
  geom_smooth() +
  geom_point(size = 3, color = swatch()[3])+
  #geom_col()+
  scale_x_continuous(breaks = sort(unique(all_caps$year_month)),
                   labels = c(2:12,1:2))+
  labs(x = "Mes", y = "Prevalencia de viruela")
dev.off()

# look at precip over the year
pdf("output_plots/rain_2023.pdf", width = 8)
rain %>%
  ggplot(aes(x = year_month, y = precip)) +
  geom_smooth(se=F, color = swatch()[4]) +
  geom_point(size = 3, color = swatch()[4])+
  scale_x_continuous(breaks = sort(unique(all_caps$year_month)),
                     labels = c(2:12,1:2))+
  labs(x = "Mes", y = "Precipitatción (mm)") +
theme(panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA))
dev.off()
pdf("output_plots/rain_2023.pdf", width = 8)
prain<- rain %>%
  ggplot(aes(x = year_month, y = precip)) +
  geom_col(fill = swatch()[4])+
  scale_x_continuous(breaks = sort(unique(all_caps$year_month)),
                     labels = c(2:12,1:2))+
  labs(x = "Mes", y = "Precipitatción (mm)") +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))
dev.off()
pdf("output_plots/rain_prevalence.pdf",height = 9)
pprev/prain
dev.off()
# recapture rate ----------------------------------------------------------


# hay que sacar algunos "bands" malos, ej NA, PETx, FULx MYI NO ANILLADO,
#captures %>% filter(ncaps > 1) %>%
#  dplyr::select(band, species, date, ncaps) %>%
#  arrange(band)

cap_clean <- captures %>%
  filter(!band == "NO ANILLO") %>%
  filter(!is.na(band)) %>%
  filter(species %in% c("PAR", "FOR", "FUL", "MAG", "SCA", "GAMO",
                        "CRA", "MYI"))
head(cap_clean)
# Quedar solo las primeras capturas de cada pajaro
cap_clean <- cap_clean %>% distinct(band, .keep_all = T) %>%
  mutate(pox_iur = as.factor(pox_iur)) %>%
  filter(!is.na(pox_iur))

# Una figura de recapturas
pdf("output_plots/recaptures_2023.pdf", height = 7, width = 10)
cap_clean %>% ggplot(aes(x = ncaps, fill = pox_iur)) +
  geom_bar() +
  labs(x = "Número de capturas", y = "Cuenta", fill = "Estado de infección") +
  scale_fill_manual(values =  swatch()[c(4,3,2)],
  labels = c("Infectado", "Recuperado", "No infectado"))
dev.off()

# La prueba estadistica
glm(ncaps~ relevel(pox_iur, ref = "U") , family = "quasipoisson", data = cap_clean) %>% summary

cap_clean %>% filter(species == "GAMO") %>%
  glm(ncaps ~ relevel(pox_iur, ref = "U"), family = "quasipoisson", data = . ) %>% summary


# create a figure to show what happens to mockingbirds over time
gamos <- captures %>% filter(species == "GAMO") #%>%
 #distinct(band, .keep_all = T) %>%
  #filter(ncaps > 1) %>%



pdf("output_plots/capture_n_gamos.pdf", width = 10)
gamos %>%
  ggplot(aes(x = capturen, color = pox_iur,
             group = as.factor(band), y = reorder(as.factor(band), ncaps, decreasing = T))) +
  geom_line(aes(color = pox_iur), linewidth = 1) +
  geom_point(size = 2) +
  theme(axis.text.y = element_blank()) +
  labs(x = "Número de captura",y = "Individuo", color = ("Estado de pox"))+
 scale_color_manual(values = swatch()[c(5,8,6)],
  labels = c("Infectado", "Recuperado", "No infectado"))
dev.off()

gamos %>%
distinct(band, .keep_all = T) %>% group_by(pox_iur) %>% summarize(mean(ncaps))

head(captures)
library(lmerTest)
glmer(ncaps ~ pox_iur + (1|band), family = "poisson", data = captures) %>%
  summary

captures %>% filter(species == "GAMO") %>%
  glmer(ncaps ~ pox_iur + (1|band), family = "poisson", data = .) %>%
  summary
captures %>%
  filter(capturen == 1) %>%
  filter (species == "GAMO" ) %>%
  group_by (pox_iur) %>%
  summarize (caps = plotrix::std.error (ncaps))
