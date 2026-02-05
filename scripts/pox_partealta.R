# packages
library(ggplot2) # para figuras
library(ggthemr) # para figuras
library(tidyr)
library(lme4)
library(lmerTest)

#install.packages("multcomp")
#install.packages("emmeans")# instalar estos paquete para hacer post-hoc tests
#install.packages("multcompView")

library(multcompView)
library(emmeans)
library(multcomp)
library(scales)

# for multinomial regression and chisq tests
devtools::install_github("ebbertd/chisq.posthoc.test")
library(chisq.posthoc.test)
library(mclogit)
library(dplyr) # para limpieza de datos


# set the theme

ggthemr('light', layout = "clean", text_size = 26)

# clean data  -------------------------------------------------------------


# data
datos <- read.csv("raw_data/base_de_datos_aves2022.csv")
data23 <- read.csv("raw_data/Campo2023_base_datos.csv")
data23 <- data23 %>% mutate(Net_no = as.character(Net_no))
datos <- bind_rows(datos, data23)
datos <- datos %>% filter(Species !="")

datos$Pox_IUR <- as.factor(datos$Pox_IUR)

table(datos$Species)
table(datos$Pox_IUR)
table(datos$type_of_site)
table(datos$zone)
# editar numbres y sitios
datos <- datos %>% mutate(type_of_site = replace(type_of_site, type_of_site == "bosque", "Bosque"))
#datos<- datos[datos$site != "El Barranco",]
datos <- datos %>% filter(site!= "El Barranco")
#UNIR FOR
datos <- datos %>% mutate(Species = replace(Species, Species == "FOR", "FOR "))
datos <- datos %>% mutate(Pox_IUR = replace(Pox_IUR, Pox_IUR == "u", "U"))

datos<- datos %>% mutate(zone=replace(zone, zone== "occidente", "Occidente"))
datos<- datos %>% mutate(zone=replace(zone, zone== "santa Rosa", "Santa Rosa"))
datos<- datos %>% mutate(Pox_IUR = replace(Pox_IUR, Pox_IUR == "u", "U"))
#sacar u despues de pasar a U
datos<- datos[datos$Pox_IUR != "u",] # sacar este hibrido
datos<- datos %>% filter(Pox_IUR!= "u") #sacar el hibrido de idioma dplyr
#ej, si quieres sacar los ANI, los Gallareta, etc. puedes hacerlo de la misma forma
# datos<- datos[datos$Species!= "FOR",]
# datos <- datos %>% filter(Species!= "FOR")
#SACAR u
datos <- datos[datos$Pox_IUR != "u",] # sacar este hibrido
datos<- datos[datos$Pox_IUR != "u",]
datos <- datos %>% filter(Pox_IUR!= "u") %>% droplevels # sacar la opcion de U
table(datos$Pox_IUR)
#sacar N en poxiur
datos<- datos[datos$Pox_IUR!= "N",]
datos <- datos %>% filter(Pox_IUR!= "N")
table(datos$Pox_IUR)
head(datos)
# crear columns para infected, recovered, uninfected
datos <- datos %>% mutate (infected = case_when(Pox_IUR == "I" ~ 1,
                                                TRUE ~ 0),
                           recovered = case_when(Pox_IUR == "R" ~ 1,
                                               TRUE ~ 0),
                           uninfected = case_when(Pox_IUR == "U" ~ 1,
                                                  TRUE ~0),
                           sign_pox = case_when(Pox_IUR == "I" ~ 1,
                                                Pox_IUR == "R" ~ 1,
                                                Pox_IUR == "U" ~ 0))

datos_clean <- datos %>% filter(!grepl("FUL", Band)) %>%
  filter(!grepl("PET", Band)) %>%
  filter(!grepl("MYI", Band)) %>%
  filter(!grepl("PAR", Band)) %>%
  filter(!Band == "NO ANILLADO") %>%
  filter(!is.na(Band))
datos <- datos %>% mutate(type_of_site = replace(type_of_site, type_of_site == "bosque", "Bosque"))
write.csv(datos_clean, "formatted_data/formatted_capture_data.csv")

#Diversidad de especies  -------------------------------------------------------------------
#install.packages("iNEXT")
library(iNEXT)
# Necesitamos un matrix de especies por tipo de sitio
sitioxspp <- table(datos$Species, datos$type_of_site) %>% as.data.frame.matrix
head(sitioxspp)
abun.raw <- iNEXT(sitioxspp, datatype="abundance")
ggiNEXT(abun.raw, type=1, se=TRUE, facet.var="None", color.var="Assemblage", grey=FALSE)
head(datos)
table(datos$type_of_site, datos$Species)
head(abun.raw)
table(datos$Species, datos$type_of_site)






#por especie
head(datos)
#PREVALENCIA POR ESPECIE
#limpiar los datos de especie
datos<- datos[datos$Species != "ANI",]
datos <- datos %>% filter(Species!= "ANI")  #sacar el hibrido de idioma dplyr
datos <- datos[datos$Species != "GALLARETA",] # sacar este hibrido
datos <- datos %>% filter(Species != "GALLARETA") #sacar el hibrido de idioma dplyr

datos <- datos[datos$Species != "MAR",] # sacar este hibrido
datos <- datos %>% filter(Species != "MAR") #sacar el hibrido de idioma dplyr
datos <- datos[datos$Species != "SCA",] # sacar este hibrido
datos <- datos %>% filter(Species != "SCA") #sacar el hibrido de idioma dplyr

datos <- datos[datos$Species != "MELA",] # sacar este hibrido
datos <- datos %>% filter(Species != "MELA") #sacar el hibrido de idioma dplyr
head(datos)
#CREAR COLUMNA codigos por nombres comunes
datos<-datos %>% mutate (Nombre_comun = case_when(Species == "CRA" ~ "pinzon vegetariano",
                                            Species == "FOR " ~ "pinzon mediano de tierra",
                                            Species =="FUL" ~ "pinzon pequeño de tierra",
                                            Species== "GAMO" ~ "cucuve",
                                            Species== "MYI" ~ "papamosca",
                                            Species=="OLI" ~ "pinzon cantor",
                                            Species=="PAL" ~ "pinzon carpintero",
                                            Species =="PAR" ~ "pinzon pequeño de arbol",
                                            Species=="PET" ~ "canario maria",
                                            Species=="PSI" ~ "pinzon grande de arbol",
                                            TRUE ~ Species))
datos<-datos %>% mutate (Nombre_comun = case_when(Species == "CRA" ~ "Vegetarian finch",
                                                  Species == "FOR " ~ "Med. ground finch",
                                                  Species =="FUL" ~ "Sm. ground finch",
                                                  Species== "GAMO" ~ "G. mockingbird",
                                                  Species== "MYI" ~ "G. flycatcher",
                                                  Species=="OLI" ~ "Warbler finch",
                                                  Species=="PAL" ~ "Woodpecker finch",
                                                  Species =="PAR" ~ "Sm. tree finch",
                                                  Species=="PET" ~ "Yellow warbler",
                                                  Species=="PSI" ~ "Lg. tree finch",
                                                  Species=="MAG" ~ "Lg. ground finch",
                                                  TRUE ~ Species))



# Prevalence by species ------------------------------------------


total_poresp <- table(datos$Nombre_comun, datos$Pox_IUR) %>% as.data.frame.matrix # hay que convertirlo en data frame
total_poresp$Infected <- 100 * total_poresp$I / (total_poresp$I + total_poresp$R + total_poresp$U)
total_poresp$Recovered <- 100 * total_poresp$R / (total_poresp$I + total_poresp$R + total_poresp$U)
total_poresp$Nombre_comun <- rownames(total_poresp)

total_porespIR <- total_poresp %>% select(Infected, Recovered, Nombre_comun) %>%
                  pivot_longer(cols = c(Infected, Recovered), names_to = "Status")
pdf("output_plots/ir_prevalence.pdf", width = 10, height = 7)
total_porespIR %>%
  ggplot(aes(y = value, x = Nombre_comun, fill = Status)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=30)) +
  labs (x = "", y = "Percent of captures")
dev.off()

total_intensidad <- table(datos$Nombre_comun, datos$Pox_scale) %>% as.data.frame.matrix
total_intensidad$total <- table(datos$Nombre_comun)
total_intensidad <- total_intensidad %>% mutate(perc_a = 100 * A/total,
                                                perc_b = 100 * B/total,
                                                perc_c = 100 * C/total,
                                                perc_d = 100 * D/total) %>%
  tibble::rownames_to_column(var = "Common_name") %>%
  dplyr::select(Common_name, perc_a, perc_b, perc_c, perc_d) %>%
  pivot_longer(cols = c(perc_a, perc_b, perc_c, perc_d)) %>% as.data.frame %>%
  mutate(name = as.factor(name)) %>%
  filter(Common_name != "Lg. ground finch", Common_name != "Lg. tree finch")



total_intensidad %>% ggplot(aes(x = Common_name, y = value, fill = forcats::fct_rev(name)))+
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=22)) +
  labs (x = "", y = "Prevalence (%)", fill = "Infection severity") +
  scale_fill_discrete(labels=c('4', "3", "2", "1"))

filter(datos, Species == "MAG")
head(datos)
datos %>% filter(year == 2023) %>% nrow
datos %>% filter(year == 2023) %>% filter (Recaptured == 1)%>% nrow


# barplot(total_poresp$prevalencia,
#         names.arg = rownames(total_poresp),
#                col = palette() ,
#         xlab = "Especies",
#         ylab = "Prevalencia (%)",
#         main = "Prevalencia de POX",
#         axes = names(total_poresp), las=2)
pdf("output_plots/species_prevalence.pdf", height = 7, width = 10)
total_poresp %>%
  ggplot(aes(x = Nombre_comun, y = prevalencia)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=22)) +
  labs (x = "", y = "Prevalence (%)")
dev.off()

pdf("output_plots/speciesabundprev.pdf", height = 7, width = 10)
datos %>%
  ggplot(aes(x = Nombre_comun, fill = Pox_IUR)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=22)) +
  labs (x = "", y = "Captures", fill = "Pox status") +
  scale_fill_discrete(labels=c('Infected', 'Recovered', 'Uninfected'))
dev.off()

datos %>%
  ggplot(aes(x = Nombre_comun, fill = Pox_IUR)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=22)) +
  labs (x = "", y = "Captures", fill = "Pox status") +
  scale_fill_discrete(labels=c('Infected', 'Recovered', 'Uninfected'))


#prevalencia
total_poresp %>%
  ggplot(aes(x = Nombre_comun, y = prevalencia ))+
  geom_col()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14))+
           labs (
             x = "",
             y = "Prevalencia (%)")


# Figuras
# Mirando numeros de infectados

# pregunta de prevalencia vs. numero de infectados.
#
# Por ejemplo, si capturamos 100 FUL, y 20 tienen pox, que es la prevalencia?
# Prevalencia = el porcentaje (o la proporcion) que tienen pox
# en este caso = 20 / 100 = 20%
#
# Ej: si capturamos 400 FOR, y 30 tiene pox, que es la prevalencia? = 7.5%
# Cual especie tiene mas pox? seria FUL porque tiene mas prevalencia (la
# prevalencia es mas alta) en FUL.
#
# Ej: si capturamos 50 CRA y 20 tiene pox, que es la prevalencia? 40%
#
# listo con prevalencia
#ANALISIS DE DATOS DE PREVALENCIA figura basica:



# Prevalence over time  ---------------------------------------------------


table(datos$Pox_IUR)
datos$month_year <- paste(datos$month, datos$year, sep = "/")
total_pormes  <- table(datos$month_year, datos$Pox_IUR) %>% as.data.frame.matrix #hay que convertirlo en data frame
total_pormes$prevalencia <- 100 * total_pormes$I / (total_pormes$I + total_pormes$R + total_pormes$U)
total_pormes$Month <-c("Oct. '22", "Feb '23", "Mar '22", "Mar '23", "Apr '22", "Apr '23", "May '22", "May '23", "Jun '22", "Jul '22", "Aug '22", "Sep '22")
total_pormes$monthcode <- rownames(total_pormes)

#ordenar los meses
total_pormes$Month <- factor(total_pormes$Month, levels = c("Mar '22", "Apr '22", "May '22", "Jun '22", "Jul '22", "Aug '22", "Sep '22","Oct. '22", "Feb '23", "Mar '23", "Apr '23", "May '23"))
pdf("output_plots/prevalence_permonth.pdf", height = 7, width = 10)
total_pormes %>%
  ggplot(aes(x = Month, y = prevalencia ))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=26))+
  labs (
    x = "",
    y = "Prevalence (%)")
dev.off()
total_pormes %>%
  ggplot(aes (x = prevalencia, y = rain))+
  geom_point()
total_pormes

# Prevalence by habitat ---------------------------------------------------


#sitio no hay diferencia
# metadata<- datos %>%  select(site,type_of_site) %>% unique
# table(metadata$type_of_site)
# head(datos)
# table(datos$Pox_IUR)
total_porsitios  <- table(datos$type_of_site, datos$Pox_IUR) %>% as.data.frame.matrix #hay que convertirlo en data frame
total_porsitios$prevalencia <- 100 * total_porsitios$I / (total_porsitios$I + total_porsitios$R + total_porsitios$U)
rownames(total_porsitios) # nombres de las filas
total_porsitios$Habitat <- c("Farm","Arid/Urban", "Scalesia Forest", "Pasture")

# barplot(total_porsitios$prevalencia,
#         names.arg = rownames(total_porsitios),
#         col = c("red","green","blue"),
#
#         xlab = "Sitio",
#         ylab = "Prevalencia (%)",
#         main = " Prevalencia de POX ",
# )
pdf("output_plots/pox_site.pdf", height = 7, width = 7)
total_porsitios %>%
  ggplot(aes(x = Habitat, y = prevalencia))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs (
    x = "",
    y = "Prevalence (%)") +
theme(text = element_text(size = 28), axis.text.x = element_text(size = 28))
dev.off()
pdf("output_plots/pox_site.pdf", height = 7, width = 7)
datos %>%
  ggplot(aes(x = type_of_site, fill = Pox_IUR)) +
  geom_bar(position = "fill") +
  theme(text = element_text(size = 22)) +
  scale_y_continuous(labels=percent_format()) +
  labs (y = "Proportion", x = "", fill = "Pox status") +
  scale_fill_discrete(labels=c('Infected', 'Recovered', 'Uninfected')) +
  scale_x_discrete(labels = c("Agricultural", "Forest", "Pasture"))
dev.off()

# stat test


table (datos$Pox_IUR, datos$type_of_site) %>% chisq.posthoc.test()


mblogit(Pox_IUR ~ type_of_site, random = ~1|site, data = datos) %>% summary
mblogit(Pox_IUR ~ Species, random = ~1|site, data = datos) %>% summary



# Recaptura  --------------------------------------------------------------
# anadir una columna con el numero de capturas (ncaps)
datos <- datos %>% group_by(Band, Species) %>% mutate(ncaps = n()) %>% as.data.frame
filter(datos, Band == "SM405") %>% dplyr::select (Band, date, site)
# hay que sacar algunos "bands" malos, ej NA, PETx, FULx MYI NO ANILLADO,
datos %>% filter(ncaps > 1) %>% dplyr::select(Band, Species, date, ncaps) %>% arrange(Band)

datos_clean <- datos %>% filter(!grepl("FUL", Band)) %>%
  filter(!grepl("PET", Band)) %>%
  filter(!grepl("MYI", Band)) %>%
  filter(!grepl("PAR", Band)) %>%
  filter(!Band == "NO ANILLADO") %>%
  filter(!is.na(Band))
#averiguar que todos que queden tienen anillos verdades
datos_clean %>% dplyr::filter(ncaps > 1) %>% dplyr::select(Band, Species, date, ncaps) %>% arrange(Band)

# Quedar solo las primeras capturas de cada pajaro
datos_clean <- datos_clean %>% distinct(Band, .keep_all = T) %>%
              mutate(Pox_IUR = as.factor(Pox_IUR))


# Una figura de recapturas
pdf("output_plots/recaptures.pdf", height = 7, width = 10)
datos_clean %>% ggplot(aes(x = ncaps, fill = Pox_IUR)) +
  geom_bar() +
  labs(x = "Number of individual captures", y = "Count", fill = "Pox status") +
  scale_fill_discrete(labels = c("Infected", "Recovered", "Uninfected")) +
  theme(text = element_text(size=22))
dev.off()
# La prueba estadistica
glm(ncaps~ relevel(Pox_IUR, ref = "U") , family = "quasipoisson", data = datos_clean) %>% summary


# Estadistica  ------------------------------------------------------------

head(datos)
# 1. ?varia la prevalencia entre espeices?
# Usaremos un "generalized linear mixed model" para probar si la probabilidad
# de ser infectado depende de la especie. Los (1|site) se llaman
# "random effects." Esto controla por variantes en sitio
# teniamos repeticiones
mod1 <- glmer(infected ~ Species + (1|site), data = datos, family = "binomial" )
summary(mod1) #ver el resultado de este modelo

# hacer comparaciones entre todos los pares de especies
# la letra nos dice cuales grupos son distintos. Por ej, GAMO es la unica especie
# en el grupo "C," lo que quiere decir es que GAMO tiene una major probabilidad
# de ser infectado que todas las otras especies.
# Al contrario, PET es en el grupo "A", lo que quiere decir es que la probabilidad
# de inffeccion en esta especie es menor que las otras especies.
emmeans(mod1, specs =  pairwise ~ Species ) %>%
  cld(.,
      Letters = letters,
      alpha = 0.05)

# ahora probar diferencias entre tipos de habitat
mod2 <- glmer(infected ~ type_of_site + (1|site), data = datos, family = "binomial")
summary(mod2) # no hay muchas differencias entre tipos de estudios
emmeans(mod2, specs = pairwise ~ type_of_site) %>%
  cld(., Letters = letters)

# probar si la taza de pox varia entre meses
mod3 <- glm(infected ~as.factor(month), data = datos, family = "binomial" )
emmeans(mod3, specs = pairwise ~ month) %>%
  cld(., Letters = letters) # Abril tenia mas pox, septiembre tenia poco pox, pero aparte todos los meses no tienen diferencias significativas



# emigration --------------------------------------------------------------

head(datos_clean)
transitions <- datos_clean %>% filter(ncaps > 1) %>% group_by(Band) %>% summarize(nsites =n_distinct(site))
emigrants <- transitions %>% filter(nsites > 1) %>% pull (Band)
filter(datos, Band %in% emigrants) %>% filter(year == 2022) %>% select (date, zone, site, Band, Species, ncaps) %>% arrange(Band)
