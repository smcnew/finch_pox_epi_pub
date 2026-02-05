# Instalacion de paquetes: devtools y ggthemr para figuras, solo hay que
# hacerlo una vez
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("devtools")
# devtools::install_github('Mikata-Project/ggthemr')

#kiman baya es ganadera
#los manzaillos debe ser bosque
# packages
library(lme4)
library(lmerTest) # para estadistica
library(ggplot2) # para figuras
library(ggthemr) # para figuras
library(tidyr) # para limpieza de datos
library(dplyr) # para limpieza de datos
# hoja de datos
install.packages("usethis")

datos <- read.csv("MISTNETTING_BIRD_DL_FINAL_2022.csv")
datos$Pox_IUR <- as.factor(datos$Pox_IUR) # debe ser una "categoria" no "caracter"


# Quiero una columna que es binary, si o no infectado,
# Tambien quiero una columa que es binary, si o no recuperado
datos <- datos %>% mutate(pox = case_when(Pox_IUR == "I" ~ "1",
                                 Pox_IUR == "R" ~ "0",
                                 Pox_IUR == "U" ~ "0")) %>%
                   mutate(pox = as.factor(pox)) %>%
                   mutate(recuperado = case_when(Pox_IUR == "I" ~ "0",
                                                  Pox_IUR == "R" ~ "1",
                                                  Pox_IUR == "U" ~ "0")) %>%
                  mutate(recuperado = as.factor(recuperado))

datos <- datos %>% mutate(type_of_site = case_when(Site == "Quiman Valle" ~ "Ganadera",
                                                   TRUE ~ type_of_site ))

# limpiar los datos de especie
# Sacar species que tienen menos de 10 registros
datos <- datos %>% group_by(Species) %>%
  filter(n() > 10)

# Limpieza de datos -------------------------------------------------------

# forma de sacar especies particulares
# datos <- datos %>% filter(Species != "SCA/FOR OR SCA/FUL?") %>%
#                  filter(Species != "MAR") # no se que es mar
table(datos$Species)


# Fijar pox
#datos <- datos %>% mutate(Pox_IUR = replace(Pox_IUR, Pox_IUR == "u", "U"))
##ya lo limpie en la hoja de excel
table(datos$Pox_IUR)

# Por las dudas, averiguar las medidas
str(datos$Mass) # esto nos dice que es numerico, bien
summary(datos$Mass) # fiajmos que las masas son entre 7 - 89g, no hay errores obvios

summary(datos$Tarsus) # alguien tiene un tarsus = 62 centimetros... error?
filter(datos, Tarsus > 60) # mirar este grande, es un FUL. Sacamos, porque no se que debe tener como tarsus
datos <- filter(datos, Tarsus < 60) # quedemos solo los casos con tarsus < 60

summary(datos$Wing) #aca alguien tiene un Wing 123 mm. Error?
filter(datos, Wing > 100) # no son errores, todas las aves con Wing > 100 son cucuves, gallaretas, cuquillos

# mirar sitios, arreglar Shiss to Schiss y Ramon
table(datos$Site)
datos <- datos %>% mutate(Site = replace(Site, Site == "Pieter Shiss", "Pieter Schiss"))
datos <- datos %>% mutate(Site = replace(Site, Site == "Miguel Ram\xf3n", "Miguel Ramon"))


# Recaptura  --------------------------------------------------------------

# anadir una columna con el total de capturas por individuo (ncaps)
# anadir una columna con el numero de captura por cada instancia
# anadir una columna para ver cuantas distinctas estados de pox tienen. Por
# ejemplo, si se captura siempre en U, seria "1", si se captura U y despues I, seria "2"

datos <- datos %>% group_by(Species) %>%
  filter(n() > 10)


datos <- datos %>% group_by(Band, Species) %>%
  mutate(ncaps = n()) %>%
  mutate(capturen = 1:n()) %>%
  mutate(n_estados = n_distinct(Pox_IUR)) %>% # numero de estados de pox, ej. si se capturan I y despues U, seria 2
  as.data.frame
write.csv(datos, "data_with_capturen.csv", row.names = F)
table(datos$Site)

# Un ejemplo: mirar un pajaro capturado muchas veces
filter(datos, Band == "SM405") %>% select (Band, Date, Site, ncaps, capturen, n_estados, Pox_IUR, Notes)

# hay que sacar algunos "bands" malos, ej NA, PETx, FULx MYI NO ANILLADO,
datos %>% filter(ncaps > 1) %>% select(Band, Species, Date, ncaps) %>% arrange(Band)

datos_clean <- datos %>% filter(!grepl("FUL", Band)) %>%
                         filter(!grepl("PET", Band)) %>%
                         filter(!grepl("MYI", Band)) %>%
                         filter(!grepl("ANI", Band)) %>%
                         filter(!Band == "NO ANILLADO") %>%
                         filter(!Band == "PAR") %>%
                         filter(!is.na(Band))
#averiguar que todos que queden tienen anillos verdades
datos_clean %>% filter(ncaps > 1) %>%
  select(Band, Species, Date, ncaps, Pox_IUR, capturen, n_estados,) %>%
  arrange(Band)

filter(datos_clean, Pox_IUR == "N")
datos_clean <- datos_clean %>%
  mutate(Pox_IUR = replace(Pox_IUR, Pox_IUR == "N", "U")) %>%
  mutate(Pox_IUR = replace(Pox_IUR, Pox_IUR == "u", "U"))

datos_clean <- datos_clean %>% droplevels()

# Ver los cambios en estado de infeccion
datos_clean %>% filter(n_estados >1 ) %>%
  ggplot(aes(x = as.numeric(capturen), y = as.numeric(Pox_IUR), color = Band)) +
  geom_point() + geom_line() +
  scale_y_continuous(name = "Estado de pox", limits = c(1,3),
                     breaks = c(1,2,3), labels = c("Infectado", "Recuperado", "No inf."))+
  labs(x = "Numero de captura")

# Quedar solo las primeras capturas de cada pajaro
datos_clean <- datos_clean %>% distinct(Band, .keep_all = T) %>%
  drop_na(Pox_IUR) %>%
  filter(ncaps < 4)
head(datos_clean)

# Una figura de recapturas
pdf("grafico_recaps.pdf")
datos_clean %>% ggplot(aes(x = ncaps, fill = Pox_IUR)) +
  geom_bar() +
  scale_fill_discrete(labels = c('Infectado','Recuperado','No infectado')) +
  labs(x = "Numero de capturas", y = "Numero de individuos", fill = "Pox status")
dev.off()

  datos_clean <- datos_clean %>% filter(is.na(Pox_IUR))
# La prueba estadistica
glm(ncaps ~ Pox_IUR , family = "quasipoisson", data = datos_clean) %>% summary


# Un ejemplo para entender LM
# pregunta: los FOR pesan mas que los FUL?
datosf <- filter(datos_clean, Species == "CRA"|  Species == "FOR" | Species == "FUL") %>% droplevels
str(datosf$Tarsus)
boxplot(Tarsus ~ Species, data = datosf)
lm(Tarsus ~ Species, data = datosf) %>% summary

# En comparacion con los infectados, lo no infectados tienen
# una mayor probabilidad de ser recapturados varias veces.
# No sabemos que esta pasando con los que no recapturamos. Puede
# ser que no cayeron de nuevo porque estaban enfermos y estaban escondidos en un arbol
# Puede ser que no cayeron porque se mudaron a otro territoreo.
# Sin embargo, estos datos nos dan un piste muy importante y relevante
# que *algo* esta pasando con los infectados, y podemos suponer que parte
# de la razon porque no los encontramos de nuevo es porque se murieron.
#




# Estadistica -------------------------------------------------------------

# Hay variacion entre habitats?
head(datos)

# Queremos que la referencia sea Bosque (lo natural)
datos$type_of_site <- relevel(as.factor(datos$type_of_site), ref = "Bosque")

# Modelo basico: la probabilidad de infeccion depende del tipo de sitio?
# Usamos un "generalized linear mixed model" que tambien tiene site como
# "random effect." Necesitamos el random effect porque visitamos varios sitios
# varias veces, y puede haber variacion entre los sitios que afecta los resultados
# aparte de su categorizacion.

glmer(pox ~ type_of_site + (1|Site), family = "binomial", data = datos) %>% summary

# no hay diferencias significativas

# Varia la prevalencia de pox entre especies?
# Necesitamos eligir una referencia de comparacion. Eligo FUL porque es lo mas commun
datos$Species <- relevel(as.factor(datos$Species), ref = "FUL")
glmer(pox ~ Species + (1|Site), family = "binomial", data = datos) %>% summary


#Esto nos dice: en comparision con FUL, GAMOS tienen una probabilidad mayor
#de infeccion, y los PET tienen una probabilidad menor de infeccion.

# Que tal recuperados? hay variacion en la proporcion de recuperados
# La mayoria de especies tienen menos recuperados que los FUL, la diferencia
# es significativa en OLI, PAL, y PET
glmer(recuperado ~ Species + (1|Site), family = "binomial", data = datos) %>% summary


# Varia pox trans tiempo? No hay un trend linear, pero necesitamos otro modelo

# glmer(pox ~ day_of_year + (1|Site), family = "binomial", data = datos) %>% summary

# Otra idea, usar un GAM (generalized additive model). Este modelo nos ayuda
# a comparar tendencias no-lineares
install.packages("mgcv")
library(mgcv)

mod1 <- gam(pox ~ s(day_of_year), data = datos, family = "binomial")

summary(mod1) # algo de variacion pero p = 0.07 asi que no es significativa

plot(mod1, xlab ="Dia del año", ylab = "Probabilidad de pox")


# Figuras -----------------------------------------------------------------
# Usando ggplot y ggthemr
ggthemr("solarized") # eligo colores a traves de ggthemr https://github.com/Mikata-Project/ggthemr

# Un ejemplo de una figura de totales de infecciones
datos_pajaros %>%
  ggplot(aes(x = Species, fill = Pox_IUR)) +
  geom_bar() +
  labs(
    x = "Especies",
    y = "Numero",
    fill = "Pox",
    title = "Total de capturas"
  ) +
  scale_fill_discrete(labels=c('Infectado', 'Recuperado', 'No infectados'))

# Proporciones
datos_pajaros %>%
  ggplot(aes(x = Species, fill = Pox_IUR)) +
  geom_bar(position = "fill") +
  labs(
    x = "Especies",
    y = "Proporcion",
    fill = "Pox",
    title = "Proporcion de infeccion"
  )

# Prevalencia
prevalencia <- table(datos_pajaros$Species, datos_pajaros$Pox_IUR ) %>% as.data.frame.matrix %>%
  mutate(prev = 100* I / (I + R +U))
prevalencia$Species <- rownames(prevalencia) # anadir una columna con especie

prevalencia %>% ggplot(aes(x = Species, y = prev)) +
  geom_bar(stat = "identity", fill = "grey31") #cambiar color si quieres


# Diversidad de especies  -------------------------------------------------------------------
#install.packages("iNEXT")
library(iNEXT)
# Necesitamos un matrix de especies por tipo de sitio
sitioxspp <- table(datos$Species, datos$type_of_site) %>% as.data.frame.matrix

head(sitioxspp)
abun.raw <- iNEXT(sitioxspp, datatype="abundance")
ggiNEXT(abun.raw, type=1, se=TRUE, facet.var="None", color.var="Assemblage", grey=FALSE)

#
# Otras cosas, no organizadas ---------------------------------------------



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

# Una figura basica:
#
head(datos)
total_poresp <- table(datos$Species, datos$Pox_IUR) %>% as.data.frame.matrix # hay que convertirlo en data frame

total_poresp$prevalencia <- 100 * total_poresp$I / (total_poresp$I + total_poresp$R + total_poresp$U)
head(total_poresp)
rownames(total_poresp) # nombres de las filas
barplot(total_poresp$prevalencia,
        names.arg = rownames(total_poresp))


# Diversidad por sitio
head(datos)
metadata <- datos %>% select(Site, type_of_site) %>% unique
table(metadata$type_of_site)

head(datos)
library("khroma")
bright <- colour("bright")
plot_scheme(bright(7), colours = TRUE, names = TRUE, size = 0.9)
bright(3)
table(datos$Species, datos$type_of_site)

# un ejemplo
table(datos$Pox_IUR, datos$Species) %>%
  barplot(col = bright(3),
          ylab = "Numero de capturas")


# Estadistica
# Modelo linear (linear model) y general linear model

head(datos)
plot(Mass ~ Wing, data = datos)
abline(h = 40)
abline(lm(Mass ~ Wing, data = datos))

# un modelo linear basico
# se usa cuando la variable dependiente es continuo y un numero (ej. peso/masa)
lm(Mass ~ Wing, data = datos) %>% summary

head(datos)

# generalized linear model
# se usa cuando la variable dependiente es si/no, binary (?) (ej. infectado, no infectado)
# as.factor = categoria en vez de caracter, que es palabras
glm(as.factor(Pox_IUR) ~ Temperature,
    data = datos[datos$Pox_IUR != "R",], family = "binomial") %>% summary

# p value = 0.08, no es significativo. P < 0,05, ahi decimos que hay una
# relacion significtiva. Sin embargo, hay una tendencia para ver mas
# pox (?) cuando hace mas frio


# pca ---------------------------------------------------------------------

pcatest <- datos %>% filter(Species %in% c("FOR")) %>% as.data.frame %>% select(Tarsus, Wing, Beak_length) %>%
  na.omit %>% prcomp

summary(pcatest)
pcatest
