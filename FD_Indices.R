#Analisis estructural y de diversidad funcional

##Carga de paquetes
# Load necessary libraries
{
library(readxl)
library(dplyr)
library(writexl)
library(ggplot2)
library(stats)
library(ggmosaic)
library(data.table)
library(jsonlite)
library(curl)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(FD)
library(AICcmodavg)
library(nlme)
library(GA)
library(gawdis)
library(lme4)
library(pbkrtest)
library(lmerTest)
library(emmeans)
library(readr)
library(gow)
library(multcompView)
}
paquetes=c("ggplot2", "fitdistrplus","cowplot", "qqplotr", "ggeffects", "GGally" , "broom", "doBy",
           "corrplot", "pROC","DHARMa", "multcomp", "broom.mixed", "glmmTMB","gamlss.dist") 
lapply(paquetes,FUN=require,character.only=T)

options(digits=3) 

  
setwd("C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Proyecto_EspeciesPiedemonte-Ciudad/AnalisisDatos")

##CARGA DE DATOS
modelo = read.csv("ModeloFiltros_Diversidad.csv", header=T,  sep = ",", stringsAsFactors = T)
traits = read.csv("species_traits.csv", header=T, sep = ",", row.names = 1)
View(traits)
View(modelo)

traits$Grime=as.factor(traits$Grime)
traits$strata=as.factor(traits$strata)
traits$Disp=as.factor(traits$Disp)
traits$Poll=as.factor(traits$Poll)

cover = read.csv("species_cover.csv", header=T, sep = ",")
View(cover)

#exploración de objetos (rasgos)
str(traits)
plot(traits$SLA)
summary(traits$SLA)
plot(traits$SDM)
summary(traits$SDM)
plot(traits$Hmax)
summary(traits$Hmax)
row.names(traits[29,])
which.min(traits$Hmax)
?which.min

#Chequear estratos
strata=as.data.frame(unique(traits$strata))
View(strata)
Arbustos = length(traits$strata[traits$strata=="Shrub"])
Herbaceas = length(traits$strata[traits$strata=="Herbaceous"])
Arboles = length(traits$strata[traits$strata=="Tree"])
Total = length(traits$strata)
View(Arbustos)

#porcentaje de estratos
Arbustos/Total*100
Herbaceas/Total*100
Arboles/Total*100

#Valores faltantes
sum(is.na(traits[,2:7]))
(31*100)/792 #3.91% de Rasgos faltantes
(98*100)/792 #12.4% de rasgos tomados de bilbiografía
(156*100)/792 #19.7% medidos individualmente

100-(3.91+12.4+19.7) #76.6 con TRY

#correlacion entre caracteres numericos
cacacteres.cor = cor(traits [,c(2:4)],use = "complete.obs")
cacacteres.cor
#plot de correlacion
x=traits$SLA
y=traits$Hmax
plot(x,y, xlab = "SDM", ylab = "H", pch=19, col = "blue")
lm_model = lm(y ~ x)
abline(lm_model, lwd = 3)

#matriz solo con los caracteres numericos
traits.numeric = traits[,c(2:4)]
view(traits.numeric)

#matriz completa (sin el origen nativo/exótico)
traits.complete = traits[,c(1:7)]
View(traits.complete)

#Cobertura
#dataframe cobertura de especies en los muestreos
str(cover)
View(cover)
#Transformar a una matriz para poder incluirla en la función dbFD
cover.matrix = as.matrix(cover[,2:132])
View(cover.matrix)
#agregar la columna de plot
rownames(cover.matrix)=cover$plot
View(cover.matrix)
print(cover.matrix)

#chequear que las especies sean las mismas en ambos set de datos
identical(sort(colnames(cover.matrix)), sort(rownames(traits.complete)))

#Si da False, Identificar cuales son las especies que no coinciden en nombres
col_names = colnames(cover.matrix)
row_names = rownames(traits.complete)
# Species in col_names but not in row_names
diff1 = setdiff(col_names, row_names)
print(diff1)
# Species in row_names but not in col_names
diff2 = setdiff(row_names, col_names)
print(diff2)

#Corregir los nombres que salen erroneos (nombre erroneo <- nombre correcto)
colnames(cover.matrix)[colnames(cover.matrix) == "Echinochloa_crus.galli"] <- "Echinochloa_crus-galli"
colnames(cover.matrix)[colnames(cover.matrix) == "Veronica_anagallis.aquatica"] <- "Veronica_anagallis-aquatica"

#Volver a verificar
identical(sort(colnames(cover.matrix)), sort(rownames(traits.complete)))

#Hacer que también coincidan el orden de las columnas y tablas 
cover.matrix <- cover.matrix[, rownames(traits.complete), drop = FALSE]

#Bases de datos finales a usar en los cálculos de FDIV
View(cover.matrix)
View(traits.complete)
View(traits.numeric)

#Calculo de Indices de Diversidad y Estructura funcional, prestar atencion al argumento corrección

#Rasgos Numericos (FD.numeric)
FD.numeric=dbFD(traits.numeric,cover.matrix,w.abun = TRUE, stand.x = TRUE,
                 calc.CWM = TRUE, calc.FDiv = TRUE, corr="cailliez")
FD.numeric = as.data.frame(FD.numeric)
View(FD.numeric)

#Agregado de columna de plot
FD.numeric$plot = cover$plot
FD.numeric=as_tibble(FD.numeric)

#Eliminar valores NA para que las estadisticas corran bien (se elimina el plot "urbano9")
FD.numeric = na.omit(FD.numeric)
View(FD.numeric)

#Eliminación del plot "urbano10" por exceso de outliers
FD.numeric = FD.numeric[FD.numeric$plot != "Urbano10", ]
View(FD.numeric)

################################################################################

#Rasgos Numericos y Categoricos (FD.Complete)

#La CWM de rasgos Categóricos se cálcula como el rasgo que más se repite en el sitio
Disp=traits.complete$Disp = as.factor(traits.complete$Disp)
traits.complete$Poll = as.factor(traits.complete$Poll)
traits.complete$strata = as.factor(traits.complete$strata)
traits.complete$Grime = as.factor(traits.complete$Grime)
distance_matrix=gowdis(traits.complete)
view(distance_matrix)
View(traits.complete)
str(traits.complete)

FD.complete <- dbFD(traits.complete, cover.matrix, w.abun = TRUE, stand.x = TRUE,
                    calc.CWM = T, calc.FDiv = TRUE, corr = "cailliez")
FD.complete = as.data.frame(FD.complete)
View(FD.complete)

#Agregado de columnas de plot, sitio y mosaico
FD.complete$plot = cover$plot
FD.complete=as_tibble(FD.complete)
View(FD.complete)

#Guardar tabla
write.csv(FD.complete, "FD.csv", row.names = FALSE)

#######################Proporción de exoticas###################################
exoticover = read.csv("data_piedemonte-urbano.csv", header=T, sep = ",")
View(exoticover)

library(dplyr)

exotic_cover_by_plot = exoticover %>%
  filter(Origin == "exotica") %>%
  group_by(plot) %>%
  summarise(ExoticRelCover = sum(StrataCover)) %>%
  ungroup()
View(exotic_cover_by_plot)

####################Proporción de ornamentales##################################

ornamentalcover = read.csv("data_piedemonte-urbano.csv", header=T, sep = ",")
View(ornamentalcover)

library(dplyr)

ornamental_cover_by_plot = ornamentalcover %>%
  filter(Ornamental == "ornamental") %>%
  group_by(plot) %>%
  summarise(ExoticRelCover = sum(StrataCover)) %>%
  ungroup()
View(ornamental_cover_by_plot)


##########################Diversidad Taxonómica################################

#Dataframe con covertura relativa de especies x plot
View(cover.matrix)
library(vegan)

# Compute species richness
richness = rowSums(cover.matrix > 0) # Count nonzero entries per site

# Compute Shannon Index
shannon = diversity(cover.matrix, index = "shannon")

# Compute Simpson Index (D = dominance, 1-D = diversity)
simpson = diversity(cover.matrix, index = "simpson")

# Compute Pielou's Evenness (Shannon divided by max possible diversity)
pielou = shannon / log(richness)

# Compute Hill Numbers (Effective number of species)
hill_1 = exp(shannon) # Hill number of order 1
hill_2 = 1 / (1 - simpson) # Hill number of order 2

# Combine all indices into a data frame
taxonomic_indices <- data.frame(
  Site = rownames(cover.matrix),
  Richness = richness,
  Shannon = shannon,
  Simpson = simpson,
  Evenness = pielou,
  Hill_1 = hill_1,
  Hill_2 = hill_2
)

print(taxonomic_indices)
View(taxonomic_indices)

write.csv(taxonomic_indices, "taxonomic_indices.csv", row.names = FALSE)

#Analisis de cantidad de especies por sitio ripario#
#agregar columna rip_modif a cover
cover_hab_trans = left_join(cover, modelo %>% select(plot, rip_modif), by = "plot")
View(cover_hab_trans)
unique(cover_hab_trans$rip_modif)
cover_hab_trans$rip_modif <- trimws(tolower(cover_hab_trans$rip_modif))
cover_hab_trans$rip_modif <- factor(cover_hab_trans$rip_modif,
                                    levels = c("null", "medium", "high"),
                                    labels = c("Null", "Medium", "High"))
species_cols <- setdiff(names(cover_hab_trans), c("plot", "rip_modif"))

# Convert to presence/absence
cover_hab_trans[species_cols] <- cover_hab_trans[species_cols] > 0

# Subset by rip_modif
null_species   <- colSums(cover_hab_trans[cover_hab_trans$rip_modif == "Null", species_cols]) > 0
medium_species <- colSums(cover_hab_trans[cover_hab_trans$rip_modif == "Medium", species_cols]) > 0
high_species   <- colSums(cover_hab_trans[cover_hab_trans$rip_modif == "High", species_cols]) > 0

# Unique and shared species
only_null   <- null_species & !medium_species & !high_species
only_medium <- medium_species & !null_species & !high_species
only_high   <- high_species & !null_species & !medium_species
all_three   <- null_species & medium_species & high_species

# Print counts
cat("Species only in 'Null':", sum(only_null), "\n")
cat("Species only in 'Medium':", sum(only_medium), "\n")
cat("Species only in 'High':", sum(only_high), "\n")
cat("Species in all three categories:", sum(all_three), "\n")

null_only_species   <- names(only_null)[only_null]
null_only_species
medium_only_species <- names(only_medium)[only_medium]
medium_only_species
high_only_species   <- names(only_high)[only_high]
high_only_species
shared_species      <- names(all_three)[all_three]
shared_species
