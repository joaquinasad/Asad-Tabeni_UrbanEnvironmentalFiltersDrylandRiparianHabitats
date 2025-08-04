####Relaciones Filtros-Diversidad####

# Paquetes requeridos
paquetes <- c(
  "ggplot2", "ggmcmc", "agricolae", "fitdistrplus", "MASS", "cowplot", "qqplotr", 
  "ggeffects", "GGally", "broom", "doBy", "corrplot", "DHARMa", "pROC", "multcomp", 
  "multcompView", "car", "broom.mixed", "glmmTMB", "gamlss.dist", "bayesplot", 
  "reshape2", "gridExtra", "brms", "emmeans", "DirichletReg", "readxl", "tidyr", 
  "dplyr", "writexl", "stats", "ggrip_modif", "data.table", "jsonlite", "curl", 
  "tidyverse", "vegan", "FD", "AICcmodavg", "nlme", "GA", "gawdis", "lme4", 
  "pbkrtest", "lmerTest", "lavaan", "readr"
)
# Cargar paquetes
invisible(lapply(paquetes, require, character.only = TRUE))
options(digits = 3)

setwd("C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Proyecto_EspeciesPiedemonte-Ciudad/AnalisisDatos")

##CARGA DE DATOS
modelo = read.csv("ModeloFiltros_Diversidad.csv", header=T,  sep = ",", stringsAsFactors = T)
View(modelo)

modelo$rip_modif=as.factor(modelo$rip_modif)
summary(modelo$rip_modif)
sd(modelo$Richness)

#correlacion de indices de diversidad funcional
FD.cor = cor(modelo [,c(22:26)],use = "complete.obs")
FD.cor

citation("FD")
citation("lavaan")
citation("vegan")
citation("lme4")
citation("glmmTMB")
citation("DirichletReg")
citation("MASS")


#plots rapidos de relación lineal entre variables
#X=filtro
x=modelo$imprevious_area
#y= patron de diversidad o estructura funcional
y=modelo$mean_ST
{ 
plot(x,y, xlab = "", ylab = "", pch=19, col = "darkgrey")
lm_model = lm(y ~ x)
abline(lm_model, lwd = 5)
x_pred <- seq(min(x), max(x), length.out = 1000) 
pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
polygon(c(x_pred, rev(x_pred)), 
        c(pred[, 2], rev(pred[, 3])), 
        col = rgb(1, 0, 0, alpha = 0.1), border = NA)
r2 <- summary(lm_model)$r.squared
mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")
}

###Distribución de las variables respuesta (Solo modificar la variable post $)
normal=fitdist(modelo$CWM.SLA,"norm")
poisson=fitdist(modelo$CWM.SLA,"pois")
lognormal=fitdist(modelo$CWM.SLA,"lnorm") 
negbinom=fitdist(modelo$CWM.SLA,"nbinom")

cdf.sp=cdfcomp(list(normal, poisson,lognormal, negbinom), main="", legendtext =c("Normal", "Poisson","Lognormal","BinNeg"),fitcol = c("orange", "blue","red", "green"), plotstyle ="ggplot")+
  geom_line(size=1.2)+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16), 
        legend.position = c(0.70,0.25),
        legend.text=element_text(size=16))
qq.sp=qqcomp(list(normal, poisson,lognormal, negbinom), main="",fitcol = c("orange","blue","red", "green"), plotstyle 	="ggplot")+
  geom_line(size=1.2)+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position ="none")
#graficos de distribuciónes
plot_grid(cdf.sp, qq.sp, ncol = 2)
gofstat(list(normal, poisson, lognormal, negbinom))$aic

#logaritmizar
modelo$log_SLA = log(modelo$CWM.SLA)

hist(modelo$CWM.SLA)
hist(modelo$log_SLA)

####### OBJETIVO 1 - ESTRUCTURA FUNCIONAL Y FILTROS AMBIENTALES URBANOS ########

############### oNE-wAY ANALYSIS OF VARIANCE ###################################

### FILTRO: TRANSFORMACIÓN DEL HÁBITAT
#Variable explicativa: Nivel de transformación del hábitat

###Rasgos categóricos
################# Rasgos Categóricos - Transformación Riparia ##################
data = read.csv("data_piedemonte-urbano.csv", header=T, sep = ",", stringsAsFactors = T)
View(data)
View(modelo)
colnames(data)

data$CSR=as.factor(data$CSR)

#Ver rasgos
summary(data$CSR)

C_strat_cover = data %>%
  filter(CSR == "C") %>%
  group_by(plot) %>%
  summarise(C = sum(StrataCover)) %>%
  ungroup()
View(C_strat_cover)

CR_strat_cover = data %>%
  filter(CSR == "CR") %>%
  group_by(plot) %>%
  summarise(CR = sum(StrataCover)) %>%
  ungroup()
View(CR_strat_cover)

CS_strat_cover = data %>%
  filter(CSR == "CS") %>%
  group_by(plot) %>%
  summarise(CS = sum(StrataCover)) %>%
  ungroup()
View(CS_strat_cover)

R_strat_cover = data %>%
  filter(CSR == "R") %>%
  group_by(plot) %>%
  summarise(R = sum(StrataCover)) %>%
  ungroup()
View(R_strat_cover)

RS_strat_cover = data %>%
  filter(CSR == "RS") %>%
  group_by(plot) %>%
  summarise(RS = sum(StrataCover)) %>%
  ungroup()
View(RS_strat_cover)

S_strat_cover = data %>%
  filter(CSR == "S") %>%
  group_by(plot) %>%
  summarise(S = sum(StrataCover)) %>%
  ungroup()
View(S_strat_cover)

write_xlsx(grime_relcover, "grime.xlsx")

#Unir tablas
df_list = list(C_strat_cover, CR_strat_cover, CS_strat_cover, R_strat_cover, RS_strat_cover, S_strat_cover)
grime_relcover = Reduce(function(x, y) merge(x, y, by = "plot", all = TRUE), df_list)
View(grime_relcover)
write_xlsx(grime_relcover, "grime.xlsx")

#agregar columna rip_modif
grime_relcover = left_join(grime_relcover, modelo %>% select(plot, rip_modif), by = "plot")
#Cambiar formato de tabla para mejor analisis
grime_relcover = grime_relcover %>%
  pivot_longer(cols = c(C, CR, CS, R, RS, S),
               names_to = "Strategy",
               values_to = "RelativeCover") %>%
  filter(!is.na(RelativeCover))
#Transformar porcentage en relativo
grime_relcover$RelativeCover = grime_relcover$RelativeCover * 0.01
View(grime_relcover)

grimeplot = ggplot(grime_relcover, aes(x = rip_modif, y = RelativeCover, fill = Strategy)) +
  geom_bar(stat = "identity", position = "fill", color = "transparent") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "R" = "#fdda24",
    "CR" = "#89a04d",
    "C" = "#3c58a0",
    "CS" = "#a656a6",
    "RS" = "#e41a1c",
    "S" = "#8b4513"
  )) +
  labs(x = "", y = "Species Cover", fill = "Ecological Strategy 
System of Grime") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", size=14),
    axis.title = element_text(color = "black", size=14),
    legend.text = element_text(color = "black", size=10),
    legend.title = element_text(color = "black", size=10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )
print(grimeplot)

#otro plot de boxplots
grimeplot2=ggplot(data = grime_relcover, aes(x = rip_modif, 
                                             y = RelativeCover, 
                                             color = Strategy)) + 
  geom_boxplot(outlier.shape = NA, size = 1) + 
  facet_grid() + 
  theme_bw() +
  xlab("Riparian modification") +
  ylab("Species Cover (%)") +
  scale_color_manual(values = c("#3c58a0", "#89a04d", "#a656a6", "#fdda24", "#e41a1c", "#8b4513")) +  
  scale_fill_manual(values = c("#3c58a0", "#89a04d", "#a656a6", "#fdda24", "#e41a1c", "#8b4513")) +  
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 9, face = "bold")) +
  labs(color = "Grime CSR 
       Strategy", fill = "")
print(grimeplot2)

#análisis estadistico GLM con variable respuesta proporciones derivada de conteos
#Distribución beta: variable respuesta entre 0 y 1
grime_relcover$Strategy=as.factor(grime_relcover$Strategy)
summary(grime_relcover)

#Modelo
grime_wide = grime_relcover %>%
  tidyr::pivot_wider(names_from = Strategy, values_from = RelativeCover, values_fill = 0)
View(grime_wide)

# Seleccionar solo las columnas de estrategias
Y = as.matrix(grime_wide[, c("C", "CR", "CS", "R", "RS", "S")])

# Crear objeto Dirichlet
Y_dirichlet = DR_data(Y)
View(Y_dirichlet)

# Agregar columna de rip_modif al objeto de modelado
grime_wide$rip_modif <- as.factor(grime_wide$rip_modif)

# Modelo: proporciones ~ grado de modificación ribereña
m_grime = DirichReg(Y_dirichlet ~ rip_modif, data = grime_wide)
summary(m_grime)
#Ejemplo S. Cambio no significativo entre high a medimum y cambio significativo entre high y null

########################## Dispersal Syndrome ##################################
###Calcular proporciones
data$Disp=as.factor(data$Disp)
summary(data$Disp)

Anemocoria_strat_cover = data %>%
  filter(Disp == "Anemochory") %>%
  group_by(plot) %>%
  summarise(Anemochory = sum(StrataCover)) %>%
  ungroup()
View(Anemocoria_strat_cover)

Autocoria_strat_cover = data %>%
  filter(Disp == "Autochory") %>%
  group_by(plot) %>%
  summarise(Autochory = sum(StrataCover)) %>%
  ungroup()
View(Autocoria_strat_cover)

Endozoocoria_strat_cover = data %>%
  filter(Disp == "Endozoochory") %>%
  group_by(plot) %>%
  summarise(Endozoochory = sum(StrataCover)) %>%
  ungroup()
View(Endozoocoria_strat_cover)

Epizoocoria_strat_cover = data %>%
  filter(Disp == "Epizoochory") %>%
  group_by(plot) %>%
  summarise(Epizoochory = sum(StrataCover)) %>%
  ungroup()
View(Epizoocoria_strat_cover)

Hidrocoria_strat_cover = data %>%
  filter(Disp == "Hidrochory") %>%
  group_by(plot) %>%
  summarise(Hidrochory = sum(StrataCover)) %>%
  ungroup()
View(Hidrocoria_strat_cover)

#Unir tablas
df_list = list(Anemocoria_strat_cover, Autocoria_strat_cover, Hidrocoria_strat_cover, Endozoocoria_strat_cover, Epizoocoria_strat_cover)
disp_relcover = Reduce(function(x, y) merge(x, y, by = "plot", all = TRUE), df_list)
#agregar columna rip_modif
disp_relcover = left_join(disp_relcover, modelo %>% select(plot, rip_modif), by = "plot")
View(disp_relcover)

write_xlsx(disp_relcover, "disp.xlsx")

#Cambiar formato de tabla para mejor analisis
disp_relcover = disp_relcover %>%
  pivot_longer(cols = c(Anemochory, Autochory, Hidrochory, Endozoochory, Epizoochory),
               names_to = "Strategy",
               values_to = "RelativeCover") %>%
  filter(!is.na(RelativeCover))
#Transformar porcentage en relativo
disp_relcover$RelativeCover = disp_relcover$RelativeCover * 0.01
View(disp_relcover)
summary(disp_relcover)

ggplot(disp_relcover, aes(x = rip_modif, y = RelativeCover, fill = Strategy)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(x = "modif", y = "relcover") +
  theme_minimal()

dispplot = ggplot(disp_relcover, aes(x = rip_modif, y = RelativeCover, fill = Strategy)) +
  geom_bar(stat = "identity", position = "fill", color = "transparent") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "Anemochory" = "#fdda24",
    "Autochory" = "#89a04d",
    "Hidrochory" = "#3c58a0",
    "Endozoochory" = "#a656a6",
    "Epizoochory" = "#e41a1c",
    "NA" = "#8b4513"
  )) +
  labs(x = "", y = "Species Cover", fill = "Dispersal Strategy") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", size=14),
    axis.title = element_text(color = "black", size=14),
    legend.text = element_text(color = "black", size=10),
    legend.title = element_text(color = "black", size=10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )
print(dispplot)

#análisis estadistico GLM con variable respuesta proporciones derivada de conteos
#Distribución beta: variable respuesta entre 0 y 1
disp_relcover$Strategy=as.factor(disp_relcover$Strategy)
summary(disp_relcover)

#Modelo
disp_wide = disp_relcover %>%
  tidyr::pivot_wider(names_from = Strategy, values_from = RelativeCover, values_fill = 0)
View(disp_wide)
summary(disp_wide)

# Seleccionar solo las columnas de estrategias
Y = as.matrix(disp_wide[, c("Anemochory", "Autochory", "Hidrochory", "Endozoochory", "Epizoochory")])

# Crear objeto Dirichlet
Y_dirichlet = DR_data(Y)
View(Y_dirichlet)

# Agregar columna de rip_modif al objeto de modelado
disp_wide$rip_modif <- as.factor(disp_wide$rip_modif)

# Modelo: proporciones ~ grado de modificación ribereña
m_disp = DirichReg(Y_dirichlet ~ rip_modif, data = disp_wide)
summary(m_disp)

###Pollination Syndrome###
data$Poll=as.factor(data$Poll)
summary(data$Poll)

Anemogamia_strat_cover = data %>%
  filter(Poll == "Anemophily") %>%
  group_by(plot) %>%
  summarise(Anemophily = sum(StrataCover)) %>%
  ungroup()
View(Anemogamia_strat_cover)

Zoogamia_strat_cover = data %>%
  filter(Poll == "Zoophily") %>%
  group_by(plot) %>%
  summarise(Zoophily = sum(StrataCover)) %>%
  ungroup()
View(Zoogamia_strat_cover)

Autogamia_strat_cover = data %>%
  filter(Poll == "Autogamy") %>%
  group_by(plot) %>%
  summarise(Autogamy = sum(StrataCover)) %>%
  ungroup()
View(Autogamia_strat_cover)

#Unir tablas
df_list = list(Anemogamia_strat_cover, Zoogamia_strat_cover, Autogamia_strat_cover)
poll_relcover = Reduce(function(x, y) merge(x, y, by = "plot", all = TRUE), df_list)
#agregar columna rip_modif
poll_relcover = left_join(poll_relcover, modelo %>% select(plot, rip_modif), by = "plot")
View(poll_relcover)

write_xlsx(poll_relcover, "poll.xlsx")

#Cambiar formato de tabla para mejor analisis
poll_relcover = poll_relcover %>%
  pivot_longer(cols = c(Anemophily, Zoophily, Autogamy),
               names_to = "Strategy",
               values_to = "RelativeCover") %>%
  filter(!is.na(RelativeCover))
#Transformar porcentage en relativo
poll_relcover$RelativeCover = poll_relcover$RelativeCover * 0.01
View(poll_relcover)

pollplot = ggplot(poll_relcover, aes(x = rip_modif, y = RelativeCover, fill = Strategy)) +
  geom_bar(stat = "identity", position = "fill", color = "transparent") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "Anemophily" = "#fdda24",
    "Autogamy" = "#89a04d",
    "NA" = "#3c58a0",
    "Zoophily" = "#a656a6",
    "NA" = "#e41a1c",
    "NA" = "#8b4513"
  )) +
  labs(x = "", y = "Species Cover", fill = "Pollination Strategy") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", size=14),
    axis.title = element_text(color = "black", size=14),
    legend.text = element_text(color = "black", size=10),
    legend.title = element_text(color = "black", size=10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )
print(pollplot)

#análisis estadistico GLM con variable respuesta proporciones derivada de conteos
#Distribución beta: variable respuesta entre 0 y 1
poll_relcover$Strategy=as.factor(poll_relcover$Strategy)
summary(poll_relcover)

#Modelo
poll_wide = poll_relcover %>%
  tidyr::pivot_wider(names_from = Strategy, values_from = RelativeCover, values_fill = 0)
View(poll_wide)

# Seleccionar solo las columnas de estrategias
Y = as.matrix(poll_wide[, c("Anemophily", "Zoophily", "Autogamy")])

# Crear objeto Dirichlet
Y_dirichlet = DR_data(Y)
View(Y_dirichlet)

# Agregar columna de rip_modif al objeto de modelado
poll_wide$rip_modif <- as.factor(poll_wide$rip_modif)

# Modelo: proporciones ~ grado de modificación ribereña
m_poll = DirichReg(Y_dirichlet ~ rip_modif, data = poll_wide)
summary(m_poll)

###Unir plots rasgos categoricos
#Combinación de los 3 boxplots de estructura funcional
combined_plot = plot_grid(grimeplot, dispplot, pollplot, 
                          ncol = 3, align = "v", axis = "tb")

plotcategorical = plot_grid(
  combined_plot,
  ggdraw() + draw_label("", fontface = 'bold', size = 10),
  ncol = 1,
  rel_heights = c(1, 0.09) # adjust this to control spacing
)
print(plotcategorical)

### Rasgos Numericos

### Hmax: Distribucón Normal
summary(modelo[,c("CWM.Hmax", "rip_modif")])
#modelo
CWM.HmaxRip=lm(CWM.Hmax~rip_modif, data=modelo)
boxplot(CWM.Hmax~rip_modif, data=modelo) 
summary(CWM.HmaxRip)
#Intercept: grupo de referencia (rip_modif:high)
#estimate: diferenecia entre la media con el grupo de referencia
anova(CWM.HmaxRip)
#F<0.005→ diferencia singnificativa entre grupos
#test a posteriori (p<0.05)
CWM.HmaxRip.post=glht(model=CWM.HmaxRip, linfct=mcp(rip_modif="Tukey"))
summary(CWM.HmaxRip.post, test=adjusted("fdr"))

#Analisis de residuos
res.CWM.HmaxRip=augment(CWM.HmaxRip)
CWM.HmaxRip.res=ggplot(data=res.CWM.HmaxRip, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

CWM.HmaxRip.res.vs.expl=ggplot(data=res.CWM.HmaxRip, aes(x=rip_modif, y=.resid))+
  geom_violin(alpha=0.9)+
  geom_jitter (size=2, width=0.1)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="rip_modif", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(),
        axis.text = element_text(size=16))

CWM.HmaxRip.Cook=ggplot(data=res.CWM.HmaxRip, aes(x=1:nrow(res.CWM.HmaxRip),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

CWM.HmaxRip.qq=ggplot(res.CWM.HmaxRip, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(CWM.HmaxRip.res, CWM.HmaxRip.res.vs.expl, CWM.HmaxRip.qq, CWM.HmaxRip.Cook, ncol=2) 

#Residuales ajustan bastante bien

#Presentación de resultados
#CI
pred.CWM.HmaxRip=ggpredict(CWM.HmaxRip, terms = c("rip_modif"))
pred.CWM.HmaxRip

#letras para diferencia de medias
em = emmeans(CWM.HmaxRip, ~ rip_modif)
em
letras = cld(em, Letters = letters, adjust = "Tukey") 

letras = letras %>%
  mutate(rip_modif = as.factor(rip_modif)) %>%
  select(rip_modif, .group)
medias = modelo %>%
  group_by(rip_modif) %>%
  summarise(CWM.Hmax = max(CWM.Hmax)) %>%
  left_join(letras, by = "rip_modif")

# Boxplot 
boxplotCWM.Hmax=ggplot(modelo, aes(x = rip_modif, y = CWM.Hmax, fill = rip_modif)) +
  geom_boxplot(color = "black", outlier.shape = NA) +  
  geom_jitter(width = 0.1, size = 1, shape = 21, color = "black", fill = "black", stroke = 0.3) +  
  scale_fill_grey() +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(color = "black", size=14),
    axis.title = element_text(color = "black", size=14)
  ) +
  labs(x = "", y = "Hmax (m)") +
  geom_text(data = medias, aes(x = rip_modif, y = CWM.Hmax + 0.001, label = .group),
            color = "black", size = 6)
print(boxplotCWM.Hmax)

###SLA: Distribución binomial negativa

library(MASS)    
library(emmeans) 
library(DHARMa)  
library(ggplot2)
library(dplyr)
library(ggeffects)
library(cowplot)

summary(modelo[, c("CWM.SLA", "rip_modif")])
#modelo
CWM.SLA.nb <- glm.nb(CWM.SLA ~rip_modif, data = modelo)
summary(CWM.SLA.nb)
#Hay mucha sobredispersión
anova(CWM.SLA.nb, test = "Chisq")

CWM.SLA.nb.emm <- emmeans(CWM.SLA.nb, ~ rip_modif)
pairs(CWM.SLA.nb.emm, adjust = "tukey")

# Letters for groups
letras <- cld(CWM.SLA.nb.emm, Letters = letters, adjust = "tukey") 
letras <- letras %>%
  mutate(rip_modif = as.factor(rip_modif)) %>%
  select(rip_modif, .group)

#Model residual diagnostics
simres <- simulateResiduals(fittedModel = CWM.SLA.nb)
plot(simres)  # general plots

res.CWM.SLA <- augment(CWM.SLA.nb)  

CWM.SLA.nb.res <- ggplot(data = res.CWM.SLA, aes(x = .fitted, y = .resid)) +
  geom_point(size = 2) +
  theme_bw() +
  geom_hline(yintercept = 0) + 
  labs(x = "Fitted values", y = "Standardized residuals") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

CWM.SLA.nb.res.vs.expl <- ggplot(data = res.CWM.SLA, aes(x = rip_modif, y = .resid)) +
  geom_violin(alpha = 0.9) +
  geom_jitter(size = 2, width = 0.1) +
  theme_bw() +
  geom_hline(yintercept = 0) + 
  labs(x = "Riparian modification", y = "Standardized residuals") +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 16))

CWM.SLA.nb.Cook <- ggplot(data = res.CWM.SLA, aes(x = 1:nrow(res.CWM.SLA), y = .cooksd)) +
  geom_linerange(aes(ymin = 0, ymax = .cooksd)) +
  theme_bw() +
  labs(x = "Data point", y = "Cook's distance") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

CWM.SLA.nb.qq <- ggplot(res.CWM.SLA, aes(sample = .std.resid)) +
  stat_qq() + stat_qq_line() +
  theme_bw() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

plot_grid(CWM.SLA.nb.res, CWM.SLA.nb.res.vs.expl, CWM.SLA.nb.qq, CWM.SLA.nb.Cook, ncol = 2)

#Resultados
pred.CWM.SLA.nb <- ggpredict(CWM.SLA.nb, terms = c("rip_modif"))
pred.CWM.SLA.nb
medias <- modelo %>%
  group_by(rip_modif) %>%
  summarise(CWM.SLA = max(CWM.SLA)) %>%
  left_join(letras, by = "rip_modif")

boxplotCWM.SLA <- ggplot(modelo, aes(x = rip_modif, y = CWM.SLA, fill = rip_modif)) +
  geom_boxplot(color = "black", outlier.shape = NA) +  
  geom_jitter(width = 0.1, size = 1, shape = 21, color = "black", fill = "black", stroke = 0.3) +  
  scale_fill_grey() +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(color = "black", size=14),
    axis.title = element_text(color = "black", size = 14)
  ) +
  labs(x = "", y = "SLA (mg/mm²)") +
  geom_text(data = medias, aes(x = rip_modif, y = CWM.SLA + 1, label = .group),
            color = "black", size = 6)
print(boxplotCWM.SLA)

###SDM: Distribución bionmial negativa
summary(modelo[, c("CWM.SDM", "rip_modif")])
#modelo
CWM.SDM.nb <- glm.nb(CWM.SDM ~rip_modif, data = modelo)
summary(CWM.SDM.nb)
anova(CWM.SDM.nb, test = "Chisq")

CWM.SDM.nb.emm <- emmeans(CWM.SDM.nb, ~ rip_modif)
pairs(CWM.SDM.nb.emm, adjust = "tukey")

# Letters for groups
letras <- cld(CWM.SDM.nb.emm, Letters = letters, adjust = "tukey") 
letras <- letras %>%
  mutate(rip_modif = as.factor(rip_modif)) %>%
  select(rip_modif, .group)

simres <- simulateResiduals(fittedModel = CWM.SDM.nb)
plot(simres)

res.CWM.SDM <- augment(CWM.SDM.nb)  

CWM.SDM.nb.res <- ggplot(data = res.CWM.SDM, aes(x = .fitted, y = .resid)) +
  geom_point(size = 2) +
  theme_bw() +
  geom_hline(yintercept = 0) + 
  labs(x = "Fitted values", y = "Standardized residuals") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

CWM.SDM.nb.res.vs.expl <- ggplot(data = res.CWM.SDM, aes(x = rip_modif, y = .resid)) +
  geom_violin(alpha = 0.9) +
  geom_jitter(size = 2, width = 0.1) +
  theme_bw() +
  geom_hline(yintercept = 0) + 
  labs(x = "Riparian modification", y = "Standardized residuals") +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 16))

CWM.SDM.nb.Cook <- ggplot(data = res.CWM.SDM, aes(x = 1:nrow(res.CWM.SDM), y = .cooksd)) +
  geom_linerange(aes(ymin = 0, ymax = .cooksd)) +
  theme_bw() +
  labs(x = "Data point", y = "Cook's distance") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

CWM.SDM.nb.qq <- ggplot(res.CWM.SDM, aes(sample = .std.resid)) +
  stat_qq() + stat_qq_line() +
  theme_bw() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

plot_grid(CWM.SDM.nb.res, CWM.SDM.nb.res.vs.expl, CWM.SDM.nb.qq, CWM.SDM.nb.Cook, ncol = 2)

#Resultados
pred.CWM.SDM.nb <- ggpredict(CWM.SDM.nb, terms = c("rip_modif"))
pred.CWM.SDM.nb
medias <- modelo %>%
  group_by(rip_modif) %>%
  summarise(CWM.SDM = max(CWM.SDM)) %>%
  left_join(letras, by = "rip_modif")

boxplotCWM.SDM <- ggplot(modelo, aes(x = rip_modif, y = CWM.SDM, fill = rip_modif)) +
  geom_boxplot(color = "black", outlier.shape = NA) +  
  geom_jitter(width = 0.1, size = 1, shape = 21, color = "black", fill = "black", stroke = 0.3) +  
  scale_fill_grey() +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(color = "black", size=14),
    axis.title = element_text(color = "black", size=14)
  ) +
  labs(x = "", y = "SDM (mg)") +
  geom_text(data = medias, aes(x = rip_modif, y = CWM.SDM + 1, label = .group),
            color = "black", size = 6)
print(boxplotCWM.SDM)

#Combinación de los 3 boxplots de estructura funcional
library(cowplot)
combined_plots <- plot_grid(boxplotCWM.SLA, boxplotCWM.SDM, boxplotCWM.Hmax,
                            ncol = 3, align = "v")
x_label <- ggdraw() + 
  draw_label("", fontface = 'bold', size = 12, x = 0.5, hjust = 0.5)
boxplotstructure <- plot_grid(combined_plots, x_label, ncol = 1, rel_heights = c(1, 0.07))
print(boxplotstructure)


### FILTRO: FRAGMENTACIÓN DEL HÁBITAT

#Variable explicativa: Índice de división del paisaje (LDI)

mCWM.Hmax=lm(CWM.Hmax ~ fragmentation, data=modelo)   
summary(mCWM.Hmax)
anova(mCWM.Hmax, test="F")
#significativo

#Analisis de residuos
res.mCWM.Hmax=augment(mCWM.Hmax)

mCWM.Hmax.res=ggplot(data=res.mCWM.Hmax, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+ theme_bw()+ geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mCWM.Hmax.res.vs.expl=ggplot(data=res.mCWM.Hmax, aes(x=fragmentation, y=.resid))+
  geom_point(size=2)+theme_bw()+geom_hline(yintercept = 0)+ 
  labs(x="fragmentation", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(), axis.text = element_text(size=16))

mCWM.Hmax.qq=ggplot(res.mCWM.Hmax, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mCWM.Hmax.Cook=ggplot(data=res.mCWM.Hmax, aes(x=1:nrow(res.mCWM.Hmax),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(mCWM.Hmax.res, mCWM.Hmax.res.vs.expl, mCWM.Hmax.qq, mCWM.Hmax.Cook, ncol=2)
#Ajustan mas o menos bien

#Gráfico
pred.CWM.Hmax <- ggpredict(mCWM.Hmax, terms = c("fragmentation"))
plotHmax_frag <- plot(pred.CWM.Hmax) +
  labs(x = "Fragmentation (LDI)", y = "Hmax (m)") +
  geom_line(size = 1.2) +
  theme_minimal(base_family = "") + 
  scale_y_continuous(breaks = seq(1, 2.5, 0.3)) +
  theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    plot.title = element_blank(),
    panel.grid = element_blank(),          # elimina grilla
    panel.border = element_rect(colour = "black", fill = NA, size = 1), # marco negro
    axis.line = element_line(color = "black")
  )
print(plotHmax_frag)

###SLA##########################################################################
m.nb.SLA=glm.nb(CWM.SLA ~ fragmentation, data = modelo)
anova(m.nb.SLA, test="Chi")
summary(m.nb.SLA)

#R² McFadden
1-(m.nb.SLA$deviance/m.nb.SLA$null.deviance) 
100-((28.572*100)/119.199)

#Residuos
# Analisis de residuos de m.nb
resid.m.nb.SLA=augment(m.nb.SLA) 
res.m.nb.SLA=simulateResiduals(fittedModel=m.nb.SLA,n=1e3, integerResponse=T)
resid.m.nb.SLA$.std.resid=residuals(res.m.nb.SLA, quantileFunction = qnorm) # convierte RQR a dist Normal
# qqplot
qq.m.nb.SLA=ggplot(data=resid.m.nb.SLA, mapping=aes(sample = .std.resid)) + 
  theme_bw()+
  stat_qq()+
  theme_bw()+
  stat_qq_line()+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
# Residuales vs fitted
res.fit.m.nb.SLA= ggplot(data=resid.m.nb.SLA, aes(x=.fitted,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuales")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
# Residuales vs variable explicativa 
res.FRAG.m.nb=ggplot(data=resid.m.nb.SLA, aes(x=fragmentation,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="FRAG",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
Cook.m.nb.SLA=ggplot(data=resid.m.nb.SLA, aes(x=1:nrow(resid.m.nb.SLA),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Distancia Cook")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
plot_grid(res.fit.m.nb.SLA,qq.m.nb.SLA,res.FRAG.m.nb,Cook.m.nb.SLA, ncol=2)
#Los residuales NO ajustan MUY bien

#Curva predicha
pred.m.nb.SLA <- ggpredict(m.nb.SLA, terms = c("fragmentation"))

# Plot
plotSLA_frag <- plot(pred.m.nb.SLA) + 
  theme_bw() +  
  labs(x = "Fragmentation (LDI)", y = "SLA (mm²/mg)") +
  geom_line(size=1.2)+
  scale_y_continuous(breaks = seq(10, 24, 4)) +
  theme(
    panel.grid = element_blank(),           # Remove grid
    axis.title = element_text(size = 15, color = "black"), # Black axis titles
    axis.text = element_text(size = 15, color = "black"),  # Black axis numbers
    plot.title = element_blank(),
    legend.position = c(0.85, 0.85),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Black frame
  ) +
  geom_line(linewidth = 1.2)  # Make the regression line thicker

print(plotSLA_frag)

###SDM##########################################################################
m.nb.SDM=glm.nb(CWM.SDM ~ fragmentation, data = modelo)
anova(m.nb.SDM, test="Chi")
summary(m.nb.SDM)

#R² McFadden
1-(m.nb.SDM$deviance/m.nb.SDM$null.deviance) #o
100-((38.296*100)/84.587)

#Residuos
# Analisis de residuos de m.nb
resid.m.nb.SDM=augment(m.nb.SDM) 
res.m.nb.SDM=simulateResiduals(fittedModel=m.nb.SDM,n=1e3, integerResponse=T)
resid.m.nb.SDM$.std.resid=residuals(res.m.nb.SDM, quantileFunction = qnorm) # convierte RQR a dist Normal
# qqplot
qq.m.nb.SDM=ggplot(data=resid.m.nb.SDM, mapping=aes(sample = .std.resid)) + 
  theme_bw()+
  stat_qq()+
  theme_bw()+
  stat_qq_line()+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
# Residuales vs fitted
res.fit.m.nb.SDM= ggplot(data=resid.m.nb.SDM, aes(x=.fitted,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuales")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
# Residuales vs variable explicativa 
res.FRAG.m.nb=ggplot(data=resid.m.nb.SDM, aes(x=fragmentation,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="FRAG",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
Cook.m.nb.SDM=ggplot(data=resid.m.nb.SDM, aes(x=1:nrow(resid.m.nb.SDM),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Distancia Cook")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
plot_grid(res.fit.m.nb.SDM,qq.m.nb.SDM,res.FRAG.m.nb,Cook.m.nb.SDM, ncol=2)
#Los residuales NO ajustan MUY bien

#Curva predicha
pred.m.nb.SDM <- ggpredict(m.nb.SDM, terms = c("fragmentation"))

# Plot
plotSDM_frag <- plot(pred.m.nb.SDM) + 
  theme_bw() +  
  labs(x = "Fragmentation (LDI)", y = "SDM (mg)") +
  geom_line(size=1.2)+
  scale_y_continuous(breaks = seq(0, 12, 3)) +
  theme(
    panel.grid = element_blank(),           # Remove grid
    axis.title = element_text(size = 15, color = "black"), # Black axis titles
    axis.text = element_text(size = 15, color = "black"),  # Black axis numbers
    plot.title = element_blank(),
    legend.position = c(0.85, 0.85),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Black frame
  ) +
  geom_line(linewidth = 1.2)  # Make the regression line thicker

print(plotSDM_frag)

#Union de gráficos
plotSLA_frag_clean <- plotSLA_frag + theme(axis.title.x = element_blank())
plotSDM_frag_clean <- plotSDM_frag + theme(axis.title.x = element_blank())
plotHmax_frag_clean <- plotHmax_frag + theme(axis.title.x = element_blank())

combined_plots <- plot_grid(plotSLA_frag_clean, plotSDM_frag_clean,
                            plotHmax_frag_clean,
                            ncol = 3, align = "v")
x_label <- ggdraw() + 
  draw_label("", fontface = 'bold', size = 12, x = 0.5, hjust = 0.5)
final_plotfragstruct <- plot_grid(combined_plots, x_label, ncol = 1, rel_heights = c(1, 0.07))
print(final_plotfragstruct)

#Figura 2

plot_grid(plotcategorical, boxplotstructure, final_plotfragstruct, nrow=3)

############################# Human interventions ###########################

#Variables explicativaS: cobertura de especies ornamentales y exóticas

###########################Exotic and Ornamental Proportions####################
data$Origin=as.factor(data$Origin)

exotic_strat_cover = data %>%
  filter(Origin == "Exotic") %>%
  group_by(plot) %>%
  summarise(Exotic = sum(StrataCover)) %>%
  ungroup()
View(exotic_strat_cover)

Native_strat_cover = data %>%
  filter(Origin == "Native") %>%
  group_by(plot) %>%
  summarise(Native = sum(StrataCover)) %>%
  ungroup()
View(Native_strat_cover)

#Unir tablas
df_list = list(exotic_strat_cover, Native_strat_cover)
origin_relcover = Reduce(function(x, y) merge(x, y, by = "plot", all = TRUE), df_list)
View(origin_relcover)
#agregar columna rip_modif
origin_relcover = left_join(origin_relcover, modelo %>% select(plot, rip_modif), by = "plot")
#Cambiar formato de tabla para mejor analisis
origin_relcover = origin_relcover %>%
  pivot_longer(cols = c(Native, Exotic),
               names_to = "Origin",
               values_to = "RelativeCover") %>%
  filter(!is.na(RelativeCover))
#Transformar porcentage en relativo
origin_relcover$RelativeCover = origin_relcover$RelativeCover * 0.01

origineplot1 = ggplot(origin_relcover, aes(x = rip_modif, y = RelativeCover, fill = Origin)) +
  geom_bar(stat = "identity", position = "fill", color = "transparent") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "Native" = "#fdda24",
    "Exotic" = "#89a04d",
    "C" = "#3c58a0",
    "CS" = "#a656a6",
    "RS" = "#e41a1c",
    "S" = "#8b4513"
  )) +
  labs(x = "", y = "Species Cover", fill = "Origin") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", size=12),
    axis.title = element_text(color = "black", size=15),
    legend.text = element_text(color = "black", size=12),
    legend.title = element_text(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )
print(origineplot1)

##ornamental

data$Ornamental=as.factor(data$Ornamental)
View(data)
Ornamental_strat_cover = data %>%
  filter(Ornamental == "Ornamental") %>%
  group_by(plot) %>%
  summarise(Ornamental = sum(StrataCover)) %>%
  ungroup()
View(ornamental_strat_cover)

NoOrnamental_strat_cover = data %>%
  filter(Ornamental == "No Ornamental") %>%
  group_by(plot) %>%
  summarise(No_ornamental = sum(StrataCover)) %>%
  ungroup()
View(NoOrnamental_strat_cover)

#Unir tablas
df_list = list(Ornamental_strat_cover, NoOrnamental_strat_cover)
Ornamental_relcover = Reduce(function(x, y) merge(x, y, by = "plot", all = TRUE), df_list)
View(Ornamental_relcover)
#agregar columna rip_modif
Ornamental_relcover = left_join(Ornamental_relcover, modelo %>% select(plot, rip_modif), by = "plot")
#Cambiar formato de tabla para mejor analisis
Ornamental_relcover = Ornamental_relcover %>%
  pivot_longer(cols = c(Ornamental, No_ornamental),
               names_to = "Ornamental",
               values_to = "RelativeCover") %>%
  filter(!is.na(RelativeCover))
#Transformar porcentage en relativo
Ornamental_relcover$RelativeCover = Ornamental_relcover$RelativeCover * 0.01

Ornamentaleplot1 = ggplot(Ornamental_relcover, aes(x = rip_modif, y = RelativeCover, fill = Ornamental)) +
  geom_bar(stat = "identity", position = "fill", color = "transparent") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "No Ornamental" = "#fdda24",
    "Ornamental" = "#89a04d",
    "C" = "#3c58a0",
    "CS" = "#a656a6",
    "RS" = "#e41a1c",
    "S" = "#8b4513"
  )) +
  labs(x = "", y = "Species Cover", fill = "Ornamental") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", size=12),
    axis.title = element_text(color = "black", size=15),
    legend.text = element_text(color = "black", size=12),
    legend.title = element_text(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )
print(Ornamentaleplot1)

#Modelos lineales multiples → Distribuciones normales

###HMax#########################################################################
#Modelo
modelo_interv=lm(CWM.Hmax~ornamental_abund+exoticsp_abund, data=modelo)
summary(modelo_interv)

#Intercept: Media de la variable respuesta cuando la abundancia de especies exoticas y ornamentales es 0
#El cambio de 1 unidad de abundancia de exoticas disminuye Hmax en 0.677 m

#Visualizacion de resultados
tidy(modelo_interv, conf.int=T)
ggcoef(modelo_interv, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))
#estimate: pendiente de las variables explicativas

#Analisis de residuales
res.modelo_interv=augment(modelo_interv)
res.modelo_interv

modelo_interv.res=ggplot(data=res.modelo_interv, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.ornamental=ggplot(data=res.modelo_interv, aes(x=ornamental_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="ornamental abund", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.exotic=ggplot(data=res.modelo_interv, aes(x=exoticsp_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="exoticsp_abund", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.Cook=ggplot(data=res.modelo_interv, aes(x=1:nrow(res.modelo_interv),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.qq=ggplot(res.modelo_interv, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_interv.res, modelo_interv.res.vs.ornamental, modelo_interv.res.vs.exotic, modelo_interv.qq,modelo_interv.Cook, ncol=3)

#Residuales se ajustan mas o menos

#Curvas condicionales (3 dimensiones)

pred.modelo_inerv=ggpredict(modelo_interv, terms = c("exoticsp_abund", "ornamental_abund"))
pred.modelo_inerv

# Curva condicional de predicciones 
plot(pred.modelo_inerv) +
  labs(x="exoticsp_abund", y="CWM.Hmax")+
  geom_point(data=modelo, aes(x=exoticsp_abund, y=CWM.Hmax), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))

#GLM → Distribuciones binomiales negativas (SLA y SDM)
summary(modelo)

###SLA##########################################################################
View(modelo)
m_interv_SLA=glm.nb(CWM.SLA~ornamental_abund+exoticsp_abund, data=modelo)
summary(m_interv_SLA)

#Intercepto: El valor de CWM.SLA cuando la abundancia de ornamentales y exoticas es 0
#Radio entre el Residual Deviance y los grados de libertad
42.282/40 #= 1.06 → no hay sobredispersión

#el cambio de una unidad de la cobertura de especies exoticas es mucho más fuerte en la SLA que las ornamentales

#pseudo R² McFadden
R2.m_interv_SLA=1-(m_interv_SLA$deviance/m_interv_SLA$null.deviance)
R2.m_interv_SLA #0.309
anova(m_interv_SLA, test="F")

#Visualización de resultados
tidy(m_interv_SLA, conf.int=T)
ggcoef(m_interv_SLA, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))

#Residuales
res.m_interv_SLA=augment(m_interv_SLA)
res.m_interv_SLA

#Generamos por simulación los residuales dun-smith
resid.m_interv_SLA=simulateResiduals(fittedModel=m_interv_SLA,n=1e3, integerResponse=T, refit=F,plot=F)

res.m_interv_SLA$.std.resid=residuals(resid.m_interv_SLA, quantileFunction = qnorm)

#QQ
qq.m_interv_SLA=ggplot(data=res.m_interv_SLA, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

# Residuales vs fitted
res.fit.m_interv_SLA= ggplot(data=res.m_interv_SLA, aes(x=.fitted,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuales")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

# Residuales vs variable explicativa 1
res.exotic.m_interv_SLA=ggplot(data=res.m_interv_SLA, aes(x=exoticsp_abund,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="exotic",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 2
res.ornamental.m_interv_SLA=ggplot(data=res.m_interv_SLA, aes(x=ornamental_abund,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="ornamental",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

Cook.m_interv_SLA=ggplot(data=res.m_interv_SLA, aes(x=1:nrow(res.m_interv_SLA),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Distancia Cook")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

plot_grid(res.fit.m_interv_SLA, qq.m_interv_SLA, res.exotic.m_interv_SLA, res.ornamental.m_interv_SLA, Cook.m_interv_SLA, ncol=2)

#Residuales mas o menos con distribución normal

pred=ggpredict(m_interv_SLA, terms=c("exoticsp_abund"))
# Visualizacion del modelo final 
plot(pred)+
  theme_bw()+
  labs(x="exoticsp_abund", y="SLA")+ 
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

###SDM##########################################################################
m_interv_SDM=glm.nb(CWM.SDM~ornamental_abund+exoticsp_abund, data=modelo)
summary(m_interv_SDM)

#Intercepto: El valor de CWM.SDM cuando la abundancia de ornamentales y exoticas es 0
#Radio entre el Residual Deviance y los grados de libertad
42.094/40 #= 1.05 → no hay sobredispersión

#el cambio de una unidad de la cobertura de especies exoticas es más fuerte en la SDM que las ornamentales, sin embargo, ambas influencian la SDM

#pseudo R² McFadden
R2.m_interv_SDM=1-(m_interv_SDM$deviance/m_interv_SDM$null.deviance)
R2.m_interv_SDM #0.358
anova(m_interv_SDM, test="F")

#Visualización de resultados
tidy(m_interv_SDM, conf.int=T)
ggcoef(m_interv_SDM, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))

#Residuales
res.m_interv_SDM=augment(m_interv_SDM)
res.m_interv_SDM

#Generamos por simulación los residuales dun-smith
resid.m_interv_SDM=simulateResiduals(fittedModel=m_interv_SDM,n=1e3, integerResponse=T, refit=F,plot=F)

res.m_interv_SDM$.std.resid=residuals(resid.m_interv_SDM, quantileFunction = qnorm)

#QQ
qq.m_interv_SDM=ggplot(data=res.m_interv_SDM, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

# Residuales vs fitted
res.fit.m_interv_SDM= ggplot(data=res.m_interv_SDM, aes(x=.fitted,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuales")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

# Residuales vs variable explicativa 1
res.exotic.m_interv_SDM=ggplot(data=res.m_interv_SDM, aes(x=exoticsp_abund,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="exotic",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 2
res.ornamental.m_interv_SDM=ggplot(data=res.m_interv_SDM, aes(x=ornamental_abund,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="ornamental",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

Cook.m_interv_SDM=ggplot(data=res.m_interv_SDM, aes(x=1:nrow(res.m_interv_SDM),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Distancia Cook")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

plot_grid(res.fit.m_interv_SDM, qq.m_interv_SDM, res.exotic.m_interv_SDM, res.ornamental.m_interv_SDM, Cook.m_interv_SDM, ncol=2)

#Residuales mas o menos con distribución normal

pred=ggpredict(m_interv_SDM, terms=c("exoticsp_abund"))
# Visualizacion del modelo final 
plot(pred)+
  theme_bw()+
  labs(x="exoticsp_abund", y="SDM")+ 
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))


######################### URBAN LOCAL ENVIRONMENT ##############################
##CARGA DE DATOS
modelo = read.csv("ModeloFiltros_Diversidad.csv", header=T,  sep = ",", stringsAsFactors = T)
View(modelo)
#Variables explicativaS: N, P, pH, CEA, COS, mean_ST, shading

#Modelos lineales multiples → Distribuciones normales
summary(modelo)

#Logaritmizar variables de suelo para reducir el orden de magnitud (skewness)
modelo$N_log = log(modelo$N)
modelo$P_log = log(modelo$P)
modelo$COS_log = log(modelo$COS)
modelo$CEA_log = log(modelo$CEA)

#Agregar la variable "shading" creando un índice proporcional a las horas de luz 
modelo$shading=1-(modelo$light_hours/max(modelo$light_hours))
View(modelo)

#eliminar plots 9 y 10 urbano (faltan los datos de suelos)
modelo = modelo %>%
  filter(!plot %in% c("U9", "U10"))

corr1 = modelo[, c("N_log", "P_log", "COS_log", "CEA_log", "pH", "mean_ST", "light_hours")]
my_cor <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  cor_test <- cor.test(x, y, use = "complete.obs")
  r <- cor_test$estimate  # Coeficiente de correlación
  p_value <- cor_test$p.value  # Valor p
  r2 <- round(r^2, 2)  # R²
  label <- paste0("R² = ", r2, "\np = ", format(p_value, digits = 3))
  ggplot(data, mapping) +
    annotate("text",
             x = mean(range(x, na.rm = TRUE)),
             y = mean(range(y, na.rm = TRUE)),
             label = label, size = 5) +
    theme_void()
}
ggpairs(corr1,
        upper = list(continuous = my_cor),
        lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.5)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)))

#Se elimina COS del modelo por alta relación con N y media relación con pH

#Analisis de variables explicativas
boxplot(pH~rip_modif, data = modelo)
boxplot(CEA~rip_modif, data=modelo)
colnames(modelo)

#Modelos lineales multiples → Distribuciones normales

###HMax#########################################################################
#Estandarizamos las variables explocativas
data_env.CWM.Hmax=as.data.frame(cbind(CWM.Hmax=modelo$CWM.Hmax, scale(modelo[,c("N_log", "P_log", "CEA_log", 
                                                                    "pH", "mean_ST", "shading")], center=T, scale=T)))
View(data_env.CWM.Hmax)

#Modelo
modelo_env=lm(CWM.Hmax~N_log+P_log+CEA_log+pH+mean_ST+shading,data=data_env.CWM.Hmax)
summary(modelo_env) 

#intercept=  en Hmax sitio con N, P, CEA, pH, ST, shading promedio
#Slope(estimate)= el grado de variación de Hmax cuando la variable explicativa aumenta una unidad

#Visualizacion de resultados
tidy(modelo_env, conf.int=T)
ggcoef(modelo_env, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))

#Analisis de residuales
res.modelo_env=augment(modelo_env)
res.modelo_env

modelo_env.res=ggplot(data=res.modelo_env, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.N=ggplot(data=res.modelo_env, aes(x=N_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="N", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_env.res.vs.P=ggplot(data=res.modelo_env, aes(x=P_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="P", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.pH=ggplot(data=res.modelo_env, aes(x=pH, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="pH", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.CEA=ggplot(data=res.modelo_env, aes(x=CEA_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="CEA", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.shade=ggplot(data=res.modelo_env, aes(x=light_hours, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="shade", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.temp=ggplot(data=res.modelo_env, aes(x=mean_ST, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="temp", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.Cook=ggplot(data=res.modelo_env, aes(x=1:nrow(res.modelo_env),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.qq=ggplot(res.modelo_env, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_env.res,modelo_env.res.vs.N,modelo_env.res.vs.P,modelo_env.res.vs.pH,modelo_env.res.vs.CEA,modelo_env.res.vs.shade,modelo_env.res.vs.temp,modelo_env.qq,modelo_env.Cook, ncol=3)

#Los residuales se ajustan decentemente

#Curvas condicionales (7 dimensiones)

pred.modelo_env=ggpredict(modelo_env, terms = c("mean_ST"))
pred.modelo_env

# Curva condicional de predicciones 
plot(pred.modelo_env) +
  labs(x="temp", y="CWM.Hmax")+
  geom_point(data=data_env.CWM.Hmax, aes(x=mean_ST, y=CWM.Hmax), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))


#GLM → Distribuciones binomiales negativas (SLA y SDM)


##SLA##########################################################################
#Estandarizamos las variables explicativas
data_env.CWM.SLA=as.data.frame(cbind(CWM.SLA=modelo$CWM.SLA, scale(modelo[,c("N_log", "P_log", "CEA_log", 
                                                                             "pH", "mean_ST", "shading")], center=T, scale=T)))
colnames(data_env.CWM.SLA)
View(data_env.CWM.SLA)

m_env_SLA=glm.nb(CWM.SLA~N_log+P_log+CEA_log+pH+mean_ST+shading, data=data_env.CWM.SLA)
summary(m_env_SLA)

#intercept =  es SLA sitio con N, P, CEA, pH, ST, light_hours promedio
#Slope(estimate)= el grado de variación de SLA cuando la variable explicativa aumenta una unidad

#ratio de resudal deviance y DF
40.888/34 #=1.2 aceptable no hay sobredispersión

#pseudo R² McFadden
R2.m_env_SLA=1-(m_env_SLA$deviance/m_env_SLA$null.deviance)
R2.m_env_SLA #0.609
anova(m_env_SLA, test="F")

#Visualización de resultados
tidy(m_env_SLA, conf.int=T)
ggcoef(m_env_SLA, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))

#Residuales
res.m_env_SLA=augment(m_env_SLA)
res.m_env_SLA

#Generamos por simulación los residuales dun-smith
resid.m_env_SLA=simulateResiduals(fittedModel=m_env_SLA,n=1e3, integerResponse=T, refit=F,plot=F)

res.m_env_SLA$.std.resid=residuals(resid.m_env_SLA, quantileFunction = qnorm)

#QQ
qq.m_env_SLA=ggplot(data=res.m_env_SLA, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

# Residuales vs fitted
res.fit.m_env_SLA= ggplot(data=res.m_env_SLA, aes(x=.fitted,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuales")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

# Residuales vs variable explicativa 1
res.N.m_env_SLA=ggplot(data=res.m_env_SLA, aes(x=N_log,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="N",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 2
res.P.m_env_SLA=ggplot(data=res.m_env_SLA, aes(x=P_log,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="P",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 3
res.CEA.m_env_SLA=ggplot(data=res.m_env_SLA, aes(x=CEA_log,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="CEA",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 4
res.pH.m_env_SLA=ggplot(data=res.m_env_SLA, aes(x=pH,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="pH",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 5
res.shade.m_env_SLA=ggplot(data=res.m_env_SLA, aes(x=light_hours,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="shade",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 6
res.temp.m_env_SLA=ggplot(data=res.m_env_SLA, aes(x=mean_ST,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="temp",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

Cook.m_env_SLA=ggplot(data=res.m_env_SLA, aes(x=1:nrow(res.m_env_SLA),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Distancia Cook")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

plot_grid(res.fit.m_env_SLA, qq.m_env_SLA, res.N.m_env_SLA, res.P.m_env_SLA, res.CEA.m_env_SLA, res.pH.m_env_SLA, res.shade.m_env_SLA, res.temp.m_env_SLA, Cook.m_env_SLA, ncol=3)

#Residuales mas o menos con distribución normal

pred=ggpredict(m_env_SLA, terms=c("CEA_log"))
# Visualizacion del modelo final 
plot(pred)+
  theme_bw()+
  labs(x="CEA", y="SLA")+ 
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

###SDM##########################################################################
data_env.CWM.SDM=as.data.frame(cbind(CWM.SDM=modelo$CWM.SDM, scale(modelo[,c("N_log", "P_log", "CEA_log", 
                                                                             "pH", "mean_ST", "shading")], center=T, scale=T)))
View(data_env.CWM.SDM)

m_env_SDM=glm.nb(CWM.SDM~N_log+P_log+CEA_log+pH+mean_ST+shading, data=data_env.CWM.SDM)
summary(m_env_SDM)

#intercept =  es SDM sitio con N, P, CEA, pH, ST, light_hours promedio
#Slope(estimate)= el grado de variación de SDM cuando la variable explicativa aumenta una unidad

#ratio de resudal deviance y DF
36.488/34 #=1.07 perfecto! no hay sobredispersión

#pseudo R² McFadden
R2.m_env_SDM=1-(m_env_SDM$deviance/m_env_SDM$null.deviance)
R2.m_env_SDM #0.601
anova(m_env_SDM, test="F")

#Visualización de resultados
tidy(m_env_SDM, conf.int=T)
ggcoef(m_env_SDM, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))

#Residuales
res.m_env_SDM=augment(m_env_SDM)
res.m_env_SDM

#Generamos por simulación los residuales dun-smith
resid.m_env_SDM=simulateResiduals(fittedModel=m_env_SDM,n=1e3, integerResponse=T, refit=F,plot=F)

res.m_env_SDM$.std.resid=residuals(resid.m_env_SDM, quantileFunction = qnorm)

#QQ
qq.m_env_SDM=ggplot(data=res.m_env_SDM, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

# Residuales vs fitted
res.fit.m_env_SDM= ggplot(data=res.m_env_SDM, aes(x=.fitted,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuales")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

# Residuales vs variable explicativa 1
res.N.m_env_SDM=ggplot(data=res.m_env_SDM, aes(x=N_log,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="N",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 2
res.P.m_env_SDM=ggplot(data=res.m_env_SDM, aes(x=P_log,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="P",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 3
res.CEA.m_env_SDM=ggplot(data=res.m_env_SDM, aes(x=CEA_log,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="CEA",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 4
res.pH.m_env_SDM=ggplot(data=res.m_env_SDM, aes(x=pH,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="pH",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 5
res.shade.m_env_SDM=ggplot(data=res.m_env_SDM, aes(x=light_hours,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="shade",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

# Residuales vs variable explicativa 6
res.temp.m_env_SDM=ggplot(data=res.m_env_SDM, aes(x=mean_ST,y=.std.resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="temp",y="Residuales")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18))

Cook.m_env_SDM=ggplot(data=res.m_env_SDM, aes(x=1:nrow(res.m_env_SDM),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Distancia Cook")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

plot_grid(res.fit.m_env_SDM, qq.m_env_SDM, res.N.m_env_SDM, res.P.m_env_SDM, res.CEA.m_env_SDM, res.pH.m_env_SDM, res.shade.m_env_SDM, res.temp.m_env_SDM, Cook.m_env_SDM, ncol=3)

#Residuales mas o menos con distribución normal

pred=ggpredict(m_env_SDM, terms=c("pH"))
# Visualizacion del modelo final 
plot(pred)+
  theme_bw()+
  labs(x="pH", y="SDM")+ 
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

####### OBJETIVO 2 - DIVERSIDAD FUNCIONAL Y FILTROS AMBIENTALES URBANOS ########

############### oNE-wAY ANALYSIS OF VARIANCE ###################################

##CARGA DE DATOS
modelo = read.csv("ModeloFiltros_Diversidad.csv", header=T,  sep = ",", stringsAsFactors = T)
View(modelo)

### FILTRO: TRANSFORMACIÓN DEL HÁBITAT
#Variable explicativa: Nivel de transformación del hábitat

###FRic####
summary(modelo[,c("FRic", "rip_modif")])
#modelo
FRicRip=lm(FRic~rip_modif, data=modelo)
boxplot(FRic~rip_modif, data=modelo) 
summary(FRicRip)
#Intercept: grupo de referencia (rip_modif:high)
#estimate: diferenecia entre la media con el grupo de referencia
anova(FRicRip)
#F<0.001→ diferencia singnificativa entre grupos
#test a posteriori (p<0.05)
FRicRip.post=glht(model=FRicRip, linfct=mcp(rip_modif="Tukey"))
summary(FRicRip.post, test=adjusted("fdr"))

#Analisis de residuos
res.FRicRip=augment(FRicRip)
FRicRip.res=ggplot(data=res.FRicRip, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

FRicRip.res.vs.expl=ggplot(data=res.FRicRip, aes(x=rip_modif, y=.resid))+
  geom_violin(alpha=0.9)+
  geom_jitter (size=2, width=0.1)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="rip_modif", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(),
        axis.text = element_text(size=16))

FRicRip.Cook=ggplot(data=res.FRicRip, aes(x=1:nrow(res.FRicRip),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

FRicRip.qq=ggplot(res.FRicRip, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(FRicRip.res, FRicRip.res.vs.expl, FRicRip.qq, FRicRip.Cook, ncol=2) 

#Residuales ajustan bastante bien

#Presentación de resultados
#CI
pred.FRicRip=ggpredict(FRicRip, terms = c("rip_modif"))
pred.FRicRip

#letras para diferencia de medias
em = emmeans(FRicRip, ~ rip_modif)
em
letras = cld(em, Letters = letters, adjust = "Tukey") 

letras = letras %>%
  mutate(rip_modif = as.factor(rip_modif)) %>%
  select(rip_modif, .group)
medias = modelo %>%
  group_by(rip_modif) %>%
  summarise(FRic = max(FRic)) %>%
  left_join(letras, by = "rip_modif")

# Boxplot 
boxplotFRic=ggplot(modelo, aes(x = rip_modif, y = FRic, fill = rip_modif)) +
  geom_boxplot(color = "black", outlier.shape = NA) +  
  geom_jitter(width = 0.1, size = 1, shape = 21, color = "black", fill = "black", stroke = 0.3) +  
  scale_fill_grey() +  
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15)
  ) +
  labs(x = "", y = "FRic") +
  geom_text(data = medias, aes(x = rip_modif, y = FRic + 0.08, label = .group),
            color = "black", size = 6)
print(boxplotFRic)

###FEve###
summary(modelo[,c("FEve", "rip_modif")])
#modelo
FEveRip=lm(FEve~rip_modif, data=modelo)
boxplot(FEve~rip_modif, data=modelo) 
summary(FEveRip)
#Intercept: grupo de referencia (rip_modif:high)
#estimate: diferenecia entre la media con el grupo de referencia
anova(FEveRip)
#F<0.005→ diferencia singnificativa entre grupos
#test a posteriori (p<0.05)
FEveRip.post=glht(model=FEveRip, linfct=mcp(rip_modif="Tukey"))
summary(FEveRip.post, test=adjusted("fdr"))

#Analisis de residuos
res.FEveRip=augment(FEveRip)
FEveRip.res=ggplot(data=res.FEveRip, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

FEveRip.res.vs.expl=ggplot(data=res.FEveRip, aes(x=rip_modif, y=.resid))+
  geom_violin(alpha=0.9)+
  geom_jitter (size=2, width=0.1)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="rip_modif", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(),
        axis.text = element_text(size=16))

FEveRip.Cook=ggplot(data=res.FEveRip, aes(x=1:nrow(res.FEveRip),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

FEveRip.qq=ggplot(res.FEveRip, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(FEveRip.res, FEveRip.res.vs.expl, FEveRip.qq, FEveRip.Cook, ncol=2) 

#Residuales ajustan bastante bien

#Presentación de resultados
#CI
pred.FEveRip=ggpredict(FEveRip, terms = c("rip_modif"))
pred.FEveRip

#letras para diferencia de medias
em = emmeans(FEveRip, ~ rip_modif)
em
letras = cld(em, Letters = letters, adjust = "Tukey") 

letras = letras %>%
  mutate(rip_modif = as.factor(rip_modif)) %>%
  select(rip_modif, .group)
medias = modelo %>%
  group_by(rip_modif) %>%
  summarise(FEve = max(FEve)) %>%
  left_join(letras, by = "rip_modif")

# Boxplot 
boxplotFEve=ggplot(modelo, aes(x = rip_modif, y = FEve, fill = rip_modif)) +
  geom_boxplot(color = "black", outlier.shape = NA) +  
  geom_jitter(width = 0.1, size = 1, shape = 21, color = "black", fill = "black", stroke = 0.3) +  
  scale_fill_grey()+
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(color = "black", size=15),
    axis.title = element_text(color = "black", size=15)
  ) +
  labs(x = "", y = "FEve") +
  geom_text(data = medias, aes(x = rip_modif, y = FEve + 0.05, label = .group),
            color = "black", size = 6)
print(boxplotFEve)

###FDiv###
summary(modelo[,c("FDiv", "rip_modif")])
#modelo
FDivRip=lm(FDiv~rip_modif, data=modelo)
boxplot(FDiv~rip_modif, data=modelo) 
summary(FDivRip)
#Intercept: grupo de referencia (rip_modif:high)
#estimate: diferenecia entre la media con el grupo de referencia
anova(FDivRip)
#F<0.005→ diferencia singnificativa entre grupos
#test a posteriori (p<0.05)
FDivRip.post=glht(model=FDivRip, linfct=mcp(rip_modif="Tukey"))
summary(FDivRip.post, test=adjusted("fdr"))

#Analisis de residuos
res.FDivRip=augment(FDivRip)
FDivRip.res=ggplot(data=res.FDivRip, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

FDivRip.res.vs.expl=ggplot(data=res.FDivRip, aes(x=rip_modif, y=.resid))+
  geom_violin(alpha=0.9)+
  geom_jitter (size=2, width=0.1)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="rip_modif", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(),
        axis.text = element_text(size=16))

FDivRip.Cook=ggplot(data=res.FDivRip, aes(x=1:nrow(res.FDivRip),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

FDivRip.qq=ggplot(res.FDivRip, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(FDivRip.res, FDivRip.res.vs.expl, FDivRip.qq, FDivRip.Cook, ncol=2) 

#Residuales ajustan bastante bien

#Presentación de resultados
#CI
pred.FDivRip=ggpredict(FDivRip, terms = c("rip_modif"))
pred.FDivRip

#letras para diferencia de medias
em = emmeans(FDivRip, ~ rip_modif)
em
letras = cld(em, Letters = letters, adjust = "Tukey") 

letras = letras %>%
  mutate(rip_modif = as.factor(rip_modif)) %>%
  select(rip_modif, .group)
medias = modelo %>%
  group_by(rip_modif) %>%
  summarise(FDiv = max(FDiv)) %>%
  left_join(letras, by = "rip_modif")

# Boxplot 
boxplotFDiv=ggplot(modelo, aes(x = rip_modif, y = FDiv, fill = rip_modif)) +
  geom_boxplot(color = "black", outlier.shape = NA) +  
  geom_jitter(width = 0.1, size = 1, shape = 21, color = "black", fill = "black", stroke = 0.3) +  
  scale_fill_grey() +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15)
  ) +
  labs(x = "", y = "FDiv") +
  geom_text(data = medias, aes(x = rip_modif, y = FDiv + 0.02, label = .group),
            color = "black", size = 6)
print(boxplotFDiv)

###RaoQ###
summary(modelo[,c("RaoQ", "rip_modif")])
#modelo
RaoQRip=lm(RaoQ~rip_modif, data=modelo)
boxplot(RaoQ~rip_modif, data=modelo) 
summary(RaoQRip)
#Intercept: grupo de referencia (rip_modif:high)
#estimate: diferenecia entre la media con el grupo de referencia
anova(RaoQRip)
#F<0.005→ diferencia singnificativa entre grupos
#test a posteriori (p<0.05)
RaoQRip.post=glht(model=RaoQRip, linfct=mcp(rip_modif="Tukey"))
summary(RaoQRip.post, test=adjusted("fdr"))

#Analisis de residuos
res.RaoQRip=augment(RaoQRip)
RaoQRip.res=ggplot(data=res.RaoQRip, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

RaoQRip.res.vs.expl=ggplot(data=res.RaoQRip, aes(x=rip_modif, y=.resid))+
  geom_violin(alpha=0.9)+
  geom_jitter (size=2, width=0.1)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="rip_modif", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(),
        axis.text = element_text(size=16))

RaoQRip.Cook=ggplot(data=res.RaoQRip, aes(x=1:nrow(res.RaoQRip),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

RaoQRip.qq=ggplot(res.RaoQRip, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(RaoQRip.res, RaoQRip.res.vs.expl, RaoQRip.qq, RaoQRip.Cook, ncol=2) 

#Residuales ajustan bastante bien

#Presentación de resultados
#CI
pred.RaoQRip=ggpredict(RaoQRip, terms = c("rip_modif"))
pred.RaoQRip

#letras para diferencia de medias
em = emmeans(RaoQRip, ~ rip_modif)
em
letras = cld(em, Letters = letters, adjust = "Tukey") 

letras = letras %>%
  mutate(rip_modif = as.factor(rip_modif)) %>%
  select(rip_modif, .group)
medias = modelo %>%
  group_by(rip_modif) %>%
  summarise(RaoQ = max(RaoQ)) %>%
  left_join(letras, by = "rip_modif")

# Boxplot 
boxplotRaoQ=ggplot(modelo, aes(x = rip_modif, y = RaoQ, fill = rip_modif)) +
  geom_boxplot(color = "black", outlier.shape = NA) +  
  geom_jitter(width = 0.1, size = 1, shape = 21, color = "black", fill = "black", stroke = 0.3) +  
  scale_fill_grey()+
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(color = "black", size=15),
    axis.title = element_text(color = "black", size=15)
  ) +
  labs(x = "", y = "RaoQ") +
  geom_text(data = medias, aes(x = rip_modif, y = RaoQ + 0.01, label = .group),
            color = "black", size = 6)
print(boxplotRaoQ)

plot_grid(boxplotFRic, boxplotFEve, boxplotFDiv,boxplotRaoQ, ncol=4 )

#Agregado de 1 solo titulo de eje x
combined_plot = plot_grid(boxplotFRic, boxplotFEve, boxplotFDiv, boxplotRaoQ, 
                          ncol = 4, align = "v", axis = "tb")

boxplotFdiversity = plot_grid(
  combined_plot,
  ggdraw() + draw_label("", fontface = 'bold', size = 15),
  ncol = 1,
  rel_heights = c(1, 0.09) # adjust this to control spacing
)

print(boxplotFdiversity)

###HABITAT FRAGMENTATION###
mFRic=lm(FRic ~ fragmentation, data=modelo)   
summary(mFRic)
anova(mFRic, test="F")

#Analisis de residuos
res.mFRic=augment(mFRic)

mFRic.res=ggplot(data=res.mFRic, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+ theme_bw()+ geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mFRic.res.vs.expl=ggplot(data=res.mFRic, aes(x=fragmentation, y=.resid))+
  geom_point(size=2)+theme_bw()+geom_hline(yintercept = 0)+ 
  labs(x="fragmentation", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(), axis.text = element_text(size=16))

mFRic.qq=ggplot(res.mFRic, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mFRic.Cook=ggplot(data=res.mFRic, aes(x=1:nrow(res.mFRic),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(mFRic.res, mFRic.res.vs.expl, mFRic.qq, mFRic.Cook, ncol=2)
#Ajustan mas o menos bien

#Curva predicha por el modelo
pred.FRic <- ggpredict(mFRic, terms = c("fragmentation"))
plotFRic_frag <- plot(pred.FRic) +
  labs(x = "Fragmentation (LDI)", y = "FRic") +
  geom_line(size = 1.2) +
  scale_y_continuous(breaks = seq(1, 1.5, 0.1)) +
  theme_minimal(base_family = "") + 
  theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    plot.title = element_blank(),
    panel.grid = element_blank(),          # elimina grilla
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black")
  )
print(plotFRic_frag)

###FEve#########################################################################
#MODELO
mFEve=lm(FEve ~ fragmentation, data=modelo)   
summary(mFEve)
anova(mFEve, test="F")

#No significativo

#Analisis de residuos
res.mFEve=augment(mFEve)

mFEve.res=ggplot(data=res.mFEve, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+ theme_bw()+ geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mFEve.res.vs.expl=ggplot(data=res.mFEve, aes(x=fragmentation, y=.resid))+
  geom_point(size=2)+theme_bw()+geom_hline(yintercept = 0)+ 
  labs(x="fragmentation", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(), axis.text = element_text(size=16))

mFEve.qq=ggplot(res.mFEve, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mFEve.Cook=ggplot(data=res.mFEve, aes(x=1:nrow(res.mFEve),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(mFEve.res, mFEve.res.vs.expl, mFEve.qq, mFEve.Cook, ncol=2)
#Ajustan mas o menos bien

#Gráfico
pred.FEve <- ggpredict(mFEve, terms = c("fragmentation"))
plotFEve_frag <- plot(pred.FEve) +
  labs(x = "Fragmentation (LDI)", y = "FEve") +
  geom_line(size = 1.2) +
  theme_minimal(base_family = "") + 
  theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    plot.title = element_blank(),
    panel.grid = element_blank(),          # elimina grilla
    panel.border = element_rect(colour = "black", fill = NA, size = 1), # marco negro
    axis.line = element_line(color = "black")
  )
print(plotFEve_frag)

###FDiv#########################################################################
#MODELO
mFDiv=lm(FDiv ~ fragmentation, data=modelo)   
summary(mFDiv)
anova(mFDiv, test="F")
#No significativo

#Analisis de residuos
res.mFDiv=augment(mFDiv)

mFDiv.res=ggplot(data=res.mFDiv, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+ theme_bw()+ geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mFDiv.res.vs.expl=ggplot(data=res.mFDiv, aes(x=fragmentation, y=.resid))+
  geom_point(size=2)+theme_bw()+geom_hline(yintercept = 0)+ 
  labs(x="fragmentation", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(), axis.text = element_text(size=16))

mFDiv.qq=ggplot(res.mFDiv, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mFDiv.Cook=ggplot(data=res.mFDiv, aes(x=1:nrow(res.mFDiv),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(mFDiv.res, mFDiv.res.vs.expl, mFDiv.qq, mFDiv.Cook, ncol=2)
#Ajustan mas o menos bien

#Gráfico
pred.FDiv <- ggpredict(mFDiv, terms = c("fragmentation"))
plotFDiv_frag <- plot(pred.FDiv) +
  labs(x = "Fragmentation (LDI)", y = "FDiv") +
  geom_line(size = 1.2) +
  scale_y_continuous(breaks = seq(0.7, 0.9, 0.02)) +
  theme_minimal(base_family = "") + 
  theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    plot.title = element_blank(),
    panel.grid = element_blank(),          # elimina grilla
    panel.border = element_rect(colour = "black", fill = NA, size = 1), # marco negro
    axis.line = element_line(color = "black")
  )
print(plotFDiv_frag)

###RaoQ#########################################################################
mRaoQ=lm(RaoQ ~ fragmentation, data=modelo)   
summary(mRaoQ)
anova(mRaoQ, test="F")

#Residuales
res.mRaoQ=augment(mRaoQ)

mRaoQ.res=ggplot(data=res.mRaoQ, aes(x=.fitted, y=.resid))+
  geom_point(size=2)+ theme_bw()+ geom_hline(yintercept = 0)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mRaoQ.res.vs.expl=ggplot(data=res.mRaoQ, aes(x=fragmentation, y=.resid))+
  geom_point(size=2)+theme_bw()+geom_hline(yintercept = 0)+ 
  labs(x="fragmentation", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y = element_blank(), axis.text = element_text(size=16))

mRaoQ.qq=ggplot(res.mRaoQ, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

mRaoQ.Cook=ggplot(data=res.mRaoQ, aes(x=1:nrow(res.mRaoQ),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+labs(x="dato", y="Distancia de Cook")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(mRaoQ.res, mRaoQ.res.vs.expl, mRaoQ.qq, mRaoQ.Cook, ncol=2)

pred.RaoQ <- ggpredict(mRaoQ, terms = c("fragmentation"))
plotRaoQ_frag <- plot(pred.RaoQ) +
  labs(x = "Fragmentation (LDI)", y = "RaoQ") +
  geom_line(size = 1.2) +
  theme_minimal(base_family = "") + 
  scale_y_continuous(breaks = seq(0.055, 0.08, 0.005)) +
  theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    plot.title = element_blank(),
    panel.grid = element_blank(),          # elimina grilla
    panel.border = element_rect(colour = "black", fill = NA, size = 1), # marco negro
    axis.line = element_line(color = "black")
  )
print(plotRaoQ_frag)

#Union de gráficos
plotFRic_frag_clean <- plotFRic_frag + theme(axis.title.x = element_blank())
plotFEve_frag_clean <- plotFEve_frag + theme(axis.title.x = element_blank())
plotFDiv_frag_clean <- plotFDiv_frag + theme(axis.title.x = element_blank())
plotRaoQ_frag_clean <- plotRaoQ_frag + theme(axis.title.x = element_blank())

combined_plots <- plot_grid(plotFRic_frag_clean, plotFEve_frag_clean,
                            plotFDiv_frag_clean, plotRaoQ_frag_clean,
                            ncol = 4, align = "v")
x_label <- ggdraw() + 
  draw_label("", fontface = 'bold', size = 16, x = 0.5, hjust = 0.5)
final_plotfragdiv <- plot_grid(combined_plots, x_label, ncol = 1, rel_heights = c(1, 0.07))

print(final_plotfragdiv)

#Figura 3
x11()
plot_grid(boxplotFdiversity, final_plotfragdiv, nrow=2)

#################### Human Interventions ######################################
###FRic########################################################################
colnames(modelo)
#Modelo
modelo_interv=lm(FRic~ornamental_abund+exoticsp_abund, data=modelo)
summary(modelo_interv)

#Intercept: Media de la variable respuesta cuando la abundancia de especies exoticas y ornamentales es 0
#El cambio de 1 unidad de abundancia de ornamentales aumenta FRic en 1.162

1.162/0.229 #especies ornamentales 5% más importante

#Visualizacion de resultados
tidy(modelo_interv, conf.int=T)
ggcoef(modelo_interv, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))
#estimate: pendiente de las variables explicativas

#Analisis de residuales
res.modelo_interv=augment(modelo_interv)
res.modelo_interv

modelo_interv.res=ggplot(data=res.modelo_interv, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.ornamental=ggplot(data=res.modelo_interv, aes(x=ornamental_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="ornamental abund", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.exotic=ggplot(data=res.modelo_interv, aes(x=exoticsp_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="exoticsp_abund", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.Cook=ggplot(data=res.modelo_interv, aes(x=1:nrow(res.modelo_interv),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.qq=ggplot(res.modelo_interv, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_interv.res, modelo_interv.res.vs.ornamental, modelo_interv.res.vs.exotic, modelo_interv.qq,modelo_interv.Cook, ncol=3)

#Residuales se ajustan "decentemente"

#Curvas condicionales (3 dimensiones)

pred.modelo_inerv=ggpredict(modelo_interv, terms = c("ornamental_abund", "exoticsp_abund"))
pred.modelo_inerv

# Curva condicional de predicciones 
plot(pred.modelo_inerv) +
  labs(x="ornamental_abund", y="FRic")+
  geom_point(data=modelo, aes(x=ornamental_abund, y=FRic), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))


###FEve#########################################################################
#Modelo
modelo_interv=lm(FEve~ornamental_abund+exoticsp_abund, data=modelo)
summary(modelo_interv)

#Intercept: Media de la variable respuesta cuando la abundancia de especies exoticas y ornamentales es 0
#El cambio de 1 unidad de abundancia de ornamentales disminuye FEve en 0.3445

0.3445/0.0569 #especies ornamentales 6% más importante

#Visualizacion de resultados
tidy(modelo_interv, conf.int=T)
ggcoef(modelo_interv, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))
#estimate: pendiente de las variables explicativas

#Analisis de residuales
res.modelo_interv=augment(modelo_interv)
res.modelo_interv

modelo_interv.res=ggplot(data=res.modelo_interv, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.ornamental=ggplot(data=res.modelo_interv, aes(x=ornamental_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="ornamental abund", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.exotic=ggplot(data=res.modelo_interv, aes(x=exoticsp_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="exoticsp_abund", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.Cook=ggplot(data=res.modelo_interv, aes(x=1:nrow(res.modelo_interv),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.qq=ggplot(res.modelo_interv, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_interv.res, modelo_interv.res.vs.ornamental, modelo_interv.res.vs.exotic, modelo_interv.qq,modelo_interv.Cook, ncol=3)

#Residuales se ajustan "casi perfecto"

#Curvas condicionales (3 dimensiones)

pred.modelo_inerv=ggpredict(modelo_interv, terms = c("ornamental_abund", "exoticsp_abund"))
pred.modelo_inerv

# Curva condicional de predicciones 
plot(pred.modelo_inerv) +
  labs(x="ornamental_abund", y="FEve")+
  geom_point(data=modelo, aes(x=ornamental_abund, y=FEve), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))

###FDIV#########################################################################
#Modelo
modelo_interv=lm(FDiv~ornamental_abund+exoticsp_abund, data=modelo)
summary(modelo_interv)

#Intercept: Media de la variable respuesta cuando la abundancia de especies exoticas y ornamentales es 0

#Visualizacion de resultados
tidy(modelo_interv, conf.int=T)
ggcoef(modelo_interv, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))
#estimate: pendiente de las variables explicativas

#Analisis de residuales
res.modelo_interv=augment(modelo_interv)
res.modelo_interv

modelo_interv.res=ggplot(data=res.modelo_interv, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.ornamental=ggplot(data=res.modelo_interv, aes(x=ornamental_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="ornamental abund", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.exotic=ggplot(data=res.modelo_interv, aes(x=exoticsp_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="exoticsp_abund", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.Cook=ggplot(data=res.modelo_interv, aes(x=1:nrow(res.modelo_interv),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.qq=ggplot(res.modelo_interv, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_interv.res, modelo_interv.res.vs.ornamental, modelo_interv.res.vs.exotic, modelo_interv.qq,modelo_interv.Cook, ncol=3)

#Residuales se ajustan "mas o menos aceptable"

#Curvas condicionales (3 dimensiones)

pred.modelo_inerv=ggpredict(modelo_interv, terms = c("ornamental_abund", "exoticsp_abund"))
pred.modelo_inerv

# Curva condicional de predicciones 
plot(pred.modelo_inerv) +
  labs(x="ornamental_abund", y="FDiv")+
  geom_point(data=modelo, aes(x=ornamental_abund, y=FDiv), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))

###RaoQ#########################################################################
#Modelo
modelo_interv=lm(RaoQ~ornamental_abund+exoticsp_abund, data=modelo)
summary(modelo_interv)

#Intercept: Media de la variable respuesta cuando la abundancia de especies exoticas y ornamentales es 0

#Visualizacion de resultados
tidy(modelo_interv, conf.int=T)
ggcoef(modelo_interv, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))
#estimate: pendiente de las variables explicativas

#Analisis de residuales
res.modelo_interv=augment(modelo_interv)
res.modelo_interv

modelo_interv.res=ggplot(data=res.modelo_interv, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.ornamental=ggplot(data=res.modelo_interv, aes(x=ornamental_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="ornamental abund", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_interv.res.vs.exotic=ggplot(data=res.modelo_interv, aes(x=exoticsp_abund, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="exoticsp_abund", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.Cook=ggplot(data=res.modelo_interv, aes(x=1:nrow(res.modelo_interv),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_interv.qq=ggplot(res.modelo_interv, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_interv.res, modelo_interv.res.vs.ornamental, modelo_interv.res.vs.exotic, modelo_interv.qq,modelo_interv.Cook, ncol=3)

#Residuales se ajustan mas o menos bien

#Curvas condicionales (3 dimensiones)

pred.modelo_inerv=ggpredict(modelo_interv, terms = c("ornamental_abund", "exoticsp_abund"))
pred.modelo_inerv

# Curva condicional de predicciones 
plot(pred.modelo_inerv) +
  labs(x="ornamental_abund", y="RaoQ")+
  geom_point(data=modelo, aes(x=ornamental_abund, y=RaoQ), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))

######################## Urban Local Environment ###############################
##CARGA DE DATOS
modelo = read.csv("ModeloFiltros_Diversidad.csv", header=T,  sep = ",", stringsAsFactors = T)
View(modelo)
#Variables explicativaS: N, P, pH, CEA, COS, mean_ST, light_hours

#Modelos lineales multiples → Distribuciones normales
summary(modelo)

#Logaritmizar variables de suelo para reducir el orden de magnitud (skewness)
modelo$N_log = log(modelo$N)
modelo$P_log = log(modelo$P)
modelo$COS_log = log(modelo$COS)
modelo$CEA_log = log(modelo$CEA)

#Agregar la variable "shading" creando un índice proporcional a las horas de luz 
modelo$shading=1-(modelo$light_hours/max(modelo$light_hours))
View(modelo)

#eliminar plots 9 y 10 urbano (faltan los datos de suelos)
modelo = modelo %>%
  filter(!plot %in% c("U9", "U10"))

corr1 = modelo[, c("N_log", "P_log", "COS_log", "CEA_log", "pH", "mean_ST", "shading")]
my_cor <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  cor_test <- cor.test(x, y, use = "complete.obs")
  r <- cor_test$estimate  # Coeficiente de correlación
  p_value <- cor_test$p.value  # Valor p
  r2 <- round(r^2, 2)  # R²
  label <- paste0("R² = ", r2, "\np = ", format(p_value, digits = 3))
  ggplot(data, mapping) +
    annotate("text",
             x = mean(range(x, na.rm = TRUE)),
             y = mean(range(y, na.rm = TRUE)),
             label = label, size = 5) +
    theme_void()
}
ggpairs(corr1,
        upper = list(continuous = my_cor),
        lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.5)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.5)))

#Se elimina COS del modelo por alta relación con N y media relación con pH
#Estandarizamos las variables explicativas
data_env.FRic=as.data.frame(cbind(FRic=modelo$FRic, scale(modelo[,c("N_log", "P_log", "CEA_log", 
                                                                    "pH", "mean_ST", "shading")], center=T, scale=T)))
modelo_env=lm(FRic~N_log+P_log+CEA_log+pH+mean_ST+shading,data=data_env.FRic)
summary(modelo_env) 

#intercept= FRic en sitio con N, P, CEA, pH, ST, shading promedio
#Slope(estimate)= el grado de variación de FRic cuando la variable explicativa aumenta una unidad

#Visualizacion de resultados
tidy(modelo_env, conf.int=T)
ggcoef(modelo_env, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=15))

#Analisis de residuales
res.modelo_env=augment(modelo_env)
res.modelo_env

modelo_env.res=ggplot(data=res.modelo_env, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.N=ggplot(data=res.modelo_env, aes(x=N_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="N", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_env.res.vs.P=ggplot(data=res.modelo_env, aes(x=P_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="P", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.pH=ggplot(data=res.modelo_env, aes(x=pH, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="pH", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.CEA=ggplot(data=res.modelo_env, aes(x=CEA_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="CEA", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.shade=ggplot(data=res.modelo_env, aes(x=shading, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="shade", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.temp=ggplot(data=res.modelo_env, aes(x=mean_ST, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="temp", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.Cook=ggplot(data=res.modelo_env, aes(x=1:nrow(res.modelo_env),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.qq=ggplot(res.modelo_env, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_env.res,modelo_env.res.vs.N,modelo_env.res.vs.P,modelo_env.res.vs.pH,modelo_env.res.vs.CEA,modelo_env.res.vs.shade,modelo_env.res.vs.temp,modelo_env.qq,modelo_env.Cook, ncol=3)

#Los residuales se ajustan decentemente

#Curvas condicionales (7 dimensiones)

pred.modelo_env=ggpredict(modelo_env, terms = c("N_log", "pH"))
pred.modelo_env

# Curva condicional de predicciones 
plot(pred.modelo_env) +
  labs(x="pH", y="FRic")+
  geom_point(data=data_env.FRic, aes(x=pH, y=FRic), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))

###FEve#########################################################################
#Estandarizamos las variables explicativas
data_env.FEve=as.data.frame(cbind(FEve=modelo$FEve, scale(modelo[,c("N_log", "P_log", "CEA_log", 
                                                                    "pH", "mean_ST", "shading")], center=T, scale=T)))
View(data_env.FEve)

#Modelo
modelo_env=lm(FEve~N_log+P_log+CEA_log+pH+mean_ST+shading,data=data_env.FEve)
summary(modelo_env) 

#intercept= FRic en sitio con N, P, CEA, pH, ST, shading promedio
#Slope(estimate)= el grado de variación de FRic cuando la variable explicativa aumenta una unidad

#Visualizacion de resultados
tidy(modelo_env, conf.int=T)
ggcoef(modelo_env, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))

#Analisis de residuales
res.modelo_env=augment(modelo_env)
res.modelo_env

modelo_env.res=ggplot(data=res.modelo_env, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.N=ggplot(data=res.modelo_env, aes(x=N_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="N", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_env.res.vs.P=ggplot(data=res.modelo_env, aes(x=P_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="P", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.pH=ggplot(data=res.modelo_env, aes(x=pH, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="pH", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.CEA=ggplot(data=res.modelo_env, aes(x=CEA_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="CEA", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.shade=ggplot(data=res.modelo_env, aes(x=shading, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="shade", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.temp=ggplot(data=res.modelo_env, aes(x=mean_ST, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="temp", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.Cook=ggplot(data=res.modelo_env, aes(x=1:nrow(res.modelo_env),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.qq=ggplot(res.modelo_env, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_env.res,modelo_env.res.vs.N,modelo_env.res.vs.P,modelo_env.res.vs.pH,modelo_env.res.vs.CEA,modelo_env.res.vs.shade,modelo_env.res.vs.temp,modelo_env.qq,modelo_env.Cook, ncol=3)

#Los residuales se ajustan decentemente

#Curvas condicionales (7 dimensiones)

pred.modelo_env=ggpredict(modelo_env, terms = c("N_log", "pH"))
pred.modelo_env

# Curva condicional de predicciones 
plot(pred.modelo_env) +
  labs(x="pH", y="FEve")+
  geom_point(data=data_env.FEve, aes(x=pH, y=FEve), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))

###FDiv#########################################################################
#Estandarizamos las variables explocativas
data_env.FDiv=as.data.frame(cbind(FDiv=modelo$FDiv, scale(modelo[,c("N_log", "P_log", "CEA_log", 
                                                                    "pH", "mean_ST", "shading")], center=T, scale=T)))
View(data_env.FDiv)

#Modelo
modelo_env=lm(FDiv~N_log+P_log+CEA_log+pH+mean_ST+shading,data=data_env.FDiv)
summary(modelo_env) 

#intercept= FRic en sitio con N, P, CEA, pH, ST, shading promedio
#Slope(estimate)= el grado de variación de FRic cuando la variable explicativa aumenta una unidad

#Visualizacion de resultados
tidy(modelo_env, conf.int=T)
ggcoef(modelo_env, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))

#Analisis de residuales
res.modelo_env=augment(modelo_env)
res.modelo_env

modelo_env.res=ggplot(data=res.modelo_env, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.N=ggplot(data=res.modelo_env, aes(x=N_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="N", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_env.res.vs.P=ggplot(data=res.modelo_env, aes(x=P_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="P", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.pH=ggplot(data=res.modelo_env, aes(x=pH, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="pH", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.CEA=ggplot(data=res.modelo_env, aes(x=CEA_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="CEA", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.shade=ggplot(data=res.modelo_env, aes(x=shading, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="shade", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.temp=ggplot(data=res.modelo_env, aes(x=mean_ST, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="temp", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.Cook=ggplot(data=res.modelo_env, aes(x=1:nrow(res.modelo_env),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.qq=ggplot(res.modelo_env, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_env.res,modelo_env.res.vs.N,modelo_env.res.vs.P,modelo_env.res.vs.pH,modelo_env.res.vs.CEA,modelo_env.res.vs.shade,modelo_env.res.vs.temp,modelo_env.qq,modelo_env.Cook, ncol=3)

#Los residuales se ajustan decentemente

#Curvas condicionales (7 dimensiones)

pred.modelo_env=ggpredict(modelo_env, terms = c("N_log", "pH"))
pred.modelo_env

# Curva condicional de predicciones 
plot(pred.modelo_env) +
  labs(x="pH", y="FDiv")+
  geom_point(data=data_env.FDiv, aes(x=pH, y=FDiv), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))

###RaoQ#########################################################################
#Estandarizamos las variables explocativas
data_env.RaoQ=as.data.frame(cbind(RaoQ=modelo$RaoQ, scale(modelo[,c("N_log", "P_log", "CEA_log", 
                                                                    "pH", "mean_ST", "shading")], center=T, scale=T)))
View(data_env.RaoQ)

#Modelo
modelo_env=lm(RaoQ~N_log+P_log+CEA_log+pH+mean_ST+shading,data=data_env.RaoQ)
summary(modelo_env) 

#intercept= FRic en sitio con N, P, CEA, pH, ST, shading promedio
#Slope(estimate)= el grado de variación de FRic cuando la variable explicativa aumenta una unidad

#Visualizacion de resultados
tidy(modelo_env, conf.int=T)
ggcoef(modelo_env, exclude_intercept=T,errorbar_height=.2,color = "black",sort = "ascending")+ 
  geom_point(aes(size=10),col="black")+ 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18))

#Analisis de residuales
res.modelo_env=augment(modelo_env)
res.modelo_env

modelo_env.res=ggplot(data=res.modelo_env, aes(x=.fitted, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_smooth(se=F)+ 
  labs(x="fitted values", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.N=ggplot(data=res.modelo_env, aes(x=N_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="N", y="stnd residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))

modelo_env.res.vs.P=ggplot(data=res.modelo_env, aes(x=P_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="P", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.pH=ggplot(data=res.modelo_env, aes(x=pH, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="pH", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.CEA=ggplot(data=res.modelo_env, aes(x=CEA_log, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="CEA", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.shade=ggplot(data=res.modelo_env, aes(x=shading, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="shade", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.res.vs.temp=ggplot(data=res.modelo_env, aes(x=mean_ST, y=.resid))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept = 0)+ 
  labs(x="temp", y="stnd residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.Cook=ggplot(data=res.modelo_env, aes(x=1:nrow(res.modelo_env),y=.cooksd))+
  geom_linerange(aes(ymin=0, ymax=.cooksd))+
  theme_bw()+
  labs(x="data point", y="Cook distance")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

modelo_env.qq=ggplot(res.modelo_env, mapping=aes(sample = .std.resid)) + 
  stat_qq()+stat_qq_line()+
  theme_bw()+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), axis.text = element_text(size=16))

plot_grid(modelo_env.res,modelo_env.res.vs.N,modelo_env.res.vs.P,modelo_env.res.vs.pH,modelo_env.res.vs.CEA,modelo_env.res.vs.shade,modelo_env.res.vs.temp,modelo_env.qq,modelo_env.Cook, ncol=3)

#Los residuales se ajustan decentemente

#Curvas condicionales (7 dimensiones)

pred.modelo_env=ggpredict(modelo_env, terms = c("CEA_log", "pH"))
pred.modelo_env

# Curva condicional de predicciones 
plot(pred.modelo_env) +
  labs(x="CEA", y="RaoQ")+
  geom_point(data=data_env.RaoQ, aes(x=CEA_log, y=RaoQ), size=2, inherit.aes = F)+
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD", "Mean", "+1 SD"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = c(0.9,0.85),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        panel.grid.major=element_line(color="black"))

################################ Figura S1 #####################################
colnames(modelo)
View(modelo)

#N
x = modelo$imprevious_area
y = modelo$N_log
plot(x, y, xlab = "", ylab = "LogTN (ppm)", pch = 19, col = "darkgrey")
lm_model = lm(y ~ x)
abline(lm_model, lwd = 5)
x_pred <- seq(min(x), max(x), length.out = 1000)
pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
polygon(c(x_pred, rev(x_pred)),
        c(pred[, 2], rev(pred[, 3])),
        col = rgb(1, 0, 0, alpha = 0.1), border = NA)
r2 <- summary(lm_model)$r.squared
mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")

# Record the plot
N_plot <- recordPlot()
print(N_plot)

#P
  x=modelo$imprevious_area
  y=modelo$P_log
  plot(x,y, xlab = "", ylab = "LogTP (ppm)", pch=19, col = "darkgrey")
  lm_model = lm(y ~ x)
  abline(lm_model, lwd = 5)
  x_pred <- seq(min(x), max(x), length.out = 1000) 
  pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
  polygon(c(x_pred, rev(x_pred)), 
          c(pred[, 2], rev(pred[, 3])), 
          col = rgb(1, 0, 0, alpha = 0.1), border = NA)
  r2 <- summary(lm_model)$r.squared
  mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")

# Record the plot
P_plot <- recordPlot()
print(P_plot)

#EC
  x=modelo$imprevious_area
  y=modelo$CEA_log
  plot(x,y, xlab = "", ylab = "LogEC (mS/cm)", pch=19, col = "darkgrey")
  lm_model = lm(y ~ x)
  abline(lm_model, lwd = 5)
  x_pred <- seq(min(x), max(x), length.out = 1000) 
  pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
  polygon(c(x_pred, rev(x_pred)), 
          c(pred[, 2], rev(pred[, 3])), 
          col = rgb(1, 0, 0, alpha = 0.1), border = NA)
  r2 <- summary(lm_model)$r.squared
  mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")

  # Record the plot
  CEA_plot <- recordPlot()
  print(CEA_plot)
  
#pH
  x=modelo$imprevious_area
  y=modelo$pH
  plot(x,y, xlab = "", ylab = "pH", pch=19, col = "darkgrey")
  lm_model = lm(y ~ x)
  abline(lm_model, lwd = 5)
  x_pred <- seq(min(x), max(x), length.out = 1000) 
  pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
  polygon(c(x_pred, rev(x_pred)), 
          c(pred[, 2], rev(pred[, 3])), 
          col = rgb(1, 0, 0, alpha = 0.1), border = NA)
  r2 <- summary(lm_model)$r.squared
  mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")

  # Record the plot
  pH_plot <- recordPlot()
  print(pH_plot)  

#UHI
  x=modelo$imprevious_area
  y=modelo$mean_ST
  plot(x,y, xlab = "", ylab = "MNLST of January (UHI - °C)", pch=19, col = "darkgrey")
  lm_model = lm(y ~ x)
  abline(lm_model, lwd = 5)
  x_pred <- seq(min(x), max(x), length.out = 1000) 
  pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
  polygon(c(x_pred, rev(x_pred)), 
          c(pred[, 2], rev(pred[, 3])), 
          col = rgb(1, 0, 0, alpha = 0.1), border = NA)
  r2 <- summary(lm_model)$r.squared
  mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")

  # Record the plot
  UHI_plot <- recordPlot()
  print(UHI_plot)  

#Shading
  x=modelo$imprevious_area
  y=modelo$shading
  plot(x,y, xlab = "", ylab = "Shading Index", pch=19, col = "darkgrey")
  lm_model = lm(y ~ x)
  abline(lm_model, lwd = 5)
  x_pred <- seq(min(x), max(x), length.out = 1000) 
  pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
  polygon(c(x_pred, rev(x_pred)), 
          c(pred[, 2], rev(pred[, 3])), 
          col = rgb(1, 0, 0, alpha = 0.1), border = NA)
  r2 <- summary(lm_model)$r.squared
  mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")

  # Record the plot
  shading_plot <- recordPlot()
  print(shading_plot)

print(N_plot)
print(P_plot)
print(CEA_plot)
print(pH_plot)
print(UHI_plot)
print(shading_plot)

View(modelo)

#Fragmentation
x=modelo$imprevious_area
y=modelo$fragmentation
plot(x,y, xlab = "", ylab = "LDI", pch=19, col = "darkgrey")
lm_model = lm(y ~ x)
abline(lm_model, lwd = 5)
x_pred <- seq(min(x), max(x), length.out = 1000) 
pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
polygon(c(x_pred, rev(x_pred)), 
        c(pred[, 2], rev(pred[, 3])), 
        col = rgb(1, 0, 0, alpha = 0.1), border = NA)
r2 <- summary(lm_model)$r.squared
mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")

# Record the plot
fragment_plot <- recordPlot()
print(fragment_plot)

#Exotic
x=modelo$imprevious_area
y=modelo$exoticsp_abund
plot(x,y, xlab = "", ylab = "Exotic Species Relative Cover", pch=19, col = "darkgrey")
lm_model = lm(y ~ x)
abline(lm_model, lwd = 5)
x_pred <- seq(min(x), max(x), length.out = 1000) 
pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
polygon(c(x_pred, rev(x_pred)), 
        c(pred[, 2], rev(pred[, 3])), 
        col = rgb(1, 0, 0, alpha = 0.1), border = NA)
r2 <- summary(lm_model)$r.squared
mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")

# Record the plot
exotic_plot <- recordPlot()
print(exotic_plot)

#Ornamental
x=modelo$imprevious_area
y=modelo$ornamental_abund
plot(x,y, xlab = "", ylab = "Ornamental Species Relative Cover", pch=19, col = "darkgrey")
lm_model = lm(y ~ x)
abline(lm_model, lwd = 5)
x_pred <- seq(min(x), max(x), length.out = 1000) 
pred <- predict(lm_model, newdata = data.frame(x = x_pred), interval = "confidence")
polygon(c(x_pred, rev(x_pred)), 
        c(pred[, 2], rev(pred[, 3])), 
        col = rgb(1, 0, 0, alpha = 0.1), border = NA)
r2 <- summary(lm_model)$r.squared
mtext(bquote(R^2 == .(round(r2, 3))), side = 3, line = 1, cex = 1.2, col = "black")

# Record the plot
exotic_plot <- recordPlot()
print(exotic_plot)

boxplot(imprevious_area ~ rip_modif, data=modelo, xlab= "Riparian Habitat Transformation", ylab = "Relative Imprevious Area (500m buffer)")
