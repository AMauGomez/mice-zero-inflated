# =============================================================================
# 03_model_selection.R
# Original analysis code from the thesis (UNAM, 2025). Kept as it was run;
# only absolute paths were replaced by relative ones and this header added.
# Comments are in Spanish, as written during the project.
#
# Model building for the imputation models: listwise deletion, stepwise selection by AIC for the
# logistic (zero part), Poisson and negative binomial models, 10-fold cross-validation, dispersion
# checks, and the final ZIP / ZINB models (glmmTMB) used to specify the MICE predictor matrices.
# Variable candidates came from the correlogram and PCA in 01 plus subject-matter judgement.
#
# Not runnable from this repo: the SISVER microdata are not public.
# Data files are expected under data/ (not included). See ../README.md.
# =============================================================================

library(tidyverse)
library(mice)
library(ggmice)
library(readxl)
library(naniar)
library(MASS)
library(GGally)
library(overdisp)
library(AER)
library(performance)
library(DHARMa)
library(glmmTMB)
library("pscl")
library(fitdistrplus)
library(lme4)
library(boot)
library(stargazer)
library(psych)
library(ROCR) 

datos<-read.csv("data/base_covid_procesada.csv")

datos_1<-datos[,1:25]
datos_2<-datos[,26:50]

x11()
vis_1<-vis_miss(datos_1, warn_large_data=FALSE)
vis_1
x11()
vis_2<-vis_miss(datos_2, warn_large_data=FALSE)
vis_2

###SELECCIONAMOS SOLO LAS VARIABLES NUMERICAS PARA INCLUIRLOS EN EL MICE####
datos<-datos[,4:50]

datos<-datos[, -2]
datos<-datos[, - (5:6)]

###ANALISIS DE LA DISTRIBUCIÓN DE LA VARIABLE DE INTERÉS####

p<-ggplot(datos, aes(x=DIASSINTOMAS)) + 
  geom_histogram(color="black", fill="lightblue", bins=100)+
  xlim(-1, 25)+
  geom_histogram(aes(y=..density..), colour="green", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  ylab("Densidad") + xlab("Días con Síntomas") + 
  theme(plot.title = element_text(size=22))+
  theme_gray(base_size = 13)+
  ggtitle("Histograma de los Días con Síntomas en Pacientes Sospechosos de SARS-COV-2")


x11()
p

#Hacer tabla o mencionar  los siguientes datos#3
b<-as.data.frame(datos$DIASSINTOMAS)

stargazer(b, title="Análisis Descriptivo de los Dás con Síntomas",type="latex")

###Contamos los 0####

ceros<-sum(datos$DIASSINTOMAS==0)
ceros
porcentaje_ceros<-(ceros/nrow(datos))*100
porcentaje_ceros

### Revisamos a que distribución se ajustan mejor BINOMIAL O POISSON####

### PRUEBA PARA POISSON ####
poi<-fitdist(datos$DIASSINTOMAS, "pois")

poi

##Hacemos una prueba de hipotesis#

poi_fit<-gofstat(poi)
poi_fit$chisqpvalue

#Se rechaza H_0, por lo tanto no se distribuye Poison 

x11()
plot(poi)

###PRUEBA PARA BINOM NEGATIVA####

bino<-fitdist(datos$DIASSINTOMAS, "nbinom")

summary(bino)

##Hacemos una prueba de hipotesis

bino_fit<-gofstat(bino)
bino_fit$chisqpvalue


#Se rechaza H_0, por lo cual no es binomial##

x11()
plot(bino)


##Comparamos los cuantiles de ambas distribuciones

comp_qq <- qqcomp(list(poi, bino), 
                  xlegend = "bottomright",
                  xlab = "Cuantiles Teóricos",
                  ylab = "Cuantiles Empíricos",
                  main = "Comparación de Cuantiles",
                  fitlty = 1,
                  fitcol = c("orange","blue"),
                  plotstyle = "ggplot",
                  addlegend = TRUE)
x11()
comp_qq


comp_den<-denscomp(list(poi, bino), 
                   xlegend = "bottomright",
                   fitlty = 1,
                   fitcol = c("orange","blue"),
                   xlim = c(0,25),
                   xlab = "Datos",
                   ylab = "Densidad",
                   plotstyle = "ggplot",
                   main = "Comparación de Densidades",
                   addlegend = TRUE)
x11()
comp_den


###SELECCIONAMOS SOLO LAS VARIABLES NUMERICAS PARA INCLUIRLOS EN EL MICE####

datos<-read.csv("data/base_covid_procesada.csv")

cero_totales<-sum(if_else(datos$DIASSINTOMAS==0,1,0))
cero_totales
datos<-datos[,4:50]

datos<-datos[, -2]
datos<-datos[, - (5:6)]

##visualizamos los datos##

datos_1<-datos[,1:22]
datos_2<-datos[,22:44]

x11()
vis_1<-vis_miss(datos_1, warn_large_data=FALSE)
vis_1
x11()
vis_2<-vis_miss(datos_2, warn_large_data=FALSE)
vis_2






####MUY IMPORTANTE HACEMOS LISTWISE DELETION##
datos<-na.omit(datos)

describe(datos$DIASSINTOMAS)

var(datos$DIASSINTOMAS)
mean(datos$DIASSINTOMAS)
sd(datos$DIASSINTOMAS)

ceros<-sum(datos$DIASSINTOMAS==0)
ceros
porcentaje_ceros<-(ceros/nrow(datos))*100
porcentaje_ceros


####REGRESION LOGISTICA STEPWISE PARA RESDEFIN ####
m_log_resde_com <- glm(RESDEFIN~., data = datos, family = "binomial")


summary(m_log_com)
m_log_resde_com_step<-step(m_log_resde_com,
                      direction = "both", steps = 1000, k = 2)

#verifiacamos significancia#

summary(m_log_resde_com_step)


###MODELO DEFINITIVO RESDEFIN STEPWISE SIGNIFICATIVAS####
#nuevo modelo selección del stepwise basado en logica y p values##


m_log_resde_com_step<-glm(RESDEFIN ~ SEXO + INTUBADO + DIGCLINE + EDAD + 
                            DIASSINTOMAS + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + 
                            DISNEA + DIARREA + MIALGIAS + DOLABDO + ANOSMIA + 
                            DISGEUSIA + EPOC + INMUSUPR + VIH.SIDA + 
                            OTRACON + ENFCARDI + OBESIDAD + INSRENCR + TABAQUIS , data = datos, family = "binomial")
summary(m_log_resde_com_step)

vif(m_log_resde_com_step)
#predicciones#

datos_1<-datos
predicciones_resdefin<-predict(m_log_resde_com_step, datos_1, type="response")
predicciones_resdefin<-as.data.frame(predicciones_resdefin)

predicciones_resdefin_0_1<-cbind(datos_1$RESDEFIN, predicciones_resdefin)
predicciones_resdefin_0_1<-as.data.frame(predicciones_resdefin_0_1)

predicciones_resdefin_0_1<-predicciones_resdefin_0_1 %>%   mutate(pred_0_1= if_else(predicciones_resdefin_0_1$predicciones_resdefin>.70,1,0))

#creamos la matriz de confusión##

confux<-table(predicciones_resdefin_0_1$`datos_1$RESDEFIN`, predicciones_resdefin_0_1$pred_0_1)
confux

library(caret) 
resumen_confux<-confusionMatrix(confux)
resumen_confux
###Comparación de los erroes del CV para RESDEFIN####
cv_log_com_step<-cv.glm(datos, m_log_com_step , K=10 )
cv_log_com_step$delta

cv_zip_visual<-cv(m_log_com_step,  datos,  criterion = mse)
cv_zip_visual





###Modelos POISSON para el DIASSINTOMAS BINARIO####

###Revisar sobre dispersión###

###Cambiamos a la variable diassintomas a 0 y 1####

datos$DIASSINTOMAS<-if_else(datos$DIASSINTOMAS==0, 1, 0)
datos$DIASSINTOMAS<-as.integer(datos$DIASSINTOMAS)

##Creamos un modelos logistico para sintomas sí  por pasos

cero_logis_com<- glm(datos$DIASSINTOMAS ~. , family=binomial(link='logit'),data=datos)
summary(cero_logis_com)

cero_logis_com_step<-step(cero_logis_com, direction = "both", k=2)

summary(cero_logis_com_step)

###MODELO DEFINITIVO DIASINTOMAS LOGISTICA STEPWISE SIGNIFICATIVAS####

cero_logis_com_step<-glm(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
                           EDAD + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + DISNEA + 
                           DIARREA + CEFALEA + MIALGIAS + ARTRAL + ATAEDOGE + 
                           RINORREA  + CONJUN + CIANOSIS + INISUBIS + DISGEUSIA + 
                           EPOC + INMUSUPR + ENFCARDI + OBESIDAD + INSRENCR + TXANTIVI + 
                           VACUNADO + RESDEFIN + ANTIPIRETICOS , family=binomial(link='logit'),data=datos)

vif(cero_logis_com_step)

##creamos las predicciones##

#predicciones#

datos_1<-datos
predicciones_1<-predict(cero_logis_com_step, datos_1, type="response")
predicciones_1<-as.data.frame(predicciones_1)

datos_com_pred_ceros<-cbind(datos_1$DIASSINTOMAS, predicciones_1)
datos_com_pred_ceros<-as.data.frame(datos_com_pred_ceros)

datos_com_pred_ceros<-datos_com_pred_ceros %>% mutate(pred_0_1=if_else(datos_com_pred_ceros$predicciones_1>.80,1,0))

#creamos la matriz de confusión##

confux_dias<-table(datos_com_pred_ceros$`datos_1$DIASSINTOMAS`, datos_com_pred_ceros$pred_0_1 )
confux_dias

resumen_confux_dias<-confusionMatrix(confux_dias)
resumen_confux_dias





###Modelos POISSON para el DIASSINTOMAS NUMERICO####

#USAREMOS LAS VARIABLES DEL MODELO ANTERIOR COMO INICIALES COMPLETOS A STEPWISE##
datos<-read.csv("data/base_covid_procesada.csv")

datos<-datos[,4:50]

datos<-datos[, -2]
datos<-datos[, - (5:6)]
datos<-na.omit(datos)

##Cortamos los diasisntomas a 30##
datos<-datos %>% filter(DIASSINTOMAS<30)


m_pois_com<-glm(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
                  EDAD + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + DISNEA + 
                  DIARREA + CEFALEA + MIALGIAS + ARTRAL + ATAEDOGE + RINORREA + 
                  CONJUN + CIANOSIS + INISUBIS + DISGEUSIA + EPOC + INMUSUPR + 
                  ENFCARDI + OBESIDAD + INSRENCR + TXANTIVI + VACUNADO + RESDEFIN + 
                  ANTIPIRETICOS, family = poisson(), data=datos)
summary(m_pois_com)


m_pois_com_step<-step(m_pois_com,
                      direction = "both", steps = 1000, k = 2)


summary(m_pois_com_step)


###MODELO DEFINITIVO DIASINTOMAS NUMERICO STEPWISE SIGNIFICATIVAS####

m_pois_com_step<-glm(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
                  EDAD + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + DISNEA + 
                  DIARREA + CEFALEA + MIALGIAS + ARTRAL + ATAEDOGE + 
                  RINORREA  + CIANOSIS + INISUBIS + DISGEUSIA + 
                  EPOC + INMUSUPR + ENFCARDI + OBESIDAD + INSRENCR + TXANTIVI + 
                  VACUNADO + RESDEFIN + ANTIPIRETICOS, family = poisson(), data=datos)


vif(m_pois_com_step)

##evaluación del modelo step##

cv_m_pois_com_step<-cv.glm(datos, m_pois_com_step, K=10 )
cv_m_pois_com_step$delta


cv_m_pois_com_step_1<-cv(m_pois_com_step,  datos,  criterion = mse)
cv_m_pois_com_step_1


##creamos las predicciones##

datos_1<-datos
predicciones_2<-predict(m_pois_com, datos_1, type="response")
predicciones_2<-as.data.frame(predicciones_2)

datos_pois_pred<-cbind(datos_1$DIASSINTOMAS, predicciones_2)
datos_pois_pred<-as.data.frame(datos_pois_pred)
datos_pois_pred$predicciones_2<-round(datos_pois_pred$predicciones_2)


##Calculamos el mse#

cv_m_pois_com_step_1<-cv(m_pois_com_step,  datos,  criterion = mse)
cv_m_pois_com_step_1





###Modelos BINOM NEG para el DIASSINTOMAS NUMERICO####
##USAREMOS STEPWISE##

datos<-read.csv("data/base_covid_procesada.csv")

datos<-datos[,4:50]

datos<-datos[, -2]
datos<-datos[, - (5:6)]
datos<-na.omit(datos)

##Cortamos los diasisntomas a 30##
datos<-datos %>% filter(DIASSINTOMAS<30)

#CREAMOS LOS MODELOS PARA STEPWISE TOMANDO COMO REFERENCIA LAS VARIABLES SELCCIOANDAS EN POISSON COMO EL COMPLETO#
m_binom_com<-glm.nb(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
                      EDAD + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + DISNEA + 
                      DIARREA + CEFALEA + MIALGIAS + ARTRAL + ATAEDOGE + 
                      RINORREA  + CIANOSIS + INISUBIS + DISGEUSIA + 
                      EPOC + INMUSUPR + ENFCARDI + OBESIDAD + INSRENCR + TXANTIVI + 
                      VACUNADO + RESDEFIN + ANTIPIRETICOS, data=datos)

m_binom_com_step<-step(m_binom_com,
                       direction = "both", steps = 1000, k = 2)
summary(m_binom_com_step)


###MODELO DEFINITIVO DIASINTOMAS NEGBINOM NUMERICO STEPWISE SIGNIFICATIVAS####

m_binom_com_step<-glm.nb(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
                           EDAD + DIAGPROB_IRAG + FIEBRE + TOS + DISNEA + 
                           DIARREA + ARTRAL + ATAEDOGE + RINORREA + 
                           CIANOSIS + INISUBIS + DISGEUSIA + EPOC + INMUSUPR + 
                           ENFCARDI + OBESIDAD + INSRENCR + RESDEFIN + 
                           ANTIPIRETICOS, data=datos)

summary(m_binom_com_step)
vif(m_binom_com_step)
##evaluación del modelo step definitivo##

###Comparación de los erroes del CV BINOM ####

cv_binom_com_step<-cv.glm(datos, m_binom_com_step, K=10 )
cv_binom_com_step$delta
## con la libreria cv##

cv_binom_com_step_1<-cv(m_binom_com_step,  datos,  criterion = mse)
cv_binom_com_step_1


##NO EXISTE DIFERENCIA ENTRE LOS ERRORES, SIN EMBARGO, EL NEG BINOM TIENE MENOS VARIABLES#

###EVALUACIÓN DE SOBREDISPERSIÓN POISSON Y BINO NEG####
###DISPERSIONES POISSON y BINOM PERFORMANCE####

dis_pois<-check_overdispersion(m_pois_com_step)
dis_pois

dis_binom<-check_overdispersion(m_binom_com_step)
dis_binom

##Notamos que tenemos sobre en el Poisson y sub en el BINOM##



###Modelos  CERO INFLADO POISSON para el DIASSINTOMAS NUMERICO####
##USAMOS LAS REGRESIONES LOGISTICAS ANTES ENCONTRADAS para modelar los CEROS##
##y LA REGRESIÓN POISON PARA CONTEOS
zip_def<-glmmTMB(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
                   EDAD + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + DISNEA + 
                   DIARREA + CEFALEA + MIALGIAS + ARTRAL + ATAEDOGE + 
                   RINORREA  + CIANOSIS + INISUBIS + DISGEUSIA + 
                   EPOC + INMUSUPR + ENFCARDI + OBESIDAD + INSRENCR + TXANTIVI + 
                   VACUNADO + RESDEFIN + ANTIPIRETICOS,
                    zi=~ SEXO + INTUBADO + DIGCLINE + 
                   EDAD + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + DISNEA + 
                   DIARREA + CEFALEA + MIALGIAS + ARTRAL + ATAEDOGE + 
                   RINORREA  + CONJUN + CIANOSIS + INISUBIS + DISGEUSIA + 
                   EPOC + INMUSUPR + ENFCARDI + OBESIDAD + INSRENCR + TXANTIVI + 
                   VACUNADO + RESDEFIN + ANTIPIRETICOS,
                    family=poisson(), data=datos)

summary(zip_def)

###MODELO DEFINITIVO DIASINTOMAS ZIP SIGNIFICATIVAS####

zip_def<-glmmTMB(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
                   EDAD +  DISNEA + 
                   DIARREA + CEFALEA + ARTRAL +
                   RINORREA  + CIANOSIS + INISUBIS + DISGEUSIA + 
                   EPOC + ENFCARDI + OBESIDAD + INSRENCR + 
                   VACUNADO + RESDEFIN + ANTIPIRETICOS,
                 zi=~ INTUBADO + DIGCLINE + 
                   EDAD + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + DISNEA + 
                   DIARREA + CEFALEA + MIALGIAS + ATAEDOGE + 
                   RINORREA  + CONJUN + CIANOSIS + INISUBIS + 
                   INMUSUPR + ENFCARDI + INSRENCR + TXANTIVI + 
                   VACUNADO + RESDEFIN + ANTIPIRETICOS,
                 family=poisson(), data=datos)

summary(zip_def)

#EVALUAMOS EL MODELO##
cv_zip_def<-cv(zip_def,  datos,  criterion = mse)
cv_zip_def

####PRUEBAS DE CERO INFLADO POISSON####
#simulamos los residuales del modelos glm#

residu_zip_def<-simulateResiduals(zip_def)
x11()
testZeroInflation(residu_zip_def)
check_zeroinflation(residu_zip_def)
#Notamos que es muy bueno ya que no esta inflando ceros#


###REVISAMOS SI LOS MODELOS ZIP MEJORAN LA SOBREDISPERSIÓN####
##Hay que usar la simulaciones de residuos#
check_overdispersion(residu_zip_def)

##PERSISTE LA SOBREDISPERSIÓN PERO SOLUCIONAMOS LO CERO INFLADO##

### Modelo CERO INFLADO BINOMNEG para el DIASSINTOMAS NUMERICO####
##PARTIMOS DEL MODELO ZIP#


zinb_def<-glmmTMB(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
                    EDAD +  DISNEA + 
                    DIARREA + CEFALEA + ARTRAL +
                    RINORREA  + CIANOSIS + INISUBIS + DISGEUSIA + 
                    EPOC + ENFCARDI + OBESIDAD + INSRENCR + 
                    VACUNADO + RESDEFIN + ANTIPIRETICOS,
                  zi=~ INTUBADO + DIGCLINE + 
                    EDAD + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + DISNEA + 
                    DIARREA + CEFALEA + MIALGIAS + ATAEDOGE + 
                    RINORREA  + 
                    INMUSUPR + ENFCARDI + TXANTIVI + 
                    VACUNADO + RESDEFIN,
                     family=nbinom2(), data=datos)

summary(zinb_def)

#EVALUAMOS EL MODELO##
cv_zinb_def<-cv(zinb_def,  datos,  criterion = mse)
cv_zinb_def

####PRUEBAS DE CERO INFLADO POISSON####
#simulamos los residuales del modelos glm#

residu_zinb_def<-simulateResiduals(zinb_def)
x11()
testZeroInflation(residu_zinb_def)
check_zeroinflation(residu_zinb_def)

#Notamos que es muy bueno ya que no esta inflando ceros#


###REVISAMOS SI LOS MODELOS ZINB MEJORAN LA SOBREDISPERSIÓN####
##Hay que usar la simulaciones de residuos#
check_overdispersion(residu_zinb_def)

#NOTAMOS QUE MEJORA YA QUE NO HAY SOBREDISPERSIÓN NI ESTA CERO INFLADO ES EL MODELO IDEAL##
###El MODELO ZIBN ES EL MEJOR###
