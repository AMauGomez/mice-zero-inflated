# =============================================================================
# 06_evaluation_metrics.R
# Original analysis code from the thesis (UNAM, 2025). Kept as it was run;
# only absolute paths were replaced by relative ones and this header added.
# Comments are in Spanish, as written during the project.
#
# Compares imputed vs. true (masked) values: MSE, RMSE, MAE, exact matches, and zeros recovered,
# for each method and missingness scenario (Table 5.12 of the thesis).
#
# Not runnable from this repo: the SISVER microdata are not public.
# Data files are expected under data/ (not included). See ../README.md.
# =============================================================================

library(tidyverse)
library(mice)
library(lattice)
library(pan)
library(ggmice)
library(readxl)
library(naniar)
library(MASS)
library(psych)
library(dplyr)
library(GGally)
library(overdisp)
library(AER)
library(performance)
library("DHARMa")
library("glmmTMB")
library("pscl")
library("gamlss")
library(countimp)
library(caTools)
library(broom.mixed)
library(lattice)
library(xtable)
library(stargazer)
library("MLmetrics")

##CARGAMOS LAS BASES IMPUTADAS##

datos_pmm_10<-read.csv("data/imputations/imputacion_COMPLETA_pmm_10.csv")
datos_pmm_10<-datos_pmm_10 %>% filter(datos_pmm_10$split==TRUE)


datos_zip_10<-read.csv("data/imputations/imputacion_COMPLETA_zip_10.csv")
datos_zip_10<-datos_zip_10 %>% filter(datos_zip_10$split==TRUE)

datos_zinb_10<-read.csv("data/imputations/imputacion_COMPLETA_zinb_10.csv")
datos_zinb_10<-datos_zinb_10 %>% filter(datos_zinb_10$split==TRUE)

datos_pmm_5<-read.csv("data/imputations/imputacion_COMPLETA_pmm_5.csv")
datos_pmm_5<-datos_pmm_5 %>% filter(datos_pmm_5$split==TRUE)

datos_zip_5<-read.csv("data/imputations/imputacion_COMPLETA__zip_5.csv")
datos_zip_5<-datos_zip_5 %>% filter(datos_zip_5$split==TRUE)

datos_zinb_5<-read.csv("data/imputations/imputacion_COMPLETA__zinb_5.csv")
datos_zinb_5<-datos_zinb_5 %>% filter(datos_zinb_5$split==TRUE)




###MSE IMPUT 10##

mse_pmm_10<-MSE(datos_pmm_10$DIASSINTOMAS, datos_pmm_10$DIASSINTOMAS.VERDADEROS)
mse_pmm_10


mse_zip_10<-MSE(datos_zip_10$DIASSINTOMAS, datos_zip_10$DIASSINTOMAS.VERDADEROS)
mse_zip_10

mse_zinb_10<-MSE(datos_zinb_10$DIASSINTOMAS, datos_zinb_10$DIASSINTOMAS.VERDADEROS)
mse_zinb_10  
  
mse_pmm_10
mse_zip_10
mse_zinb_10


rmse_pmm_10<-sqrt(mse_pmm_10)
rmse_zip_10<-sqrt(mse_zip_10)
rmse_zinb_10<-sqrt(mse_zinb_10)

rmse_pmm_10
rmse_zip_10
rmse_zinb_10



##MAE##

mae_pmm_10<-MAE(datos_pmm_10$DIASSINTOMAS, datos_pmm_10$DIASSINTOMAS.VERDADEROS)

mae_zip_10<-MAE(datos_zip_10$DIASSINTOMAS, datos_zip_10$DIASSINTOMAS.VERDADEROS)

mae_zinb_10<-MAE(datos_zinb_10$DIASSINTOMAS, datos_zinb_10$DIASSINTOMAS.VERDADEROS)

mae_pmm_10
mae_zip_10
mae_zinb_10


###MSE INPUT 5##

mse_pmm_5<-MSE(datos_pmm_5$DIASSINTOMAS, datos_pmm_5$DIASSINTOMAS.VERDADEROS)
mse_pmm_5


mse_zip_5<-MSE(datos_zip_5$DIASSINTOMAS, datos_zip_5$DIASSINTOMAS.VERDADEROS)
mse_zip_5

mse_zinb_5<-MSE(datos_zinb_5$DIASSINTOMAS, datos_zinb_5$DIASSINTOMAS.VERDADEROS)
mse_zinb_5  

mse_pmm_5
mse_zip_5
mse_zinb_5


rmse_pmm_5<-sqrt(mse_pmm_5)
rmse_zip_5<-sqrt(mse_zip_5)
rmse_zinb_5<-sqrt(mse_zinb_5)

rmse_pmm_5
rmse_zip_5
rmse_zinb_5



##MAE##

mae_pmm_5<-MAE(datos_pmm_5$DIASSINTOMAS, datos_pmm_5$DIASSINTOMAS.VERDADEROS)

mae_zip_5<-MAE(datos_zip_5$DIASSINTOMAS, datos_zip_5$DIASSINTOMAS.VERDADEROS)

mae_zinb_5<-MAE(datos_zinb_5$DIASSINTOMAS, datos_zinb_5$DIASSINTOMAS.VERDADEROS)

mae_pmm_5
mae_zip_5
mae_zinb_5

###PORCENTAJE DE ACERTADOS##

pre_pmm_5<-sum(if_else(datos_pmm_5$DIASSINTOMAS==datos_pmm_5$DIASSINTOMAS.VERDADEROS, 1,0))
pre_pmm_10<-sum(if_else(datos_pmm_10$DIASSINTOMAS==datos_pmm_10$DIASSINTOMAS.VERDADEROS, 1,0))


pre_zip_5<-sum(if_else(datos_zip_5$DIASSINTOMAS==datos_zip_5$DIASSINTOMAS.VERDADEROS,1,0))
pre_zip_10<-sum(if_else(datos_zip_10$DIASSINTOMAS==datos_zip_10$DIASSINTOMAS.VERDADEROS,1,0))

pre_zinb_5<-sum(if_else(datos_zinb_5$DIASSINTOMAS==datos_zinb_5$DIASSINTOMAS.VERDADEROS,1,0))
pre_zinb_10<-sum(if_else(datos_zinb_10$DIASSINTOMAS==datos_zinb_10$DIASSINTOMAS.VERDADEROS,1,0))


pre_pmm_5
pre_pmm_10

pre_zip_5
pre_zip_10


pre_zinb_5
pre_zinb_10




##Ceros acertados##

ceros_to<-sum(if_else(datos_zinb_10$DIASSINTOMAS.VERDADEROS==0,1,0))
ceros_to
##PMM 10 y 5##

datos_pmm_10<-datos_pmm_10 %>%  mutate(cero=if_else(datos_pmm_10$DIASSINTOMAS.VERDADEROS==0,1,0))

##filtramos los ceros verdaderos
datos_pmm_10_0<-datos_pmm_10 %>% filter(datos_pmm_10$cero==TRUE)
a<-dim(datos_pmm_10_0)
a  
datos_pmm_10_0<-sum(if_else(datos_pmm_10_0$DIASSINTOMAS==datos_pmm_10_0$DIASSINTOMAS.VERDADEROS,1,0))
datos_pmm_10_0
  


datos_pmm_5<-datos_pmm_5 %>%  mutate(cero=if_else(datos_pmm_5$DIASSINTOMAS.VERDADEROS==0,1,0))

##filtramos los ceros verdaderos
datos_pmm_5_0<-datos_pmm_5 %>% filter(datos_pmm_5$cero==TRUE)
a<-dim(datos_pmm_5_0)
a  
datos_pmm_5_0<-sum(if_else(datos_pmm_5_0$DIASSINTOMAS==datos_pmm_5_0$DIASSINTOMAS.VERDADEROS,1,0))
datos_pmm_5_0



##ZIP 5 Y 10 ###


datos_zip_10<-datos_zip_10 %>%  mutate(cero=if_else(datos_zip_10$DIASSINTOMAS.VERDADEROS==0,1,0))

#filtramos los ceros verdaderos
datos_zip_10_0<-datos_zip_10 %>% filter(datos_zip_10$cero==TRUE)
a<-dim(datos_zip_10_0)
a  
datos_zip_10_0<-sum(if_else(datos_zip_10_0$DIASSINTOMAS==datos_zip_10_0$DIASSINTOMAS.VERDADEROS,1,0))

datos_zip_10_0


datos_zip_5<-datos_zip_5 %>%  mutate(cero=if_else(datos_zip_5$DIASSINTOMAS.VERDADEROS==0,1,0))

datos_zip_5_0<-datos_zip_5 %>% filter(datos_zip_5$cero==TRUE)
a<-dim(datos_zip_5_0)
a  
datos_zip_5_0<-sum(if_else(datos_zip_5_0$DIASSINTOMAS==datos_zip_5_0$DIASSINTOMAS.VERDADEROS,1,0))
datos_zip_5_0




##ZINB 5 Y 10 ###


datos_zinb_10<-datos_zinb_10 %>%  mutate(cero=if_else(datos_zinb_10$DIASSINTOMAS.VERDADEROS==0,1,0))

#filtramos los ceros verdaderos
datos_zinb_10_0<-datos_zinb_10 %>% filter(datos_zinb_10$cero==TRUE)
a<-dim(datos_zinb_10_0)
a  
datos_zinb_10_0<-sum(if_else(datos_zinb_10_0$DIASSINTOMAS==datos_zinb_10_0$DIASSINTOMAS.VERDADEROS,1,0))

datos_zinb_10_0


datos_zinb_5<-datos_zinb_5 %>%  mutate(cero=if_else(datos_zinb_5$DIASSINTOMAS.VERDADEROS==0,1,0))

datos_zinb_5_0<-datos_zinb_5 %>% filter(datos_zinb_5$cero==TRUE)
a<-dim(datos_zinb_5_0)
a  
datos_zinb_5_0<-sum(if_else(datos_zinb_5_0$DIASSINTOMAS==datos_zinb_5_0$DIASSINTOMAS.VERDADEROS,1,0))
datos_zinb_5_0

####EXPORTACION DE TABLAS##

file.choose()


tiempos<-read.csv("data/Tiempos mice.csv")


xtable(tiempos, type="latex")


file.choose()
errores<-read.csv( "data/RESUMEN IMPUTACIONES.csv")


xtable(errores, type="latex")

#setwd("data")
#write.csv(entidad_casos_sexo, file = "entidad_casos_sexo.csv", row.names = FALSE)



summary(datos$DIASSINTOMAS)
