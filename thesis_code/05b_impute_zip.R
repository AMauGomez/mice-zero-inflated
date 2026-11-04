# =============================================================================
# 05b_impute_zip.R
# Original analysis code from the thesis (UNAM, 2025). Kept as it was run;
# only absolute paths were replaced by relative ones and this header added.
# Comments are in Spanish, as written during the project.
#
# Imputes days with symptoms with the zero-inflated Poisson method of countimp (method = 'zip').
# The predictor matrix sets which covariates enter the count part and which the zero part.
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


###CRONOMETRO###

tic <- function () { 
  now <- proc.time()
  function () { 
    proc.time() - now 
  }
}

###IMPUTACION PMM 5, 10,10##

###CARGAMOS LOS DATOS ####
set.seed(1005)

datos<-read.csv("data/datos_imp_5.csv")


###SELECCIONAMOS SOLO LAS VARIABLES NUMERICAS PARA INCLUIRLOS EN EL MICE####
datos<-datos[,-(45)]

###IMPUTACION PMMM CON MATRIZ DE PRED HECHA A MANO###
##CARGAMOS LA MATRIZ DE PREDICCION##

pred_matrix_zip<-read.csv("predictor_matrices/predictor_matrix_zip.csv")
##CAMBIAMOS EL FORMATO##
rownames(pred_matrix_zip)<-pred_matrix_zip$X

pred_matrix_zip<-pred_matrix_zip[,2:45]
pred_matrix_zip<-as.matrix(pred_matrix_zip)


##hacemos la imp inicial##
set.seed(105)
ini_zip<-countimp(datos, maxit=0, print = FALSE)

###NO IMPUTAMOS DIAS SINTOMAS ESO SE HARA POR SEPARADO DEBIDO A LA RAM##

meth_ini_zip<-ini_zip$method
meth_ini_zip
meth_ini_zip[5]<-"zip"
meth_ini_zip
##imputación con ZIP##

toc <- tic()

zip_5<-countimp(datos, method=meth_ini_zip, pred= pred_matrix_zip, maxit=5, m=5)

toc()

###IMPUTACION COMPLETA###

imputacion_zip_com_5<-complete(zip_5)

##contamos LOS NA###

missing<-apply(X = is.na(imputacion_zip_com_5), MARGIN = 2, FUN = sum)
missing


missing_percent<-apply(X = is.na(imputacion_zip_com_5), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


mean(is.na(imputacion_zip_com_5))
sum(is.na(datos))


##Exportamos la imputación PMM##

# setwd("data")
write.csv(imputacion_zip_com_5, file = "imputacion_COMPLETA__zip_5.csv")


x11()
plot_trace(zip_5, c("DIASSINTOMAS"))

x11()
densityplot(zip_5, ~ DIASSINTOMAS)


##POOLED PHASE###
fit_zip_5<-with(zip_5, glmmTMB(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
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
                                 family=poisson()))

pool(fit_zip_5)


##IMPUTACIONES AGRUPADAS###
stacked <- complete(zip_5)

fit<- glmmTMB(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
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
              family=poisson(), data=stacked)
fit









###IMPUTACION 10 % 10, 15, 15 PMM##

datos_10<-read.csv("data/datos_imp_10.csv")

datos_10<-datos_10[,-(45)]

##hacemos la imp inicial##

set.seed(110)
ini_zip_10<-mice(datos_10, maxit=0, print = FALSE)

###NO IMPUTAMOS DIAS SINTOMAS ESO SE HARA POR SEPARADO DEBIDO A LA RAM##

meth_ini_zip_10<-ini_zip_10$method
meth_ini_zip_10
meth_ini_zip_10[5]<-"zip"
meth_ini_zip_10
##imputación con PMM##

toc <- tic()

pmm_10<-countimp(datos_10, method=meth_ini_zip_10, pred= pred_matrix_zip, maxit=10, m=10)


toc()
zip_10<-pmm_10
###IMPUTACION COMPLETA###

imputacion_zip_com_10<-complete(zip_10)

##contamos LOS NA###

missing<-apply(X = is.na(imputacion_zip_com_10), MARGIN = 2, FUN = sum)
missing


missing_percent<-apply(X = is.na(imputacion_zip_com_10), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


mean(is.na(imputacion_zip_com_10))
sum(is.na(datos_10))


##Exportamos la imputación PMM##

# setwd("data")
write.csv(imputacion_zip_com_10, file = "imputacion_COMPLETA_zip_10.csv")



x11()
plot_trace(zip_10, c("DIASSINTOMAS"))

x11()
densityplot(zip_10, ~ DIASSINTOMAS)



##POOLED PHASE###
fit_zip_10<-with(zip_10, glmmTMB(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
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
                                 family=poisson()))

pool(fit_zip_10)


##IMPUTACIONES AGRUPADAS###
stacked <- complete(zip_10)

fit<- glmmTMB(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
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
               family=poisson(), data=stacked)
coef(fit)

fit


###IMPUTACION 15 % 15, 15, 15 PMM##

##CARGAMOS LOS DATOS##

set.seed(120)

datos_20<-read.csv("data/datos_imp_15.csv")

datos_20<-datos_20[,-(45)]
##hacemos la imp inicial##
ini_zip_20<-countimp(datos_20, maxit=0, print = FALSE)

###INI METHOD##

meth_ini_zip_20<-ini_zip_20$method
meth_ini_zip_20
meth_ini_zip_20[5]<-"zip"
meth_ini_zip_20
##imputación con PMM##

toc <- tic()

zip_20<-countimp(datos_20, method=meth_ini_zip_20, pred= pred_matrix_zip, maxit=15, m=15)

toc()

###IMPUTACION COMPLETA###

imputacion_zip_com_20<-complete(zip_20)

##contamos LOS NA###

missing<-apply(X = is.na(imputacion_zip_com_20), MARGIN = 2, FUN = sum)
missing


missing_percent<-apply(X = is.na(imputacion_zip_com_20), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


mean(is.na(imputacion_zip_com_20))
sum(is.na(datos_20))


##Exportamos la imputación PMM##

# setwd("data")
write.csv(imputacion_zip_com_20, file = "imputacion_COMPLETA_zip_15.csv")


x11()
plot_trace(zip_20, c("DIASSINTOMAS"))

x11()
densityplot(zip_20, ~ DIASSINTOMAS)




##POOLED PHASE###






