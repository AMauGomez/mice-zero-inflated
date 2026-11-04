# =============================================================================
# 05a_impute_pmm.R
# Original analysis code from the thesis (UNAM, 2025). Kept as it was run;
# only absolute paths were replaced by relative ones and this header added.
# Comments are in Spanish, as written during the project.
#
# Imputes days with symptoms with Predictive Mean Matching (mice). m = iterations = 5 for the 5%
# scenario and 10 for the 10% scenario. Trace plots and densities for convergence.
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

pred_matrix_pmm<-read.csv("predictor_matrices/predictor_matrix_pmm.csv")
##CAMBIAMOS EL FORMATO##
rownames(pred_matrix_pmm)<-pred_matrix_pmm$X

pred_matrix_pmm<-pred_matrix_pmm[,2:45]
pred_matrix_pmm<-as.matrix(pred_matrix_pmm)


##hacemos la imp inicial##
set.seed(110)
ini_pmm<-mice(datos, maxit=0, print = FALSE)

###DEFINIMOS EL METODO DE IMPUTACION##

meth_ini_pmm<-ini_pmm$method
meth_ini_pmm
meth_ini_pmm[5]<-"pmm"
meth_ini_pmm
##imputación con PMM##

toc <- tic()

pmm_5<-mice(datos, method=meth_ini_pmm, pred= pred_matrix_pmm, maxit=5, m=5, set.seed(105))

toc()

###IMPUTACION COMPLETA###

imputacion_pmm_com_5<-complete(pmm_5)

##contamos LOS NA###

missing<-apply(X = is.na(imputacion_pmm_com_5), MARGIN = 2, FUN = sum)
missing


missing_percent<-apply(X = is.na(imputacion_pmm_com_5), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


mean(is.na(imputacion_pmm_com_5))
sum(is.na(datos))


##Exportamos la imputación PMM##

# setwd("data")
write.csv(imputacion_pmm_com_5, file = "imputacion_COMPLETA_pmm_5.csv")


x11()
plot_trace(pmm_5, c("DIASSINTOMAS"))

x11()
densityplot(pmm_5, ~ DIASSINTOMAS)




###IMPUTACION 10 % 10, 15, 15 PMM##

datos_10<-read.csv("data/datos_imp_10.csv")

datos_10<-datos_10[,-(45)]

##hacemos la imp inicial##
set.seed(110)
ini_pmm<-mice(datos_10, maxit=0, print = FALSE)

###NO IMPUTAMOS DIAS SINTOMAS ESO SE HARA POR SEPARADO DEBIDO A LA RAM##

meth_ini_pmm<-ini_pmm$method
meth_ini_pmm
meth_ini_pmm[5]<-"pmm"
meth_ini_pmm
##imputación con PMM##

toc <- tic()

pmm_10<-mice(datos_10, method=meth_ini_pmm, pred= pred_matrix_pmm, maxit=10, m=10, set.seed(110))

toc()

###IMPUTACION COMPLETA###

imputacion_pmm_com_10<-complete(pmm_10)

##contamos LOS NA###

missing<-apply(X = is.na(imputacion_pmm_com_10), MARGIN = 2, FUN = sum)
missing


missing_percent<-apply(X = is.na(imputacion_pmm_com_10), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


mean(is.na(imputacion_pmm_com_10))
sum(is.na(datos))


##Exportamos la imputación PMM##

# setwd("data")
write.csv(imputacion_pmm_com_10, file = "imputacion_COMPLETA_pmm_10.csv")


x11()
plot_trace(pmm_10, c("DIASSINTOMAS"))

x11()
densityplot(pmm_10, ~ DIASSINTOMAS,  from=0,to=25)

pmm_11<-pmm_10
pmm_11<-filter(pmm_11$data$DIASSINTOMAS<30)

pmm_11<-pmm_11 %>% filter(DIASSINTOMAS<30)


summary(pmm_11)
x11()
densityplot(pmm_11, ~ DIASSINTOMAS)

