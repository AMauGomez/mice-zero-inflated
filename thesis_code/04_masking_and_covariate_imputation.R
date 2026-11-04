# =============================================================================
# 04_masking_and_covariate_imputation.R
# Original analysis code from the thesis (UNAM, 2025). Kept as it was run;
# only absolute paths were replaced by relative ones and this header added.
# Comments are in Spanish, as written during the project.
#
# Creates the evaluation sets: masks 5%, 10% (and 15%) of the observed days with symptoms
# completely at random (caTools::sample.split), keeps the true values, and imputes the missing
# values of all other covariates with PMM using a hand-built predictor matrix.
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

###
toc <- tic()

toc()

###CARGAMOS LOS DATOS ####
set.seed(999)

datos<-read.csv("data/base_covid_procesada.csv")


###SELECCIONAMOS SOLO LAS VARIABLES NUMERICAS PARA INCLUIRLOS EN EL MICE####
datos<-datos[,4:50]

datos<-datos[, -2]
datos<-datos[, - (5:6)]

###CREAMOS 5% DE  DATOS FALTANTES AL AZAR####

split=sample.split(datos$DIASSINTOMAS, SplitRatio = .05)

split<-as.vector(split)

datos<-datos %>% mutate(split=split)

####Cambiamos a NA la muestra que se eligio#

datos$DIASSINTOMAS<-ifelse(datos$split=="TRUE", "NA", datos$DIASSINTOMAS)

datos$DIASSINTOMAS<-as.numeric(datos$DIASSINTOMAS)

##Quitamos la colunma split pero la guaradamos##

datos_1<-datos[, 1:44]

##CONTAMOS LOS NA por columna##

missing<-apply(X = is.na(datos_1), MARGIN = 2, FUN = sum)
missing

missing_percent<-apply(X = is.na(datos_1), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)

###IMPUTAMOS LAS DEMAS VARIABLES POR PMM EXEPTO DIASSINTOMAS##

###IMPUTACION PMMM CON MATRIZ DE PRED HECHA A MANO###
##CARGAMOS LA MATRIZ DE PREDICCION##


pred_matrix_pmm<-read.csv("predictor_matrices/predictor_matrix_pmm.csv")

##CAMBIAMOS EL FORMATO##
rownames(pred_matrix_pmm)<-pred_matrix_pmm$X

pred_matrix_pmm<-pred_matrix_pmm[,2:45]
pred_matrix_pmm<-as.matrix(pred_matrix_pmm)



##Visualizamos la matriz de pred#

x11()

plot_pred(pred_matrix_pmm)


##hacemos la imp inicial##
ini_pmm<-mice(datos_1, maxit=0, print = FALSE)

###NO IMPUTAMOS DIAS SINTOMAS ESO SE HARA POR SEPARADO DEBIDO A LA RAM##

meth_ini_pmm<-ini_pmm$method
meth_ini_pmm
meth_ini_pmm[5]<-""
meth_ini_pmm
##imputación con PMM##

toc <- tic()

pmm<-mice(datos_1, method=meth_ini_pmm, pred= pred_matrix_pmm, maxit=12, m=10)

toc()

###IMPUTACION COMPLETA###

imputacion_pmm_def<-complete(pmm)

##contamos LOS NA###

missing<-apply(X = is.na(imputacion_pmm_def), MARGIN = 2, FUN = sum)
missing


missing_percent<-apply(X = is.na(imputacion_pmm_def), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


mean(is.na(imputacion_pmm_def))
sum(is.na(datos_1))


##Exportamos la imputación PMM##

#setwd("data")
#write.csv(imputacion_pmm_def, file = "imputacion_pmm_def_5_12_10.csv")



###PLOTEAMOS LAS VARIABLE SCON MAYOR FALTA DE DATOS##
#X11()
#densityplot(pmm, ~ INTUBADO | .imp)

x11()
plot_trace(pmm, c("INTUBADO"))

x11()
densityplot(pmm, ~ INTUBADO)



##ANOSMIA##
x11()
plot_trace(pmm, c("ANOSMIA"))

x11()
densityplot(pmm, ~ ANOSMIA)


###DISGEUSIA##

x11()
plot_trace(pmm, c("DISGEUSIA"))

x11()
densityplot(pmm, ~ DISGEUSIA)

###TXCROBIA##


x11()
plot_trace(pmm, c("TXCROBIA"))

x11()
densityplot(pmm, ~ TXCROBIA)


###RESDEFIN##

x11()
plot_trace(pmm, c("RESDEFIN"))

x11()
densityplot(pmm, ~ RESDEFIN)












##CREAMOS LA FALTA DE DATOS EN 3 DISTINTOS PORCENTAJES EN LA BASE YA IMPUTADA##

##USAMOS LA BASE YA IMPUTADA QUE LLAMAREMOS DATOS_IMP##

datos_imp<-read.csv("data/datos imputados listo para la variable sintomas.csv")


missing_percent<-apply(X = is.na(datos_imp), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


sum(is.na(datos_imp))



###IMPUTACION PMM AL 5% DE FALTA DE DATOS CON M=10 Y MAXIT=10 ###
##DIAS SINTOMAS##

set.seed(12345)

###CREAMOS 5% DE  DATOS FALTANTES AL AZAR####

datos_imp<-read.csv("data/datos imputados listo para la variable sintomas.csv")

split=sample.split(datos_imp$DIASSINTOMAS, SplitRatio = .05)

split<-as.vector(split)

datos_imp_5<-datos_imp %>% mutate(split=split)


####Cambiamos a NA la muestra que se eligio#

datos_imp_5$DIASSINTOMAS<-ifelse(datos_imp_5$split=="TRUE", "NA", datos_imp_5$DIASSINTOMAS)

datos_imp_5$DIASSINTOMAS<-as.numeric(datos_imp_5$DIASSINTOMAS)


##CONTAMOS LOS NA por columna##

missing<-apply(X = is.na(datos_imp_5), MARGIN = 2, FUN = sum)
missing

missing_percent<-apply(X = is.na(datos_imp_5), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


###GUARDAMOS LA BASE CON %5 DE FALTA DE DATOS##
# setwd("data")
write.csv(datos_imp_5, file = "datos_imp_5.csv")






###IMPUTACION PMM AL 10% DE FALTA DE DATOS CON M=10 Y MAXIT=10 ###
##DIAS SINTOMAS##

set.seed(10)

###CREAMOS 5% DE  DATOS FALTANTES AL AZAR####

datos_imp<-read.csv("data/datos imputados listo para la variable sintomas.csv")

split=sample.split(datos_imp$DIASSINTOMAS, SplitRatio = .10)

split<-as.vector(split)

datos_imp_10<-datos_imp %>% mutate(split=split)


####Cambiamos a NA la muestra que se eligio#

datos_imp_10$DIASSINTOMAS<-ifelse(datos_imp_10$split=="TRUE", "NA", datos_imp_10$DIASSINTOMAS)

datos_imp_10$DIASSINTOMAS<-as.numeric(datos_imp_10$DIASSINTOMAS)


##CONTAMOS LOS NA por columna##

missing<-apply(X = is.na(datos_imp_10), MARGIN = 2, FUN = sum)
missing

missing_percent<-apply(X = is.na(datos_imp_10), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


###GUARDAMOS LA BASE CON %5 DE FALTA DE DATOS##
# setwd("data")
write.csv(datos_imp_10, file = "datos_imp_10.csv")



set.seed(15)

###CREAMOS 15% DE  DATOS FALTANTES AL AZAR####

datos_imp<-read.csv("data/datos imputados listo para la variable sintomas.csv")

split=sample.split(datos_imp$DIASSINTOMAS, SplitRatio = .15)

split<-as.vector(split)

datos_imp_15<-datos_imp %>% mutate(split=split)


####Cambiamos a NA la muestra que se eligio#

datos_imp_15$DIASSINTOMAS<-ifelse(datos_imp_15$split=="TRUE", "NA", datos_imp_15$DIASSINTOMAS)

datos_imp_15$DIASSINTOMAS<-as.numeric(datos_imp_15$DIASSINTOMAS)


##CONTAMOS LOS NA por columna##

missing<-apply(X = is.na(datos_imp_15), MARGIN = 2, FUN = sum)
missing

missing_percent<-apply(X = is.na(datos_imp_15), MARGIN = 2, FUN = mean)

round(missing_percent*100, 2)


###GUARDAMOS LA BASE CON %5 DE FALTA DE DATOS##
# setwd("data")
write.csv(datos_imp_15, file = "datos_imp_15.csv")





###IMPUTACION PMM AL 5% DE FALTA DE DATOS CON M=10 Y MAXIT=10 ###
##DIAS SINTOMAS##


##Quitamos la colunma split pero la guaradamos##
##QUITMAOS LA NUMERACION SI ESTA DUPLICADA##

datos_imp_5<-datos[, 1:44]

###DEFINIMOS LA IMPUTACION PMM##

meth_ini_pmm<-ini_pmm$method
meth_ini_pmm
meth_ini_pmm[5]<-"pmm"
meth_ini_pmm
##imputación con PMM##

toc <- tic()

pmm<-mice(datos_imp, method=meth_ini_pmm, pred= pred_matrix_pmm, maxit=10, m=10, set.seed(105))

toc()











###APLICAMOS EL MODELO MICE CERO INFLADO POISSON####

##CARGAMOS LA MATRIZ DE PREDICCION##

pred_matrix<-read.csv("predictor_matrices/predictor_matrix_zip.csv")

##CAMBIAMOS EL FORMATO##
rownames(pred_matrix)<-pred_matrix$X

pred_matrix<-pred_matrix[,2:45]
pred_matrix<-as.matrix(pred_matrix)

pred_matrix_zip<-pred_matrix

##Visualizamos la matriz de pred#


x11()
plot_pred(pred_matrix_zip)


##hacemos la imp inicial##
ini_zip<-countimp(datos_1, maxit=0, print = FALSE)

meth_ini_zip<-ini_zip$method
meth_ini_zip

##declaramos el metodo paras dias sintomas##

meth_ini_zip[5]<-"zip"

meth_ini_zip

####IMPUTACIÓN ZIP CON LA MATRIZ DE PRED YA DEFINIDA##


##imputación con zip##

toc <- tic()

zip<-countimp(datos, method=meth_ini_zip, pred= pred_matrix_zip, maxit=1, m=1)
#datos_1
toc()

summary(zip)


#ploteamos

x11()

plot(zip, c("DIASSINTOMAS"))


x11()
plot_trace(zip, c("DIASSINTOMAS"))


imputacion_zip<-complete(zip)

##Exportamos la imputación ZIP##

#setwd("data")
#write.csv(imputacion_zip, file = "imputacion_zip_5_3_3.csv")



##Pooled phase###


#fit_zip<-with(zip, glmmTMB(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
#                        EDAD +  DISNEA + 
 #                       DIARREA + CEFALEA + ARTRAL +
  #                      RINORREA  + CIANOSIS + INISUBIS + DISGEUSIA + 
   #                    VACUNADO + RESDEFIN + ANTIPIRETICOS,
    #                  zi=~ INTUBADO + DIGCLINE + 
     #                   EDAD + DIAGPROB_IRAG + FIEBRE + TOS + ODINOGIA + DISNEA + 
      #                  DIARREA + CEFALEA + MIALGIAS + ATAEDOGE + 
       #                 RINORREA  + CONJUN + CIANOSIS + INISUBIS + 
        #                INMUSUPR + ENFCARDI + INSRENCR + TXANTIVI + 
         #               VACUNADO + RESDEFIN + ANTIPIRETICOS,
          #            family=poisson()))


#pool(fit_zip)


#xtable(a, type="latex" )

#stargazer(pool(fit_zip), type="latex" )

###Graficamos las densidades##

3x11()
#densityplot(zip$imp)

###IMPUTACION PMMM CON MATRIZ DE PRED HECHA A MANO###


##CARGAMOS LA MATRIZ DE PREDICCION##


pred_matrix_pmm<-read.csv("predictor_matrices/predictor_matrix_pmm.csv")

##CAMBIAMOS EL FORMATO##
rownames(pred_matrix_pmm)<-pred_matrix_pmm$X

pred_matrix_pmm<-pred_matrix_pmm[,2:45]
pred_matrix_pmm<-as.matrix(pred_matrix_pmm)



##Visualizamos la matriz de pred#


x11()
plot_pred(pred_matrix_pmm)


##hacemos la imp inicial##
ini_pmm<-mice(datos_1, maxit=0, print = FALSE)

###NO IMPUTAMOS DIAS SINTOMAS ESO SE HARA POR SEPARADO DEBIDO A LA RAM##

meth_ini_pmm<-ini_pmm$method
meth_ini_pmm
meth_ini_pmm[5]<-""

##imputación con PMM##

toc <- tic()

pmm<-mice(datos_1, method=meth_ini_pmm, pred= pred_matrix_pmm, maxit=3, m=3)

toc()

summary(pmm)


#ploteamos

x11()

plot(pmm, c("DIASSINTOMAS"))


x11()
plot_trace(pmm, c("DIASSINTOMAS"))


imputacion_pmm<-complete(pmm)



##Pooled phase###


#fit_pmm<-with(pmm, glmmTMB(DIASSINTOMAS ~ SEXO + INTUBADO + DIGCLINE + 
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




print(pool(fit_pmm))


##Graficamos las imputaciones###
x11()
densityplot(pmm)




