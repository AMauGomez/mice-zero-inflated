# =============================================================================
# 00_data_processing.R
# Original analysis code from the thesis (UNAM, 2025). Kept as it was run;
# only absolute paths were replaced by relative ones and this header added.
# Comments are in Spanish, as written during the project.
#
# Data processing of the SISVER extract: patient ID, recoding of 'SE IGNORA' to NA, variable-by-
# variable cleaning (sex, intubation, clinical diagnosis, age, days with symptoms, lab result,
# antipyretics), conversion of yes/no items to binary, and the correlation matrix used for
# variable selection.
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


### Datos originales 

#set.seed(123)
#datos<-read.csv("data/BASE_COMPLETA (1) cs.csv")

#datos<- as.data.frame(datos)

#datos[datos== ''] <- NA
#datos[datos== 'SE IGNORA'] <- NA

#datos_1<-datos[,1:30]
#datos_2<-datos[, 31:60]
#datos_3<-datos[,61:90]
#datos_4<-datos[, 91:113]

#vis_miss(datos_1, warn_large_data=FALSE)
#vis_miss(datos_2, warn_large_data=FALSE)
#vis_miss(datos_3, warn_large_data=FALSE)
#vis_miss(datos_4, warn_large_data=FALSE)

#summary(datos)


##PROCESAMIENTO DE BASE DATOS DE 90 VARIABLES (base_reducida)####

datos<-read.csv("data/Base reducida.csv")
datos<- as.data.frame(datos)

###Camabiamos "se ignora" a "NA"
datos[datos== ''] <- NA
datos[datos== 'SE IGNORA'] <- NA

#resumen<-describe(datos)
#resumen

####PROCESAMIENTO VARIABLE POR VARIABLE#####

###CREAMOS EL ID###
datos <- datos %>% mutate(ID = row_number(), .before=CVEENTUNI)

###CVEENTIUNI####
datos$CVEENTUNI<-as.numeric(datos$CVEENTUNI)


###SEXO####
#Contamos HOMBRES Y MUJERES EN LA MUESTRA

hombres_1_0<-ifelse(datos$SEXO=="MASCULINO", 1,0)
hombres_total<-sum(hombres_1_0)

hombres_total


mujeres_0_1<-ifelse(datos$SEXO=="FEMENINO", 1, 0)
mujeres_total=sum(mujeres_0_1)
mujeres_total


###POBLACIÓN TOTAL####

poblacion<-hombres_total+mujeres_total
poblacion

###Cambio de variables a binario o factores#####


datos$SEXO<-if_else(datos$SEXO=='MASCULINO',1,0)

hombres<-sum(datos$SEXO==1)
hombres


###EVOLUCIÓN  (Falta ordenarlos para hacer logit ordinal)#####

datos$EVOLUCI=factor(datos$EVOLUCI)


###INTUBADO#####

sum(is.na(datos$INTUBADO))

datos$INTUBADO<-if_else(datos$INTUBADO=="NO",0, 1)
datos$INTUBADO<-as.numeric(datos$INTUBADO)


###DIGCLINE#####

datos$DIGCLINE<-if_else(datos$DIGCLINE=="NO", 0, 1)

###EDAD####
datos$EDAD<-as.numeric(datos$EDAD)

###DIASSINTOAMAS#####
datos$DIASSINTOMAS<-as.numeric(datos$DIASSINTOMAS)

###DIAGPROB####
## cambiamos la  variable de 2 factores a una binaria
#donde 1 es IRAG 
datos$DIAGPROB<-factor(datos$DIAGPROB)

datos$DIAGPROB<-if_else(datos$DIAGPROB=="INFECCION RESPIRATORIA AGUDA GRAVE (IRAG)", 1, 0)
datos<-rename(datos, DIAGPROB_IRAG=DIAGPROB)

###CONVERSION DE VARIABLES DICOTOMICAS 13-48 A VARIABLES BINARIAS ####

datos_sub<-datos[, 13:48]


variables_si_no<-function(x){
  result<-ifelse(x=="NO", 0, 1)
  return(result)
}


datos_sub_0_1<-apply(datos_sub, 2, variables_si_no)
datos_sub_0_1<-as.data.frame(datos_sub_0_1)


datos[,13:48]<-datos_sub_0_1

###RESDEFIN####

datos$RESDEFIN<-factor(datos$RESDEFIN)

#Contamos las observaciones por factor

conteo_RESDEFIN<-table(datos$RESDEFIN)
conteo_RESDEFIN

#A las categorías de RESDEFIN

#RECHAZADA
#NO RECIBIDA
#NO ADECUADO
#NO AMPLIFICO
#NO SUBTIPIFICADO

#le pondremos NA como ejercicio para Usar el MICE 

datos[datos== 'RECHAZADA'] <- NA
datos[datos== 'NO RECIBIDA'] <- NA
datos[datos== 'NO ADECUADO'] <- NA
datos[datos== 'NO AMPLIFICO'] <- NA
datos[datos== 'NO SUBTIPIFICADO'] <- NA

#VERIFICAMOS LOS CONTEOS PARA VER QUE SIGA IGUAL 
conteo_RESDEFIN_NA<-table(datos$RESDEFIN)
conteo_RESDEFIN_NA

#A la variable RESDEFIN la convertimos en binaria con 1 para SARSCOV2 que es lo 
#que nos interesa para hacer una matriz de correlación

datos$RESDEFIN<-if_else(datos$RESDEFIN=="SARS-CoV-2", 1, 0)

### ANTIPIRETICOS####

datos$ANTIPIRETICOS<-if_else(datos$ANTIPIRETICOS=="NO", 0, 1)

###Creamos una matriz de corrrelaciones para las variables 0-1 DE DATOS_1####

SEXO<-datos[, "SEXO"]

INTUBADO<-datos[, "INTUBADO"]

DIGCLINE<-datos[, "DIGCLINE"]


variables_binarias<- cbind(SEXO,INTUBADO, DIGCLINE, datos[,12:50])

cor_matrix_1<-cor(variables_binarias, use="pairwise.complete.obs")

cor_matrix_1<-as.data.frame(cor_matrix_1)


###EXPORTAMOS LA BASE Y LA MATRIZ DE CORRELACIONES####

#setwd("data")
#write.csv(datos, file = "base_covid_procesada.csv", row.names = FALSE)


#write.csv(cor_matrix_1, file = "matriz_de_correlación_procesada.csv", row.names = FALSE)









