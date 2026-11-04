# =============================================================================
# 01_exploratory_analysis.R
# Original analysis code from the thesis (UNAM, 2025). Kept as it was run;
# only absolute paths were replaced by relative ones and this header added.
# Comments are in Spanish, as written during the project.
#
# Exploratory analysis of the SISVER extract: recoding of 'SE IGNORA' and invalid lab results to NA,
# missing-data patterns, epidemic curves and descriptive statistics of days with symptoms.
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
library(psych)
library("xtable")
library(ggplot2)
library("ggpmisc")

datos_VIS_MIS<-read.csv("data/BASE_COMPLETA vis mis.csv")

datos_VIS_MIS<- as.data.frame(datos_VIS_MIS)

datos_VIS_MIS[datos_VIS_MIS== ''] <- NA
datos_VIS_MIS[datos_VIS_MIS== 'SE IGNORA'] <- NA

#A las categorías de RESDEFIN

#RECHAZADA
#NO RECIBIDA
#NO ADECUADO
#NO AMPLIFICO
#NO SUBTIPIFICADO

#le pondremos NA como ejercicio para Usar el MICE 

datos_VIS_MIS[datos_VIS_MIS== 'RECHAZADA'] <- NA
datos_VIS_MIS[datos_VIS_MIS== 'NO RECIBIDA'] <- NA
datos_VIS_MIS[datos_VIS_MIS== 'NO ADECUADO'] <- NA
datos_VIS_MIS[datos_VIS_MIS== 'NO AMPLIFICO'] <- NA
datos_VIS_MIS[datos_VIS_MIS== 'NO SUBTIPIFICADO'] <- NA


#datos_base_reducida<-read.csv("data/Base reducida.csv")

#Dividimos la base para observarla de mejor forma#

datos_VIS_MIS_1<-datos_VIS_MIS[,1:30]

datos_VIS_MIS_2<-datos_VIS_MIS[,31:61]


##Visualizamos la base_vis_mis###
#Graficamos los datos##

#x11()
#vis_mis_1<-vis_miss(datos_VIS_MIS_1, warn_large_data=FALSE)
#vis_mis_1
#x11()
#vis_mis_2<-vis_miss(datos_VIS_MIS_2, warn_large_data=FALSE)
#vis_mis_2


##contamos los NA##
sum(is.na(datos_VIS_MIS))

###Algunos datos de analisis descriptivo####

##entidad y edad##

entidad_casos_edad<-describeBy(datos_VIS_MIS$EDAD, datos_VIS_MIS$ENTIDAD, mat = TRUE) 
entidad_casos_edad<-as.data.frame(entidad_casos_edad)

# setwd("data")
#write.csv(entidad_casos_edad, file = "entidad_casos_edad.csv", row.names = FALSE)


###entidad y sexo##

datos_VIS_MIS$SEXO<-ifelse(datos_VIS_MIS$SEXO=="MASCULINO", 1,0)

entidad_casos_sexo<-describeBy(datos_VIS_MIS$SEXO, datos_VIS_MIS$ENTIDAD, mat=TRUE)
entidad_casos_sexo


# setwd("data")
#write.csv(entidad_casos_sexo, file = "entidad_casos_sexo.csv", row.names = FALSE)


###entidad USMI##
#remove(datos_VIS_MIS$`if_else(datos_VIS_MIS$ORIGEN == "USMI", 1, 0)`)

datos_VIS_MIS<-mutate(datos_VIS_MIS, USMI=if_else(datos_VIS_MIS$ORIGEN=="USMI", 1, 0))

entidad_USMI<-describeBy(datos_VIS_MIS$USMI, datos_VIS_MIS$ENTIDAD, mat=TRUE)
entidad_USMI<-as.data.frame(entidad_USMI)
entidad_USMI

#write.csv(entidad_USMI, file = "entidad_USMI.csv", row.names = FALSE)

##ENTIDAD HOSPITALIZADO##

datos_VIS_MIS<-mutate(datos_VIS_MIS, hospitalizado=if_else(datos_VIS_MIS$TIPACIEN=="HOSPITALIZADO", 1,0))
entidad_hospi<-describeBy(datos_VIS_MIS$hospitalizado, datos_VIS_MIS$ENTIDAD, mat=TRUE)
entidad_hospi<-as.data.frame(entidad_hospi)
entidad_hospi

#write.csv(entidad_hospi, file = "entidad_hospi.csv", row.names = FALSE)

##entidad y resdefin###

datos_VIS_MIS$RESDEFIN<-ifelse(datos_VIS_MIS$RESDEFIN=="SARS-CoV-2", 1,0)

entidad_resdefin<-describeBy(datos_VIS_MIS$RESDEFIN, datos_VIS_MIS$ENTIDAD, mat = TRUE)

entidad_resdefin<-as.data.frame(entidad_resdefin)


# setwd("data")
#write.csv(entidad_resdefin, file = "entidad_resdefin.csv", row.names = FALSE)

##entidad y evolución# defunciones##

datos_VIS_MIS<-mutate(datos_VIS_MIS, defuncion=if_else(datos_VIS_MIS$EVOLUCI=="DEFUNCION", 1,0))

entidad_defuncion<-describeBy(datos_VIS_MIS$defuncion, datos_VIS_MIS$ENTIDAD, mat=TRUE)
entidad_defuncion<-as.data.frame(entidad_defuncion)

total_defu<-sum(datos_VIS_MIS$defuncion)
total_defu
##defunciones po covid##

datos_VIS_MIS<-mutate(datos_VIS_MIS, defu_covid=if_else(datos_VIS_MIS$RESDEFIN==1 & datos_VIS_MIS$defuncion==1, 1, 0 ))
entidad_defuncion_covid<-describeBy(datos_VIS_MIS$defu_covid, datos_VIS_MIS$ENTIDAD, mat=TRUE)

entidad_defuncion_covid<-as.data.frame(entidad_defuncion_covid)
entidad_defuncion_covid
# setwd("data")
#write.csv(entidad_defuncion_covid, file = "entidad_defuncion_covid.csv", row.names = FALSE)

total_defu_COVID<-sum(datos_VIS_MIS$defu_covid, na.rm =TRUE)
total_defu_COVID

table(is.na(datos_VIS_MIS$defu_covid))


##entidad y diagnostico de neumonia##

datos_VIS_MIS<-mutate(datos_VIS_MIS, neumo=if_else(datos_VIS_MIS$DIGCLINE=="SI", 1, 0))

entidad_neumonia<-describeBy(datos_VIS_MIS$neumo, datos_VIS_MIS$ENTIDAD, mat=TRUE)
entidad_neumonia<-as.data.frame(entidad_neumonia)
entidad_neumonia

#write.csv(entidad_neumonia, file = "entidad_neumonia.csv", row.names = FALSE)

##hacer una estadistica de los díassintomas##

entidad_diassintomas<-describeBy(datos_VIS_MIS$DIASSINTOMAS, datos_VIS_MIS$ENTIDAD, mat=TRUE)
entidad_diassintomas<-as.data.frame(entidad_diassintomas)
entidad_diassintomas

#write.csv(entidad_diassintomas, file = "entidad_diassintomas.csv", row.names = FALSE)

##entidad IRAG##

datos_VIS_MIS<-mutate(datos_VIS_MIS, IRAG=if_else(datos_VIS_MIS$DIAGPROB=="INFECCION RESPIRATORIA AGUDA GRAVE (IRAG)", 1, 0))

entidad_irag<-describeBy(datos_VIS_MIS$IRAG, datos_VIS_MIS$ENTIDAD, mat=TRUE)
entidad_irag<-as.data.frame(entidad_irag)
entidad_irag


#write.csv(entidad_irag, file = "entidad_irag.csv", row.names = FALSE)

##Cargamos las tablas ordenadas##
#tabla_1<-read.csv("data/entidad_casos_edad - copia_1.csv")


#tabla_2<-read.csv("data/entidad_casos_edad - copia - 2.csv")
#tabla_2<-tabla_2[1:33,]

#tabla_1_latex<-xtable(tabla_1, type="latex")

#tabla_1_latex
#tabla_2_latex<-xtable(tabla_2, type="latex")
#tabla_2_latex


##hacer boxplot###

##Tipo de paciente##
#correr unicamente la base de datos para la siguiente sección##

##Serie de tiempo##
# Datos con fechas y variables. La columna 'date' (fecha) es de clase "Date"

fecha_ingreso<-as.Date(datos_VIS_MIS$FECINGRE)

fecha_ini_sintom<-as.Date(datos_VIS_MIS$FECINISI)

RESDEFIN<-datos_VIS_MIS$RESDEFIN

##Serie de casos confirmados##

datos_serie<-datos_VIS_MIS%>% mutate(Fecha=as.Date(fecha_ingreso)) %>% group_by(Fecha) %>% summarise("Casos Confirmados de SARS-COV-2 por Día"=sum(RESDEFIN, na.rm = TRUE))

serie_1<-ggplot(datos_serie, aes(x = Fecha, y = `Casos Confirmados de SARS-COV-2 por Día`)) +
  geom_line()+
  geom_smooth(se = T)+
  scale_x_date(date_breaks = "2 weeks",   limit = c(as.Date("2020-01-01"), as.Date("2020-06-29"))) +
  stat_peaks(geom = "point", span = 15, color = "red", size = 4) +
  stat_peaks(geom = "label", span = 15, color = "steelblue3", angle = -25,
             hjust = 0, x.label.fmt = "%Y-%m-%d") +
  stat_peaks(geom = "rug", span = 15, color = "red", sides = "b")+
  theme(plot.title = element_text(size=22))+
  theme_gray(base_size = 13)+
  ggtitle("Serie de Tiempo Casos Confirmados SARS-COV-2")

x11()

serie_1

##serie casos acumulados##

datos_serie<-datos_serie %>% mutate("Casos Confirmados Acumulados de SARS-COV-2 por Día"=cumsum(datos_serie$`Casos Confirmados de SARS-COV-2 por Día`))

serie_acumulados<-ggplot( datos_serie, aes(x=Fecha, y=`Casos Confirmados Acumulados de SARS-COV-2 por Día`)) +
  geom_line()+
  geom_smooth(se = T)+
  scale_x_date(date_breaks = "2 weeks",   limit = c(as.Date("2020-01-01"), as.Date("2020-06-29"))) +
  stat_peaks(geom = "point", span = 19, color = "red", size = 4) +
  stat_peaks(geom = "rug", span = 19, color = "red", sides = "b")+
  theme(plot.title = element_text(size=22))+
  theme_gray(base_size = 13)+
  ggtitle("Serie de Tiempo Casos Confirmados Acumulados SARS-COV-2")

x11()

serie_acumulados



##Graficas descriptivas##
datos<-datos_VIS_MIS

datos$OCUPACIO<-as.factor(datos$OCUPACIO)

Ocupación<-table(datos$OCUPACIO)
Ocupación<-as.data.frame(Ocupación)

Ocupación$Var1<-if_else(Ocupación$Var1=="COMERCIANTES DE MERCADOS FIJOS O AMBULANTES", "COMERCIANTES", Ocupación$Var1)
Ocupación$Var1<-if_else(Ocupación$Var1=="GERENTES O PROPIETARIOS DE EMPRESAS O NEGOCIOS", "EMPRESARIOS", Ocupación$Var1)
Ocupación$Var1<-if_else(Ocupación$Var1=="JUBILADO / PENSIONADO", "JUBILADO",Ocupación$Var1)
Ocupación$Var1<-if_else(Ocupación$Var1=="OTROS TRABAJADORES DE LA SALUD", "OTROS SALUD", Ocupación$Var1 )
Ocupación$Var1<-if_else(Ocupación$Var1=="OTROS PROFESIONISTAS", "PROFESIONISTAS", Ocupación$Var1 )

# Basic Barplot
x11()

rm(Ocupacion_barplot)

par(mar=c(8,4,4,4))


Ocupacion_barplot <- barplot(Ocupación$Freq, border=F , names.arg=Ocupación$Var1 , cex.names = .85,
                  las=2 , 
                  col=c("steelblue", "lightgreen", "grey", "lightblue"), 
                  ylim=c(0,135500) , 
                  main="Ocupación de Pacientes Sospechosos a SARS-COV-2" )

#text(Ocupacion_barplot, Ocupación$Freq , paste("n: ", Ocupación$Freq, sep="") , cex=.8) 

##boxplot de sintomas#

datos<-read.csv("data/Base reducida.csv")


###CONVERSION DE VARIABLES DICOTOMICAS 13-48 A VARIABLES BINARIAS ####

datos_sub<-datos[, 12:47]


variables_si_no<-function(x){
  result<-ifelse(x=="NO", 0, 1)
  return(result)
}

datos_sub_0_1<-apply(datos_sub, 2, variables_si_no)
datos_sub_0_1<-as.data.frame(datos_sub_0_1)
datos_sub_0_1<-na.omit(datos_sub_0_1)


promedio_sintomas<-colMeans(datos_sub_0_1)

promedio_sintomas<-100*promedio_sintomas

promedio_sintomas<-as.data.frame(promedio_sintomas)

nombres_var_sin<-as.vector(colnames(datos_sub_0_1))

rm(sintomas_barplot)

sintomas_barplot <- barplot(promedio_sintomas$promedio_sintomas, border=F , names.arg=nombres_var_sin , cex.names = 1,
                             las=2 , 
                             col=c("steelblue", "lightgreen", "grey", "lightblue"), 
                             ylim=c(0,80) , 
                             main="Porcentaje de Pacientes con Síntomas Relativos a SARS-COV-2")


#restringir los dias con sintomas#

