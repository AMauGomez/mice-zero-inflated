# =============================================================================
# 02_distribution_and_dispersion.R
# Original analysis code from the thesis (UNAM, 2025). Kept as it was run;
# only absolute paths were replaced by relative ones and this header added.
# Comments are in Spanish, as written during the project.
#
# Which distribution fits days with symptoms? Chi-square goodness-of-fit for Poisson and negative
# binomial, density and quantile comparisons, overdispersion tests (AER, performance) and
# zero-inflation tests with simulated residuals (DHARMa) for Poisson, NB, ZIP and ZINB.
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
library("DHARMa")
library("glmmTMB")
library("pscl")
library(fitdistrplus)
library("lme4")
library(boot)

datos<-read.csv("data/base_covid_procesada.csv")

###SELECCIONAMOS SOLO LAS VARIABLES NUMERICAS PARA INCLUIRLOS EN EL MICE####
datos<-datos[,4:50]

datos<-datos[, -2]
datos<-datos[, - (5:6)]

###Histograma varaible DIASSINTOMAS####

x11()
h_1<-hist(datos$DIASSINTOMAS, xlim = c(0,25), ylim = c(0,100000), breaks = 10000)

summary(datos$DIASSINTOMAS)

mean(datos$DIASSINTOMAS)
var(datos$DIASSINTOMAS)

###Contamos los 0####

ceros<-sum(datos$DIASSINTOMAS==0)

porcentaje_ceros<-(ceros/nrow(datos))*100
porcentaje_ceros

### Revisamos a que distribución se ajustan mejor BINOMIAL O POISSON####

### PRUEBA PARA POISSON ####
poi<-fitdist(datos$DIASSINTOMAS, "pois")

summary(poi)

##Hacemos una prueba de hipotesis

poi_fit<-gofstat(poi)
poi_fit$chisqpvalue

#Se rechaza H_0, por lo tanto no se distribuye Poison 

x11()
plot(poi)

##Notamos que las densidades no se asemejan###

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
               xlab = "Cuantiles teóricos",
               ylab = "Cuantiles empíricos",
               fitlty = 1,
               fitcol = c("red1","blue"),
               plotstyle = "ggplot",
               addlegend = TRUE)
x11()
comp_qq


comp_den<-denscomp(list(poi, bino), 
         xlegend = "bottomright",
         fitlty = 1,
         fitcol = c("orange","blue"),
         plotstyle = "ggplot",
         addlegend = TRUE)
x11()
comp_den

#Aunque con nuestos datos se acercan a la densidad de la dist Binomial Negativa
#podemos notar que el QQ-plot no 

###Creamos los modelos poisson para identificar sobredispersión####
#Ya que sabemos que los datos son de recuento 
##La selección de variables se hizo de acuerdo al PCA, COrrelación y Experiencia

m_pois_visual<-glm(DIASSINTOMAS~ANTIPIRETICOS + FIEBRE +ARTRAL+TOS + MIALGIAS +CALOFRIOS + DISNEA+ TXCROBIA+DIAGPROB_IRAG + RESDEFIN+ DIGCLINE + INISUBIS + ATAEDOGE + DOTORACI + ODINOGIA, family = poisson(), data=datos)

summary(m_pois_visual)


m_pois_com<-glm(DIASSINTOMAS~., family = poisson(), data=datos)
summary(m_pois_com)


#step(m_pois_com, direction = "both")


###DISPERSIONES  POISSON AER####
dis_visual<-dispersiontest(m_pois_visual)
dis_visual
dis_com<-dispersiontest(m_pois_com)
dis_com

###DISPERSIONES POISSON PERFORMANCE####

dis_visual_2<-check_overdispersion(m_pois_visual)
dis_visual_2
dis_com_2<-check_overdispersion(m_pois_com)
dis_com_2


###EVALUAMOS EL MODELO POISSON###
#Gráficas
x11()
plot(m_pois_visual)
x11()
plot(m_pois_com)

##Validación cruzada

cv_pois_vis<-cv.glm(datos, m_pois_visual, k=10)




###Creamos los modelos binomial negativo para identificar sobredispersión####

m_binom_visual<-glm.nb(DIASSINTOMAS~ANTIPIRETICOS + FIEBRE +ARTRAL+TOS + MIALGIAS +CALOFRIOS + DISNEA+ TXCROBIA+DIAGPROB_IRAG + RESDEFIN+ DIGCLINE + INISUBIS + ATAEDOGE + DOTORACI + ODINOGIA, data=datos)

summary(m_binom_visual)


m_binom_com<-glm.nb(DIASSINTOMAS~., data=datos)
summary(m_binom_com)

###DISPERSIONES  BINOM CON  AER NO SE PUEDEN HACRE####

###DISPERSIONES BINOM PERFORMANCE####

dis_visual_2_nb<-check_overdispersion(m_binom_visual)
dis_visual_2_nb
dis_com_2_nb<-check_overdispersion(m_binom_com)
dis_com_2_nb



###Creamos los modelos para verificar si esta cero inflado####
###CERO INFLADO POISSON####

zip_visual<-glmmTMB(DIASSINTOMAS~ANTIPIRETICOS + FIEBRE +ARTRAL+TOS + MIALGIAS +CALOFRIOS + DISNEA+
                      TXCROBIA+DIAGPROB_IRAG + RESDEFIN+ DIGCLINE + INISUBIS + ATAEDOGE + DOTORACI + ODINOGIA,
                    zi=~ANTIPIRETICOS + FIEBRE +ARTRAL+TOS + MIALGIAS +CALOFRIOS + DISNEA+
                      TXCROBIA+DIAGPROB_IRAG + RESDEFIN+ DIGCLINE + INISUBIS + ATAEDOGE + DOTORACI + ODINOGIA,
                    family=poisson(), data=datos)

summary(zip_visual)


##HACER EL COMPLETO CON LA LIBRERIA PSIC##
#zip_com<-glmmTMB(DIASSINTOMAS~ ANTIPIRETICOS + FIEBRE +ARTRAL+TOS , zi=~ ANTIPIRETICOS + FIEBRE +ARTRAL+TOS , family=poisson(), data = datos_3)


# nbin <- glmmTMB(Bats ~ Avg.Biomass + Num.Orders + Avg.Understory + Avg.Midstory + 
#Avg.Canopy.Cover + Perc.Dec.Dom + Avg.Snags + Avg.Bat.Date + Avg.Bat.Temp +
# Bat.Elevation + Bat.Moon + Bat.Water.Feat + Avg.Biomass + Num.Orders + 
#Avg.Bat.Temp*Avg.Bat.Date + Avg.Biomass*Year + Year + Habitat + 
#offset(log(Survey.Nights)) + (1|Site), 
#data = insect.data, 
#ziformula = ~0, 
#family = nbinom2)

#summary(zip_com)

####PRUEBAS DE CERO INFLADO POISSON####

##DHARMA#3
#simulamos los residuales del modelos glm#

residu_pois_vis<-simulateResiduals(m_pois_visual)

x11()
testZeroInflation(residu_pois_vis)

##El modelo esta inflado por ceros


#simulamos los residuales del ZIP#

residu_zip_vis<-simulateResiduals(zip_visual)
x11()
testZeroInflation(residu_zip_vis)

#Notamos que es muy bueno ya que no esta inflando ceros#


###PERFORMANCE##

check_zeroinflation(m_pois_visual)

##Los modelos cero inflados hacer la prueba con las simulaciones
check_zeroinflation(residu_zip_vis)

###REVISAMOS SI LOS MODELOS ZIP MEJORAN LA SOBREDISPERSIÓN####
##Hay que usar la simulaciones de residuos#
check_overdispersion(residu_zip_vis)

#NOTAMOS QUE MEJORA LA SOBREDISPERSIÓN DE dis_visual, pero sigue existiendo sobrediapersión##

###Creamos los modelos binomiasles para verificar si esta cero inflado####



zinb_visual<-glmmTMB(DIASSINTOMAS~ANTIPIRETICOS + FIEBRE +ARTRAL+TOS + MIALGIAS +CALOFRIOS + DISNEA+
                       TXCROBIA+DIAGPROB_IRAG + RESDEFIN+ DIGCLINE + INISUBIS + ATAEDOGE + DOTORACI + ODINOGIA,
                     zi=~ANTIPIRETICOS + FIEBRE +ARTRAL+TOS + MIALGIAS +CALOFRIOS + DISNEA+
                       TXCROBIA+DIAGPROB_IRAG + RESDEFIN+ DIGCLINE + INISUBIS + ATAEDOGE + DOTORACI + ODINOGIA,
                     family=nbinom2(), data=datos)

summary(zinb_visual)


###CERO INFLADO BINOMIAL####


##DHARMA#3
#simulamos los residuales del modelos glm#

residu_nb_vis<-simulateResiduals(m_binom_visual)
x11()
testZeroInflation(residu_nb_vis)

##El modelo nb de glm esta inflado por ceros tiene radio=1.24 aceptable ##


##simulamos los residuales del ZINB#

residu_zinb_vis<-simulateResiduals(zinb_visual)

x11()
testZeroInflation(residu_zinb_vis)

#Notamos que es muy bueno ya que no esta inflando ceros#


###PERFORMANCE##
#con modelo glm#
check_zeroinflation(m_binom_visual)

##Los modelos cero inflados hacer la prueba con las simulaciones
check_zeroinflation(residu_zinb_vis)

###REVISAMOS SI LOS MODELOS ZINB MEJORAN LA SOBREDISPERSIÓN####
##Hay que usar la simulaciones de residuos#
check_overdispersion(residu_zinb_vis)

#NOTAMOS QUE MEJORA YA QUE NO HAY SOBREDISPERSIÓN NI ESTA CERO INFLADO ES EL MODELO IDEAL##
###El MODELO ZIBN ES EL MEJOR####


























