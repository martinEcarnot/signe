library(MASS)
library(FactoMineR)
library(signal)
library(plyr)
library(caret)
library(dplyr)
library(prospectr)
library(sampling)
library(rnirs)
library(ggplot2)
library(plotly)
library("cowplot")
library("gridGraphics")
library(MetStaT)

rm(list = ls())


source('Script_R_2020/adj_asd.R')
source('Script_R_2020/SIGNE_load.R')
source('Script_R_2020/SIGNE_maha0.R')
source("Script_R_2020/sp2dfclo.R")

##Fonctions utiles :
#plotsp, fonction d'affichage du package rnirs.


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
#set.seed(1)

brb3="C:/Users/avitvale/Documents/Test/globalmatrix"
load(file=brb3)
sp=globalmatrix

L= 2*(1:996)
sp2=sp[,L]
sp=sp2[sample(1:5354,3500),] #On choisit aléatoirement 2/3 des données (sert à réduire le temps de calcul. On suppose que ca n'impact pas les résultats).

# Data Filter

### FIXATION DES PARAMETRES UTILISES:
## Nombre de repetitions de la boucle de PLSDA:
repet= 2
## Parametres du Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
p=2 #2
n=11#11 #Faire avec un n plus gros ?
m=1 #1
## Nombre de VL max autorisees
ncmax=35 #35
## Nombre de groupes de CV
k=2 #2

## PLSDA ##


### Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
sp_pre=adj_asd(sp,c(276,676)) #Retrouvé empiriquement, ne correspondait pas à ce qui etait indiqué precedemment.

## SNV
sp_pre=t(scale(t(sp_pre)))

## Derivation Savitsky Golay
sp=savitzkyGolay(sp_pre, m = m, p = p, w = n)


sp=sp[,-(665:695)] #Coupure du spectre autour des résidus de sauts de detecteurs. Ne semble pas encore effacer completement les sauts.
sp=sp[,-(250:280)]




class=as.factor(substr(rownames(sp),9,9)) #9,9
classclo=as.factor(substr(rownames(sp),9,13)) #9,13

## Variable qui mesure le nombre de classes
c=length(levels(class))
cclo=length(levels(classclo))



# Création des jeux de calibration/ validation      (Enjeu d'éviter que des spectres trop semblables soient en calibration et en validation)
# On créé un facteur datclone qui groupe un clone à 1 date
datclone=substr(rownames(sp),1,13)
ndc=length(unique(datclone))
# On créé un facteur souche qui groupe les 6 spectres de chaque souche
numsp=as.numeric(substr(rownames(sp),15,16))
souche=cut(numsp, breaks = c(0,6,12,18),labels=c("s1","s2","s3"))
endroit=as.numeric(substr(rownames(sp),15,16))%%3

rownames(sp)=paste(rownames(sp),souche,endroit)

sp=sp2dfclo(sp,class,classclo)
sp=cbind(sp,datclone,souche)
str(sp)




sp$y1=as.numeric(as.factor(sp$y1))
sp$annee=as.numeric(as.factor(substr(rownames(sp),1,4)))
sp$parcelle=as.numeric(as.factor(substr(rownames(sp),18,18)))
sp$position=as.numeric(as.factor(substr(rownames(sp),23,23)))
sp$jour=as.numeric(as.factor(substr(rownames(sp),1,8)))
sp$mois=as.numeric(as.factor(substr(rownames(sp),5,6)))



## Choix ici du sous-ensemble auquel on applique ASCA (si on le fait sur toutes les données, sauter cette étape)
sp=sp[which(sp$parcelle==4),] #G
sp=sp[which(sp$y1==1),]
#




splevels=data.frame(cepage=as.numeric(sp$y1), annee=as.numeric(sp$annee), parcelle=as.numeric(sp$parcelle), position=as.numeric(sp$position), souche=as.numeric(sp$souche) )


fact=matrix(cbind(sp[,1],sp[,2],sp[,5],sp[,6],sp[,7],sp[,8],sp[,10]), ncol=7) # 1:cépage. 2:clone. 3:souche. 4:année. 5:parcelle. 6:position. 7:mois.
#fact=matrix(cbind(sp[,1],sp[,6]), ncol=2)
#fact=as.matrix(sp[,c(1,6)])
res=ASCA.Calculate(sp$x, fact)

ASCA.DoPermutationTest(res)
