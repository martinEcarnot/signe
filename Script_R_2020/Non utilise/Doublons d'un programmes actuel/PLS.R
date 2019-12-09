library(pls)
library(MASS)
library(FactoMineR)
library(signal)
library(plyr)
library(dplyr)
library(caret)

rm(list = ls())

source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha0.R')


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
# set.seed(1)


#"C:/Users/Noémie/Desktop/SFE/Caracteristiques_agro/globalmatrix"
brb="C:/Users/avitvale/Documents/matrice toutes dates/Globalmatrix/globalmatrix"
load(file=brb)

## Filtrage des spectres aberrants
# globalmatrix=globalmatrix[globalmatrix[,500]>0.6,]
# globalmatrix=globalmatrix[globalmatrix[,1]<0.2,]
# globalmatrix=globalmatrix[globalmatrix[,2000]<0.25,]

# Data Filter
# Select dates
dates=list(
  #"20180619T"
  #,"20180619P"
  #,"20180619N"
  # ,"20180627T"
  # ,"20180627P"
  # ,"20180627N"
  # ,"20180704P"
  # ,"20180709T"
  # ,"20180709P"
  # ,"20180709N"
  # ,"20180710P"
  # ,"20180816T"
  # ,"20180816P"
  # ,"20180816N"
  # ,"20180817T"
  # ,"20180817P"
  # ,"20180817N"
  # ,"20170524P"
  # ,"20170529P"
  # ,"20170606P"
  # ,"20170612P"
  # ,"20170619P"
  # ,"20170626P"
  # ,"20170703P"
  # ,"20170710P"
  # ,"20170717P"
  # ,"20170724P"
  # ,"20170731P"
  # ,"20180823P"
  # ,"20180823T"
  # ,"20180823N"
  # ,"20180731P"
  # ,"20180731T"
  # ,"20180731N"
  # ,"20180810P"
  # ,"20180810T"
  # ,"20180810N"
  # ,"20180724P"
  #"20180724T"
  "2010809T"
  # ,"20180724N"
  # "20180816X"
  # ,"20180816E"
  # ,"20180710X"
  # ,"20180710E"
  # "20180816A"
  # ,"20180816E"
  # ,"20180731A"
  # ,"20180731E"
  
  
)

Esca <- read.table("C:/Users/Noémie/Desktop/Test_poids_ok/Esca.csv",
                         header=TRUE, sep=";",dec=".",row.names=1, check.names=FALSE,
                         fileEncoding="latin1")
### FIXATION DES PARAMETRES UTILISES:
## Nombre de repetitions de la boucle de PLSDA:
repet=4
## Parametres du Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
p=2
n=11
m=1
## Nombre de VL max autorisees
ncmax=10
## Taille de l'echantillon de validation (1/v):
#v=3
## Nombre de groupes de CV
k=5
##Variable pour la moyenne des spectres
z=6
moy= 220

mean1_final=matrix(nrow = 1, ncol = ncol(globalmatrix))
mean2_final=matrix(nrow = 220, ncol = ncol(globalmatrix))

mean1_final=colMeans(globalmatrix[(1:z),])

for(j in 1:moy) {
mean2=colMeans(globalmatrix[((z+1):(z+6)),])
z=z+6
#print(mean2)
#print(z)

mean2_final[j,]=mean2
}

sp=rbind(mean1_final,mean2_final)

### Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
sp=adj_asd(sp,c(602,1402))
## Reduction des variables (extremites bruitees)
sp=sp[,seq(51,ncol(sp)-30,1)]
## SNV
sp=t(scale(t(sp)))
## Derivation Savitsky Golay
sp=t(apply(sp,1,sgolayfilt,p=p,n=n,m=m))

##Reunion des spectres et des valeurs de reference
rownames(sp)=rownames(Esca)

spok=data.frame(Esca)
spok$sp = sp

x=substr(rownames(spok),1,9)=="MPG 10/33"
spok=spok[(x==FALSE),]

RMSEP_final=matrix(nrow = repet, ncol = ncmax)
# creation de la matrice des dv et perok maximaux
maxi_final=matrix(nrow= repet, ncol = 2)


#for(j in 1:repet) {
# creation des jeux d'apprentissage et validation
flds <- createFolds(1:nrow(spok), k = k)
predm0=as.data.frame(matrix(nrow = nrow(spok), ncol = ncmax))

# Boucle sur les groupes de CV
# for (i in 1:k) {
  id_val=sort(unlist(flds[1])) #identifiants du jeu de validation sous forme de liste
  spok_val=spok[id_val,]        # matrice du jeu de validation
  spok_cal=spok[-id_val,]      #matrice du jeu de calibration composée de tout ce qui n'est pas en validation

## PLS and application to have loadings and scores
 rpls=pls:: plsr(Degré~sp, ncomp = 10, data = spok_cal, validation ="CV")
RMSEP=pls::RMSEP(rpls, estimate="all", newdata = spok_val)
print(RMSEP)
R2=pls::R2(rpls,estimate="all", newdata =spok_val)
print(R2)

#RMSEP_final[j,]=RMSEP

#}

 #plot(RMSEP, legendpos = "topright")
plot(R2)
plot(rpls, ncomp =2, asp = 1, line = TRUE)
predplot(rpls, ncomp =2, asp = 1, line = TRUE, newdata = spok_val)
print(mvrValstats(rpls,estimate = "test",newdata=spok_val, ncomp = 2))
print(mvrValstats(rpls,estimate = "train",newdata=spok_val, ncomp = 2))
 # predm0=  predict(rpls,ncomp = , newdata = spok_val)
 # plot(predm0)




