library(pls)
library(MASS)
library(mixOmics)
library(FactoMineR)
library(signal)
library(plyr)
library(caret)
library(dplyr)
library(nirs)
library(MetStaT)

rm(list = ls())

#source('C:/Users/Noémie/Desktop/SFE/Script_R/adj_asd.R')
#source('C:/Users/Noémie/Desktop/SFE/Script_R/SIGNE_load.R')
#source('C:/Users/Noémie/Desktop/SFE/Script_R/SIGNE_maha.R')
#source('C:/Users/Noémie/Desktop/SFE/Script_R/SIGNE_maha0.R')

source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha0.R')


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
# set.seed(1)

brb="C:/Users/avitvale/Documents/matrice toutes dates/Globalmatrix/globalmatrix"
load(file=brb)


# Data Filter
# Select dates
dates=list(
  #"20180619T"
  "20180619P"
  #,"20180619N"
  # ,"20180627T"
   ,"20180627P"
  # ,"20180627N"
   ,"20180704P"
  # ,"20180709T"
   ,"20180709P"
  # ,"20180709N"
  # ,"20180710P"
  # ,"20180816T"
#   ,"20180816P"
  # ,"20180816N"
  # ,"20180817T"
#   ,"20180817P"
  # ,"20180817N"
#   ,"20170524P"
 #  ,"20170529P"
 #  ,"20170606P"
  # ,"20170612P"
   ,"20170619P"
   ,"20170626P"
   ,"20170703P"
   ,"20170710P"
  # ,"20170717P"
  # ,"20170724P"
  # ,"20170731P"
   #,"20180823P"
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
  #"2010809T"
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

iok=substr(rownames(globalmatrix),1,9) %in% dates
globalmatrix=globalmatrix[iok,]

# #Parametres du Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
# p=2
# n=11
# m=1
# # nombre de DV max autorisees
# ncmax=15
# #Taille de l'echantillon de validation (1/v):
# #v=3
# # Nombre de groupes de CV
# k=6

## LDA ##
sp=globalmatrix

## Pretraitements
# # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651??me l.o.) et 1800 (1451))
sp=adj_asd(sp,c(602,1402))
# # Reduction des variables (extremites bruitees)
sp=sp[,seq(+1,ncol(sp)-30,1)]
# # SNV
sp=t(scale(t(sp)))
# # # Derivation Savitsky Golay
# sp=t(apply(sp,1,sgolayfilt,p=p,n=n,m=m))


titre=rownames(sp)
# ##Pour avoir les bons noms de lignes apres moyenne
# f=(substr(titre,15,16))== "01" | (substr(titre,15,16))== "04"| (substr(titre,15,16))== "07"| (substr(titre,15,16))== "10"| (substr(titre,15,16))== "13"| (substr(titre,15,16))== "16"
# spf=globalmatrix[f==TRUE,]


# ##Variable pour la moyenne des spectres
# z=3
# moy= 1187
# 
# mean1_final=matrix(nrow = 1, ncol = ncol(globalmatrix))
# mean2_final=matrix(nrow = 1187, ncol = ncol(globalmatrix))
# 
# mean1_final=colMeans(globalmatrix[(1:z),])
# 
# for(j in 1:moy) {
#   mean2=colMeans(globalmatrix[((z+1):(z+3)),])
#   z=z+3
#   #print(mean2)
#   #print(z)
#   
#   mean2_final[j,]=mean2
# }
# 
# sp=rbind(mean1_final,mean2_final)

# rownames(sp)=rownames(spf)

# titre2=rownames(sp)
###Rajouter nom cepage
##Filtre
ca=(substr(titre,11,13))== "015" |  (substr(titre,11,13))== "169" |  (substr(titre,11,13))== "685"
cb=(substr(titre,11,13))== "222" |  (substr(titre,11,13))== "509" |  (substr(titre,11,13))== "787"
cc=(substr(titre,11,13))== "471" |  (substr(titre,11,13))== "525" |  (substr(titre,11,13))== "747" |  (substr(titre,11,13))== "877"
##Separation en 3 matrices
spa=sp[ca==TRUE,]
spb=sp[cb==TRUE,]
spc=sp[cc==TRUE,]
##Rajoute le nom du cepage au nom de la ligne
rownames(spa)=paste(rownames(spa),"C",sep = "-")
rownames(spb)=paste(rownames(spb),"G",sep = "-")
rownames(spc)=paste(rownames(spc),"S",sep = "-")
##Recombine les 3 matrices pour reformer sp
sp_Cep=rbind(spa,spb,spc)

titre2=rownames(sp_Cep)
###Rajouter numero de souche
s1=(substr(titre2,15,16))== "01" | (substr(titre2,15,16))== "02"| (substr(titre2,15,16))== "03"| (substr(titre2,15,16))== "04"| (substr(titre2,15,16))== "05"| (substr(titre2,15,16))== "06"
s2=(substr(titre2,15,16))== "07" | (substr(titre2,15,16))== "08"| (substr(titre2,15,16))== "09"| (substr(titre2,15,16))== "10"| (substr(titre2,15,16))== "11"| (substr(titre2,15,16))== "12"
s3=(substr(titre2,15,16))== "13" | (substr(titre2,15,16))== "14"| (substr(titre2,15,16))== "15"| (substr(titre2,15,16))== "16"| (substr(titre2,15,16))== "17"| (substr(titre2,15,16))== "18"
##Separation en 3 matrices
sp1=sp_Cep[s1==TRUE,]
sp2=sp_Cep[s2==TRUE,]
sp3=sp_Cep[s3==TRUE,]
##Rajoute le numero de la souche au nom de la ligne
rownames(sp1)=paste(rownames(sp1),"1",sep = "-")
rownames(sp2)=paste(rownames(sp2),"2",sep = "-")
rownames(sp3)=paste(rownames(sp3),"3",sep = "-")
##Recombine les 3 matrices pour reformer sp
sp_sou=rbind(sp1,sp2,sp3)

titre3=rownames(sp_sou)
###Rajouter si jeune ou vieille feuille
J=(substr(titre3,15,16))== "01" | (substr(titre3,15,16))== "02"| (substr(titre3,15,16))== "03"| (substr(titre3,15,16))== "07"| (substr(titre3,15,16))== "08"| (substr(titre3,15,16))== "09"| (substr(titre3,15,16))== "13"| (substr(titre3,15,16))== "14"| (substr(titre3,15,16))== "15"
V=(substr(titre3,15,16))== "04" | (substr(titre3,15,16))== "05"| (substr(titre3,15,16))== "06"| (substr(titre3,15,16))== "10"| (substr(titre3,15,16))== "11"| (substr(titre3,15,16))== "12"| (substr(titre3,15,16))== "16"| (substr(titre3,15,16))== "17"| (substr(titre3,15,16))== "18"
##Separation en 2 matrices
spJ=sp_sou[J==TRUE,]
spV=sp_sou[V==TRUE,]
##Rajoute le numero de la souche au nom de la ligne
rownames(spJ)=paste(rownames(spJ),"J",sep = "-")
rownames(spV)=paste(rownames(spV),"V",sep = "-")
##Recombine les 2 matrices pour reformer sp
sp_ok=rbind(spJ,spV)

# titre4=rownames(sp_age)
# ###Rajouter nom du lieu
# E=(substr(titre4,1,4))== "2017" | (substr(titre4,1,8))== "20180619"| (substr(titre4,1,8))== "20180627"| (substr(titre4,1,8))== "20180704"| (substr(titre4,1,8))== "20180709"| (substr(titre4,1,8))== "20180816"
# X=(substr(titre4,1,8))== "20180710"| (substr(titre4,1,8))== "20180817"
# B=(substr(titre4,1,8))== "20180724" | (substr(titre4,1,8))== "20180810"
# A=(substr(titre4,1,8))== "20180731" | (substr(titre4,1,8))== "20180823"
# ##Separation en 4 matrices
# spE=sp_age[E==TRUE,]
# spX=sp_age[X==TRUE,]
# spB=sp_age[B==TRUE,]
# spA=sp_age[A==TRUE,]
# ##Rajoute le du lieu au nom de la ligne
# rownames(spE)=paste(rownames(spE),"E",sep = "-")
# rownames(spX)=paste(rownames(spX),"X",sep = "-")
# rownames(spB)=paste(rownames(spB),"B",sep = "-")
# rownames(spA)=paste(rownames(spA),"A",sep = "-")
# ##Recombine les 4 matrices pour reformer sp
# sp_ok=rbind(spE,spX,spB,spA)


###Classes
#Années
Annee=as.factor(substr(rownames(sp_ok),1,4))
c1=length(levels(Annee))
#rownames(Annee)=rownames(spf)
#Date
Date=as.factor(substr(rownames(sp_ok),1,8))
c2=length(levels(Date))
#Clone
Clone=as.factor(substr(rownames(sp_ok),11,13))
c3=length(levels(Clone))
#Cepage
Cepage=as.factor(substr(rownames(sp_ok),18,18))
c4=length(levels(Cepage))
#Souche
Souche=as.factor(substr(rownames(sp_ok),20,20))
c5=length(levels(Souche))
#Age
Age=as.factor(substr(rownames(sp_ok),22,22))
c6=length(levels(Age))
#Lieu
# Lieu=as.factor(substr(rownames(sp_ok),24,24))
# c7=length(levels(Lieu))

Facteurs=cbind(Annee, Date, Cepage,Clone,Cepage,Souche,Age)#Lieu)

##Base de donnees avec facteurs -> Facteurs
# Facteurs=data.frame(as.numeric(Annee), as.numeric(Date), as.numeric(Clone), as.numeric(Cepage), as.numeric(Souche),as.numeric(Age), as.numeric(Lieu))
# Facteurs$sp_ok=sp_ok
rownames(Facteurs)=rownames(sp_ok)

##ASCA
Asca= ASCA.Calculate(sp_ok, Facteurs, equation.elements = "3,4,34", scaling = FALSE,
               only.means.matrix = FALSE, use.previous.asca = NULL)

