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
sp_pre=adj_asd(sp,c(552,1352)) #Retrouvé empiriquement, ne correspondait pas à ce qui etait indiqué precedemment.

## SNV
sp_pre=t(scale(t(sp_pre)))

## Derivation Savitsky Golay
sp=savitzkyGolay(sp_pre, m = m, p = p, w = n)


sp=sp[,-(1330:1390)] #Coupure du spectre autour des résidus de sauts de detecteurs. Ne semble pas encore effacer completement les sauts.
sp=sp[,-(500:560)]


#Changement de l'ordre des spectres, par exemple en les classant par date/par cepage/par clone... Classement fin en en faisant plusieurs à la suite
#Changer l'ordre est notamment utile pour les fonctions d'affichage.

#Laisser le premier pour un affichage clone/cepage, laisser les trois pour un affichage clone/cepage.
L=unique(substr(rownames(sp),9,13))
L2=sort(L)
sp2=sp

for (i in 1:length(L2)){
  N=which(substr(rownames(sp),9,13)==L2[i])
  sp2=rbind(sp2,sp[N,])
}

sp3=sp2[(length(sp[,1])+1):length(sp2[,1]),]
rownames(sp3)=substr(rownames(sp3),1,18)

sp=sp3




L=unique(substr(rownames(sp),1,8))
L2=sort(L)
sp2=sp

for (i in 1:length(L2)){
  N=which(substr(rownames(sp),1,8)==L2[i])
  sp2=rbind(sp2,sp[N,])
}

sp3=sp2[(length(sp[,1])+1):length(sp2[,1]),]
rownames(sp3)=substr(rownames(sp3),1,18)

sp=sp3


L=unique(substr(rownames(sp),9,9))
L2=sort(L)
sp2=sp

for (i in 1:length(L2)){
  N=which(substr(rownames(sp),9,9)==L2[i])
  sp2=rbind(sp2,sp[N,])
}

sp3=sp2[(length(sp[,1])+1):length(sp2[,1]),]
rownames(sp3)=substr(rownames(sp3),1,18)

sp=sp3



##Permet de voir si, avec des trucs aleatoires, ca ferait des prédictions "illusoires"
# L=c("C 015", "C 169", "C 685", "G 222", "G 509", "G 787", "S 471", "S 525", "S 747", "S 877")
# alea=L[sample(1:10,length(sp[,1]),replace = T) ]
#
# rownames(sp)=paste(rownames(sp),alea)



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




perok_finalm0=matrix(nrow = repet, ncol = ncmax)
perok_finalm=matrix(nrow = repet, ncol = ncmax)

perok_finalm0C=matrix(nrow = repet, ncol = ncmax)
perok_finalm0G=matrix(nrow = repet, ncol = ncmax)
perok_finalm0S=matrix(nrow = repet, ncol = ncmax)
perok_finalmC=matrix(nrow = repet, ncol = ncmax)
perok_finalmG=matrix(nrow = repet, ncol = ncmax)
perok_finalmS=matrix(nrow = repet, ncol = ncmax)



###separation validation calibration PLSDA###
#set.seed(1) # fixe le tirage aleatoire

###Le resultat de cette boucle est touché par la surcalibration : il peut être utilisé pour comparer differentes methodes entre elles, mais le pourcentage de reussite en lui-meme doit etre consideré avec prudence.
for(j in 1:repet) {
  # On selectionne le jeu de validation de maniere à ce que tous les datclone soient représentés et 1 souche sur les 3 tirée random
  m=mstage(sp,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,ndc)), method=list("srswor","srswor"))
  spval=getdata(sp,m)[[2]]
  #

  idval=which(rownames(sp)  %in%  rownames(spval))
  idvalC=which((rownames(sp)  %in%  rownames(spval)) & sp$y1=="C")
  idvalG=which((rownames(sp)  %in%  rownames(spval)) & sp$y1=="G")
  idvalS=which((rownames(sp)  %in%  rownames(spval)) & sp$y1=="S")
  #
  # ##On selectionne les spectres ayant ces numeros dans le jeu de validation, les autres vont dans le jeu de calibration
  spval=sp[idval,]
  spcal=sp[-idval,]

  classval=sp$y1[idval]
  classcal=sp$y1[-idval]

  classvalclo=sp$y2[idval]
  classcalclo=sp$y2[-idval]


  spvalC=sp[idvalC,]
  spcalC1=sp[-idvalC,]
  spcalC=spcalC1[which(spcalC1$y1=="C"),]

  classvalC=sp$y2[idvalC]
  classvalC=droplevels(classvalC)

  classcalC1=sp$y2[-idvalC]
  classcalC=classcalC1[which(classcalC1=="C 015" | classcalC1=="C 169" | classcalC1=="C 685")]
  classcalC=droplevels(classcalC)



  spvalG=sp[idvalG,]
  spcalG1=sp[-idvalG,]
  spcalG=spcalG1[which(spcalG1$y1=="G"),]

  classvalG=sp$y2[idvalG]
  classvalG=droplevels(classvalG)

  classcalG1=sp$y2[-idvalG]
  classcalG=classcalG1[which(classcalG1=="G 222" | classcalG1=="G 509" | classcalG1=="G 787")]
  classcalG=droplevels(classcalG)



  spvalS=sp[idvalS,]
  spcalS1=sp[-idvalS,]
  spcalS=spcalS1[which(spcalS1$y1=="S"),]

  classvalS=sp$y2[idvalS]
  classvalS=droplevels(classvalS)

  classcalS1=sp$y2[-idvalS]
  classcalS=classcalS1[which(classcalS1=="S 471" | classcalS1=="S 525" | classcalS1=="S 747" | classcalS1=="S 877")]
  classcalS=droplevels(classcalS)

  ### Creation des jeux d'apprentissage et validation
  predm=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))
  predm0=as.data.frame(matrix(nrow = length(classcal), ncol = ncmax))
  spcaldef=spcal # spcal deflaté du(des) groupe(s) de CV déjà validés

  predmC=as.data.frame(matrix(nrow = length(classvalC), ncol = ncmax))
  predm0C=as.data.frame(matrix(nrow = length(classcalC), ncol = ncmax))
  spcaldefC=spcalC

  predmG=as.data.frame(matrix(nrow = length(classvalG), ncol = ncmax))
  predm0G=as.data.frame(matrix(nrow = length(classcalG), ncol = ncmax))
  spcaldefG=spcalG

  predmS=as.data.frame(matrix(nrow = length(classvalS), ncol = ncmax))
  predm0S=as.data.frame(matrix(nrow = length(classcalS), ncol = ncmax))
  spcaldefS=spcalS

## Boucle CV
  for (i in 1:k) {
    print(i)
    m=mstage(spcaldef,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,ndc)))
    spvalCV=getdata(spcaldef,m)[[2]]

    idvalCV =which(rownames(spcal)  %in%  rownames(spvalCV))
    idvalCVC =which(rownames(spcalC)  %in%  rownames(spvalCV))
    idvalCVG =which(rownames(spcalG)  %in%  rownames(spvalCV))
    idvalCVS =which(rownames(spcalS)  %in%  rownames(spvalCV))

    spcaldef=spcaldef[-(which(rownames(spcaldef)  %in%  rownames(spvalCV))),]
    spcaldefC=spcaldefC[-(which(rownames(spcaldefC)  %in%  rownames(spvalCV))),]
    spcaldefG=spcaldefG[-(which(rownames(spcaldefG)  %in%  rownames(spvalCV))),]
    spcaldefS=spcaldefS[-(which(rownames(spcaldefS)  %in%  rownames(spvalCV))),]

    spvalCV=spcal[idvalCV,]
    classvalCV=classcal[idvalCV]  #identifiants des classes du jeu de validation
    spcalCV=spcal[-idvalCV,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    classcalCV=classcal[-idvalCV] #identifiants des classes du jeu de calibration


    spvalCVC=spcalC[idvalCVC,]
    classvalCVC=classcalC[idvalCVC]
    spcalCVC=spcalC[-idvalCVC,]
    classcalCVC=classcalC[-idvalCVC]


    spvalCVG=spcalG[idvalCVG,]
    classvalCVG=classcalG[idvalCVG]
    spcalCVG=spcalG[-idvalCVG,]
    classcalCVG=classcalG[-idvalCVG]


    spvalCVS=spcalS[idvalCVS,]
    classvalCVS=classcalS[idvalCVS]
    spcalCVS=spcalS[-idvalCVS,]
    classcalCVS=classcalS[-idvalCVS]


    # ## PLSDA and application to have loadings and scores
    rplsda=caret::plsda(spcalCV$x, classcalCV,ncomp=ncmax)
    sccalCV=rplsda$scores
    spvalCV_c=scale(spvalCV$x,center=rplsda$Xmeans,scale = F)
    scvalCV=spvalCV_c%*%rplsda$projection    #score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

    rplsdaC=caret::plsda(spcalCVC$x, classcalCVC, ncomp=ncmax)
    sccalCVC=rplsdaC$scores
    spvalCVC_c=scale(spvalCVC$x,center=rplsdaC$Xmeans,scale = F)
    scvalCVC=spvalCVC_c%*%rplsdaC$projection

    rplsdaG=caret::plsda(spcalCVG$x, classcalCVG, ncomp=ncmax)
    sccalCVG=rplsdaG$scores
    spvalCVG_c=scale(spvalCVG$x,center=rplsdaG$Xmeans,scale = F)
    scvalCVG=spvalCVG_c%*%rplsdaG$projection

    rplsdaS=caret::plsda(spcalCVS$x, classcalCVS, ncomp=ncmax)
    sccalCVS=rplsdaS$scores
    spvalCVS_c=scale(spvalCVS$x,center=rplsdaS$Xmeans,scale = F)
    scvalCVS=spvalCVS_c%*%rplsdaS$projection



    for (ii in 2:ncmax) {
      ## Validation
      predm0[idvalCV,ii]=SIGNE_maha0(sccalCV[,1:ii], classcalCV, scvalCV[,1:ii])$class
      predm0C[idvalCVC,ii]=SIGNE_maha0(sccalCVC[,1:ii], classcalCVC, scvalCVC[,1:ii])$class
      predm0G[idvalCVG,ii]=SIGNE_maha0(sccalCVG[,1:ii], classcalCVG, scvalCVG[,1:ii])$class
      predm0S[idvalCVS,ii]=SIGNE_maha0(sccalCVS[,1:ii], classcalCVS, scvalCVS[,1:ii])$class
    }
  }


  ## Table de contingence CV
  tsm0=lapply(as.list(predm0), classcal, FUN = table)
  tsm0C=lapply(as.list(predm0C), classcalC, FUN = table)
  tsm0G=lapply(as.list(predm0G), classcalG, FUN = table)
  tsm0S=lapply(as.list(predm0S), classcalS, FUN = table)
  ## Matrice mauvais classements par clone CV
  diagsm0=lapply(tsm0, FUN = diag)
  diagsm0C=lapply(tsm0C, FUN = diag)
  diagsm0G=lapply(tsm0G, FUN = diag)
  diagsm0S=lapply(tsm0S, FUN = diag)
  ## Pourcentage de bien classes CV
  perokm0 =100*unlist(lapply(diagsm0, FUN = sum))/length(classcal)
  perokm0C =100*unlist(lapply(diagsm0C, FUN = sum))/length(classcalC)
  perokm0G =100*unlist(lapply(diagsm0G, FUN = sum))/length(classcalG)
  perokm0S =100*unlist(lapply(diagsm0S, FUN = sum))/length(classcalS)
  ### Enregistrement des matrices de resultat final CV
  ##Remplissage de la matrice des perok finale
  perok_finalm0[j,]=perokm0
  perok_finalm0C[j,]=perokm0C
  perok_finalm0G[j,]=perokm0G
  perok_finalm0S[j,]=perokm0S

  # ## PLSDA sur le jeu de validation
  rplsda=caret::plsda(spcal$x, classcal,ncomp=ncmax)
  sccal=rplsda$scores
  spval_c=scale(spval$x,center=rplsda$Xmeans,scale = F)
  scval=spval_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
  for (ii in 2:ncmax) {predm[,ii]=SIGNE_maha0(sccal[,1:ii], classcal, scval[,1:ii])$class}
  tsm=lapply(as.list(predm), classval, FUN = table)
  diagsm=lapply(tsm, FUN = diag)
  perokm =100*unlist(lapply(diagsm, FUN = sum))/length(idval)
  perok_finalm[j,]=perokm


  rplsdaC=caret::plsda(spcalC$x, classcalC,ncomp=ncmax)
  sccalC=rplsdaC$scores
  spvalC_c=scale(spvalC$x,center=rplsdaC$Xmeans,scale = F)
  scvalC=spvalC_c%*%rplsdaC$projection


  rplsdaG=caret::plsda(spcalG$x, classcalG,ncomp=ncmax)
  sccalG=rplsdaG$scores
  spvalG_c=scale(spvalG$x,center=rplsdaG$Xmeans,scale = F)
  scvalG=spvalG_c%*%rplsdaG$projection


  rplsdaS=caret::plsda(spcalS$x, classcalS,ncomp=ncmax)
  sccalS=rplsdaS$scores
  spvalS_c=scale(spvalS$x,center=rplsdaS$Xmeans,scale = F)
  scvalS=spvalS_c%*%rplsdaS$projection

  for (ii in 2:ncmax) {
    predmC[,ii]=SIGNE_maha0(sccalC[,1:ii], classcalC, scvalC[,1:ii])$class
    predmG[,ii]=SIGNE_maha0(sccalG[,1:ii], classcalG, scvalG[,1:ii])$class
    predmS[,ii]=SIGNE_maha0(sccalS[,1:ii], classcalS, scvalS[,1:ii])$class
    }

  tsmC=lapply(as.list(predmC), classvalC, FUN = table)
  diagsmC=lapply(tsmC, FUN = diag)
  perokmC =100*unlist(lapply(diagsmC, FUN = sum))/length(idvalC)
  perok_finalmC[j,]=perokmC

  tsmG=lapply(as.list(predmG), classvalG, FUN = table)
  diagsmG=lapply(tsmG, FUN = diag)
  perokmG =100*unlist(lapply(diagsmG, FUN = sum))/length(idvalG)
  perok_finalmG[j,]=perokmG

  tsmS=lapply(as.list(predmS), classvalS, FUN = table)
  diagsmS=lapply(tsmS, FUN = diag)
  perokmS =100*unlist(lapply(diagsmS, FUN = sum))/length(idvalS)
  perok_finalmS[j,]=perokmS
}

plot(colMeans(perok_finalm0), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalm), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
perok_finalm

plot(colMeans(perok_finalm0C), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalmC), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)

plot(colMeans(perok_finalm0G), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalmG), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)

plot(colMeans(perok_finalm0S), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalmS), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)


#stop()









###### Knn

idval=which(substr(rownames(sp),1,4)=="2017")
spval=sp[idval,]



spval=sp[idval,]
spcal=sp[-idval,]
classval=class[idval]
classcal=class[-idval]
classvalclo=classclo[idval]
classcalclo=classclo[-idval]


#predmFK=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))

### knnplsda lwplsda  ####
predmFK=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))

rplsdaL=lwplsdalm(Xr=spcal$x, Yr=as.character(classcal), Xu=spval$x, Yu=as.character(classval), diss="mahalanobis", ncomp=27, ncompdis=27, h=1.5, k=1300, print=T)
# predmFK= rplsdaL$fit$y1[rplsdaL$fit$ncomp==25&rplsdaL$fit$k==400&rplsdaL$fit$ncompdis==27&rplsdaL$fit$h==1.5]
# tsmK=table(predmFK, classval)
# diagsmK=diag(tsmK)
# perokmK=100*sum(diagsmK)/length(idval)
# perokmK
# tsmK
for (ii in 2:ncmax) {predmFK[,ii]=rplsdaL$fit$y1[rplsdaL$fit$ncomp==ii&rplsdaL$fit$k==1300&rplsdaL$fit$ncompdis==27&rplsdaL$fit$h==1.5]}
tsm=lapply(as.list(predmFK), classval, FUN = table)
diagsm=lapply(tsm, FUN = diag)
perokm =100*unlist(lapply(diagsm, FUN = sum))/length(idval)

plot(perokm, xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
perokm
tsm[25]


erreurs=spval[-which(predmFK==classval),]
bons=spval[which(predmFK==classval),]




rownames(erreurs)
length(which(substr(rownames(spval),9,9)=="S"))/length(rownames(spval))
length(which(substr(rownames(erreurs),9,9)=="S"))/length(rownames(erreurs))
length(which(substr(rownames(bons),9,9)=="S"))/length(rownames(bons))

length(which(substr(rownames(erreurs),9,9)=="C"))/length(which(substr(rownames(spval),9,9)=="C"))


unique(rownames(sp[which(substr(rownames(sp),1,8)=="20170612"),]))

rownames(sp[which(substr(rownames(sp),1,4)=="2017"),])
unique(substr(rownames(sp),1,8))

unique(substr(rownames(erreurs),9,13))

length(which(substr(rownames(sp),1,8)=="20180724"))
rownames(sp[which(substr(rownames(sp),18,18)=="g"),])

length(which(substr(rownames(sp),1,4)=="2017"))

names(rplsdaL)
head(rplsdaL$y)
head(rplsdaL$fit)
head(rplsdaL$r)
z <- err(rplsdaL, ~ ncomp + k + ncompdis)

u <- z
u$group <- paste("k=", u$k, ", ncompdis=", u$ncompdis, sep = "")
#plotmse(u, group = "group")

ggplot(data = u,aes(x=ncomp,y=errp,group = group,color =group))+ geom_line()

z[z$errp == min(z$errp), ]
plotmse(z, nam = "errp", group =)



