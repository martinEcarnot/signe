# Fonction ls.str()


library(MASS)
# library(mixOmics)
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

rm(list = ls())

source('Script_R_2020/adj_asd.R')
source('Script_R_2020/SIGNE_load.R')
# source('C:/Users/No?mie/Desktop/SFE/Script_R/SIGNE_maha.R')
source('Script_R_2020/SIGNE_maha0.R')
source("Script_R_2020/sp2dfclo.R")

# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
#set.seed(1)

# brb3="~/Documents/NICOLAS/Stage de fin annee au Grau du roi/globalmatrixN1"
brb3="C:/Users/avitvale/Documents/Test/globalmatrix"
load(file=brb3)
globalmatrixN1=globalmatrix


# Data Filter



# iok=substr(rownames(globalmatrixN1),1,9) %in% dates
sp=globalmatrixN1 #[iok,]


L=unique(substr(rownames(sp),9,13))
L2=sort(L)

# sp2=as.data.frame(matrix(ncol = 5))
# names(sp2)=names(sp)
sp2=sp

for (i in 1:length(L2)){
  N=which(substr(rownames(sp),9,13)==L2[i])
  sp2=rbind(sp2,sp[N,])
}

sp3=sp2[(length(sp[,1])+1):length(sp2[,1]),]
rownames(sp3)=substr(rownames(sp3),1,18)

sp=sp3



#####sp=sp[which(substr(rownames(sp),18,18)=="G"),]


#
# date=as.factor(substr(rownames(sp2),5,8))
# DATE=c("0703","0704","0710","0709","0702")
# #DATE="0702"
# sp3=sp2[which(date  %in%  DATE),]
#
# parc=as.factor(substr(rownames(sp3),18,18))
# PARC="G"
# sp=sp3[which(parc %in% PARC),]
#
# annee=as.factor(substr(rownames(sp4),1,4))
# ANNEE=c("2017","2018")
# sp=sp4[which(annee  %in%  ANNEE),]

# annee=as.factor(substr(rownames(sp4),1,4))
# ANNEE="2019"
# sp=sp4[which(annee %in% ANNEE),]

#length(which(substr(rownames(sp),9,9)=="C"))

#unique(rownames(sp3)[which(substr(rownames(sp3),1,8)=="20180710")])


#aC= substr(rownames(sp),1,4)=="2017"
#sp =sp[(aC==TRUE),]

## Creation de la matrice de classes
# class=as.factor(substr(rownames(sp),11,13))

##Permet de voir si, avec des trucs aléatoires, ca ferait des prédictions "illusoires"
# L=c("C 015", "C 169", "C 685", "G 222", "G 509", "G 787", "S 471", "S 525", "S 747", "S 877")
# alea=L[sample(1:10,length(sp[,1]),replace = T) ]
#
# rownames(sp)=paste(rownames(sp),alea)



class=as.factor(substr(rownames(sp),9,9)) #9,9
classclo=as.factor(substr(rownames(sp),9,13)) #9,13

## Variable qui mesure le nombre de classes
c=length(levels(class))
cclo=length(levels(classclo))






# Création des jeux de calibration/ validation
# On créé un facteur datclone qui groupe un clone à 1 date
datclone=substr(rownames(sp),1,13)
ndc=length(unique(datclone))
# On créé un facteur souche qui groupe les 6 spectres de chaque souche
numsp=as.numeric(substr(rownames(sp),15,16))
souche=cut(numsp, breaks = c(0,6,12,18),labels=c("s1","s2","s3"))  # paste(datclone,cut(numsp, breaks = c(0,6,12,18),labels=c("s1","s2","s3")))

sp=sp2dfclo(sp,class,classclo)
sp=cbind(sp,datclone,souche)   # mutate(sp,datclone=substr(titre,1,13), souche=substr(souche,15,16))

##Mais... On se sert pas du y ajouté dans le data.frame plus tard, je crois. Bizarre.

# Le tirage sera fait plus loin dans la boucle
### FIXATION DES PARAMETRES UTILISES:
## Nombre de repetitions de la boucle de PLSDA:
repet= 2
## Parametres du Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
p=2
n=11
m=1
## Nombre de VL max autorisees
ncmax=35
## Nombre de groupes de CV
k=2

## PLSDA ##

### Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
#sp_pre=adj_asd(sp$x,c(602,1402))
sp_pre=sp$x

## Reduction des variables (extremites bruitees)
# sp=sp[,seq(51,ncol(sp)-30,1)]
## Coupure du spectre a 1000nm
#spx=sp[,seq(+1,601,1)]
#matplot(t(spx),pch = ".",xlab = "Longueurs d'ondes (nm)", ylab = "Transflectance")
## SNV
sp_pre=t(scale(t(sp_pre)))

## Derivation Savitsky Golay
sp$x=savitzkyGolay(sp_pre, m = m, p = p, w = n)



perok_finalm0=matrix(nrow = repet, ncol = ncmax)
perok_finalm=matrix(nrow = repet, ncol = ncmax)
#perok_final=matrix(nrow = repet, ncol = ncmax)

#perok_finalm0clo=matrix(nrow = repet, ncol = ncmax)
#perok_finalmclo=matrix(nrow = repet, ncol = ncmax)
perok_finalm0C=matrix(nrow = repet, ncol = ncmax)
perok_finalm0G=matrix(nrow = repet, ncol = ncmax)
perok_finalm0S=matrix(nrow = repet, ncol = ncmax)
perok_finalmC=matrix(nrow = repet, ncol = ncmax)
perok_finalmG=matrix(nrow = repet, ncol = ncmax)
perok_finalmS=matrix(nrow = repet, ncol = ncmax)
#perok_finalclo=matrix(nrow = repet, ncol = ncmax)

#sp=sp[-idval,]

###s?paration validation calibration PLSDA###
#set.seed(1) # fixe le tirage aleatoire
for(j in 1:repet) {
  # On selectionne le jeu de validation de manière à ce que tous les datclone soient représentés et 1 souche sur les 3 tirée random
  m=mstage(sp,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,ndc)))
  spval=getdata(sp,m)[[2]]
  #

  idval=which(rownames(sp)  %in%  rownames(spval))
  idvalC=which((rownames(sp)  %in%  rownames(spval)) & sp$y1=="C")
  idvalG=which((rownames(sp)  %in%  rownames(spval)) & sp$y1=="G")
  idvalS=which((rownames(sp)  %in%  rownames(spval)) & sp$y1=="S")
  #
  # ##On selectionne les spectres ayant ces num?ros dans le jeu de validation, les autres vont dans le jeu de calibration
  spval=sp[idval,]
  spcal=sp[-idval,]

  classval=class[idval]
  classcal=class[-idval]

  classvalclo=classclo[idval]
  classcalclo=classclo[-idval]


  spvalC=sp[idvalC,]
  spcalC1=sp[-idvalC,]
  spcalC=spcalC1[which(spcalC1$y1=="C"),]

  classvalC=classclo[idvalC]
  classvalC=droplevels(classvalC)

  classcalC1=classclo[-idvalC]
  classcalC=classcalC1[which(classcalC1=="C 015" | classcalC1=="C 169" | classcalC1=="C 685")]
  classcalC=droplevels(classcalC)



  spvalG=sp[idvalG,]
  spcalG1=sp[-idvalG,]
  spcalG=spcalG1[which(spcalG1$y1=="G"),]

  classvalG=classclo[idvalG]
  classvalG=droplevels(classvalG)

  classcalG1=classclo[-idvalG]
  classcalG=classcalG1[which(classcalG1=="G 222" | classcalG1=="G 509" | classcalG1=="G 787")]
  classcalG=droplevels(classcalG)



  spvalS=sp[idvalS,]
  spcalS1=sp[-idvalS,]
  spcalS=spcalS1[which(spcalS1$y1=="S"),]

  classvalS=classclo[idvalS]
  classvalS=droplevels(classvalS)

  classcalS1=classclo[-idvalS]
  classcalS=classcalS1[which(classcalS1=="S 471" | classcalS1=="S 525" | classcalS1=="S 747" | classcalS1=="S 877")]
  classcalS=droplevels(classcalS)

  #Ca, ca doit être bon

  # ## Creation des jeux d'apprentissage et validation
  predm=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))
  predm0=as.data.frame(matrix(nrow = length(classcal), ncol = ncmax))
  spcaldef=spcal # spcal deflaté du(des) groupe(s) de CV déjà validés

  predmC=as.data.frame(matrix(nrow = length(classvalC), ncol = ncmax))
  predm0C=as.data.frame(matrix(nrow = length(classcalC), ncol = ncmax))
  spcaldefC=spcalC # spcal deflaté du(des) groupe(s) de CV déjà validés

  predmG=as.data.frame(matrix(nrow = length(classvalG), ncol = ncmax))
  predm0G=as.data.frame(matrix(nrow = length(classcalG), ncol = ncmax))
  spcaldefG=spcalG # spcal deflaté du(des) groupe(s) de CV déjà validés

  predmS=as.data.frame(matrix(nrow = length(classvalS), ncol = ncmax))
  predm0S=as.data.frame(matrix(nrow = length(classcalS), ncol = ncmax))
  spcaldefS=spcalS # spcal deflaté du(des) groupe(s) de CV déjà validés


  # spcal=sp
  # spcaldef=spcal
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

    # # En mettant une autre année en validation
    # idvalCV =which(substr(rownames(spcal),1,4)  %in%  '2018')
    spvalCV=spcal[idvalCV,]
    classvalCV=classcal[idvalCV]  #identifiants des classes du jeu de validation
    spcalCV=spcal[-idvalCV,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    classcalCV=classcal[-idvalCV] #identifiants des classes du jeu de calibration


    spvalCVC=spcalC[idvalCVC,]
    classvalCVC=classcalC[idvalCVC]  #identifiants des classes du jeu de validation
    spcalCVC=spcalC[-idvalCVC,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    classcalCVC=classcalC[-idvalCVC] #identifiants des classes du jeu de calibration
#    classcalCVC=droplevels(classcalCVC)


    spvalCVG=spcalG[idvalCVG,]
    classvalCVG=classcalG[idvalCVG]  #identifiants des classes du jeu de validation
    spcalCVG=spcalG[-idvalCVG,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    classcalCVG=classcalG[-idvalCVG] #identifiants des classes du jeu de calibration
#    classcalCVG=droplevels(classcalCVG)


    spvalCVS=spcalS[idvalCVS,]
    classvalCVS=classcalS[idvalCVS]  #identifiants des classes du jeu de validation
    spcalCVS=spcalS[-idvalCVS,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    classcalCVS=classcalS[-idvalCVS] #identifiants des classes du jeu de calibration
#    classcalCVS=droplevels(classcalCVS)


    # ## PLSDA and application to have loadings and scores
    rplsda=caret::plsda(spcalCV$x, classcalCV,ncomp=ncmax)
    sccalCV=rplsda$scores
    spvalCV_c=scale(spvalCV$x,center=rplsda$Xmeans,scale = F)
    scvalCV=spvalCV_c%*%rplsda$projection    #score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

    rplsdaC=caret::plsda(spcalCVC$x, classcalCVC, ncomp=ncmax)
    sccalCVC=rplsdaC$scores
    spvalCVC_c=scale(spvalCVC$x,center=rplsdaC$Xmeans,scale = F)
    scvalCVC=spvalCVC_c%*%rplsdaC$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

    rplsdaG=caret::plsda(spcalCVG$x, classcalCVG, ncomp=ncmax)
    sccalCVG=rplsdaG$scores
    spvalCVG_c=scale(spvalCVG$x,center=rplsdaG$Xmeans,scale = F)
    scvalCVG=spvalCVG_c%*%rplsdaG$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

    rplsdaS=caret::plsda(spcalCVS$x, classcalCVS, ncomp=ncmax)
    sccalCVS=rplsdaS$scores
    spvalCVS_c=scale(spvalCVS$x,center=rplsdaS$Xmeans,scale = F)
    scvalCVS=spvalCVS_c%*%rplsdaS$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas



    for (ii in 2:ncmax) {
      ## Validation
      predm0[idvalCV,ii]=SIGNE_maha0(sccalCV[,1:ii], classcalCV, scvalCV[,1:ii])$class
      predm0C[idvalCVC,ii]=SIGNE_maha0(sccalCVC[,1:ii], classcalCVC, scvalCVC[,1:ii])$class
      predm0G[idvalCVG,ii]=SIGNE_maha0(sccalCVG[,1:ii], classcalCVG, scvalCVG[,1:ii])$class
      predm0S[idvalCVS,ii]=SIGNE_maha0(sccalCVS[,1:ii], classcalCVS, scvalCVS[,1:ii])$class
#      predm0G[idvalCVG,ii]=SIGNE_maha0(sccalCVG[,1:ii], classcalCVG, scvalCVG[,1:ii])$classclo
      # M1= matrix(nrow= nrow(scvalCV[,1:ii]), ncol= nlevels(classcalCV))
      # M2= matrix(nrow= nrow(scvalCVG[,1:ii]), ncol= nlevels(classcalCVG))
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
  # perokm0 =100*unlist(lapply(diagsm0, FUN = sum))/length(idvalCV)
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
  scvalC=spvalC_c%*%rplsdaC$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas


  rplsdaG=caret::plsda(spcalG$x, classcalG,ncomp=ncmax)
  sccalG=rplsdaG$scores
  spvalG_c=scale(spvalG$x,center=rplsdaG$Xmeans,scale = F)
  scvalG=spvalG_c%*%rplsdaG$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas


  rplsdaS=caret::plsda(spcalS$x, classcalS,ncomp=ncmax)
  sccalS=rplsdaS$scores
  spvalS_c=scale(spvalS$x,center=rplsdaS$Xmeans,scale = F)
  scvalS=spvalS_c%*%rplsdaS$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

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


  # tsmclo=lapply(as.list(predmclo), classvalclo, FUN = table)
  # diagsmclo=lapply(tsmclo, FUN = diag)
  # perokmclo =100*unlist(lapply(diagsmclo, FUN = sum))/length(idvalclo)
  # perok_finalmclo[j,]=perokmclo

}

plot(colMeans(perok_finalm0), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalm), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)

plot(colMeans(perok_finalm0C), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalmC), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)

plot(colMeans(perok_finalm0G), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalmG), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)


plot(colMeans(perok_finalm0S), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalmS), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)


#stop()

#  C 015   C 169   C 685   G 222   G 509   G 787   S 471   S 525   S 747   S 877

#length(which(substr(rownames(spcal),20,24)=="S 877"))

#unique(substr(rownames(sp[which(substr(rownames(sp),18,18)=="g"),]),1,4))
# length(which(substr(rownames(sp[which(substr(rownames(sp),18,18)=="G"),]),9,13)=="S 877"))

# "0613" "0617" "0624" "0628" "0702" "0710" "0726" "0718" "0730"
# unique(substr(rownames(sp),1,8))
idval=which(substr(rownames(sp),1,4)=="2019" & substr(rownames(sp),18,18)=="G" & substr(rownames(sp),5,8)!="0613")
#idval=which(substr(rownames(sp),18,18)!="G")

#which(substr(rownames(sp[idval,]),9,13)=="C 015" & substr(rownames(sp[idval,]),18,18)=="G" & substr(rownames(sp[idval,]),18,18)=="G")

# idval=idval[-c(1, 96, 466, 467, 542, 1022, 1023, 1098, 1668)]
# idval=idval[-L]
# spval=sp[idval,]
# rownames(spval[c(73, 91, 109, 519, 537, 555, 1167, 1185, 1203, 1221),])
#
# which(rownames(sp %in% rownames(spval[c(73, 91, 109, 519, 537, 555, 1167, 1185, 1203, 1221),])))

# m=mstage(sp,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,ndc)))
# spval=getdata(sp,m)[[2]]
# idval=which(rownames(sp)  %in%  rownames(spval) & substr(rownames(sp),18,18)=="G")

#
# L=c(73, 91, 109, 519, 537, 555, 1167, 1185, 1203, 1221)
# L2=c(L, (2+L), (2+L),(3+L),(4+L),(5+L))
# L2
#idval=which(substr(rownames(sp),1,4)=="2019" & substr(rownames(sp),18,18)=="G")

#ncmax=300
spval=sp[idval,]
spcal=sp[-idval,]
classval=class[idval]
classcal=class[-idval]
classvalclo=classclo[idval]
classcalclo=classclo[-idval]

idcal=which(substr(rownames(sp),18,18)=="G")
classcal=classcal[which(substr(rownames(spcal),18,18)=="G")]
classcalclo=classcalclo[which(substr(rownames(spcal),18,18)=="G")]
spcal=spcal[which(substr(rownames(spcal),18,18)=="G"),]


predmF=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))
distances=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))

rplsda=caret::plsda(spcal$x, classcal,ncomp=ncmax)
sccal=rplsda$scores
spval_c=scale(spval$x,center=rplsda$Xmeans,scale = F)
scval=spval_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

for (ii in 2:ncmax) {
  predmF[,ii]=SIGNE_maha0(sccal[,1:ii], classcal, scval[,1:ii])$class
  distances[,ii]=SIGNE_maha0(sccal[,1:ii], classcal, scval[,1:ii])$dist
}
tsm=lapply(as.list(predmF), classval, FUN = table)
diagsm=lapply(tsm, FUN = diag)
perokm =100*unlist(lapply(diagsm, FUN = sum))/length(idval)

plot(perokm, xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
perokm
tsm[25]



#Pour l'étude des distances on se place, un peu arbitrairement, à 8VL.
#132B4C ou #11274A
#4CB4FF ou #4FB2FC
bleus=c("#132B4C","#4CB4FF")
Ldist=as.data.frame(matrix(nrow = length(scval[,1]), ncol = 8))
VL=18
l1=c(1,2,3,7,8,9,13,14,15)
l2=c(4,5,6,10,11,12,16,17,18)
substr(rownames(scval),15,16)
for (i in 1:length(scval[,1])){
  S=distances[VL][i,][1]+distances[VL][i,][2]+distances[VL][i,][3]
  Ldist[i,1]=min(distances[VL][i,][1]/S,distances[VL][i,][2]/S,distances[VL][i,][3]/S)

  Ldist[i,3]=substr(rownames(scval)[i],9,9)

  n=which(c("C","G","S")==Ldist[i,3])
  Ldist[i,2]=min(distances[VL][i,][1],distances[VL][i,][2],distances[VL][i,][3])
  #Ldist[i,2]=distances[VL][i,][n]

  Ldist[i,4]=as.character(predmF[VL][i,])
  Ldist[i,5]="Mal classé"
  if (Ldist[i,3]==Ldist[i,4]){
    Ldist[i,5]="Bien classé"
  }
  Ldist[i,6]=i
  Ldist[i,7]=substr(rownames(scval)[i],5,8)
  Ldist[i,8]="J"
  if (substr(rownames(scval)[i],15,16) %in% l2){
    Ldist[i,8]="V"
  }
}







bleu=rgb(0.09,0.63,1)
bleu2=rgb(0,0.43,0.8)
bleu3=rgb(0.59,0.83,1)
vert=rgb(0.10,0.94,0.36)
vert2=rgb(0.12,0.75,0.10)
vert3=rgb(0.50,0.94,0.36)
vert4=rgb(0.50,0.75,0.36)
rouge=rgb(1,0.35,0.13)
rouge2=rgb(0.8,0.35,0.13)
rouge3=rgb(1,0.55,0.33)
colo=c(rouge, rouge3, bleu, bleu3, vert, vert4)


aff <- ggplot(Ldist, aes(x=Ldist[,6], y=(Ldist[,2]),colour=paste(Ldist[,3],Ldist[,5]),date=Ldist[,7],cepage=Ldist[,3])) +
  geom_point(size=2, alpha=1) +
  scale_color_manual(values = colo)
ggplotly(aff)





idcalC=which(spcal$y1=="C")
spcalC=spcal[idcalC,]
classcalC=classcalclo[idcalC]

rplsdaC=caret::plsda(spcalC$x, classcalC,ncomp=ncmax)
sccalC=rplsdaC$scores

## On pourrait considérer qu'il y a en réalité deux grandeurs à faire varier, en n'utilisant pas le même nb de VL pour la PLSDA sur cépages et pour la PLSDA sur clones
perokmC=0
for (ii in 2:ncmax){        #NB : On met tout dans la boucle, contrairement à ce qu'on fait pour former perokm parce qu'ici length(classvalCT) est variable en f° de ii
  idvalCT=which(predmF[,25]=="C")
  spvalCT=spval[idvalCT,]
  classvalCT=classvalclo[idvalCT]

  spvalC_c=scale(spvalCT$x,center=rplsdaC$Xmeans,scale = F)
  scvalC=spvalC_c%*%rplsdaC$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
  predmC=SIGNE_maha0(sccalC[,1:ii], classcalC, scvalC[,1:ii])$class
  tsmC=table(predmC, classvalCT)
  # print(ii)
  # print(tsmC)
  # print("")
  diagsmC=diag(tsmC)
  perokmC[ii] =100*sum(diagsmC)/length(idvalCT)
}
plot(perokmC, xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)

perokmC
tsmC





idcalG=which(spcal$y1=="G")
spcalG=spcal[idcalG,]
classcalG=classcalclo[idcalG]

rplsdaG=caret::plsda(spcalG$x, classcalG,ncomp=ncmax)
sccalG=rplsdaG$scores

perokmG=0
for (ii in 2:ncmax){        #NB : On met tout dans la boucle, contrairement à ce qu'on fait pour former perokm parce qu'ici length(classvalCT) est variable en f° de ii
  idvalGT=which(predmF[,25]=="G")
  spvalGT=spval[idvalGT,]
  classvalGT=classvalclo[idvalGT]

  spvalG_c=scale(spvalGT$x,center=rplsdaG$Xmeans,scale = F)
  scvalG=spvalG_c%*%rplsdaG$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
  predmG=SIGNE_maha0(sccalG[,1:ii], classcalG, scvalG[,1:ii])$class
  classvalGT=relevel(classvalGT, "G 787")
  classvalGT=relevel(classvalGT, "G 509")
  classvalGT=relevel(classvalGT, "G 222")
  tsmG=table(predmG, classvalGT)
  # print(ii)
  # print(tsmG)
  # print("")
  diagsmG=diag(tsmG)
  perokmG[ii] =100*sum(diagsmG)/length(idvalGT)
}
plot(perokmG, xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)

perokmG
tsmG




idcalS=which(spcal$y1=="S")
spcalS=spcal[idcalS,]
classcalS=classcalclo[idcalS]

rplsdaS=caret::plsda(spcalS$x, classcalS,ncomp=ncmax)
sccalS=rplsdaS$scores

perokmS=0
for (ii in 2:ncmax){        #NB : On met tout dans la boucle, contrairement à ce qu'on fait pour former perokm parce qu'ici length(classvalCT) est variable en f° de ii
  idvalST=which(predmF[,25]=="S")
  spvalST=spval[idvalST,]
  classvalST=classvalclo[idvalST]

  spvalS_c=scale(spvalST$x,center=rplsdaS$Xmeans,scale = F)
  scvalS=spvalS_c%*%rplsdaS$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
  predmS=SIGNE_maha0(sccalS[,1:ii], classcalS, scvalS[,1:ii])$class
  classvalST=relevel(classvalST, "S 877")
  classvalST=relevel(classvalST, "S 747")
  classvalST=relevel(classvalST, "S 525")
  classvalST=relevel(classvalST, "S 471")
  tsmS=table(predmS, classvalST)
  diagsmS=diag(tsmS)
  perokmS[ii] =100*sum(diagsmS)/length(idvalST)
}
plot(perokmS, xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)

perokmS
tsmS



predmclo=c(as.character(predmC),as.character(predmG),as.character(predmS))
classvalcloT=c(as.character(classvalCT),as.character(classvalGT),as.character(classvalST))
tsmclo=table(predmclo, classvalcloT)
diagsmclo=diag(tsmclo)
perokmclo =100*sum(diagsmclo)/length(idval)
perokmclo
tsmclo



##Clones direct
predmcloD=as.data.frame(matrix(nrow = length(classvalclo), ncol = ncmax))

rplsdaD=caret::plsda(spcal$x, classcalclo,ncomp=ncmax)
sccalD=rplsdaD$scores
spval_cD=scale(spval$x,center=rplsdaD$Xmeans,scale = F)
scvalD=spval_cD%*%rplsdaD$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

for (ii in 2:ncmax) {predmcloD[,ii]=SIGNE_maha0(sccalD[,1:ii], classcalclo, scvalD[,1:ii])$class}
tsmD=lapply(as.list(predmcloD), classvalclo, FUN = table)
diagsmD=lapply(tsmD, FUN = diag)
perokmD =100*unlist(lapply(diagsmD, FUN = sum))/length(idval)

plot(perokmD, xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
perokmD
tsmD[25]


###### Knn

idval=which(substr(rownames(sp),1,4)=="2017")
spval=sp[idval,]


# m=mstage(sp,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,ndc)))
# spval=getdata(sp,m)[[2]]
# idval=which(rownames(sp)  %in%  rownames(spval))



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




# rplsdaK=knnwda(Xr=spcal$x, Yr=as.character(classcal), Xu=spval$x, Yu=as.character(classval), diss="mahalanobis", ncompdis=16, h=1.1, k=700, print=F)
# predmFK= rplsdaK$fit$y1
#
# predmFK= rplsdaL$fit$y1[rplsdaL$fit$ncomp==15&rplsdaL$fit$k==100&rplsdaL$fit$ncompdis==20]
# tsmK=table(predmFK, classval)
# diagsmK=diag(tsmK)
# perokmK=100*sum(diagsmK)/length(idval)
# perokmK
#
# TEST=c(10,20,25,30,35,40)
# predmFK=as.data.frame(matrix(nrow = length(classval), ncol = length(TEST)))
# for (ii in 1:length(TEST)) {
#   print(ii)
#   rplsdaK=knnwda(Xr=spcal$x, Yr=as.character(classcal), Xu=spval$x, Yu=as.character(classval), diss="mahalanobis", ncompdis=23, h=1.1, k=TEST[ii], print=F)
#   predmFK[,ii]= rplsdaK$fit$y1
# }
# tsmK=lapply(as.list(predmFK), classval, FUN = table)
# diagsmK=lapply(tsmK, FUN = diag)
# perokmK =100*unlist(lapply(diagsmK, FUN = sum))/length(idval)
# #VL entre 11 et 16. Ici 12.
# plot(perokmK, xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)


idvalCTK=which(predmFK=="C")
spvalCTK=spval[idvalCTK,]
classvalCTK=classvalclo[idvalCTK]

idcalCK=which(spcal$y1=="C")
spcalCK=spcal[idcalCK,]
classcalCK=classcalclo[idcalCK]

rplsdaCK=knnwda(Xr=spcalCK$x, Yr=as.character(classcalCK), Xu=spvalCTK$x, Yu=as.character(classvalCTK), diss="mahalanobis", ncompdis=20, h=1.1, k=100, print=F)

rplsdaCK=lwplsdalm(Xr=spcalCK$x, Yr=as.character(classcalCK), Xu=spvalCTK$x, Yu=as.character(classvalCTK), diss="mahalanobis", ncomp=27, ncompdis=c(18,20,22,25,27), h=c(1,1.1,1.5), k=c(500,1000,1300), print=T)


z <- err(rplsdaCK, ~ ncomp + k + ncompdis)

u <- z
u$group <- paste("k=", u$k, ", ncompdis=", u$ncompdis, sep = "")
#plotmse(u, group = "group")

ggplot(data = u,aes(x=ncomp,y=errp,group = group,color =group))+ geom_line()

z[z$errp == min(z$errp), ]
plotmse(z, nam = "errp", group =)





# rplsdaC=caret::plsda(spcalCK$x, classcalCK,ncomp=ncmax)
# sccalC=rplsdaC$scores
# spvalC_c=scale(spvalCTK$x,center=rplsdaC$Xmeans,scale = F)
# scvalC=spvalC_c%*%rplsdaC$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
# predmCK=SIGNE_maha0(sccalC[,1:25], classcalCK, scvalC[,1:25])$class

predmFK= rplsdaL$fit$y1[rplsdaL$fit$ncomp==4&rplsdaL$fit$k==1000&rplsdaL$fit$ncompdis==22]
#predmCK= rplsdaCK$fit$y1
tsmCK=table(predmCK, classvalCTK)
tsmCK
diagsmCK=diag(tsmCK)
perokmCK=100*sum(diagsmCK)/length(idvalCTK)
perokmCK


idvalGTK=which(predmFK=="G")
spvalGTK=spval[idvalGTK,]
classvalGTK=classvalclo[idvalGTK]
classvalGTK=relevel(classvalGTK, "G 787")
classvalGTK=relevel(classvalGTK, "G 509")
classvalGTK=relevel(classvalGTK, "G 222")

idcalGK=which(spcal$y1=="G")
spcalGK=spcal[idcalGK,]
classcalGK=classcalclo[idcalGK]

rplsdaGK=knnwda(Xr=spcalGK$x, Yr=as.character(classcalGK), Xu=spvalGTK$x, Yu=as.character(classvalGTK), diss="mahalanobis", ncompdis=20, h=1.1, k=100, print=F)
predmGK= rplsdaGK$fit$y1
tsmGK=table(predmGK, classvalGTK)
tsmGK
diagsmGK=diag(tsmGK)
perokmGK=100*sum(diagsmGK)/length(idvalGTK)
perokmGK


idvalSTK=which(predmFK=="S")
spvalSTK=spval[idvalSTK,]
classvalSTK=classvalclo[idvalSTK]
classvalSTK=relevel(classvalSTK, "S 877")
classvalSTK=relevel(classvalSTK, "S 747")
classvalSTK=relevel(classvalSTK, "S 525")
classvalSTK=relevel(classvalSTK, "S 471")

idcalSK=which(spcal$y1=="S")
spcalSK=spcal[idcalSK,]
classcalSK=classcalclo[idcalSK]

rplsdaSK=knnwda(Xr=spcalSK$x, Yr=as.character(classcalSK), Xu=spvalSTK$x, Yu=as.character(classvalSTK), diss="mahalanobis", ncompdis=20, h=1.1, k=10, print=F)
predmSK= rplsdaSK$fit$y1
tsmSK=table(predmSK, classvalSTK)
tsmSK
diagsmSK=diag(tsmSK)
perokmSK=100*sum(diagsmSK)/length(idvalSTK)
perokmSK


idvalcloTK=predmFK
predmcloK=c(as.character(predmCK),as.character(predmGK),as.character(predmSK))
classvalcloTK=c(as.character(classvalCTK),as.character(classvalGTK),as.character(classvalSTK))
tsmcloK=table(predmcloK, classvalcloTK)
tsmcloK
diagsmcloK=diag(tsmcloK)
perokmcloK =100*sum(diagsmcloK)/length(idvalcloTK)
perokmcloK

# ##### Lw
# #lwplsdalm
# rplsdaL=lwplsda(Xr=spcal$x, Yr=as.character(classcal), Xu=spval$x, Yu=as.character(classval), diss="euclidian", ncompdis=25, ncomp=1, h=1.1, k=1000, print=F)
# predmFL= rplsdaL$fit$y1[(1+(24*length(classval))):(25*length(classval))]
# tsmL=table(predmFL, classval)
# diagsmL=diag(tsmL)
# perokmL=100*sum(diagsmL)/length(idval)
# perokmL



print(paste0("Precision discrimination cépages : " , perokm))
print(paste0("Precision discrimination clones Cabernet-Sauvignon : " , perokmC))
print(paste0("Precision discrimination clones Gamay : " , perokmG))
print(paste0("Precision discrimination clones Syrah : " , perokmS))
print(paste0("Precision discrimination clones globale : " , perokmclo))

print(paste0("kkn : Precision discrimination cépages : " , perokmK))
print(paste0("kkn : Precision discrimination clones Cabernet-Sauvignon : " , perokmCK))
print(paste0("kkn : Precision discrimination clones Gamay : " , perokmGK))
print(paste0("kkn : Precision discrimination clones Syrah : " , perokmSK))
print(paste0("kkn : Precision discrimination clones globale : " , perokmcloK))










