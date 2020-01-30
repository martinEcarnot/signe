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

#aC= substr(rownames(sp),1,4)=="2017"
#sp =sp[(aC==TRUE),]

## Creation de la matrice de classes
# class=as.factor(substr(rownames(sp),11,13))
class=as.factor(substr(rownames(sp),9,9))
classclo=as.factor(substr(rownames(sp),11,13))

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
ncmax=25
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
perok_final=matrix(nrow = repet, ncol = ncmax)

perok_finalm0clo=matrix(nrow = repet, ncol = ncmax)
perok_finalmclo=matrix(nrow = repet, ncol = ncmax)
perok_finalmC=matrix(nrow = repet, ncol = ncmax)
perok_finalmG=matrix(nrow = repet, ncol = ncmax)
perok_finalmS=matrix(nrow = repet, ncol = ncmax)
perok_finalclo=matrix(nrow = repet, ncol = ncmax)

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

  classvalclo=class[idval]
  classcalclo=classclo[-idval]


  spvalC=sp[idvalC,]
  spcalC1=sp[-idvalC,]
  spcalC=spcalC1[which(spcalC1$y1=="C"),]

  classvalC=classclo[idvalC]
  classvalC=droplevels(classvalC)

  classcalC1=classclo[-idvalC]
  classcalC=classcalC1[which(classcalC1=="015" | classcalC1=="169" | classcalC1=="685")]
  classcalC=droplevels(classcalC)



  spvalG=sp[idvalG,]
  spcalG1=sp[-idvalG,]
  spcalG=spcalG1[which(spcalG1$y1=="G"),]

  classvalG=classclo[idvalG]
  classvalG=droplevels(classvalG)

  classcalG1=classclo[-idvalG]
  classcalG=classcalG1[which(classcalG1=="222" | classcalG1=="509" | classcalG1=="787")]
  classcalG=droplevels(classcalG)



  spvalS=sp[idvalS,]
  spcalS1=sp[-idvalS,]
  spcalS=spcalS1[which(spcalS1$y1=="S"),]

  classvalS=classclo[idvalS]
  classvalS=droplevels(classvalS)

  classcalS1=classclo[-idvalS]
  classcalS=classcalS1[which(classcalS1=="471" | classcalS1=="525" | classcalS1=="747" | classcalS1=="877")]
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

  predmclo=as.data.frame(matrix(nrow = length(classvalclo), ncol = ncmax))
  predm0clo=as.data.frame(matrix(nrow = length(classcalclo), ncol = ncmax))
  # spcal=sp
  # spcaldef=spcal
## Boucle CV
  for (i in 1:k) {
    print(i)
    m=mstage(spcaldef,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,ndc)))
    spvalCV=getdata(spcaldef,m)[[2]]

    idvalCV =which(rownames(spcal)  %in%  rownames(spvalCV))
    idvalCVC =which(rownames(spcal)  %in%  rownames(spvalCV) & spcal$y1=="C")
    idvalCVG =which(rownames(spcal)  %in%  rownames(spvalCV) & spcal$y1=="G")
    idvalCVS =which(rownames(spcal)  %in%  rownames(spvalCV) & spcal$y1=="S")

    spcaldef=spcaldef[-(which(rownames(spcaldef)  %in%  rownames(spvalCV))),]

    # # En mettant une autre année en validation
    # idvalCV =which(substr(rownames(spcal),1,4)  %in%  '2018')
    spvalCV=spcal[idvalCV,]
    classvalCV=classcal[idvalCV]  #identifiants des classes du jeu de validation
    spcalCV=spcal[-idvalCV,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    classcalCV=classcal[-idvalCV] #identifiants des classes du jeu de calibration


    spvalCVC=spcal[idvalCVC,]
    classvalCVC=classcalclo[idvalCVC]  #identifiants des classes du jeu de validation

    spcalCVC1=spcal[-idvalCVC,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    spcalCVC=spcalCVC1[which(spcalCVC1$y1=="C"),]
    classcalCVC1=classcalclo[-idvalCVC] #identifiants des classes du jeu de calibration
    classcalCVC=classcalCVC1[which(classcalCVC1=="015" | classcalCVC1=="169" | classcalCVC1=="685")]
    classcalCVC=droplevels(classcalCVC)


    spvalCVG=spcal[idvalCVG,]
    classvalCVG=classcalclo[idvalCVG]  #identifiants des classes du jeu de validation

    spcalCVG1=spcal[-idvalCVG,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    spcalCVG=spcalCVG1[which(spcalCVG1$y1=="G"),]
    classcalCVG1=classcalclo[-idvalCVG] #identifiants des classes du jeu de calibration
    classcalCVG=classcalCVG1[which(classcalCVG1=="222" | classcalCVG1=="509" | classcalCVG1=="787")]
    classcalCVG=droplevels(classcalCVG)


    spvalCVS=spcal[idvalCVS,]
    classvalCVS=classcalclo[idvalCVS]  #identifiants des classes du jeu de validation

    spcalCVS1=spcal[-idvalCVS,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    spcalCVS=spcalCVS1[which(spcalCVS1$y1=="S"),]
    classcalCVS1=classcalclo[-idvalCVS] #identifiants des classes du jeu de calibration
    classcalCVS=classcalCVS1[which(classcalCVS1=="471" | classcalCVS1=="525" | classcalCVS1=="747" | classcalCVS1=="877")]
    classcalCVS=droplevels(classcalCVS)


    # ## PLSDA and application to have loadings and scores
    rplsda=caret::plsda(spcalCV$x, classcalCV,ncomp=ncmax)
    sccalCV=rplsda$scores
    spvalCV_c=scale(spvalCV$x,center=rplsda$Xmeans,scale = F)
    scvalCV=spvalCV_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

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
      predm0clo[idvalCVC,ii]=SIGNE_maha0(sccalCVC[,1:ii], classcalCVC, scvalCVC[,1:ii])$class
      predm0clo[idvalCVG,ii]=SIGNE_maha0(sccalCVG[,1:ii], classcalCVG, scvalCVG[,1:ii])$class
      predm0clo[idvalCVS,ii]=SIGNE_maha0(sccalCVS[,1:ii], classcalCVS, scvalCVS[,1:ii])$class
      #Ca le rempli (Peut-être même pas, en fait), mais pas complètement. Est-ce normal ? Embêtant ? P-ê predm0G est mal fait
#      predm0G[idvalCVG,ii]=SIGNE_maha0(sccalCVG[,1:ii], classcalCVG, scvalCVG[,1:ii])$classclo
      # M1= matrix(nrow= nrow(scvalCV[,1:ii]), ncol= nlevels(classcalCV))
      # M2= matrix(nrow= nrow(scvalCVG[,1:ii]), ncol= nlevels(classcalCVG))
    }
  }
  ## Table de contingence CV
  tsm0=lapply(as.list(predm0), classcal, FUN = table)
  tsm0clo=lapply(as.list(predm0clo), classcalclo, FUN = table)
  ## Matrice mauvais classements par clone CV
  diagsm0=lapply(tsm0, FUN = diag)
  diagsm0clo=lapply(tsm0clo, FUN = diag)
  ## Pourcentage de bien classes CV
  perokm0 =100*unlist(lapply(diagsm0, FUN = sum))/length(classcal)
  perokm0clo =100*unlist(lapply(diagsm0clo, FUN = sum))/length(classcalclo)
  # perokm0 =100*unlist(lapply(diagsm0, FUN = sum))/length(idvalCV)
  ### Enregistrement des matrices de resultat final CV
  ##Remplissage de la matrice des perok finale
  perok_finalm0[j,]=perokm0
  perok_finalm0clo[j,]=perokm0clo

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

plot(colMeans(perok_finalm0clo), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalmclo), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)

plot(colMeans(perok_finalmC), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalmG), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalmS), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)


#stop()







###PLSDA on Maha scores
## Calibration
rplsda=caret::plsda(spcal$x, classcal,ncomp=10)
sccal=rplsda$scores

## Validation
spval_c=scale(spval$x,center=rplsda$Xmeans,scale = F)
scval=spval_c%*%rplsda$projection
resval=SIGNE_maha0(sccal[,1:10], classcal, scval[,1:10])$class

cepage=table (resval,classval)
cepage
cepage[1,1]/(cepage[1,1]+cepage[1,2]+cepage[1,3])
cepage[2,2]/(cepage[2,1]+cepage[2,2]+cepage[2,3])
cepage[3,3]/(cepage[3,1]+cepage[3,2]+cepage[3,3])

