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
source("Script_R_2020/sp2df.R")

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
rownames(sp)
## Variable qui mesure le nombre de classes
c=length(levels(class))






# Création des jeux de calibration/ validation
# On créé un facteur datclone qui groupe un clone à 1 date
datclone=substr(rownames(sp),1,13)
ndc=length(unique(datclone))
# On créé un facteur souche qui groupe les 6 spectres de chaque souche
numsp=as.numeric(substr(rownames(sp),15,16))
souche=cut(numsp, breaks = c(0,6,12,18),labels=c("s1","s2","s3"))  # paste(datclone,cut(numsp, breaks = c(0,6,12,18),labels=c("s1","s2","s3")))
sp=sp2df(sp,class)
sp=cbind(sp,datclone,souche)   # mutate(sp,datclone=substr(titre,1,13), souche=substr(souche,15,16))
# Le tirage sera fait plus loin dans la boucle
sp$datclone
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

###s?paration validation calibration PLSDA###
#set.seed(1) # fixe le tirage aleatoire
for(j in 1:repet) {

  # On selectionne le jeu de validation de manière à ce que tous les datclone soient représentés et 1 souche sur les 3 tirée random
  m=mstage(sp,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,ndc)))
  spval=getdata(sp,m)[[2]]
  #
  idval=which(rownames(sp)  %in%  rownames(spval))
  #
  # ##On selectionne les spectres ayant ces num?ros dans le jeu de validation, les autres vont dans le jeu de calibration
  spval=sp[idval,]
  spcal=sp[-idval,]
  classval=class[idval]
  classcal=class[-idval]

  # ## Creation des jeux d'apprentissage et validation
  predm=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))
  predm0=as.data.frame(matrix(nrow = length(classcal), ncol = ncmax))
  spcaldef=spcal # spcal deflaté du(des) groupe(s) de CV déjà validés

  # spcal=sp
  # spcaldef=spcal
## Boucle CV
  for (i in 1:k) {
    print(i)
    m=mstage(spcaldef,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,ndc)))
    spvalCV=getdata(spcaldef,m)[[2]]

    idvalCV =which(rownames(spcal)  %in%  rownames(spvalCV))

    spcaldef=spcaldef[-(which(rownames(spcaldef)  %in%  rownames(spvalCV))),]

    # # En mettant une autre année en validation
    # idvalCV =which(substr(rownames(spcal),1,4)  %in%  '2018')
    spvalCV=spcal[idvalCV,]
    classvalCV=classcal[idvalCV]  #identifiants des classes du jeu de validation
    spcalCV=spcal[-idvalCV,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    classcalCV=classcal[-idvalCV] #identifiants des classes du jeu de calibration

    # ## PLSDA and application to have loadings and scores
    rplsda=caret::plsda(spcalCV$x, classcalCV,ncomp=ncmax)
    sccalCV=rplsda$scores
    spvalCV_c=scale(spvalCV$x,center=rplsda$Xmeans,scale = F)
    scvalCV=spvalCV_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
    for (ii in 2:ncmax) {
      ## Validation
      predm0[idvalCV,ii]=SIGNE_maha0(sccalCV[,1:ii], classcalCV, scvalCV[,1:ii])$class
    }
    ## Package rnirs
#     for (ii in 2:ncmax) {
#     rknnwda=knnwda(spcalCV$x, classcalCV,spvalCV$x, classvalCV, ncompdis = ii, diss = "mahalanobis",k=30,h=1)
#     print(err(rknnwda))
#     predm0[idvalCV,ii]=rknnwda$fit$y1
# }



  }

  ## Table de contingence CV
  tsm0=lapply(as.list(predm0), classcal, FUN = table)
  ## Matrice mauvais classements par clone CV
  diagsm0=lapply(tsm0, FUN = diag)
  ## Pourcentage de bien classes CV
  perokm0 =100*unlist(lapply(diagsm0, FUN = sum))/length(classcal)
  # perokm0 =100*unlist(lapply(diagsm0, FUN = sum))/length(idvalCV)
  ### Enregistrement des matrices de resultat final CV
  ##Remplissage de la matrice des perok finale
  perok_finalm0[j,]=perokm0

  # browser()
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

}

plot(colMeans(perok_finalm0), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
plot(colMeans(perok_finalm), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)


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



###En fonction des clones

## Creation de la matrice de classes clones
class_clones=as.factor(substr(rownames(sp),11,13))
## Variable qui mesure le nombre de classes
c=length(levels(class_clones))

## Separation de sp_cal en 3 jeux de calibration par cepage
aC= substr(rownames(spcal),9,9)=="C"
aG= substr(rownames(spcal),9,9)=="G"
aS= substr(rownames(spcal),9,9)=="S"
#aF= substr(rownames(spcal),9,9)=="F"

spcal_C =spcal[(aC==TRUE),]
spcal_G =spcal[(aG==TRUE),]
spcal_S =spcal[(aS==TRUE),]
#spcal_F =spcal[(aF==TRUE),]

###Identifiants des matrices de calibration
##Cabernet
idcal_C=which(rownames(sp)  %in%  rownames(spcal_C))
classcal_C=droplevels(class_clones[idcal_C])

##Gamay
idcal_G=which(rownames(sp)  %in%  rownames(spcal_G))
classcal_G=droplevels(class_clones[idcal_G])
# #CF
# idcal_F=which(rownames(sp)  %in%  rownames(spcal_F))
# classcal_F=droplevels(class_clones[idcal_F])

##Syrah
idcal_S=which(rownames(sp)  %in%  rownames(spcal_S))
classcal_S=droplevels(class_clones[idcal_S])




# ###CV sur clones
# ##Formation de idcalC (Cabernet)
#   c1=idcal[1,1]
#   c2=idcal[1,2]
#   c3=idcal[1,3]
#   c4=idcal[2,1]
#   c5=idcal[2,2]
#   c6=idcal[2,3]
#
# idcalC=c(c1,c2,c3,c4,c5,c6)
# dim(idcalC)=c(3,2)
# idcalC=t(idcalC)
#
# idcal2C=matrix( ,nrow=2, ncol=3)
#
# for (i in 1:3) {
#   idcal2C[,i]=idcalC[sample(1:2),i]
# }
#
# ## Boucle CV (Cabernet)
# for (i in 1:k) {
#   commdC=paste("spvalCVC=rbind(sp",idcal2C[i,1],",sp",idcal2C[i,2],",sp",idcal2C[i,3],")",sep="")
#   eval(parse(text=commdC))
#
#   idvalCVC=which(rownames(sp_cal_C)  %in%  rownames(spvalCVC))
#
#   spvalCVC=sp_cal_C[idvalCVC,]       # matrice du jeu de validation
#   class_valCVC=droplevels(class_cal_C[idvalCVC])  #identifiants des classes du jeu de validation
#   spcalCVC=sp_cal_C[-idvalCVC,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
#   class_calCVC=droplevels(class_cal_C[-idvalCVC]) #identifiants des classes du jeu de calibration
#
#   ## PLSDA and application to have loadings and scores (Cabernet)
#   rplsdaCVC=caret::plsda(spcalCVC, class_calCVC,ncomp=ncmax)
#   sccalCVC=rplsdaCVC$scores
#   spvalCVC_c=scale(spvalCVC,center=rplsdaCVC$Xmeans,scale = F)
#   scvalCVC=spvalCVC_c%*%rplsdaCVC$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
#
#   ## Creation des jeux d'apprentissage et validation
#   iok3=substr(rownames(sp_cal_C),1,9) %in% dates
#   predm0C=as.data.frame(matrix(nrow = sum(iok3), ncol = ncmax))
#
#
#   for (ii in 2:ncmax) {
#
#     ## Validation
#     predm0C[idvalCVC,ii]=SIGNE_maha0(sccalCVC[,1:ii], class_calCVC, scvalCVC[,1:ii])$class
#   }
# }
#
#
# ## Table de contingence CV
# tsm0C=lapply(as.list(predm0C), class_cal_C, FUN = table)
#
#
# ## Matrice mauvais classements par clone CV
# diagsm0C=lapply(tsm0C, FUN = diag)
#
# ## Pourcentage de bien classes CV
# perokm0C=100*unlist(lapply(diagsm0C, FUN = sum))/length(class_valCVC)
#
# ## Pourcentage de bien classes CV
# maxiC=max(perokm0C)
# maxi.idC=which.max(perokm0C)
#
# ### Enregistrement des matrices de resultat final
# ##Remplissage de la matrice des perok finale
# perok_finalm0C[j,]=perokm0C
#
# ## Remplissage de la VL max et de son % de bons classements globaux
# maxi_finalC[j,1]=maxi.idC
# maxi_finalC[j,2]=maxiC
#
# ##Formation de idcalG (Gamay)
#   g1=idcal[1,4]
#   g2=idcal[1,5]
#   g3=idcal[1,6]
#   g4=idcal[2,4]
#   g5=idcal[2,5]
#   g6=idcal[2,6]
#
# idcalG=c(g1,g2,g3,g4,g5,g6)
# dim(idcalG)=c(3,2)
# idcalG=t(idcalG)
#
# idcal2G=matrix( ,nrow=2, ncol=3)
#
# for (i in 1:3) {
#   idcal2G[,i]=idcalG[sample(1:2),i]
# }
#
# ## Boucle CV (Gamay)
# for (i in 1:k) {
#   commdG=paste("spvalCVG=rbind(sp",idcal2G[i,1],",sp",idcal2G[i,2],",sp",idcal2G[i,3],")",sep="")
#   eval(parse(text=commdG))
#
#   idvalCVG=which(rownames(sp_cal_G)  %in%  rownames(spvalCVG))
#
#   spvalCVG=sp_cal_G[idvalCVG,]       # matrice du jeu de validation
#   class_valCVG=droplevels(class_cal_G[idvalCVG])  #identifiants des classes du jeu de validation
#   spcalCVG=sp_cal_G[-idvalCVG,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
#   class_calCVG=droplevels(class_cal_G[-idvalCVG]) #identifiants des classes du jeu de calibration
#
#   ## PLSDA and application to have loadings and scores (Gamay)
#   rplsdaCVG=caret::plsda(spcalCVG, class_calCVG,ncomp=ncmax)
#   sccalCVG=rplsdaCVG$scores
#   spvalCVG_c=scale(spvalCVG,center=rplsdaCVG$Xmeans,scale = F)
#   scvalCVG=spvalCVG_c%*%rplsdaCVG$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
#
#   ## Creation des jeux d'apprentissage et validation
#   iok4=substr(rownames(sp_cal_G),1,9) %in% dates
#   predm0G=as.data.frame(matrix(nrow = sum(iok4), ncol = ncmax))
#
#   for (ii in 2:ncmax) {
#
#     ## Validation
#     predm0G[idvalCVG,ii]=SIGNE_maha0(sccalCVG[,1:ii], class_calCVG, scvalCVG[,1:ii])$class
#   }
# }
#
# ## Table de contingence CV
# tsm0G=lapply(as.list(predm0G), class_cal_G, FUN = table)
#
# ## Matrice mauvais classements par clone CV
# diagsm0G=lapply(tsm0G, FUN = diag)
#
# ## Pourcentage de bien classes CV
# perokm0G=100*unlist(lapply(diagsm0G, FUN = sum))/length(class_valCVG)
#
# ## Pourcentage de bien classes CV
# maxiG=max(perokm0G)
# maxi.idG=which.max(perokm0G)
#
# ### Enregistrement des matrices de resultat final
# ##Remplissage de la matrice des perok finale
# perok_finalm0G[j,]=perokm0G
#
# ## Remplissage de la VL max et de son % de bons classements globaux
# maxi_finalG[j,1]=maxi.idG
# maxi_finalG[j,2]=maxiG

# ##Formation de idcalS (Syrah)
#   s1=idcal[1,7]
#   s2=idcal[1,8]
# #  s3=idcal[1,9]
# #  s4=idcal[1,10]
#   s5=idcal[2,7]
#   s6=idcal[2,8]
# #  s7=idcal[2,9]
# #  s8=idcal[2,10]
#
# idcalS= c(s1,s2,s5,s6)#c(s1,s2,s3,s4,s5,s6,s7,s8)
# dim(idcalS)=c(2,2)
# idcalS=t(idcalS)
#
# idcal2S=matrix( ,nrow=2, ncol=2)
#
# for (i in 1:2) {
#   idcal2S[,i]=idcalS[sample(1:2),i]
# }

# ## Boucle CV  (Syrah)
# for (i in 1:k) {
#   #commdS=paste("spvalCVS=rbind(sp",idcal2S[i,1],",sp",idcal2S[i,2],",sp",idcal2S[i,3],",sp",idcal2S[i,4],")",sep="")
#   commdS=paste("spvalCVS=rbind(sp",idcal2S[i,1],",sp",idcal2S[i,2],")",sep="")
#   eval(parse(text=commdS))
#
#   idvalCVS=which(rownames(sp_cal_S)  %in%  rownames(spvalCVS))
#
#   spvalCVS=sp_cal_S[idvalCVS,]       # matrice du jeu de validation
#   class_valCVS=droplevels(class_cal_S[idvalCVS])  #identifiants des classes du jeu de validation
#   spcalCVS=sp_cal_S[-idvalCVS,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
#   class_calCVS=droplevels(class_cal_S[-idvalCVS]) #identifiants des classes du jeu de calibration
#
#   ## PLSDA and application to have loadings and scores (Syrah)
#   rplsdaCVS=caret::plsda(spcalCVS, class_calCVS,ncomp=ncmax)
#   sccalCVS=rplsdaCVS$scores
#   spvalCVS_c=scale(spvalCVS,center=rplsdaCVS$Xmeans,scale = F)
#   scvalCVS=spvalCVS_c%*%rplsdaCVS$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
#
#   ## Creation des jeux d'apprentissage et validation
#    iok5=substr(rownames(sp_cal_S),1,9) %in% dates
#    predm0S=as.data.frame(matrix(nrow = sum(iok5), ncol = ncmax))
#
#   for (ii in 2:ncmax) {
#
#     ## Validation
#     predm0S[idvalCVS,ii]=SIGNE_maha0(sccalCVS[,1:ii], class_calCVS, scvalCVS[,1:ii])$class
#   }
# }

# ## Table de contingence CV
# tsm0S=lapply(as.list(predm0S), class_cal_S, FUN = table)
#
# ## Matrice mauvais classements par clone CV
# diagsm0S=lapply(tsm0S, FUN = diag)
#
# ## Pourcentage de bien classes CV
# perokm0S=100*unlist(lapply(diagsm0S, FUN = sum))/length(class_valCVS)
#
# ## Pourcentage de bien classes CV
# maxiS=max(perokm0S)
# maxi.idS=which.max(perokm0S)
#
# ### Enregistrement des matrices de resultat final
# ##Remplissage de la matrice des perok finale
# perok_finalm0S[j,]=perokm0S
#
# ## Remplissage de la VL max et de son % de bons classements globaux
# maxi_finalS[j,1]=maxi.idS
# maxi_finalS[j,2]=maxiS

###PLSDA sur clones sans CV
## Separation de sp_val en 3 jeux de validation par cepage
## Seulement avec les biens classes en cepages
test=cbind(resval,classval)
rownames(test)=rownames(sp[idval,])

z1=(substr(rownames(test),9,9))=="C"
z2=(substr(rownames(test), 9, 9))=="G"
z3=(substr(rownames(test),9,9))=="S"
#z4=(substr(rownames(test),9,9))=="F"

test1 =test[(z1==TRUE),]
test2 =test[(z2==TRUE),]
test3 =test[(z3==TRUE),]
#test4 =test[(z4==TRUE),]

test1a = test1[test1[,1]==1,]
test2a = test2[test2[,1]==2,]
test3a = test3[test3[,1]==3,]
#test4a = test4[test4[,1]==2,]

spval_C =spval[rownames(test1a),]
spval_G =spval[rownames(test2a),]
spval_S =spval[rownames(test3a),]
#spval_F =spval[rownames(test4a),]

##Avec tous
# bC= substr(rownames(sp_val),9,9)=="C"
# bG= substr(rownames(sp_val),9,9)=="G"
# bS= substr(rownames(sp_val),9,9)=="S"
#
# sp_val_C =sp_val[(bC==TRUE),]
# sp_val_G =sp_val[(bG==TRUE),]
# sp_val_S =sp_val[(bS==TRUE),]

###Selection des spectres
##Cabernet Sauvignon
idval_C=which(rownames(sp)  %in%  rownames(spval_C))
classval_C=droplevels(class_clones[idval_C])

##Gamay
idval_G=which(rownames(sp)  %in%  rownames(spval_G))
classval_G=droplevels(class_clones[idval_G])

# ##Cabernet franc
# idval_F=which(rownames(sp)  %in%  rownames(spval_F))
# classval_F= droplevels(class_clones[idval_F])

##Syrah
idval_S=which(rownames(sp)  %in%  rownames(spval_S))
classval_S=droplevels(class_clones[idval_S])



####PLSDA on Maha scores
### Calibration
## Cabernet Sauvignon
rplsdaC=caret::plsda(spcal_C$x, classcal_C,ncomp= 8)# Modifier ncmax en fonction des resultats de CV_2
sccal_C=rplsdaC$scores

##Gamay
rplsdaG=caret::plsda(spcal_G$x, classcal_G,ncomp=8)# Modifier ncmax en fonction des resultats de CV_2
sccal_G=rplsdaG$scores

##Syrah
rplsdaS=caret::plsda(spcal_S$x, classcal_S,ncomp=8)# Modifier ncmax en fonction des resultats de CV_2
sccal_S=rplsdaS$scores

##Cabernet franc
#rplsdaF=caret::plsda(spcal_F$x, classcal_F,ncomp=8)# Modifier ncmax en fonction des resultats de CV_2
#sccal_F=rplsdaF$scores

### Validation
## Cabernet Sauvignon
spval_c_C=scale(spval_C$x,center=rplsdaC$Xmeans,scale = F)
scval_C=spval_c_C%*%rplsdaC$projection
resval_C=SIGNE_maha0(sccal_C[,1:8], classcal_C, scval_C[,1:8])$class

Cabernet=table (resval_C,classval_C)
print (Cabernet)

##Gamay
spval_c_G=scale(spval_G$x,center=rplsdaG$Xmeans,scale = F)
scval_G=spval_c_G%*%rplsdaG$projection
resval_G=SIGNE_maha0(sccal_G[,1:8], classcal_G, scval_G[,1:8])$class

Gamay=table (resval_G,classval_G)
print (Gamay)

##Syrah
spval_c_S=scale(spval_S$x,center=rplsdaS$Xmeans,scale = F)
scval_S=spval_c_S%*%rplsdaS$projection
resval_S=SIGNE_maha0(sccal_S[,1:8], classcal_S, scval_S[,1:8])$class

Syrah=table (resval_S,classval_S)
print (Syrah)

# ##Cabernet franc
# spval_c_F=scale(sp_val_F,center=rplsdaF$Xmeans,scale = F)
# scval_F=spval_c_F%*%rplsdaF$projection
# resval_F=SIGNE_maha0(sccal_F[,1:8], classcal_F, scval_F[,1:8])$class
#
# CF=table (res_val_F,class_val_F)
# #print (CF)

###Calcul des pourcentages de bons classements en fonction du tirage
##Somme des clones biens classes
clones_bc= sum(diag(Cabernet))+sum(diag(Gamay))+sum(diag(Syrah))
#clones_bc= sum(diag(Cabernet))+sum(diag(Syrah))+sum(diag(CF))
#print(clones_bc)

#Somme colonne table contingence cepages
sumcolcep=apply(cepage,2,sum)

##Nombre total de clones
total=clones_bc+(sum(cepage)-sum(diag(cepage)))+(sum(Cabernet)-sum(diag(Cabernet)))+(sum(Gamay)-sum(diag(Gamay)))+(sum(Syrah)-sum(diag(Syrah)))
#total=clones_bc+(sum(cepage)-sum(diag(cepage)))+(sum(Cabernet)-sum(diag(Cabernet)))+(sum(Syrah)-sum(diag(Syrah)))+(sum(CF)-sum(diag(CF)))
total_C=sum(Cabernet)
total_G=sum(Gamay)
total_S=sum(Syrah)
#total_F=sum(CF)

total_C_cep=sum(diag(Cabernet))+(sum(Cabernet)-sum(diag(Cabernet)))+(sumcolcep[1]-cepage[1,1])
total_G_cep=sum(diag(Gamay))+(sum(Gamay)-sum(diag(Gamay)))+(sumcolcep[2]-cepage[2,2])
total_S_cep=sum(diag(Syrah))+(sum(Syrah)-sum(diag(Syrah)))+(sumcolcep[3]-cepage[3,3])
#total_F_cep=sum(diag(CF))+(sum(CF)-sum(diag(CF)))+(sumcolcep[2]-cepage[2,2])
#print(total)



## Pourcentage total de clones biens classes(prise en compte erreur cepages)
perok=100*(clones_bc/total)
#print(perok)
perok_final=perok

## Pourcentage de clones biens classes
perok_C=100*(sum(diag(Cabernet))/total_C)
perok_G=100*(sum(diag(Gamay))/total_G)
perok_S=100*(sum(diag(Syrah))/total_S)
#perok_F=100*(sum(diag(CF))/total_F)

perok_final_C[j,]=perok_C
perok_final_G[j,]=perok_G
perok_final_S[j,]=perok_S
#perok_final_F[j,]=perok_F

## Pourcentage de clones biens classes (prise en compte erreur cepages)
perok_C_cep=100*(sum(diag(Cabernet))/total_C_cep)
perok_G_cep=100*(sum(diag(Gamay))/total_G_cep)
perok_S_cep=100*(sum(diag(Syrah))/total_S_cep)
#perok_F_cep=100*(sum(diag(CF))/total_F_cep)

perok_final_C_cep[j,]=perok_C_cep
perok_final_G_cep[j,]=perok_G_cep
perok_final_S_cep[j,]=perok_S_cep
#perok_final_F_cep[j,]=perok_F_cep

##Pourcentage de cepages biens classes
perok_cepages=100*(sum(diag(cepage))/sum(cepage))
#print(perok_cepages)
perok_final_cepages[j,]=perok_cepages



#print(perok_final)
#print(perok_final_cepages)
print(perok_cepages)
#print(mean(perok_final))
print(perok_final)
#print(mean(perok_final_cepages))
#print(mean(perok_final_C))
print(perok_C)
#print(mean(perok_final_G))
print(perok_G)
#print(mean(perok_final_S))
print(perok_S)
#print(mean(perok_final_F))
#print(mean(perok_final_C_cep))
print(perok_C_cep)
#print(mean(perok_final_G_cep))
print(perok_G_cep)
#print(mean(perok_final_S_cep))
print(perok_S_cep)
#print(mean(perok_final_F_cep))

###Sorties graphiques
## Tracage de l'evolution des perok en fonction du nombre de VL utilisees (cepages)
plot(colMeans(perok_finalm0), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
legend(ncmax*2/3,15,legend=c("Maha on PLSDA scores"),
       col=c("black"), lty=1, cex=0.8)

# ## Tracage de l'evolution des perok en fonction du nombre de VL utilisees (Cabernet)
# plot(colMeans(perok_finalm0C), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
# legend(ncmax*2/3,15,legend=c("Maha on PLSDA scores"),
#        col=c("black"), lty=1, cex=0.8)
#
# ## Tracage de l'evolution des perok en fonction du nombre de VL utilisees (Gamay)
# plot(colMeans(perok_finalm0G), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
# legend(ncmax*2/3,15,legend=c("Maha on PLSDA scores"),
#        col=c("black"), lty=1, cex=0.8)
#
# ## Tracage de l'evolution des perok en fonction du nombre de VL utilisees (Syrah)
# plot(colMeans(perok_finalm0S), xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
# legend(ncmax*2/3,15,legend=c("Maha on PLSDA scores"),
#        col=c("black"), lty=1, cex=0.8)

## Clones


plot(perok_final, type="o",xlab= "Nombre de tirages", ylab = "Pourcentage de clones biens class?s",pch=19, cex=2)

##Cepages
plot(perok_final_cepages, type="o",xlab= "Nombre de tirages", ylab = "Pourcentage de c?pages biens class?s",pch=21, cex=2,bg="blue")


