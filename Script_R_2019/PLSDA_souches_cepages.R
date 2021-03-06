library(MASS)
# library(mixOmics)
library(FactoMineR)
library(signal)
library(plyr)
library(caret)
library(dplyr)
library(prospectr)
library(sampling)

# library(nirs)

rm(list = ls())

source('Script_R_2019/adj_asd.R')
source('Script_R_2019/SIGNE_load.R')
# source('C:/Users/No?mie/Desktop/SFE/Script_R/SIGNE_maha.R')
source('Script_R_2019/SIGNE_maha0.R')
source("Script_R_2019/sp2df.R")


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
#set.seed(1)

# brb3="~/Documents/NICOLAS/Stage de fin annee au Grau du roi/globalmatrixN1"
brb3="~/Documents/INRA/Projets/SIGNE/2018/GrauduRoi/PTN/globalmatrixN1"
load(file=brb3)

## Filtrage des spectres aberrants
globalmatrixN1=globalmatrixN1[globalmatrixN1[,500]>0.5,]
globalmatrixN1=globalmatrixN1[globalmatrixN1[,1]<0.2,]
globalmatrixN1=globalmatrixN1[globalmatrixN1[,2000]<0.25,]

# Data Filter
# Select dates
dates=list(
  # "20180619N",
  # "20180627N",
  # "20180619P",
  # "20180619T",
  "20180627P",
  # "20180627T",
  # "20180709N",
  # "20180709T",
  "20180709P")
  # "20180816T",
  # "20180816P")
 # "20180816N")

  # dates[1]="20180619N",
  # dates[2]="20180627N",
  # dates[3]="20180619P",
  # dates[4]="20180619T",
  # dates[5]="20180627P",
  # dates[6]="20180627T",
  # dates[7]="20180709N",
  # dates[8]="20180709T",
  # dates[9]="20180709P",
  # # dates[10]="20180709N"
  # # dates[11]="20180710P"
  # dates[10]="20180816T",
  # dates[11]="20180816P",
  # dates[12]="20180816N")
 # ,"20180817N"
  # "20170524P"
  # , "20170529P"
   # ,"20170606P"
 #   ,"20170612P"
  #   "20170619P"
 #    ,"20170626P"
  #  "20170703P"
 # ,"20170710P"
 #  ,"20170717P"
 #  ,"20170724P"
 #  ,"20170731P"
  # "20180823P"
  # ,"20180823T"
  # ,"20180823N"
 #  , "20180731P"
  # ,"20180731T"
  # ,"20180731N"
 #  ,"20180810P"
  # ,"20180810T"
  # ,"20180810N"
 #  ,"20180724P"
  # ,"20180724T"
  # ,"20180724N"
 # "20180809T"
 #  , "20180822P"
 #  , "20180730P"


iok=substr(rownames(globalmatrixN1),1,9) %in% dates
sp=globalmatrixN1[iok,]

###Ajout du cepage au nom de la ligne. !!Attention a bien commenter en fonction des cepages presents dans la BDD!!
##Filtre en fonction du cepage

titre=rownames(sp)
yC= substr(titre,11,13)== "015" |  substr(titre,11,13)== "169" |  substr(titre,11,13)== "685"
yG= substr(titre,11,13)== "222" | substr(titre,11,13)== "509" |  substr(titre,11,13)== "787"
yS= substr(titre,11,13)== "471" |  substr(titre,11,13)== "525" |  substr(titre,11,13)== "747" |  substr(titre,11,13)== "877"
#yc=(substr(titre,11,13))== "471" |  (substr(titre,11,13))== "747" |  (substr(titre,11,13))== "877" #Si on supprime le 525
#yc=(substr(titre,11,13))== "747" |  (substr(titre,11,13))== "877"# Si on compare le 877 et le 747
#yc=(substr(titre,11,13))== "471" |  (substr(titre,11,13))== "525"

##Rajoute le nom du cepage au nom de la ligne
rownames(sp)[yC]=paste(rownames(sp)[yC],"C",sep = "-")
rownames(sp)[yG]=paste(rownames(sp)[yG],"G",sep = "-")
rownames(sp)[yS]=paste(rownames(sp)[yS],"S",sep = "-")

## Creation de la matrice de classes
class=as.factor(substr(rownames(sp),18,18))
## Variable qui mesure le nombre de classes
c=length(levels(class))


# Création des jeux de calibration/ validation
# On créé un facteur datclone qui groupe un clone à 1 date
datclone=substr(titre,1,13)
ndc=length(unique(datclone))
# On créé un facteur souche qui groupe les 6 spectres de chaque souche
numsp=as.numeric(substr(titre,15,16))
souche=cut(numsp, breaks = c(0,6,12,18),labels=c("s1","s2","s3"))  # paste(datclone,cut(numsp, breaks = c(0,6,12,18),labels=c("s1","s2","s3")))
sp=sp2df(sp,class)
sp=cbind(sp,datclone,souche)   # mutate(sp,datclone=substr(titre,1,13), souche=substr(souche,15,16))
# Le tirage sera fait plus loin dans la boucle


### FIXATION DES PARAMETRES UTILISES:
## Nombre de repetitions de la boucle de PLSDA:
repet= 10
## Parametres du Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
p=2
n=11
m=1
## Nombre de VL max autorisees
ncmax=20
## Taille de l'echantillon de validation (1/v):
v=3
## Nombre de groupes de CV
k=2

## PLSDA ##

### Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
sp_pre=adj_asd(sp$x,c(602,1402))
## Reduction des variables (extremites bruitees)
# sp=sp[,seq(51,ncol(sp)-30,1)]
## Coupure du spectre a 1000nm
#spx=sp[,seq(+1,601,1)]
#matplot(t(spx),pch = ".",xlab = "Longueurs d'ondes (nm)", ylab = "Transflectance")
## SNV
sp_pre=t(scale(t(sp_pre)))
## Derivation Savitsky Golay
sp$x=savitzkyGolay(sp_pre, m = m, p = p, w = n)

## Definition des matrices de resultat final
# Creation de la matrice des perok finale
perok_final=matrix(nrow = repet, ncol = 1)
perok_final_C=matrix(nrow = repet, ncol = 1)
perok_final_G=matrix(nrow = repet, ncol = 1)
perok_final_S=matrix(nrow = repet, ncol = 1)
perok_final_F=matrix(nrow = repet, ncol = 1)
perok_final_C_cep=matrix(nrow = repet, ncol = 1)
perok_final_G_cep=matrix(nrow = repet, ncol = 1)
perok_final_S_cep=matrix(nrow = repet, ncol = 1)
perok_final_F_cep=matrix(nrow = repet, ncol = 1)
perok_final_cepages=matrix(nrow = repet, ncol = 1)
perok_finalm0=matrix(nrow = repet, ncol = ncmax)
perok_finalm0C=matrix(nrow = repet, ncol = ncmax)
perok_finalm0G=matrix(nrow = repet, ncol = ncmax)
perok_finalm0S=matrix(nrow = repet, ncol = ncmax)
perok_finalm0F=matrix(nrow = repet, ncol = ncmax)
## Creation matrice de % de mauvais classements par clone
mc=matrix(nrow = ncmax,ncol = c)

## Creation de la matrice des VL et perok maximaux
maxi_final=matrix(nrow= repet, ncol = 2)
maxi_finalC=matrix(nrow= repet, ncol = 2)
#maxi_finalG=matrix(nrow= repet, ncol = 2)
maxi_finalS=matrix(nrow= repet, ncol = 2)
maxi_finalF=matrix(nrow= repet, ncol = 2)
## Creation de la matrice de % de mauvais classements
mc_final=matrix(nrow= repet, ncol = length(levels(class)))
## Creation d'un matrice cubique pour enregistrer les tables de contingence
t_final=array(dim=c(c,c,repet))
## Noms des colonnes et des lignes
colnames(t_final)=c(basename(levels(class)))
rownames(t_final)=c(basename(levels(class)))
colnames(maxi_final)= c("maxi.id","perok max")
colnames(mc_final)= c(basename(levels(class)))

###s?paration validation calibration PLSDA###
#set.seed(1) # fixe le tirage aleatoire
for(j in 1:repet) {
#
# id2=matrix( ,nrow=1, ncol=10)
# #  for (i in 1:3) {
# #    id2[,i]=sample(id[,i],1, replace = FALSE)
#
#
#  idval=matrix( ,nrow=1, ncol=8) #matrice validation PLSDA
#  idcal=matrix( ,nrow=2, ncol=8) #matrice calibration PLSDA
#  idcal2=matrix( ,nrow=2, ncol=8) #matrice calibration PLSDA
#  # idval2=matrix( ,nrow=1, ncol=10) #idval pour CV
#  # idcal2=matrix( ,nrow=1, ncol=10) #idcal pour CV
#
#   for (i in 1:8) {
#     icol=sample(1:3,1)
#     # id2[,i]=sample(id[,i],1, replace = FALSE)
#     idval[,i]=id[icol,i]
#     idcal[,i]=id[c(1:3)[-icol],i]
#
#     }
#
#  for (i in 1:8) {
#    # icol2=sample(1:2,1)
#    idcal2[,i]=idcal[sample(1:2),i]
#    # idval2[,i]=idcal[icol2,i]
#    # idcal2[,i]=idcal[c(1:2)[-icol2],i]

     # }

 ###PLSDA cepages
#
# #`commd2=paste("sp_val=rbind(sp",idval[1,1],",sp",idval[1,2],",sp",idval[1,3],",sp",idval[1,4],",sp",idval[1,5],",sp",idval[1,6],",sp",idval[1,7],",sp",idval[1,8],",sp",idval[1,9],",sp",idval[1,10],")",sep="")
#  #commd2=paste("sp_val=rbind(sp",idval[1,1],",sp",idval[1,2],",sp",idval[1,3],",sp",idval[1,4],",sp",idval[1,5],",sp",idval[1,6],",sp",idval[1,7],",sp",idval[1,8],",sp",idval[1,9],")",sep="") # Si on supprime 1 clone
#   commd2=paste("sp_val=rbind(sp",idval[1,1],",sp",idval[1,2],",sp",idval[1,3],",sp",idval[1,4],",sp",idval[1,5],",sp",idval[1,6],",sp",idval[1,7],",sp",idval[1,8],")",sep="")
#  eval(parse(text=commd2))


 # On selectionne le jeu de validation de manière à ce que tous les datclone soient représentés et 1 souche sur les 3 tirée random
 m=mstage(sp,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,30)))
 spval=getdata(sp,m)[[2]]

 id_val=which(rownames(sp)  %in%  rownames(spval))

 ##On selectionne les spectres ayant ces num?ros dans le jeu de validation, les autres vont dans le jeu de calibration
 spval=sp[id_val,]
 spcal=sp[-id_val,]
 classval=class[id_val]
 classcal=class[-id_val]

 ## Creation des jeux d'apprentissage et validation
 predm0=as.data.frame(matrix(nrow = length(classcal), ncol = ncmax))
 spcaldef=spcal # spcal deflaté du(des) groupe(s) de CV déjà validés
## Boucle CV
   for (i in 1:k) {
    #commd=paste("sp_val=rbind(sp",idcal2[i,1],",sp",idcal2[i,2],",sp",idcal2[i,3],")",sep="")   # pour CS et G
   # commd=paste("spvalCV=rbind(sp",idcal2[i,1],",sp",idcal2[i,2],",sp",idcal2[i,3],",sp",idcal2[i,4],")",sep="") # pour S
    #commd=paste("spvalCV=rbind(sp",idcal2[i,1],",sp",idcal2[i,2],",sp",idcal2[i,3],",sp",idcal2[i,4],",sp",idcal2[i,5],",sp",idcal2[i,6],",sp",idcal2[i,7],",sp",idcal2[i,8],",sp",idcal2[i,9],",sp",idcal2[i,10],")",sep="")
   # commd=paste("spvalCV=rbind(sp",idcal2[i,1],",sp",idcal2[i,2],",sp",idcal2[i,3],",sp",idcal2[i,4],",sp",idcal2[i,5],",sp",idcal2[i,6],",sp",idcal2[i,7],",sp",idcal2[i,8],")",sep="")
   #  # commd=paste("spvalCV=rbind(sp",idcal2[i,1],",sp",idcal2[i,2],",sp",idcal2[i,3],",sp",idcal2[i,4],",sp",idcal2[i,5],",sp",idcal2[i,6],",sp",idcal2[i,7],",sp",idcal2[i,8],",sp",idcal2[i,9],")",sep="")
   #   eval(parse(text=commd))
     m=mstage(spcaldef,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,30)))
     spvalCV=getdata(spcaldef,m)[[2]]

    idvalCV =which(rownames(spcal)  %in%  rownames(spvalCV))

    spcaldef=spcaldef[-(which(rownames(spcaldef)  %in%  rownames(spvalCV))),]

    # spvalCV=sp_cal[idvalCV,]       # matrice du jeu de validation
    classvalCV=classcal[idvalCV]  #identifiants des classes du jeu de validation
    spcalCV=spcal[-idvalCV,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
    classcalCV=classcal[-idvalCV] #identifiants des classes du jeu de calibration

    ## PLSDA and application to have loadings and scores
    rplsda=caret::plsda(spcalCV$x, classcalCV,ncomp=ncmax)
    sccalCV=rplsda$scores
    spvalCV_c=scale(spvalCV$x,center=rplsda$Xmeans,scale = F)
    scvalCV=spvalCV_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

      for (ii in 2:ncmax) {

     ## Validation
     predm0[idvalCV,ii]=SIGNE_maha0(sccalCV[,1:ii], classcalCV, scvalCV[,1:ii])$class
       }
   }

 ## Table de contingence CV
  tsm0=lapply(as.list(predm0), classcal, FUN = table)

 ## Matrice mauvais classements par clone CV
  diagsm0=lapply(tsm0, FUN = diag)

 ## Pourcentage de bien classes CV
  perokm0 =100*unlist(lapply(diagsm0, FUN = sum))/length(classcal)

 ## Pourcentage de bien classes CV
  maxi=max(perokm0)
  maxi.id=which.max(perokm0)

 ### Enregistrement des matrices de resultat final
 ##Remplissage de la matrice des perok finale
 perok_finalm0[j,]=perokm0

 ## Remplissage de la VL max et de son % de bons classements globaux
 maxi_final[j,1]=maxi.id
 maxi_final[j,2]=maxi
 ## Remplissage de la matrice de mauvais classements par clone
 mc_final[j,]=mc[maxi.id,]


###PLSDA on Maha scores
## Calibration
rplsda=caret::plsda(sp_cal, class_cal,ncomp=10)
sc_cal=rplsda$scores

## Validation
sp_val_c=scale(sp_val,center=rplsda$Xmeans,scale = F)
sc_val=sp_val_c%*%rplsda$projection
res_val=SIGNE_maha0(sc_cal[,1:10], class_cal, sc_val[,1:10])$class

cepage=table (res_val,class_val)
#print (cepage)


###En fonction des clones

## Creation de la matrice de classes clones
class_clones=as.factor(substr(rownames(sp),11,13))
## Variable qui mesure le nombre de classes
c=length(levels(class_clones))

## Separation de sp_cal en 3 jeux de calibration par cepage
aC= substr(rownames(sp_cal),18,18)=="C"
#aG= substr(rownames(sp_cal),18,18)=="G"
aS= substr(rownames(sp_cal),18,18)=="S"
aF= substr(rownames(sp_cal),18,18)=="F"

sp_cal_C =sp_cal[(aC==TRUE),]
#sp_cal_G =sp_cal[(aG==TRUE),]
sp_cal_S =sp_cal[(aS==TRUE),]
sp_cal_F =sp_cal[(aF==TRUE),]

###Identifiants des matrices de calibration
##Cabernet
id_cal_C=which(rownames(sp)  %in%  rownames(sp_cal_C))
class_cal_C=droplevels(class_clones[id_cal_C])

# ##Gamay
# id_cal_G=which(rownames(sp)  %in%  rownames(sp_cal_G))
# class_cal_G=droplevels(class_clones[id_cal_G])
#CF
id_cal_F=which(rownames(sp)  %in%  rownames(sp_cal_F))
class_cal_F=droplevels(class_clones[id_cal_F])

##Syrah
id_cal_S=which(rownames(sp)  %in%  rownames(sp_cal_S))
class_cal_S=droplevels(class_clones[id_cal_S])




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
test=cbind(res_val,class_val)
rownames(test)=rownames(sp[id_val,])

z1=(substr(rownames(test),18,18))=="C"
#z2=(substr(rownames(test), 18, 18))=="G"
z3=(substr(rownames(test),18,18))=="S"
z4=(substr(rownames(test),18,18))=="F"

test1 =test[(z1==TRUE),]
#test2 =test[(z2==TRUE),]
test3 =test[(z3==TRUE),]
test4 =test[(z4==TRUE),]

test1a = test1[test1[,1]==1,]
#test2a = test2[test2[,1]==2,]
#test3a = test3[test3[,1]==3,]
test3a = test3[test3[,1]==3,]
test4a = test4[test4[,1]==2,]

sp_val_C =sp_val[rownames(test1a),]
#sp_val_G =sp_val[rownames(test2a),]
sp_val_S =sp_val[rownames(test3a),]
sp_val_F =sp_val[rownames(test4a),]

##Avec tous
# bC= substr(rownames(sp_val),18,18)=="C"
# bG= substr(rownames(sp_val),18,18)=="G"
# bS= substr(rownames(sp_val),18,18)=="S"
#
# sp_val_C =sp_val[(bC==TRUE),]
# sp_val_G =sp_val[(bG==TRUE),]
# sp_val_S =sp_val[(bS==TRUE),]

###Selection des spectres
##Cabernet Sauvignon
id_val_C=which(rownames(sp)  %in%  rownames(sp_val_C))
class_val_C=droplevels(class_clones[id_val_C])

# ##Gamay
# id_val_G=which(rownames(sp)  %in%  rownames(sp_val_G))
# class_val_G=droplevels(class_clones[id_val_G])
##Cabernet franc
id_val_F=which(rownames(sp)  %in%  rownames(sp_val_F))
class_val_F= droplevels(class_clones[id_val_F])

##Syrah
id_val_S=which(rownames(sp)  %in%  rownames(sp_val_S))
class_val_S=droplevels(class_clones[id_val_S])



####PLSDA on Maha scores
### Calibration
## Cabernet Sauvignon
rplsdaC=caret::plsda(sp_cal_C, class_cal_C,ncomp= 8)# Modifier ncmax en fonction des resultats de CV_2
sc_cal_C=rplsdaC$scores

# ##Gamay
# rplsdaG=caret::plsda(sp_cal_G, class_cal_G,ncomp=8)# Modifier ncmax en fonction des resultats de CV_2
# sc_cal_G=rplsdaG$scores

##Syrah
rplsdaS=caret::plsda(sp_cal_S, class_cal_S,ncomp=8)# Modifier ncmax en fonction des resultats de CV_2
sc_cal_S=rplsdaS$scores

##Cabernet franc
rplsdaF=caret::plsda(sp_cal_F, class_cal_F,ncomp=8)# Modifier ncmax en fonction des resultats de CV_2
sc_cal_F=rplsdaF$scores

### Validation
## Cabernet Sauvignon
sp_val_c_C=scale(sp_val_C,center=rplsdaC$Xmeans,scale = F)
sc_val_C=sp_val_c_C%*%rplsdaC$projection
res_val_C=SIGNE_maha0(sc_cal_C[,1:8], class_cal_C, sc_val_C[,1:8])$class

Cabernet=table (res_val_C,class_val_C)
#print (Cabernet)

# ##Gamay
# sp_val_c_G=scale(sp_val_G,center=rplsdaG$Xmeans,scale = F)
# sc_val_G=sp_val_c_G%*%rplsdaG$projection
# res_val_G=SIGNE_maha0(sc_cal_G[,1:8], class_cal_G, sc_val_G[,1:8])$class
#
# Gamay=table (res_val_G,class_val_G)
# #print (Gamay)

##Syrah
sp_val_c_S=scale(sp_val_S,center=rplsdaS$Xmeans,scale = F)
sc_val_S=sp_val_c_S%*%rplsdaS$projection
res_val_S=SIGNE_maha0(sc_cal_S[,1:8], class_cal_S, sc_val_S[,1:8])$class

Syrah=table (res_val_S,class_val_S)
#print (Syrah)

##Cabernet franc
sp_val_c_F=scale(sp_val_F,center=rplsdaF$Xmeans,scale = F)
sc_val_F=sp_val_c_F%*%rplsdaF$projection
res_val_F=SIGNE_maha0(sc_cal_F[,1:8], class_cal_F, sc_val_F[,1:8])$class

CF=table (res_val_F,class_val_F)
#print (CF)

###Calcul des pourcentages de bons classements en fonction du tirage
##Somme des clones biens classes
#clones_bc= sum(diag(Cabernet))+sum(diag(Gamay))+sum(diag(Syrah))
clones_bc= sum(diag(Cabernet))+sum(diag(Syrah))+sum(diag(CF))
#print(clones_bc)

#Somme colonne table contingence cepages
sumcolcep=apply(cepage,2,sum)

##Nombre total de clones
#total=clones_bc+(sum(cepage)-sum(diag(cepage)))+(sum(Cabernet)-sum(diag(Cabernet)))+(sum(Gamay)-sum(diag(Gamay)))+(sum(Syrah)-sum(diag(Syrah)))
total=clones_bc+(sum(cepage)-sum(diag(cepage)))+(sum(Cabernet)-sum(diag(Cabernet)))+(sum(Syrah)-sum(diag(Syrah)))+(sum(CF)-sum(diag(CF)))
total_C=sum(Cabernet)
#total_G=sum(Gamay)
total_S=sum(Syrah)
total_F=sum(CF)

total_C_cep=sum(diag(Cabernet))+(sum(Cabernet)-sum(diag(Cabernet)))+(sumcolcep[1]-cepage[1,1])
#total_G_cep=sum(diag(Gamay))+(sum(Gamay)-sum(diag(Gamay)))+(sumcolcep[2]-cepage[2,2])
#total_S_cep=sum(diag(Syrah))+(sum(Syrah)-sum(diag(Syrah)))+(sumcolcep[3]-cepage[3,3])
total_S_cep=sum(diag(Syrah))+(sum(Syrah)-sum(diag(Syrah)))+(sumcolcep[3]-cepage[3,3])
total_F_cep=sum(diag(CF))+(sum(CF)-sum(diag(CF)))+(sumcolcep[2]-cepage[2,2])
#print(total)

## Pourcentage total de clones biens classes(prise en compte erreur cepages)
perok=100*(clones_bc/total)
#print(perok)
perok_final[j,]=perok

## Pourcentage de clones biens classes
perok_C=100*(sum(diag(Cabernet))/total_C)
#perok_G=100*(sum(diag(Gamay))/total_G)
perok_S=100*(sum(diag(Syrah))/total_S)
perok_F=100*(sum(diag(CF))/total_F)

perok_final_C[j,]=perok_C
#perok_final_G[j,]=perok_G
perok_final_S[j,]=perok_S
perok_final_F[j,]=perok_F

## Pourcentage de clones biens classes (prise en compte erreur cepages)
perok_C_cep=100*(sum(diag(Cabernet))/total_C_cep)
#perok_G_cep=100*(sum(diag(Gamay))/total_G_cep)
perok_S_cep=100*(sum(diag(Syrah))/total_S_cep)
perok_F_cep=100*(sum(diag(CF))/total_F_cep)

perok_final_C_cep[j,]=perok_C_cep
# perok_final_G_cep[j,]=perok_G_cep
perok_final_S_cep[j,]=perok_S_cep
perok_final_F_cep[j,]=perok_F_cep

##Pourcentage de cepages biens classes
perok_cepages=100*(sum(diag(cepage))/sum(cepage))
#print(perok_cepages)
perok_final_cepages[j,]=perok_cepages

}

# print(perok_final)
# print(perok_final_cepages)
print(mean(perok_final))
print(mean(perok_final_cepages))
print(mean(perok_final_C))
#print(mean(perok_final_G))
print(mean(perok_final_S))
print(mean(perok_final_F))
print(mean(perok_final_C_cep))
#print(mean(perok_final_G_cep))
print(mean(perok_final_S_cep))
print(mean(perok_final_F_cep))

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

