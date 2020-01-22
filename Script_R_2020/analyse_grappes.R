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
library(readr)

rm(list = ls())


source('Script_R_2020/adj_asd.R')
source('Script_R_2020/SIGNE_load.R')
source('Script_R_2020/SIGNE_maha0.R')
source("Script_R_2020/sp2dfclo.R")
source("Script_R_2020/sp2df.R")
source("Script_R_2020/segmFact.R")

##Fonctions utiles :
#plotsp, fonction d'affichage du package rnirs.


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
#set.seed(1)

brb3="C:/Users/avitvale/Documents/Test/globalmatrixconservatoiregrappes"
load(file=brb3)
sp=globalmatrixconservatoiregrappes

code=substr(rownames(sp),11,13)
rep=substr(rownames(sp),15,16)
sp2=sp2dfclo(sp, code, rep)
names(sp2)=c("code", "rep", "spectre")
sp2$code=as.numeric(as.character(sp2$code))
sp2$rep=as.numeric(as.character(sp2$rep))
#
# T=sp2[1,]
# SP=sp2[1,]
# SP$rep=NA
# Truc=sp2$spectre
#
# sp3=data.frame(NA,NA,NA)
# names(sp3)=c("code","rep","spectre")
# for (i in unique(sp2$code)){
#   SP$code= i
#   SP$spectre=t(as.matrix(colMeans(sp2$spectre[which(sp2$code==i),])))
#   sp3=cbind(sp3,SP)
#   which(sp2$code==i)
# }
#
#
# plotsp(colMeans(sp2$spectre[which(sp2$code==1),]))
# plotsp(sp2$spectre[which(sp2$code==1),])


trad=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Correspondance_code_conservatoire_gamay.csv")

donnees1=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Conservatoire_2019_1.csv")

donnees2=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Conservatoire_2019_2.csv")
#donnees2=donnees2[complete.cases(donnees2),]


jonction1=left_join(trad, donnees1)
jonction1=jonction1[complete.cases(jonction1),]
jonction2=left_join(trad, donnees2)
jonction2=jonction2[complete.cases(jonction2),]

#table=left_join(sp2, jonction1)
table=left_join(sp2,jonction2)
table=table[complete.cases(table),]

# intersect(donnees1$clone, trad$clone)
# setdiff(donnees1$clone, trad$clone)
# setdiff(trad$clone, donnees1$clone)
#
# setdiff(donnees2$clone, trad$clone)
# setdiff(trad$clone, donnees2$clone)

# print(table)
# brb="C:/Users/avitvale/Documents/Test/"
# save(table, file=paste(brb,"table",sep=""))
# print(length(table))
# # write.table(globalmatrix, file=paste(brb,"globalmatrix.csv",sep=""),sep=";", quote=FALSE)
# ### END ###

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
sp_pre=adj_asd(table$spectre,c(552,1352)) #Retrouvé empiriquement, ne correspondait pas à ce qui etait indiqué precedemment.

## SNV
sp_pre=t(scale(t(sp_pre)))

## Derivation Savitsky Golay
sp=savitzkyGolay(sp_pre, m = m, p = p, w = n)


sp=sp[,-(1330:1390)] #Coupure du spectre autour des résidus de sauts de detecteurs. Ne semble pas encore effacer completement les sauts.
sp=sp[,-(500:560)]


table$spectre=sp


a=floor(length(unique(table$code))/3)
segm=list()
for (i in 1:10){
  l1=sample(3*a, a)
  l2=sample((1:(3*a))[-l1],a)
  l3=sample((1:(3*a))[-c(l1,l2)],a)
  #fitcv
  #n°rep, n0 clone      16, yiyi

  L1=vector()
  for (i in l1) {
    L1=c(L1,which(table$code==i))
      }
  L2=vector()
  for (i in l2) {
    L2=c(L2,which(table$code==i))
  }
  L3=vector()
  for (i in l3) {
    L3=c(L3,which(table$code==i))
  }

  segm_1=list(list(L1,L2,L3))
  segm=c(segm,segm_1)
}

#length(unique(table$code)) #192. 192/6=32.
names(table)
fm=fitcv(table$spectre, table$intensite_botrytis, plsr, segm, print=T, ncomp=20)

z <- mse(fm, ~ ncomp)

plotmse(z, nam = "rmsep")
plotmse(z, nam = "r2")

fm20=lapply(fm,function (x) {x[x$ncomp==1,]})
plot(fm20$y$x1,fm20$fit$x1)
hist(table$intensite_botrytis)



#Le RPD c'est l'err / écart-type des données de départ.
plot(z$rmsep/sd(table$pH))



cor(table$pH,table$pds_100baies)
cor(table$degre,table$pds_100baies)
cor(table$acidite_totale,table$pds_100baies)
cor(table$degre,table$acidite_totale)
cor(table$degre,table$pH)
cor(table$pH,table$acidite_totale)


hist(table$degre)

unique(table$code[which(table$degre>11,9)])

#spatial, moyenne, corrélation.

##Permet de voir si, avec des trucs aleatoires, ca ferait des prédictions "illusoires"
# L=c("C 015", "C 169", "C 685", "G 222", "G 509", "G 787", "S 471", "S 525", "S 747", "S 877")
# alea=L[sample(1:10,length(sp[,1]),replace = T) ]
#
# rownames(sp)=paste(rownames(sp),alea)



idval=which(substr(rownames(sp),1,4)=="2019" & substr(rownames(sp),18,18)=="G" & substr(rownames(sp),9,9)!="A" & substr(rownames(sp),5,8)!="08171" )#& substr(rownames(sp),1,8)=="20170711" )#& substr(rownames(sp),5,8)!="0702")


spval=sp[idval,]
spcal=sp[-idval,]

### Ici, choisir l'un ou l'autre selon si on veut faire l'analyse sur cepages ou sur clones.
classval=sp$y1[idval]        #Sur les cépages
classcal=sp$y1[-idval]
# classval=sp$y2[idval]         #Sur les clones
# classcal=sp$y2[-idval]


#Consitution d'un jeu de calibration à partir de ce qui n'est pas dans le jeu de validation (on empeche la meme donnée d'etre utilisée dans les deux).
idcal=                     which(((substr(rownames(sp),18,18)=="G"    | substr(rownames(sp),18,18)=="G"
) & substr(rownames(sp),9,9)!="A"    & substr(rownames(sp),1,4)!="2019" & substr(rownames(sp),1,8) %in% dates )    |    substr(rownames(sp),1,9)=="20180816G1" )
classcal=      classcal[which(((substr(rownames(spcal),18,18)=="G" | substr(rownames(spcal),18,18)=="G"
) & substr(rownames(spcal),9,9)!="A" & substr(rownames(spcal),1,4)!="2019" & substr(rownames(spcal),1,8) %in% dates ) | substr(rownames(spcal),1,9)=="20180816G1" )]
spcal=            spcal[which(((substr(rownames(spcal),18,18)=="G" | substr(rownames(spcal),18,18)=="G"
) & substr(rownames(spcal),9,9)!="A" & substr(rownames(spcal),1,4)!="2019" & substr(rownames(spcal),1,8) %in% dates ) | substr(rownames(spcal),1,9)=="20180816G1" ),]



predmF=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))
distances=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))

rplsda=caret::plsda(spcal$x, classcal,ncomp=ncmax)
sccal=rplsda$scores
spval_c=scale(spval$x,center=rplsda$Xmeans,scale = F)
scval=spval_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas


# plotsp(t(rplsda$projection)[1,])
plotsp(t(rplsda$coefficients[,1,])[1,])


for (ii in 2:ncmax) {
  predmF[,ii]=SIGNE_maha0(sccal[,1:ii], classcal, scval[,1:ii])$class
  distances[,ii]=SIGNE_maha0(sccal[,1:ii], classcal, scval[,1:ii])$dist
}

#plotsp(t(rplsda$projection)[2,], col="blue")
# plotsp(spcal$x[c(430,440,450,460,470,480,490),], col="blue")
# str(spval$x)
# rownames(spcal)[88]
#1:88 89:214 215:429 430:645 spcal


#G 2017
#B 1 6 7 16  M 20 22

#G 2019
#B 2 7 10 M 8 28

#plotsp(t(scval)[2,], col="blue")

# for (ii in 2:ncmax) {
#   predmF[,ii]=SIGNE_maha0(sccal[,-ii], classcalclo, scval[,-ii])$class
#   distances[,ii]=SIGNE_maha0(sccal[,-ii], classcalclo, scval[,-ii])$dist
#  }

# predmF[,2]=SIGNE_maha0(cbind(sccal[,1],sccal[,2],sccal[,3],sccal[,4],sccal[,5],sccal[,6],sccal[,7],sccal[,8],sccal[,9],sccal[,10]), classcal, cbind(scval[,1],scval[,2],scval[,3],scval[,4],scval[,5],scval[,6],scval[,7],scval[,8],scval[,9],scval[,10]))$class
# distances[,2]=SIGNE_maha0(cbind(sccal[,1],sccal[,2],sccal[,3],sccal[,4],sccal[,5],sccal[,6],sccal[,7],sccal[,8],sccal[,9],sccal[,10]), classcal, cbind(scval[,1],scval[,2],scval[,3],scval[,4],scval[,5],scval[,6],scval[,7],scval[,8],scval[,9],scval[,10]))$dist
#
# predmF[,2]=SIGNE_maha0(cbind(sccal[,21],sccal[,22],sccal[,20]), classcal, cbind(scval[,21],scval[,22],scval[,20]))$class
# distances[,2]=SIGNE_maha0(cbind(sccal[,21],sccal[,22],sccal[,20]), classcal, cbind(scval[,21],scval[,22],scval[,20]))$dist

#
# classval=relevel(classval, "G 787")
# classval=relevel(classval, "G 509")
# classval=relevel(classval, "G 222")
#
# classval=relevel(classval, "S 877")
# classval=relevel(classval, "S 747")
# classval=relevel(classval, "S 525")
# classval=relevel(classval, "S 471")

tsm=lapply(as.list(predmF), classval, FUN = table)
diagsm=lapply(tsm, FUN = diag)
perokm =100*unlist(lapply(diagsm, FUN = sum))/length(idval)

plot(perokm, xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
perokm
VL=23
tsm[VL]


#analyse <- manova(scval ~ substr(rownames(scval),9,9) * substr(rownames(scval),5,8) * substr(rownames(scval),11,13))
#analyse


#Pour l'étude des distances on se place, un peu arbitrairement, à 8VL.
#132B4C ou #11274A
#4CB4FF ou #4FB2FC
Ldist=as.data.frame(matrix(nrow = length(scval[,1]), ncol = 8))
#Ldist2=as.data.frame(cbind(Premier=matrix(nrow = length(scval[,1]), ncol = 4),Second=matrix(nrow = length(scval[,1]), ncol = 4),Troisieme=as.data.frame(matrix(nrow = length(scval[,1]), ncol = 8))))

l1=c(1,2,3,7,8,9,13,14,15)
l2=c(4,5,6,10,11,12,16,17,18)
###REREGARDER JEUNE/VIEUX, ET AJOUTER SELON LE POINT SUR LA FEUILLE (ET SELON LA SOUCHE ?)


numero=matrix(nrow = length(scval[,1]), ncol = 1)
clone=matrix(nrow = length(scval[,1]), ncol = 1)
date=matrix(nrow = length(scval[,1]), ncol = 1)
age=matrix(nrow = length(scval[,1]), ncol = 1)
souche2=matrix(nrow = length(scval[,1]), ncol = 1)
dista=matrix(nrow = length(scval[,1]), ncol = VL)
distamin=matrix(nrow = length(scval[,1]), ncol = VL)
diffdista=matrix(nrow = length(scval[,1]), ncol = VL)
predclone=matrix(nrow = length(scval[,1]), ncol = VL)
succes=matrix(nrow = length(scval[,1]), ncol = VL)
#Troisieme=as.data.frame(matrix(nrow = length(scval[,1]), ncol = 8))


numero= 1:length(scval[,1])
clone=   substr(rownames(scval),9,13)
cepage= substr(rownames(scval),9,9)
date= substr(rownames(scval),5,8)
souche2 = substr(rownames(scval),20,21)
endroit2 = substr(rownames(scval),23,23)


for (i in 1:length(scval[,1])){
  n=which(c("C",   "G",   "S")==cepage[i])  #"C 015",   "C 169",   "C 685" #"G 222",   "G 509",   "G 787" #"S 877", "S 747", "S 525", "S 471"
  for (j in 2:VL){
    dista[i,j]=distances[j][i,][n]             #
    distamin[i,j]=min(distances[j][i,])
    diffdista[i,j]=sort(distances[j][i,])[2]/sort(distances[j][i,])[1]
    predclone[i,j]=as.character(predmF[j][i,])     #
    succes[i,j]="Mal classé"                     #
    if (predclone[i,j]==cepage[i]){
      succes[i,j]="Bien classé"
    }
  }

  age[i]="J"
  if (as.numeric(substr(rownames(scval)[i],15,16)) %in% l2){
    age[i]="V"
  }
}

dista2=data.frame(dista=I(dista))
diffdista2=data.frame(diffdista=I(diffdista))
distamin2=data.frame(distamin=I(distamin))
predclone2=data.frame(predclone=I(predclone))
succes2=data.frame(succes=I(succes))


Ldist2=cbind(numero,   dista2, distamin2, diffdista2, clone, predclone2, date, age, succes2, souche2, endroit2)





#colour=paste(clone,succes[,VL])
aff2 <- ggplot(Ldist2, aes(x=numero, y=diffdista[,VL],colour=paste(predclone[,VL],succes[,VL]),date=date,clone=clone,predit=predclone[,VL])) + #colour=paste(predclone[,VL],succes[,VL]) #colour=paste(souche2,age)
  geom_point(size=2, alpha=1) +
  scale_color_manual(values = colo)
ggplotly(aff2)
