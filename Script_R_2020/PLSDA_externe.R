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

#Laisser le premier pour un affichage clone, laisser les trois pour un affichage cepage.
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




# L=unique(substr(rownames(sp),1,8))
# L2=sort(L)
# sp2=sp
#
# for (i in 1:length(L2)){
#   N=which(substr(rownames(sp),1,8)==L2[i])
#   sp2=rbind(sp2,sp[N,])
# }
#
# sp3=sp2[(length(sp[,1])+1):length(sp2[,1]),]
# rownames(sp3)=substr(rownames(sp3),1,18)
#
# sp=sp3
#
#
# L=unique(substr(rownames(sp),9,9))
# L2=sort(L)
# sp2=sp
#
# for (i in 1:length(L2)){
#   N=which(substr(rownames(sp),9,9)==L2[i])
#   sp2=rbind(sp2,sp[N,])
# }
#
# sp3=sp2[(length(sp[,1])+1):length(sp2[,1]),]
# rownames(sp3)=substr(rownames(sp3),1,18)
#
# sp=sp3



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


gamay=c("#4D0080", "#9684F0", "#00037D", "#2C77FC", "#006172", "#0EBAB6")
cabernet=c("#A41300", "#F0373C", "#F34307", "#F7830F", "#6E0900", "#D50B47")
syrah=c("#17371B", "#02894F", "#237709", "#34C625", "#6B8D1F", "#9DFD37", "#B77829","#F5C73B")    #"#4D3D25", "#8C5B1B")

#  C 015   C 169   C 685   G 222   G 509   G 787   S 471   S 525   S 747   S 877

# 20170524G 20170529G 20170606G 20170612G 	20170619G 20170626G 20170703G 20170710G 20170711g 20170717G 20170724G 20170731G
# 20180619G 20180627G 20180704G 20180709G 20180710g 20180724B 20180731A 20180810B 20180816G 20180817g 20180823A
# 20190613G 20190617G 20190624G 20190628g 20190702G 20190703A 20190710G 20190718G 20190723A 20190726G 20190730G 20190822B

#Sélection des dates disponibilisées ou non pour la calibration. Une date presente dans cette liste n'apparaitra par forcement en calibration, en fonction des autrs criteres.
dates= c(
  "20170524", "20170529",
  "20170606", "20170612",	"20170619", "20170626", "20170703", "20170710",
  "20170717",
  "20170724",
  "20170731",
  "20180619",
  "20180627",
  "20180704",
  "20180709",
  "20180816",
  "20190613", "20190617", "20190624", "20190702",
  "20190710",
  "20190718",
  "20190726",
  "20170711", "20180710", "20180817", "20190628")





#unique(substr(rownames(sp[which(substr(rownames(sp),18,18)=="g"),]),1,4))
# length(which(substr(rownames(sp[which(substr(rownames(sp),18,18)=="G"),]),9,13)=="S 877"))




#Consitution du jeu de validation. Les individus en faisant partie sont selectionné à partir des infos contenues dans leurs noms.
#Choix de l'annee, de la parcelle, du cepage (pour prediction clones), d'une potentielle date d'enrichissement.
#Caracteres 1 à 4 : annee. 5 à 8 : date. 9 : cepage. 9 à 13 : clone. 15 à 16 : numero mesure (de 1 à 18, normalement). 18 : parcelle. 20 à 21 : souche (prudence). 23 : emplacement sur la feuille.
#Si on veut faire sur cepage, prendre les 3 cepages differents, si on veut faire sur clone, prendre seulement l'un des cepages (via caractère n°9).
idval=which(substr(rownames(sp),1,4)=="2018" & substr(rownames(sp),18,18)=="G" & substr(rownames(sp),9,9)=="G" & substr(rownames(sp),5,8)!="08171" )#& substr(rownames(sp),1,8)=="20170711" )#& substr(rownames(sp),5,8)!="0702")


spval=sp[idval,]
spcal=sp[-idval,]

### Ici, choisir l'un ou l'autre selon si on veut faire l'analyse sur cepages ou sur clones.
# classval=sp$y1[idval]        #Sur les cépages
# classcal=sp$y1[-idval]
classval=sp$y2[idval]         #Sur les clones
classcal=sp$y2[-idval]


#Consitution d'un jeu de calibration à partir de ce qui n'est pas dans le jeu de validation (on empeche la meme donnée d'etre utilisée dans les deux).
idcal=                     which(((substr(rownames(sp),18,18)=="G"    | substr(rownames(sp),18,18)=="G"
) & substr(rownames(sp),9,9)=="G"    & substr(rownames(sp),1,4)!="2018" & substr(rownames(sp),1,8) %in% dates )    |    substr(rownames(sp),1,9)=="20180816G1" )
classcal=      classcal[which(((substr(rownames(spcal),18,18)=="G" | substr(rownames(spcal),18,18)=="G"
) & substr(rownames(spcal),9,9)=="G" & substr(rownames(spcal),1,4)!="2018" & substr(rownames(spcal),1,8) %in% dates ) | substr(rownames(spcal),1,9)=="20180816G1" )]
spcal=            spcal[which(((substr(rownames(spcal),18,18)=="G" | substr(rownames(spcal),18,18)=="G"
) & substr(rownames(spcal),9,9)=="G" & substr(rownames(spcal),1,4)!="2018" & substr(rownames(spcal),1,8) %in% dates ) | substr(rownames(spcal),1,9)=="20180816G1" ),]



predmF=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))
distances=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))

rplsda=caret::plsda(spcal$x, classcal,ncomp=ncmax)
sccal=rplsda$scores
spval_c=scale(spval$x,center=rplsda$Xmeans,scale = F)
scval=spval_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas


# plotsp(t(rplsda$projection)[1,])
plotsp(t(rplsda$coefficients[,1,])[1,]) #cepage
plotsp(t(rplsda$coefficients[,4,])[1,]) #clone Gamay


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
classval=relevel(classval, "G 787")
classval=relevel(classval, "G 509")
classval=relevel(classval, "G 222")
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
VL=3
tsm[VL]

# ############################################################################
## On veut ca. Les i, c'est chaque point
#

classvalT=as.data.frame(matrix(nrow=length(classval), ncol=ncmax))
Passe=logical(length=length(Test[,2][,1]))
PasseT=matrix(F, nrow=length(Test[,2][,1]), ncol= ncmax)
seuil=1.5
predmFT=predmF

for (j in 2:ncmax) {
  for (i in 1:length(distances[,2][,1])) {
    if (sort(distances[j][i,])[2]/sort(distances[j][i,])[1] > seuil) {
      PasseT[i,j]=T
    }
  }
  predmFT[j]=c(predmF[j][PasseT[,j],], rep(NA, nrow(predmFT)-length(predmF[j][PasseT[,j],])))
  classvalT[j]=c(classval[PasseT[,j]], rep(NA, length(classval)-length(classval[PasseT[,j]])))
}

# df <- data.frame(var=1:10)  #notice I created a data.frame vs. the vector you called
# new.col <- c(1:5)
#
# #METHOD 1
# df$new.col <- c(new.col, rep(NA, nrow(df)-length(new.col)))  #keep as integer

# mat<-matrix(NULL,nrow=1,ncol=60)
# for(i in 1:60){mat<-list(mat,matrix(FALSE,nrow=4,ncol=4))}


taille=vector()
tsmT=tsm
for (j in 2:ncmax){
  tsmT[[j]]=table(predmFT[,j][complete.cases(predmFT[,j])], classvalT[,j][complete.cases(classvalT[,j])])
  taille[j]=length(classvalT[,j][complete.cases(classvalT[,j])])/length(classvalT[,j])
}


classvalT[,j][complete.cases(classvalT[,j])]


diagsmT=lapply(tsmT, FUN = diag)

perokmT=perokm
for (j in 2:ncmax){
  perokmT[j]=100*sum(diagsmT[[j]])/length(classvalT[,j][complete.cases(classvalT[,j])])
}

plot(perokmT, xlab= "Nombre de VL", ylab = "Pourcentage de biens class?s",pch=19, cex=1.5)
perokmT
VL=23
VL=6
tsmT[VL]




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
#  n=which(c("C",   "G",   "S")==cepage[i])  #"C 015",   "C 169",   "C 685" #"G 222",   "G 509",   "G 787" #"S 877", "S 747", "S 525", "S 471"
  n=which(c("G 222",   "G 509",   "G 787")==clone[i])
  for (j in 2:VL){
    dista[i,j]=distances[j][i,][n]             #
    distamin[i,j]=min(distances[j][i,])
    diffdista[i,j]=sort(distances[j][i,])[2]/sort(distances[j][i,])[1]
    predclone[i,j]=as.character(predmF[j][i,])     #
    succes[i,j]="Mal classé"                     #
#    if (predclone[i,j]==cepage[i]){
    if (predclone[i,j]==clone[i]){
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
aff2 <- ggplot(Ldist2, aes(x=numero, y=dista[,VL],colour=paste(predclone[,VL],succes[,VL]),date=date,clone=clone,predit=predclone[,VL])) + #colour=paste(predclone[,VL],succes[,VL]) #colour=paste(souche2,age)
  geom_point(size=2, alpha=1) +
#  scale_color_manual(values = colo)
  scale_color_manual(values = gamay) + geom_hline(yintercept=1.5, linetype="dashed", color = "red")
ggplotly(aff2)
