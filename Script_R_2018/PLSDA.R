library(MASS)
library(mixOmics)
library(FactoMineR)
library(signal)
library(plyr)
library(caret)
library(dplyr)

rm(list = ls())

user="ME"
if (user=="NL") {source_dir="C:/Users/Noémie/Desktop/SFE/Script_R/"} else {source_dir="Script_R_2018"}

source(file.path(source_dir,'adj_asd.R'))
source(file.path(source_dir,'SIGNE_load.R'))
source(file.path(source_dir,'SIGNE_maha.R'))
source(file.path(source_dir,'SIGNE_maha0.R'))


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
# set.seed(1)

if (user=="NL") {
  brb="C:/Users/No?mie/Desktop/SFE/Resultats/PTN1/PLSDA/P/globalmatrix"} else {
  brb="/home/ecarnot/Documents/INRA/Projets/SIGNE/2018/GrauduRoi/PTN/globalmatrix2"}

load(file=brb)

# Modifier la date du 20180710, pour qu'elle correponde à 20180709T et N
rownames(globalmatrix)=gsub("20180710","20180709",rownames(globalmatrix))

#globalmatrix[,500]>0.6
globalmatrix=globalmatrix[globalmatrix[,500]>0.6,]
globalmatrix=globalmatrix[globalmatrix[,1]<0.2,]
globalmatrix=globalmatrix[globalmatrix[,2000]<0.25,]

# Data Filter
# Select dates
# dates=list(
  #  "20180731N"
  # ,"20180731P"
  # ,"20180731T"
  # ,"20180724N"
  # ,"20180724P"
  # ,"20180724T"
   # "20180619P"
   # ,"20180627P"
  # ,"20180704P"
   # ,"20180710P"
  # ,"20180710E"
  # ,"20180619T"
  # ,"20180619N"
   #,"20180627T"
   #,"20180627N"
   #,"20180709T"
   #,"20180709N"

# )
# iok=substr(rownames(globalmatrix),1,9) %in% dates
# globalmatrix=globalmatrix[iok,]

titre=rownames(globalmatrix)
## decommenter pour retirer les jeunes feuilles:
#w=as.numeric(substr(titre,5,6))==4 | as.numeric(substr(titre,5,6))==5 | as.numeric(substr(titre,5,6))==6 | as.numeric(substr(titre,5,6))==10 | as.numeric(substr(titre,5,6))==11 | as.numeric(substr(titre,5,6))==12 | as.numeric(substr(titre,5,6))==16 | as.numeric(substr(titre,5,6))==17 | as.numeric(substr(titre,5,6))==18
# ## decommenter pour retirer les vieilles feuilles:
# # w=as.numeric(substr(titre,5,6))==1 | as.numeric(substr(titre,5,6))==2 | as.numeric(substr(titre,5,6))==3 | as.numeric(substr(titre,5,6))==7 | as.numeric(substr(titre,5,6))==8 | as.numeric(substr(titre,5,6))==9 | as.numeric(substr(titre,5,6))==13 | as.numeric(substr(titre,5,6))==14 | as.numeric(substr(titre,5,6))==15
# ## retire les lignes correspondantes (a mettre en commentaire si pas de selection de feuilles ou cepages)
# globalmatrix=globalmatrix[(w==TRUE),]
titre=rownames(globalmatrix)
## decommenter pour retirer la syrah:
# z=as.numeric(substr(titre,1,3))==471 | as.numeric(substr(titre,1,3))==525 | as.numeric(substr(titre,1,3))==747 | as.numeric(substr(titre,1,3))==877
## decommenter pour retirer le cabernet sauvignon:
# z=as.numeric(substr(titre,1,3))==015 | as.numeric(substr(titre,1,3))==169 | as.numeric(substr(titre,1,3))==685
## decommenter pour retirer le gamay:
# z=as.numeric(substr(titre,1,3))==787 | as.numeric(substr(titre,1,3))==509 | as.numeric(substr(titre,1,3))==222
## retire les lignes correspondantes (a mettre en commentaire si pas de selection de feuilles ou cepages)
# globalmatrix=globalmatrix

## FIXATION DES PARAMETRES UTILISES:
# nombre de repetitions de la boucle de FDA:
repet= 4
#Parametres du Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
p=2
n=11
m=1
# nombre de DV max autorisees
ncmax=15
#Taille de l'echantillon de validation (1/v):
#v=3
# Nombre de groupes de CV
#k=6


# On ne retient que ceux qui sont en 3 exemplaires
titre=rownames(globalmatrix)
titrer=paste(substr(titre,1,8),substr(titre,10,16))
ok3=which(table(titrer)==3)
iok=titrer %in% names(ok3)
globalmatrix=globalmatrix[iok,]
titre=rownames(globalmatrix)


## LDA ##
sp=globalmatrix

# creation de la matrice de classes
class=as.factor(substr(rownames(sp),11,13))
# variable qui mesure le nombre de classes
c=length(levels(class))

## Pretraitements
# # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651??me l.o.) et 1800 (1451))
sp=adj_asd(sp,c(651,1451))
# # Reduction des variables (extremites bruitees)
sp=sp[,seq(+1,ncol(sp)-30,1)]
# # SNV
sp=t(scale(t(sp)))
# # Derivation Savitsky Golay
sp=t(apply(sp,1,sgolayfilt,p=p,n=n,m=m))

#brb="C:\\Users\\No?mie\\Desktop\\SFE\\Resultats\\PTN1\\PLSDA\\P\\"
#write.table(sp, file=paste(brb,"sp.csv",sep=""),sep=";", quote=FALSE)

###PLSDA###
#set.seed(1) # fixe le tirage aleatoire

# sp1 = sp[1:179, ]
# sp2 = sp[180:358, ]
# sp3 = sp[359:536, ]
# ech1 <- sample(1 : nrow(sp1),100) # echantillonne aleatoirement 100 spectres sur sp1
# ech2 <- sample(1 : nrow(sp2),100) # echantillonne aleatoirement 100 spectres sur sp2
# ech3 <- sample(1 : nrow(sp3),100) # echantillonne aleatoirement 100 spectres sur sp3
# #echantillon <- sample(seq(along = dates) , 100) # echantillonne 100 spectres par dates (marche pas!!!!!)
#
# echantillon <- rbind(ech1, ech2, ech3) #rassemble les 3 sous echantillons en 1 seul
#
# test <-



# selction des scan pris sur souche
iok=which(substr(rownames(globalmatrix),9,9)=="P")
ivalP=sort(sample(iok,length(iok)/3))
icalP=setdiff(iok,ivalP)

# On remplace P par T et N, et on sleectionne les spectre correspondants
noms_valT=gsub("P","T",titre[ivalP])
ivalT=which(titre %in% noms_valT)
noms_valN=gsub("P","N",titre[ivalP])
ivalN=which(titre %in% noms_valN)

ival=c(ivalP,ivalT,ivalN)

# On selectionne les spectres ayant ces numéros dans le jeu de validation, les autres vont dans le jeu de calibration
sp_val=sp[ival,]
sp_cal=sp[icalP,]
class_val=class[ival]
class_cal=class[icalP]

# Calibration
rplsda=caret::plsda(sp_cal, class_cal,ncomp=ncmax)
sc_cal=rplsda$scores

# Validation
sp_val_c=scale(sp_val,center=rplsda$Xmeans,scale = F)
sc_val=sp_val_c%*%rplsda$projection
res_val=SIGNE_maha0(sc_cal[,1:10], class_cal, sc_val[,1:10])$class

nval=length(res_val)
valP=res_val[1:(nval/3)]
valT=res_val[(nval/3+1):(nval*2/3)]
valN=res_val[(nval*2/3+1):nval]
