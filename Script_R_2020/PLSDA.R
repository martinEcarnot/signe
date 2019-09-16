library(MASS)
library(mixOmics)
library(FactoMineR)
library(signal)
library(plyr)
library(caret)
library(dplyr)
library(nirs)

rm(list = ls())

source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha0.R')


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
# set.seed(1)

#"C:/Users/Noémie/Desktop/SFE/Robustesse/96vsX/globalmatrix"
brb="C:/Users/avitvale/Documents/matrice toutes dates/Globalmatrix/globalmatrix"
load(file=brb)

## Filtrage des spectres aberrants
globalmatrix=globalmatrix[globalmatrix[,500]>0.6,]
globalmatrix=globalmatrix[globalmatrix[,1]<0.2,]
globalmatrix=globalmatrix[globalmatrix[,2000]<0.25,]

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
# ,"20180724T"
# ,"20180724N"
"20180816X"
 ,"20180816E"
 ,"20180710X"
 ,"20180710E"
# "20180816A"
# ,"20180816E"
# ,"20180731A"
# ,"20180731E"
  
  
)
iok=substr(rownames(globalmatrix),1,9) %in% dates
globalmatrix=globalmatrix[iok,]

### FIXATION DES PARAMETRES UTILISES:
## Nombre de repetitions de la boucle de FDA:
#repet= 4
##Parametres de Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
p=2
n=11
m=1
## Nombre de VL max autorisees
ncmax=10
##Taille de l'echantillon de validation (1/v):
#v=3
## Nombre de groupes de CV
#k=6

## PLSDA ##
sp=globalmatrix

### Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
sp=adj_asd(sp,c(602,1402))
## Reduction des variables (extremites bruitees)
sp=sp[,seq(51,ncol(sp)-30,1)]
## SNV
sp=t(scale(t(sp)))
## Derivation Savitsky Golay
sp=t(apply(sp,1,sgolayfilt,p=p,n=n,m=m))

###Ajout du cepage au nom de la ligne. !!Attention a bien commenter en fonction des cepages presents dans la BDD!!
##Filtre en fonction du cepage
titre=rownames(sp)

ya=(substr(titre,11,13))== "015" |  (substr(titre,11,13))== "169" |  (substr(titre,11,13))== "685"
yb=(substr(titre,11,13))== "222" |  (substr(titre,11,13))== "509" |  (substr(titre,11,13))== "787"
yc=(substr(titre,11,13))== "471" |  (substr(titre,11,13))== "525" |  (substr(titre,11,13))== "747" |  (substr(titre,11,13))== "877"

##Separation en 3 matrices
spa=sp[ya==TRUE,]
spb=sp[yb==TRUE,]
spc=sp[yc==TRUE,]

#Rajoute le nom du cepage au nom de la ligne
rownames(spa)=paste(rownames(spa),"C",sep = "-")
rownames(spb)=paste(rownames(spb),"G",sep = "-")
rownames(spc)=paste(rownames(spc),"S",sep = "-")

##Recombine les 3 matrices pour reformer sp
sp=rbind(spa,spb,spc)

##Creation de la matrice de classes
class=as.factor(substr(rownames(sp),18,18))
##Variable qui mesure le nombre de classes
c=length(levels(class))

###PLSDA###
#set.seed(1) # fixe le tirage aleatoire

##Selection des scan pris sur souche
iok=which(substr(rownames(sp),9,9)=="E")
ivalP=sort(sample(iok,length(iok)/3))
icalP=iok

##On remplace P par T et N, et on selectionne les spectres correspondants
noms_valT=gsub("E","X",rownames(sp)[icalP])
ivalT=which(rownames(sp) %in% noms_valT)
# noms_valN=gsub("P","N",rownames(sp)[icalP])
# ivalN=which(rownames(sp) %in% noms_valN)

ival= ivalT #c(ivalT,ivalN)
# nval=length(ival)
# res_val=as.data.frame(matrix(nrow = nval , ncol = ncmax))
# rownames(res_val)= rownames(sp[ival,])

##On selectionne les spectres ayant ces numéros dans le jeu de validation, les autres vont dans le jeu de calibration
sp_val=sp[ival,]
sp_cal=sp[icalP,]
class_val=class[ival]
class_cal=class[icalP]

###PLSDA on Maha scores
## Calibration
rplsda=caret::plsda(sp_cal, class_cal,ncomp=10)
sc_cal=rplsda$scores

## Validation
sp_val_c=scale(sp_val,center=rplsda$Xmeans,scale = F)
sc_val=sp_val_c%*%rplsda$projection
res_val=SIGNE_maha0(sc_cal[,1:10], class_cal, sc_val[,1:10])$class

# ##KNN
# #Calibration
# fm = nirs::plsda(sp_cal,class_cal,sp_val,class_val,ncomp = ncmax,typda = "daknn", k = 15) #k = nbr de voisins
# sc_cal=fm$scores

# ##Validation
# sp_val_c=scale(sp_val,center=fm$Xmeans,scale = F)
# sc_val=sp_val_c%*%fm$projection
# res_val=SIGNE_maha0(sc_cal[,1:ncmax], class_cal, sc_val[,1:ncmax])$class

# # nval=length(res_val)
# # valP=res_val[1:(nval/3)]
# # valT=res_val[(nval/3+1):(nval*2/3)]
# # valN=res_val[(nval*2/3+1):nval]
cepage=table (res_val,class_val)
print (cepage )

###En fonction des clones

## Creation de la matrice de classes clones
class_clones=as.factor(substr(rownames(sp),11,13))
## Variable qui mesure le nombre de classes
c=length(levels(class_clones))

## Separation de sp_cal en 3 jeux de calibration par cepage
aC= substr(rownames(sp_cal),18,18)=="C"
aG= substr(rownames(sp_cal),18,18)=="G"
aS= substr(rownames(sp_cal),18,18)=="S"

sp_cal_C =sp_cal[(aC==TRUE),]
sp_cal_G =sp_cal[(aG==TRUE),]
sp_cal_S =sp_cal[(aS==TRUE),]

###Identifiants des matrices de calibration
##Cabernet
id_cal_C=which(rownames(sp)  %in%  rownames(sp_cal_C))
class_cal_C=droplevels(class_clones[id_cal_C])

##Gamay
id_cal_G=which(rownames(sp)  %in%  rownames(sp_cal_G))
class_cal_G=droplevels(class_clones[id_cal_G])

##Syrah
id_cal_S=which(rownames(sp)  %in%  rownames(sp_cal_S))
class_cal_S=droplevels(class_clones[id_cal_S])

###PLSDA sur clones sans CV
## Separation de sp_val en 3 jeux de validation par cepage
## Seulement avec les biens classes en cepages
test=cbind(res_val,class_val)
rownames(test)=rownames(sp[ival,])

z1=(substr(rownames(test),18,18))=="C"
z2=(substr(rownames(test), 18, 18))=="G"
z3=(substr(rownames(test),18,18))=="S"

test1 =test[(z1==TRUE),]
test2 =test[(z2==TRUE),]
test3 =test[(z3==TRUE),]

test1a = test1[test1[,1]==1,]
test2a = test2[test2[,1]==2,]
test3a = test3[test3[,1]==3,]

sp_val_C =sp_val[rownames(test1a),]
sp_val_G =sp_val[rownames(test2a),]
sp_val_S =sp_val[rownames(test3a),]

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

##Gamay
id_val_G=which(rownames(sp)  %in%  rownames(sp_val_G))
class_val_G=droplevels(class_clones[id_val_G])


##Syrah
id_val_S=which(rownames(sp)  %in%  rownames(sp_val_S))
class_val_S=droplevels(class_clones[id_val_S])

####PLSDA on Maha scores
### Calibration
## Cabernet Sauvignon
rplsdaC=caret::plsda(sp_cal_C, class_cal_C,ncomp= 8)# Modifier ncmax en fonction des resultats de CV_2
sc_cal_C=rplsdaC$scores

##Gamay
rplsdaG=caret::plsda(sp_cal_G, class_cal_G,ncomp=8)# Modifier ncmax en fonction des resultats de CV_2
sc_cal_G=rplsdaG$scores

##Syrah
rplsdaS=caret::plsda(sp_cal_S, class_cal_S,ncomp=8)# Modifier ncmax en fonction des resultats de CV_2
sc_cal_S=rplsdaS$scores

### Validation
## Cabernet Sauvignon
sp_val_c_C=scale(sp_val_C,center=rplsdaC$Xmeans,scale = F)
sc_val_C=sp_val_c_C%*%rplsdaC$projection
res_val_C=SIGNE_maha0(sc_cal_C[,1:8], class_cal_C, sc_val_C[,1:8])$class

Cabernet=table (res_val_C,class_val_C)
print (Cabernet)

##Gamay
sp_val_c_G=scale(sp_val_G,center=rplsdaG$Xmeans,scale = F)
sc_val_G=sp_val_c_G%*%rplsdaG$projection
res_val_G=SIGNE_maha0(sc_cal_G[,1:8], class_cal_G, sc_val_G[,1:8])$class

Gamay=table (res_val_G,class_val_G)
print (Gamay)

##Syrah
sp_val_c_S=scale(sp_val_S,center=rplsdaS$Xmeans,scale = F)
sc_val_S=sp_val_c_S%*%rplsdaS$projection
res_val_S=SIGNE_maha0(sc_cal_S[,1:8], class_cal_S, sc_val_S[,1:8])$class

Syrah=table (res_val_S,class_val_S)
print (Syrah)

###Calcul des pourcentages de bons classements en fonction du tirage
##Somme des clones biens classes
clones_bc= sum(diag(Cabernet))+sum(diag(Gamay))+sum(diag(Syrah))
#clones_bc= sum(diag(Cabernet))+sum(diag(Syrah))
#print(clones_bc)

#Somme colonne table contingence cepages
sumcolcep=apply(cepage,2,sum)

##Nombre total de clones
total=clones_bc+(sum(cepage)-sum(diag(cepage)))+(sum(Cabernet)-sum(diag(Cabernet)))+(sum(Gamay)-sum(diag(Gamay)))+(sum(Syrah)-sum(diag(Syrah)))
#total=clones_bc+(sum(cepage)-sum(diag(cepage)))+(sum(Cabernet)-sum(diag(Cabernet)))+(sum(Syrah)-sum(diag(Syrah))) #sans gamay

total_C=sum(Cabernet)
total_G=sum(Gamay)
total_S=sum(Syrah)

total_C_cep=sum(diag(Cabernet))+(sum(Cabernet)-sum(diag(Cabernet)))+(sumcolcep[1]-cepage[1,1])
total_G_cep=sum(diag(Gamay))+(sum(Gamay)-sum(diag(Gamay)))+(sumcolcep[2]-cepage[2,2])
total_S_cep=sum(diag(Syrah))+(sum(Syrah)-sum(diag(Syrah)))+(sumcolcep[3]-cepage[3,3])

#print(total)

## Pourcentage total de clones biens classes(prise en compte erreur cepages)
perok=100*(clones_bc/total)

## Pourcentage de clones biens classes
perok_C=100*(sum(diag(Cabernet))/total_C)
perok_G=100*(sum(diag(Gamay))/total_G)
perok_S=100*(sum(diag(Syrah))/total_S)

## Pourcentage de clones biens classes (prise en compte erreur cepages)
perok_C_cep=100*(sum(diag(Cabernet))/total_C_cep)
perok_G_cep=100*(sum(diag(Gamay))/total_G_cep)
perok_S_cep=100*(sum(diag(Syrah))/total_S_cep)

##Pourcentage de cepages biens classes
perok_cepages=100*(sum(diag(cepage))/sum(cepage))

print(perok)
print(perok_cepages)
print(perok_C)
print(perok_G)
print(perok_S)
print(perok_C_cep)
print(perok_G_cep)
print(perok_S_cep)

# ## Clones
# plot(perok_final, type="o",xlab= "Nombre de tirages", ylab = "Pourcentage de clones biens classés",pch=19, cex=2)
# 
# ##Cepages
# plot(perok_final_cepages, type="o",xlab= "Nombre de tirages", ylab = "Pourcentage de cépages biens classés",pch=21, cex=2,bg="blue")


 

#plot
#df=data.frame(rplsda=rplsda$scores[,1:10],class)
#r=ggplot(data = df, aes(rplsda.Comp.1,rplsda.Comp.2, colour=class))
#+ geom_point()
#plot(r)

#plot(sc_cal[1:20,1],sc_val[1:20,1], col= c("red","blue") ,pch= 19)
#text(sc_cal[1:20,1],sc_val[1:20,1], substr(rownames(sp_cal),11,13) ,
#    cex=0.65, pos=1)


#plot(sc_cal[,1],sc_cal[,2], pch= 19)
#plot(sc_cal[,1],sc_cal[,2], col=dat_cal)
#+text(sc_val[,1],sc_val[,2], labels=substr(rownames(sp_val),11,12), col = "red")
#plot(sc_cal[,1],sc_cal[,2], col=as.factor(v))
#couleur <- factor(class)
#levels(couleur)<- c("red", "black", "cyan","green","magenta", "orange4", "yellow","darkgrey","blue","purple3") #colorer en fonction du type de clone
#couleur <- c("red", "black", "cyan","green","magenta", "orange4", "yellow","darkgrey","blue","purple3")
 # plot(sc_cal[,1],sc_cal[,2], pch=19, cex= 2, axes = FALSE, xlim= c(-0.20,0.30),
 #      ylim= c(-0.10,0.10), xlab = "Comp 1", ylab = "Comp 2", col=class_cal)
 # points(sc_val[,1],sc_val[,2],pch=7, cex= 2,col= class_val)
 # axis(1, pos = 0)
 # axis(2, pos = 0)
 # legend(x="bottomright", legend=c("P","N et T"), pch=c(19,7),pt.cex= 2, bty="n")
#+legend (x="bottomright", col = couleur)
#text(sc_val[,1],sc_val[,2], labels=class_val, cex=0.65, pos=1)
#text(sc_cal[,1],sc_cal[,2], labels= class_cal, cex=0.65, pos=1)

#  ##Matrice avec seulement les biens classes par cepage
# 
#  test=cbind(res_val,class_val)
#  rownames(test)=rownames(sp[ival,])
# 
#  z1=(substr(rownames(test),18,18))=="C"
#  z2=(substr(rownames(test), 18, 18))=="G"
#  z3=(substr(rownames(test),18,18))=="S"
# 
#  test1 =test[(z1==TRUE),]
#  test2 =test[(z2==TRUE),]
#  test3 =test[(z3==TRUE),]
# 
#  test1a = test1[test1[,1]==1,]
#  test2a = test2[test2[,1]==2,]
#  test3a = test3[test3[,1]==3,]
# 
#  test_cepages = rbind (test1a, test2a, test3a)
# 
#  sp_clones = sp[rownames(test_cepages),]
# 
#  nomsP1= gsub("T","P",rownames(sp_clones))
#  nomsP2= gsub("N","P",rownames(sp_clones))
# 
#  sp_clones2= sp[nomsP1,]
#  sp_clones3= sp[nomsP2,]
#  sp_clones= rbind (sp_clones,sp_clones2, sp_clones3)
# 
#  nomssp=rownames(sp_clones)
#  doublons= which(duplicated(nomssp))
#  nomssp2=nomssp[-doublons]
#  sp_clones=sp[nomssp2,]
# 
# ## Prediction apres cepages
# 
# iok2=substr(rownames(sp_clones),1,9) %in% dates
# sp_clones=sp_clones[iok2,]
# 
#  # creation de la matrice de classes
#  class=as.factor(substr(rownames(sp_clones),11,13))
#  # variable qui mesure le nombre de classes
#  c=length(levels(class))
# 
#  ###PLSDA###
#  #set.seed(1) # fixe le tirage aleatoire
# 
#  # selection des scan pris sur souche
#  iok2=which(substr(rownames(sp_clones),9,9)=="P")
#  ivalP=sort(sample(iok2,length(iok2)/3))
#  icalP=iok2
# 
#  # On remplace P par T et N, et on selectionne les spectres correspondants
#  noms_valT=gsub("P","T",rownames(sp_clones)[icalP])
#  ivalT=which(rownames(sp_clones) %in% noms_valT)
#  noms_valN=gsub("P","N",rownames(sp_clones)[icalP])
#  ivalN=which(rownames(sp_clones) %in% noms_valN)
# 
#  ival=c(ivalT,ivalN)
# 
#  # On selectionne les spectres ayant ces numéros dans le jeu de validation, les autres vont dans le jeu de calibration
#  sp_val=sp[ival,]
#  sp_cal=sp[icalP,]
#  class_val=class[ival]
#  class_cal=class[icalP]
# 
#  # Calibration
#  rplsda=caret::plsda(sp_cal, class_cal,ncomp=ncmax)
#  sc_cal=rplsda$scores
# 
#  # Validation
#  sp_val_c=scale(sp_val,center=rplsda$Xmeans,scale = F)
#  sc_val=sp_val_c%*%rplsda$projection
#  res_val=SIGNE_maha0(sc_cal[,1:ncmax], class_cal, sc_val[,1:ncmax])$class
# 
#  # nval=length(res_val)
#  # valP=res_val[1:(nval/3)]
#  # valT=res_val[(nval/3+1):(nval*2/3)]
#  # valN=res_val[(nval*2/3+1):nval]
# 
# 
#  print (table(res_val,class_val))
# 
# 
#  #plot
#  #df=data.frame(rplsda=rplsda$scores[,1:10],class)
#  #r=ggplot(data = df, aes(rplsda.Comp.1,rplsda.Comp.2, colour=class))
#  #+ geom_point()
#  #plot(r)
# 
#  #plot(sc_cal[1:20,1],sc_val[1:20,1], col= c("red","blue") ,pch= 19)
#  #text(sc_cal[1:20,1],sc_val[1:20,1], substr(rownames(sp_cal),11,13) ,
#  #    cex=0.65, pos=1)
# 
# 
#  #plot(sc_cal[,1],sc_cal[,2], pch= 19)
#  #plot(sc_cal[,1],sc_cal[,2], col=dat_cal)
#  #+text(sc_val[,1],sc_val[,2], labels=substr(rownames(sp_val),11,12), col = "red")
#  #plot(sc_cal[,1],sc_cal[,2], col=as.factor(v))
#  #couleur <- factor(class)
#  #levels(couleur)<- c("red", "black", "cyan","green","magenta", "orange4", "yellow","darkgrey","blue","purple3") #colorer en fonction du type de clone
#  #couleur <- c("red", "black", "cyan","green","magenta", "orange4", "yellow","darkgrey","blue","purple3")
 plot(sc_cal[,1],sc_cal[,2], pch=19, cex= 2, axes = FALSE, xlim= c(-0.15,0.10),
      ylim= c(-0.10,0.10), xlab = "Comp 1", ylab = "Comp 2", col=class_cal)
 points(sc_val[,1],sc_val[,2],pch=7, cex= 2,col= class_val)
 axis(1, pos = 0)
 axis(2, pos = 0)
 legend(x="bottomright", legend=c("Espiguette","Autres emplacements Espiguette"), pch=c(19,7),pt.cex= 2, bty="n")
 
 plot(sc_cal_C[,1],sc_cal_C[,2], pch=19, cex= 2, axes = FALSE, xlim= c(-0.10,0.15),
      ylim= c(-0.05,0.05), xlab = "Comp 1", ylab = "Comp 2", col=class_cal_C)
 points(sc_val_C[,1],sc_val_C[,2],pch=7, cex= 2,col= class_val_C)
 axis(1, pos = 0)
 axis(2, pos = 0)
 legend(x="bottomright", legend=c("Espiguette","Autres emplacements Espiguette"), pch=c(19,7),pt.cex= 2, bty="n")
 
 plot(sc_cal_G[,1],sc_cal_G[,2], pch=19, cex= 2, axes = FALSE, xlim= c(-0.10,0.10),
      ylim= c(-0.05,0.05), xlab = "Comp 1", ylab = "Comp 2", col=class_cal_G)
 points(sc_val_G[,1],sc_val_G[,2],pch=7, cex= 2,col= class_val_G)
 axis(1, pos = 0)
 axis(2, pos = 0)
 legend(x="bottomright", legend=c("Espiguette","Autres emplacements Espiguette"), pch=c(19,7),pt.cex= 2, bty="n")
 
 plot(sc_cal_S[,1],sc_cal_S[,2], pch=19, cex= 2, axes = FALSE, xlim= c(-0.15,0.10),
      ylim= c(-0.05,0.05), xlab = "Comp 1", ylab = "Comp 2", col=class_cal_S)
 points(sc_val_S[,1],sc_val_S[,2],pch=7, cex= 2,col= class_val_S)
 axis(1, pos = 0)
 axis(2, pos = 0)
 legend(x="bottomright", legend=c("Espiguette","Autres emplacements Espiguette"), pch=c(19,7),pt.cex= 2, bty="n")
 #plot
#  #+legend (x="bottomright", col = couleur)
#  #text(sc_val[,1],sc_val[,2], labels=class_val, cex=0.65, pos=1)
#  #text(sc_cal[,1],sc_cal[,2], labels= class_cal, cex=0.65, pos=1)
#  #c("red", "black", "cyan","green","magenta", "orange4", "yellow","darkgrey","blue","purple3")
