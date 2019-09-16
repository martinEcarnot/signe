library(MASS)
library(mixOmics)
library(FactoMineR)
library(signal)
library(plyr)
library(caret)


rm(list = ls())

source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha0.R')


## Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
#set.seed(1)

brb4="C:/Users/avitvale/Documents/matrice toutes dates/Globalmatrix/globalmatrix"
load(file=brb4)

## Filtrage des spectres aberrants
globalmatrix=globalmatrix[globalmatrix[,500]>0.6,]
globalmatrix=globalmatrix[globalmatrix[,1]<0.2,]
globalmatrix=globalmatrix[globalmatrix[,2000]<0.25,]

# Data Filter
# Select dates
dates=list(
# "20180619T"
#  "20180619P"
#  ,"20180619N"
#  ,"20180627T"
#,"20180627P"
# ,"20180627N"
 "20180704P"
# ,"20180709T"
 ,"20180709P"
# ,"20180709N"
# ,"20180710P"
# ,"20180816T"
#,"20180816P"
# ,"20180816N"
# ,"20180817T"
# ,"20180817P"
# ,"20180817N"
# "20170524P"
# ,"20170529P"
# ,"20170606P"
# ,"20170612P"
#,"20170619P"
# ,"20170626P"
,"20170703P"
 ,"20170710P"
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
  
)
iok=substr(rownames(globalmatrix),1,9) %in% dates
globalmatrix=globalmatrix[iok,]

### FIXATION DES PARAMETRES UTILISES:
## Nombre de repetitions de la boucle de PLSDA:
repet= 50
## Parametres de Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
p=2
n=11
m=1
## Nombre de VL max autorisees
ncmax=10
## Taille de l'echantillon de validation (1/v):
#v=3
## Nombre de groupes de CV (3 ou 6)
k=3

## PLSDA ##
sp=globalmatrix

## Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
sp=adj_asd(sp,c(602,1402))
# # Reduction des variables (extremites bruitees)
sp=sp[,seq(+1,ncol(sp)-30,1)]
## SNV
sp=t(scale(t(sp)))
## Derivation Savitsky Golay
sp=t(apply(sp,1,sgolayfilt,p=p,n=n,m=m))

### Ajout du cepage au nom de la ligne. !!Attention a bien commenter en fonction des cepages presents dans la BDD!!
## Filtre en fonction du cepage
titre=rownames(sp)
ya=(substr(titre,11,13))== "015" |  (substr(titre,11,13))== "169" |  (substr(titre,11,13))== "685"
yb=(substr(titre,11,13))== "222" |  (substr(titre,11,13))== "509" |  (substr(titre,11,13))== "787"
yc=(substr(titre,11,13))== "471" |  (substr(titre,11,13))== "525" |  (substr(titre,11,13))== "747" |  (substr(titre,11,13))== "877"

## Separation en 3 matrices
spa=sp[ya==TRUE,]
spb=sp[yb==TRUE,]
spc=sp[yc==TRUE,]

## Rajoute le nom du cepage au nom de la ligne
rownames(spa)=paste(rownames(spa),"C",sep = "-")
rownames(spb)=paste(rownames(spb),"G",sep = "-")
rownames(spc)=paste(rownames(spc),"S",sep = "-")

## Recombine les 3 matrices pour reformer sp
sp=rbind(spa,spb,spc)

iok2=substr(rownames(sp),1,9) %in% dates
sp=sp[iok2,]

## Creation de la matrice de classes
class=as.factor(substr(rownames(sp),11,13)) #(18,18) cepages, (11,13) clones, (9,9) mode de prelevement
## Variable qui mesure le nombre de classes
c=length(levels(class))

##Echantillonnage par "souche"
## Separation de sp en 30 matrices si on separe les souches
titre=rownames(sp)

## Cabernet Sauvignon -> spa
# x1= (substr(titre,11,16))== "015-01" |  (substr(titre,11,16))== "015-02" |  (substr(titre,11,16))== "015-03" |  (substr(titre,11,16))== "015-04" |  (substr(titre,11,16))== "015-05" |  (substr(titre,11,16))== "015-06"
# x4= (substr(titre,11,16))== "015-07" |  (substr(titre,11,16))== "015-08" |  (substr(titre,11,16))== "015-09" |  (substr(titre,11,16))== "015-10" |  (substr(titre,11,16))== "015-11" |  (substr(titre,11,16))== "015-12"
# x7= (substr(titre,11,16))== "015-13" |  (substr(titre,11,16))== "015-14" |  (substr(titre,11,16))== "015-15" |  (substr(titre,11,16))== "015-16" |  (substr(titre,11,16))== "015-17" |  (substr(titre,11,16))== "015-18"
# x2= (substr(titre,11,16))== "169-01" |  (substr(titre,11,16))== "169-02" |  (substr(titre,11,16))== "169-03" |  (substr(titre,11,16))== "169-04" |  (substr(titre,11,16))== "169-05" |  (substr(titre,11,16))== "169-06"
# x5= (substr(titre,11,16))== "169-07" |  (substr(titre,11,16))== "169-08" |  (substr(titre,11,16))== "169-09" |  (substr(titre,11,16))== "169-10" |  (substr(titre,11,16))== "169-11" |  (substr(titre,11,16))== "169-12"
# x8= (substr(titre,11,16))== "169-13" |  (substr(titre,11,16))== "169-14" |  (substr(titre,11,16))== "169-15" |  (substr(titre,11,16))== "169-16" |  (substr(titre,11,16))== "169-17" |  (substr(titre,11,16))== "169-18"
# x3= (substr(titre,11,16))== "685-01" |  (substr(titre,11,16))== "685-02" |  (substr(titre,11,16))== "685-03" |  (substr(titre,11,16))== "685-04" |  (substr(titre,11,16))== "685-05" |  (substr(titre,11,16))== "685-06"
# x6= (substr(titre,11,16))== "685-07" |  (substr(titre,11,16))== "685-08" |  (substr(titre,11,16))== "685-09" |  (substr(titre,11,16))== "685-10" |  (substr(titre,11,16))== "685-11" |  (substr(titre,11,16))== "685-12"
# x9= (substr(titre,11,16))== "685-13" |  (substr(titre,11,16))== "685-14" |  (substr(titre,11,16))== "685-15" |  (substr(titre,11,16))== "685-16" |  (substr(titre,11,16))== "685-17" |  (substr(titre,11,16))== "685-18"

###Gamay -> spb
# x1= (substr(titre,11,16))== "222-01" |  (substr(titre,11,16))== "222-02" |  (substr(titre,11,16))== "222-03" |  (substr(titre,11,16))== "222-04" |  (substr(titre,11,16))== "222-05" |  (substr(titre,11,16))== "222-06"
# x4= (substr(titre,11,16))== "222-07" |  (substr(titre,11,16))== "222-08" |  (substr(titre,11,16))== "222-09" |  (substr(titre,11,16))== "222-10" |  (substr(titre,11,16))== "222-11" |  (substr(titre,11,16))== "222-12"
# x7= (substr(titre,11,16))== "222-13" |  (substr(titre,11,16))== "222-14" |  (substr(titre,11,16))== "222-15" |  (substr(titre,11,16))== "222-16" |  (substr(titre,11,16))== "222-17" |  (substr(titre,11,16))== "222-18"
# x2= (substr(titre,11,16))== "509-01" |  (substr(titre,11,16))== "509-02" |  (substr(titre,11,16))== "509-03" |  (substr(titre,11,16))== "509-04" |  (substr(titre,11,16))== "509-05" |  (substr(titre,11,16))== "509-06"
# x5= (substr(titre,11,16))== "509-07" |  (substr(titre,11,16))== "509-08" |  (substr(titre,11,16))== "509-09" |  (substr(titre,11,16))== "509-10" |  (substr(titre,11,16))== "509-11" |  (substr(titre,11,16))== "509-12"
# x8= (substr(titre,11,16))== "509-13" |  (substr(titre,11,16))== "509-14" |  (substr(titre,11,16))== "509-15" |  (substr(titre,11,16))== "509-16" |  (substr(titre,11,16))== "509-17" |  (substr(titre,11,16))== "509-18"
# x3= (substr(titre,11,16))== "787-01" |  (substr(titre,11,16))== "787-02" |  (substr(titre,11,16))== "787-03" |  (substr(titre,11,16))== "787-04" |  (substr(titre,11,16))== "787-05" |  (substr(titre,11,16))== "787-06"
# x6= (substr(titre,11,16))== "787-07" |  (substr(titre,11,16))== "787-08" |  (substr(titre,11,16))== "787-09" |  (substr(titre,11,16))== "787-10" |  (substr(titre,11,16))== "787-11" |  (substr(titre,11,16))== "787-12"
# x9= (substr(titre,11,16))== "787-13" |  (substr(titre,11,16))== "787-14" |  (substr(titre,11,16))== "787-15" |  (substr(titre,11,16))== "787-16" |  (substr(titre,11,16))== "787-17" |  (substr(titre,11,16))== "787-18"

# ### Syrah -> spc attention 4 clones
# x1= (substr(titre,11,16))== "471-01" |  (substr(titre,11,16))== "471-02" |  (substr(titre,11,16))== "471-03" |  (substr(titre,11,16))== "471-04" |  (substr(titre,11,16))== "471-05" |  (substr(titre,11,16))== "471-06"
# x5= (substr(titre,11,16))== "471-07" |  (substr(titre,11,16))== "471-08" |  (substr(titre,11,16))== "471-09" |  (substr(titre,11,16))== "471-10" |  (substr(titre,11,16))== "471-11" |  (substr(titre,11,16))== "471-12"
# x9= (substr(titre,11,16))== "471-13" |  (substr(titre,11,16))== "471-14" |  (substr(titre,11,16))== "471-15" |  (substr(titre,11,16))== "471-16" |  (substr(titre,11,16))== "471-17" |  (substr(titre,11,16))== "471-18"
# x2= (substr(titre,11,16))== "525-01" |  (substr(titre,11,16))== "525-02" |  (substr(titre,11,16))== "525-03" |  (substr(titre,11,16))== "525-04" |  (substr(titre,11,16))== "525-05" |  (substr(titre,11,16))== "525-06"
# x6= (substr(titre,11,16))== "525-07" |  (substr(titre,11,16))== "525-08" |  (substr(titre,11,16))== "525-09" |  (substr(titre,11,16))== "525-10" |  (substr(titre,11,16))== "525-11" |  (substr(titre,11,16))== "525-12"
# x10=(substr(titre,11,16))== "525-13" |  (substr(titre,11,16))== "525-14" |  (substr(titre,11,16))== "525-15" |  (substr(titre,11,16))== "525-16" |  (substr(titre,11,16))== "525-17" |  (substr(titre,11,16))== "525-18"
# x3= (substr(titre,11,16))== "747-01" |  (substr(titre,11,16))== "747-02" |  (substr(titre,11,16))== "747-03" |  (substr(titre,11,16))== "747-04" |  (substr(titre,11,16))== "747-05" |  (substr(titre,11,16))== "747-06"
# x7= (substr(titre,11,16))== "747-07" |  (substr(titre,11,16))== "747-08" |  (substr(titre,11,16))== "747-09" |  (substr(titre,11,16))== "747-10" |  (substr(titre,11,16))== "747-11" |  (substr(titre,11,16))== "747-12"
# x11=(substr(titre,11,16))== "747-13" |  (substr(titre,11,16))== "747-14" |  (substr(titre,11,16))== "747-15" |  (substr(titre,11,16))== "747-16" |  (substr(titre,11,16))== "747-17" |  (substr(titre,11,16))== "747-18"
# x4= (substr(titre,11,16))== "877-01" |  (substr(titre,11,16))== "877-02" |  (substr(titre,11,16))== "877-03" |  (substr(titre,11,16))== "877-04" |  (substr(titre,11,16))== "877-05" |  (substr(titre,11,16))== "877-06"
# x8= (substr(titre,11,16))== "877-07" |  (substr(titre,11,16))== "877-08" |  (substr(titre,11,16))== "877-09" |  (substr(titre,11,16))== "877-10" |  (substr(titre,11,16))== "877-11" |  (substr(titre,11,16))== "877-12"
# x12=(substr(titre,11,16))== "877-13" |  (substr(titre,11,16))== "877-14" |  (substr(titre,11,16))== "877-15" |  (substr(titre,11,16))== "877-16" |  (substr(titre,11,16))== "877-17" |  (substr(titre,11,16))== "877-18"

###Tous
x1= (substr(titre,11,16))== "015-01" |  (substr(titre,11,16))== "015-02" |  (substr(titre,11,16))== "015-03" |  (substr(titre,11,16))== "015-04" |  (substr(titre,11,16))== "015-05" |  (substr(titre,11,16))== "015-06"
x11=(substr(titre,11,16))== "015-07" |  (substr(titre,11,16))== "015-08" |  (substr(titre,11,16))== "015-09" |  (substr(titre,11,16))== "015-10" |  (substr(titre,11,16))== "015-11" |  (substr(titre,11,16))== "015-12"
x21=(substr(titre,11,16))== "015-13" |  (substr(titre,11,16))== "015-14" |  (substr(titre,11,16))== "015-15" |  (substr(titre,11,16))== "015-16" |  (substr(titre,11,16))== "015-17" |  (substr(titre,11,16))== "015-18"
x2= (substr(titre,11,16))== "169-01" |  (substr(titre,11,16))== "169-02" |  (substr(titre,11,16))== "169-03" |  (substr(titre,11,16))== "169-04" |  (substr(titre,11,16))== "169-05" |  (substr(titre,11,16))== "169-06"
x12=(substr(titre,11,16))== "169-07" |  (substr(titre,11,16))== "169-08" |  (substr(titre,11,16))== "169-09" |  (substr(titre,11,16))== "169-10" |  (substr(titre,11,16))== "169-11" |  (substr(titre,11,16))== "169-12"
x22=(substr(titre,11,16))== "169-13" |  (substr(titre,11,16))== "169-14" |  (substr(titre,11,16))== "169-15" |  (substr(titre,11,16))== "169-16" |  (substr(titre,11,16))== "169-17" |  (substr(titre,11,16))== "169-18"
x3= (substr(titre,11,16))== "222-01" |  (substr(titre,11,16))== "222-02" |  (substr(titre,11,16))== "222-03" |  (substr(titre,11,16))== "222-04" |  (substr(titre,11,16))== "222-05" |  (substr(titre,11,16))== "222-06"
x13=(substr(titre,11,16))== "222-07" |  (substr(titre,11,16))== "222-08" |  (substr(titre,11,16))== "222-09" |  (substr(titre,11,16))== "222-10" |  (substr(titre,11,16))== "222-11" |  (substr(titre,11,16))== "222-12"
x23=(substr(titre,11,16))== "222-13" |  (substr(titre,11,16))== "222-14" |  (substr(titre,11,16))== "222-15" |  (substr(titre,11,16))== "222-16" |  (substr(titre,11,16))== "222-17" |  (substr(titre,11,16))== "222-18"
x4= (substr(titre,11,16))== "471-01" |  (substr(titre,11,16))== "471-02" |  (substr(titre,11,16))== "471-03" |  (substr(titre,11,16))== "471-04" |  (substr(titre,11,16))== "471-05" |  (substr(titre,11,16))== "471-06"
x14=(substr(titre,11,16))== "471-07" |  (substr(titre,11,16))== "471-08" |  (substr(titre,11,16))== "471-09" |  (substr(titre,11,16))== "471-10" |  (substr(titre,11,16))== "471-11" |  (substr(titre,11,16))== "471-12"
x24=(substr(titre,11,16))== "471-13" |  (substr(titre,11,16))== "471-14" |  (substr(titre,11,16))== "471-15" |  (substr(titre,11,16))== "471-16" |  (substr(titre,11,16))== "471-17" |  (substr(titre,11,16))== "471-18"
x5= (substr(titre,11,16))== "509-01" |  (substr(titre,11,16))== "509-02" |  (substr(titre,11,16))== "509-03" |  (substr(titre,11,16))== "509-04" |  (substr(titre,11,16))== "509-05" |  (substr(titre,11,16))== "509-06"
x15=(substr(titre,11,16))== "509-07" |  (substr(titre,11,16))== "509-08" |  (substr(titre,11,16))== "509-09" |  (substr(titre,11,16))== "509-10" |  (substr(titre,11,16))== "509-11" |  (substr(titre,11,16))== "509-12"
x25=(substr(titre,11,16))== "509-13" |  (substr(titre,11,16))== "509-14" |  (substr(titre,11,16))== "509-15" |  (substr(titre,11,16))== "509-16" |  (substr(titre,11,16))== "509-17" |  (substr(titre,11,16))== "509-18"
x6= (substr(titre,11,16))== "525-01" |  (substr(titre,11,16))== "525-02" |  (substr(titre,11,16))== "525-03" |  (substr(titre,11,16))== "525-04" |  (substr(titre,11,16))== "525-05" |  (substr(titre,11,16))== "525-06"
x16=(substr(titre,11,16))== "525-07" |  (substr(titre,11,16))== "525-08" |  (substr(titre,11,16))== "525-09" |  (substr(titre,11,16))== "525-10" |  (substr(titre,11,16))== "525-11" |  (substr(titre,11,16))== "525-12"
x26=(substr(titre,11,16))== "525-13" |  (substr(titre,11,16))== "525-14" |  (substr(titre,11,16))== "525-15" |  (substr(titre,11,16))== "525-16" |  (substr(titre,11,16))== "525-17" |  (substr(titre,11,16))== "525-18"
x7= (substr(titre,11,16))== "685-01" |  (substr(titre,11,16))== "685-02" |  (substr(titre,11,16))== "685-03" |  (substr(titre,11,16))== "685-04" |  (substr(titre,11,16))== "685-05" |  (substr(titre,11,16))== "685-06"
x17=(substr(titre,11,16))== "685-07" |  (substr(titre,11,16))== "685-08" |  (substr(titre,11,16))== "685-09" |  (substr(titre,11,16))== "685-10" |  (substr(titre,11,16))== "685-11" |  (substr(titre,11,16))== "685-12"
x27=(substr(titre,11,16))== "685-13" |  (substr(titre,11,16))== "685-14" |  (substr(titre,11,16))== "685-15" |  (substr(titre,11,16))== "685-16" |  (substr(titre,11,16))== "685-17" |  (substr(titre,11,16))== "685-18"
x8= (substr(titre,11,16))== "747-01" |  (substr(titre,11,16))== "747-02" |  (substr(titre,11,16))== "747-03" |  (substr(titre,11,16))== "747-04" |  (substr(titre,11,16))== "747-05" |  (substr(titre,11,16))== "747-06"
x18=(substr(titre,11,16))== "747-07" |  (substr(titre,11,16))== "747-08" |  (substr(titre,11,16))== "747-09" |  (substr(titre,11,16))== "747-10" |  (substr(titre,11,16))== "747-11" |  (substr(titre,11,16))== "747-12"
x28=(substr(titre,11,16))== "747-13" |  (substr(titre,11,16))== "747-14" |  (substr(titre,11,16))== "747-15" |  (substr(titre,11,16))== "747-16" |  (substr(titre,11,16))== "747-17" |  (substr(titre,11,16))== "747-18"
x9= (substr(titre,11,16))== "787-01" |  (substr(titre,11,16))== "787-02" |  (substr(titre,11,16))== "787-03" |  (substr(titre,11,16))== "787-04" |  (substr(titre,11,16))== "787-05" |  (substr(titre,11,16))== "787-06"
x19=(substr(titre,11,16))== "787-07" |  (substr(titre,11,16))== "787-08" |  (substr(titre,11,16))== "787-09" |  (substr(titre,11,16))== "787-10" |  (substr(titre,11,16))== "787-11" |  (substr(titre,11,16))== "787-12"
x29=(substr(titre,11,16))== "787-13" |  (substr(titre,11,16))== "787-14" |  (substr(titre,11,16))== "787-15" |  (substr(titre,11,16))== "787-16" |  (substr(titre,11,16))== "787-17" |  (substr(titre,11,16))== "787-18"
x10=(substr(titre,11,16))== "877-01" |  (substr(titre,11,16))== "877-02" |  (substr(titre,11,16))== "877-03" |  (substr(titre,11,16))== "877-04" |  (substr(titre,11,16))== "877-05" |  (substr(titre,11,16))== "877-06"
x20=(substr(titre,11,16))== "877-07" |  (substr(titre,11,16))== "877-08" |  (substr(titre,11,16))== "877-09" |  (substr(titre,11,16))== "877-10" |  (substr(titre,11,16))== "877-11" |  (substr(titre,11,16))== "877-12"
x30=(substr(titre,11,16))== "877-13" |  (substr(titre,11,16))== "877-14" |  (substr(titre,11,16))== "877-15" |  (substr(titre,11,16))== "877-16" |  (substr(titre,11,16))== "877-17" |  (substr(titre,11,16))== "877-18"

sp1=sp[(x1==TRUE),] #1
sp2=sp[(x2==TRUE),] #2
sp3=sp[(x3==TRUE),] #3
sp4=sp[(x4==TRUE),] #4
sp5=sp[(x5==TRUE),] #5
sp6=sp[(x6==TRUE),] #6
sp7=sp[(x7==TRUE),] #7
sp8=sp[(x8==TRUE),] #8
sp9=sp[(x9==TRUE),] #9
sp10=sp[(x10==TRUE),] #10
sp11=sp[(x11==TRUE),] #11
sp12=sp[(x12==TRUE),] #12
sp13=sp[(x13==TRUE),] #13
sp14=sp[(x14==TRUE),] #14
sp15=sp[(x15==TRUE),] #15
sp16=sp[(x16==TRUE),] #16
sp17=sp[(x17==TRUE),] #17
sp18=sp[(x18==TRUE),] #18
sp19=sp[(x19==TRUE),] #19
sp20=sp[(x20==TRUE),] #20
sp21=sp[(x21==TRUE),] #21
sp22=sp[(x22==TRUE),] #22
sp23=sp[(x23==TRUE),] #23
sp24=sp[(x24==TRUE),] #24
sp25=sp[(x25==TRUE),] #25
sp26=sp[(x26==TRUE),] #26
sp27=sp[(x27==TRUE),] #27
sp28=sp[(x28==TRUE),] #28
sp29=sp[(x29==TRUE),] #29
sp30=sp[(x30==TRUE),] #30

id=1:30
dim(id)=c(10,3)
id=t(id)

### Boucle pour effectuer plusieurs PLSDA (reduire l'impact tirage aleatoire)

## Initialisation vecteur de % de bons classements par VL
perok=vector(mode='numeric',length=ncmax)
## Creation matrice de % de mauvais classements par clone
mc=matrix(nrow = ncmax,ncol = c)

## Definition des matrices de resultat final
## Creation de la matrice des peroks finale
perok_finalm0=matrix(nrow = repet, ncol = ncmax)


## Creation de la matrice des VL et perok maximaux
maxi_final=matrix(nrow= repet, ncol = 2)
## Creation de la matrice de % de mauvais classements
mc_final=matrix(nrow= repet, ncol = length(levels(class)))
## Creation d'un matrice cubique pour enregistrer les tables de contingence
t_final=array(dim=c(c,c,repet))
## Noms des colonnes et des lignes
colnames(t_final)=c(basename(levels(class)))
rownames(t_final)=c(basename(levels(class)))
colnames(maxi_final)= c("maxi.id","perok max")
colnames(mc_final)= c(basename(levels(class)))

for(j in 1:repet) {
  
  grps=matrix( ,nrow=3, ncol=10)
  for (i in 1:10) {
    grps[,i]=id[sample(1:3),i]
  }
  
## Creation des jeux d'apprentissage et validation
#  flds <- createFolds(1:sum(iok), k = k)
  predm0=as.data.frame(matrix(nrow = sum(iok2), ncol = ncmax))
#  predPLSDAKNN=as.data.frame(matrix(nrow = sum(iok), ncol = ncmax))
  
## Boucle sur les groupes de CV
 for (i in 1:k) {
 #  commd=paste("sp_val=rbind(sp",grps[i,1],",sp",grps[i,2],",sp",grps[i,3],")",sep="")   # pour CS et G
 #  commd=paste("sp_val=rbind(sp",grps[i,1],",sp",grps[i,2],",sp",grps[i,3],",sp",grps[i,4],")",sep="") # pour S
   commd=paste("sp_val=rbind(sp",grps[i,1],",sp",grps[i,2],",sp",grps[i,3],",sp",grps[i,4],",sp",grps[i,5],",sp",grps[i,6],",sp",grps[i,7],",sp",grps[i,8],",sp",grps[i,9],",sp",grps[i,10],")",sep="")
    eval(parse(text=commd))
    
    id_val=which(rownames(sp)  %in%  rownames(sp_val))

    sp_val=sp[id_val,]       # matrice du jeu de validation 
    class_val=class[id_val]  #identifiants des classes du jeu de validation
    sp_cal=sp[-id_val,]      #matrice du jeu de calibration composée de tout ce qui n'est pas en validation
    class_cal=class[-id_val] #identifiants des classes du jeu de calibration
    
    ## PLSDA and application to have loadings and scores
    rplsda=caret::plsda(sp_cal, class_cal,ncomp=ncmax)
    sc_cal=rplsda$scores
    #sp_cal_c=scale(sp_cal,center=rplsda$Xmeans,scale = F)
    #sc_cal=sp_cal_c%*%rplsda$projection 
    sp_val_c=scale(sp_val,center=rplsda$Xmeans,scale = F)
    sc_val=sp_val_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas
    
    #rpca=prcomp(sp_cal)
    #spca_cal=rpca$x
    #spca_val=predict(rpca,sp_val)
    #fm = nirs::plsda(sp_cal,class_cal,sp_val,class_val,ncomp = ncmax,typda = "daknn", k = 15)
    
    for (ii in 2:ncmax) {
      ## Apprentissage    
      #rlda=lda(sc_cal[,1:ii], class_cal,tol=1.0e-5)
      
      ## Validation
      predm0[id_val,ii]=SIGNE_maha0(sc_cal[,1:ii], class_cal, sc_val[,1:ii])$class
      #predPLSDAKNN[id_val,ii]=
    }
  }
 
  ## Table de contingence
  tsm0=lapply(as.list(predm0), class, FUN = table)

  ## Matrice mauvais classements par clone
  #mc[i,]=(rowSums((as.matrix(t)))-diag(as.matrix(t)))/rowSums((as.matrix(t)))
  diagsm0=lapply(tsm0, FUN = diag)
  
  ## Pourcentage de bien classes
  perokm0=100*unlist(lapply(diagsm0, FUN = sum))/length(class)
  
  ## Pourcentage de bien classes
  #perok=100*unlist(lapply(diags, FUN = sum))/length(class)
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

}

## Affichage des matrices de resultat final
cat("RESULTATS FINAUX SUR", repet, "TIRAGES", "\n")
# print("Ensemble des tables de contingence:")
# print(t_final)
print("Mauvais classements par clone:")
print(mc_final)
# print("moyenne de mauvais classements par clone:")
# mc_final=mc_final[complete.cases(mc_final),]
# print(colMeans(mc_final))
# print(perok_final)
# print("Nombre de DV au max et maximum de bon classements:")
print("Pourcentages d'erreur:")
print(maxi_final)

## Tracage de l'evolution des perok en fonction du nombre de DV utilisees
plot(colMeans(perok_finalm0), xlab= "Nombre de VL", ylab = "Pourcentage de biens classés",pch=19, cex=1.5)
legend(ncmax*2/3,15,legend=c("Maha on PLSDA scores"),
       col=c("black"), lty=1, cex=0.8)
#legend(ncmax*2/3,15,legend=c("Maha on PLSDA scores", "Predict on LDA+PLSDA scores"),
#       col=c("black","blue"), lty=1, cex=0.8)

##Tracage des moyennes de mauvais classements par clone
# plot(colMeans(mc_final))
print("Moyenne de DV max et de perok:")
print(colMeans(maxi_final))

mc_final[!rowSums(!is.finite(mc_final)),]
mc_final[!is.finite(mc_final)] <- 0

mc_cep=vector(mode="logical", length = 3)
mc_cep[1]=mean(colMeans(mc_final)[1],colMeans(mc_final)[2],colMeans(mc_final)[7])
mc_cep[2]=mean(colMeans(mc_final)[3],colMeans(mc_final)[5],colMeans(mc_final)[9])
mc_cep[3]=mean(colMeans(mc_final)[4],colMeans(mc_final)[6],colMeans(mc_final)[8],colMeans(mc_final)[10])
names(mc_cep)=c("cabernet sauv","gamay","syrah")
print("pourcentage moyen d'erreurs de classement des clones par cepage:")
print(mc_cep)

############END#################
