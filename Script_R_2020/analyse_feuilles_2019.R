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
spg=globalmatrixconservatoiregrappes

code=substr(rownames(spg),11,13)
rep=substr(rownames(spg),15,16)
spg2=sp2dfclo(spg, code, rep)
names(spg2)=c("code", "rep", "spectreg")
spg2$code=as.numeric(as.character(spg2$code))
spg2$rep=as.numeric(as.character(spg2$rep))



brb3="C:/Users/avitvale/Documents/Test/globalmatrixconservatoirefeuilles2019"
load(file=brb3)
sp=globalmatrixconservatoirefeuilles2019

sp=sp[-grep("787", rownames(sp)),]
clone_m=gsub("-.*", "", rownames(sp))
clone_m=gsub(".*(N )|.*(T )", "", clone_m)
clone_m=gsub(" ", "", clone_m)
clone_m=gsub("0","",clone_m)

rep=gsub(".*-", "", rownames(sp))
rep=gsub(" B", "", rep)

sp2=sp2dfclo(sp, clone_m, rep)
names(sp2)=c("clone_m", "rep", "spectre")
sp2$rep=as.numeric(as.character(sp2$rep))



#
# test=sp2$spectre
# xm=lapply(unique(sp2$code), function(x) colMeans(sp2$spectre[which(sp2$code==x),]))
#
# xm=unlist(xm)
#
# dim(xm)=c(1992,length(xm)/1992)
#
# xm=t(xm)
# sp3=sp2df(xm,unique(code))
# names(sp3)=c("code","spectre")
# sp3$code=as.numeric(as.character(sp3$code))





trad=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Correspondance_code_conservatoire_gamay.csv")

donnees1=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Conservatoire_2019_1.csv")

donnees2=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Conservatoire_2019_2.csv")
#donnees2=donnees2[complete.cases(donnees2),]



trad2=gsub("/", "", trad$clone)
trad2=gsub(" ", "", trad2)
trad2=gsub("0", "", trad2)
trad=data.frame(trad$code, trad$clone, trad2)
names(trad)=c("code", "clone", "clone_m")

setdiff(sp2$clone_m, trad$clone_m)



sp2$clone_m=gsub("2153","ARB729",sp2$clone_m)
sp2$clone_m=gsub("CJB827","CJP827",sp2$clone_m)
sp2$clone_m=gsub("CJB827","CJP827",sp2$clone_m)
sp2$clone_m=gsub("CJB858","CJP858",sp2$clone_m)
sp2$clone_m=gsub("CJP635","CJB635",sp2$clone_m)
sp2$clone_m=gsub("CJPR2113","CJPR213",sp2$clone_m)
sp2$clone_m=gsub("CJP635","CJB635",sp2$clone_m)
sp2$clone_m=gsub("COG771","CPG771",sp2$clone_m)
sp2$clone_m[1534]="GAL243"
sp2$clone_m[grep("SMHH",sp2$clone_m)]="SMHH837"
sp2$clone_m=gsub("HPG61","MPG61",sp2$clone_m)

# rownames(sp)[grep("2 43", rownames(sp))]
# unique(D[grep("BEC",D)])
# D[1536]
# setdiff(trad2,C)[grep("CJPR", setdiff(trad2,C))]





jonction1=left_join(trad, donnees1)
jonction1=jonction1[complete.cases(jonction1),]
jonction2=left_join(trad, donnees2)
jonction2=jonction2[complete.cases(jonction2),]

#table=left_join(sp2, jonction1)
#table=left_join(sp2,jonction2)
#table=left_join(sp2,jonction1)
table_1=left_join(sp2,jonction1)
table_2=left_join(table_1,spg2)
table=left_join(table_2,jonction2)
table=table[complete.cases(table),]

table$spectre=cbind(table$spectreg, table$spectre)
# plotsp(Test[1:200,])
# plotsp(table$spectreg[1:200,])
# plotsp(table$spectre[1:200,1:200])
#
# A=table$spectre[1:200,]
#
# B=A[A[,1]<0.3,]
#
# plotsp(B)
# setdiff(tableX$clone, table$clone)
# intersect(tableX$clone, table$clone)
# tableX$clone[-which(tableX$clone %in% table$clone)]
# setdiff(table$clone, tableX$clone)

# donnees1$clone
# C=donnees1$clone
# C=C[-which(C=="clone787")]
# A=sub(" ", "", C)
# A=sub("/", "", A)
#
# B=vector()
# for (i in 1:length(unique(A))){
#   B=c(B,length(which(A==unique(A)[i])))
# }
# unique(A)[which(B!=1)]
# C[which(A=="PVM203")]

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
sp_preg=adj_asd(table$spectreg,c(552,1352))
## SNV
sp_pre=t(scale(t(sp_pre)))
sp_preg=t(scale(t(sp_preg)))

## Derivation Savitsky Golay
sp=savitzkyGolay(sp_pre, m = m, p = p, w = n)
spg=savitzkyGolay(sp_preg, m = m, p = p, w = n)


sp=sp[,-(1330:1390)] #Coupure du spectre autour des résidus de sauts de detecteurs. Ne semble pas encore effacer completement les sauts.
sp=sp[,-(500:560)]

spg=spg[,-(1330:1390)] #Coupure du spectre autour des résidus de sauts de detecteurs. Ne semble pas encore effacer completement les sauts.
spg=spg[,-(500:560)]


table$spectre=sp
table$spectreg=spg

n=4
a=floor(length(unique(table$code))/n)
segm=list()
for (i in 1:10){
  l1=sample(n*a, a)
  l2=sample((1:(n*a))[-l1],a)
  l3=sample((1:(n*a))[-c(l1,l2)],a)
  l4=sample((1:(n*a))[-c(l1,l2,l3)],a)

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
  L4=vector()
  for (i in l4) {
    L4=c(L4,which(table$code==i))
  }

  segm_1=list(list(L1,L2,L3,L4))
  segm=c(segm,segm_1)
}

#length(unique(table$code)) #192. 192/6=32.
names(table)
#length(which(table$forme=="ronde-ovoide"))
#fm=fitcv(table$spectreg, table$pds_100baies, lwplsr, segm, print=T, ncomp=15, k=1000, ncompdis=15)
fm=fitcv(table$spectre, table$acidite_totale, plsr, segm, print=T, ncomp=20)

# z <- err(fm, ~ ncomp + rep)
#
# plotmse(z, nam = "errp", group="rep")
#
# fmextrait=lapply(fm,function (x) {x[x$ncomp==3,]})
#
# table(fmextrait$y$x1,fmextrait$fit$x1)

#plotsp(sp3$spectre)

z <- mse(fm, ~ ncomp + rep)
#str(z)

#plotmse(z, nam = "rmsep", group ="rep")
plotmse(z, nam = "r2", group ="rep")

fm20=lapply(fm,function (x) {x[x$ncomp==12,]})
plot(fm20$y$x1,fm20$fit$x1)
hist(table$acidite_totale)





# ncomp <- 10
# fm <- pca(table$spectre, ncomp = ncomp)
# #fm <- pca(zXr, ncomp = ncomp, algo = pca.eigen)
# #fm <- pca(zXr, ncomp = ncomp, algo = pca.nipals)
# names(fm)
#
# A=table[133,]
# ### Scores
# which(table$spectre>0.15)
# headm(fm$Tr)
#
# comp <- c(1, 2)
# #comp <- c(2, 3)
# z <-
#   plotxy(fm$Tr[, comp])
#
# plotxy(fm$Tr[, comp], group = table$arome)
#
# plotxy(fm$Tr[, comp], label = TRUE)



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







