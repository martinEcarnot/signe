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
d=load(file=brb3)
sp=globalmatrixconservatoiregrappes
#sp=sp[,1:1100]


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



test=sp2$spectre
xm=lapply(unique(sp2$code), function(x) colMeans(sp2$spectre[which(sp2$code==x),]))

xm=unlist(xm)

dim(xm)=c(1992,length(xm)/1992)
#dim(xm)=c(1100,length(xm)/1100)

xm=t(xm)
sp3=sp2dfclo(xm,unique(code),vector(length=208))
names(sp3)=c("code","condition","spectre")
sp3$code=as.numeric(as.character(sp3$code))
# plotsp(sp3$spectre)
# plotsp(sp2$spectre)


# xm=lapply(levels(dat$ID), function(x) colMeans(dat$x[which(dat$ID==x),]))
#
# xm=unlist(xm)
#
# dim(xm)=c(2151,length(xm)/2151)
#
# ym=aggregate(dat, list(dat$ID), mean)






# trad=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Correspondance_code_conservatoire_gamay.csv")
# donnees1=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Conservatoire_2019_1.csv")
# donnees2=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Conservatoire_2019_2.csv")
###donnees2=donnees2[complete.cases(donnees2),]

# Modif MEC juillet 2022
dnew="/run/user/1000/gvfs/smb-share:server=cubebe,share=agap/daav-ge2pop/Cereales/Blé/valentin_Avit/SIGNE Documents Stage 4/Moins utile/Données"
trad=read_csv2(file = file.path(dnew,"Correspondance_code_conservatoire_gamay.csv"))
donnees1=read_csv2(file = file.path(dnew,"Conservatoire_2019_1.csv"))
donnees2=read_csv2(file = file.path(dnew,"Conservatoire_2019_2.csv"))


jonction1=left_join(trad, donnees1)
jonction1=jonction1[complete.cases(jonction1),]
jonction2=left_join(trad, donnees2)
jonction2=jonction2[complete.cases(jonction2),]

#table=left_join(sp2, jonction1)
#table=left_join(sp2,jonction2)
table_1=left_join(sp3,jonction1)
table=left_join(table_1,jonction2)
table=table[complete.cases(table),]

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
n=11 #11 #Faire avec un n plus gros ?
m=1 #1
## Nombre de VL max autorisees
ncmax=35 #35
## Nombre de groupes de CV
k=2 #2

## PLSDA ##


### Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
sp_pre=adj_asd(table$spectre,c(552,1352)) #Retrouvé empiriquement, ne correspondait pas à ce qui etait indiqué precedemment.
#sp_pre=adj_asd(table$spectre,c(552))
sp_pre=sp_pre+000.1
## SNV
sp_pre=t(scale(t(sp_pre)))
#sp_pre=log(sp_pre)


## Derivation Savitsky Golay
sp=savitzkyGolay(sp_pre, m = m, p = p, w = n)
#sp=detrend(sp_pre, method="poly")


sp=sp[,-(1330:1390)] #Coupure du spectre autour des résidus de sauts de detecteurs. Ne semble pas encore effacer completement les sauts.
sp=sp[,-(500:560)]


table$spectre=sp

n=4
a=floor(length(unique(table$code))/n)
segm=list()
for (i in 1:3){
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
#
# critere1=which(table$pds_100baies<200)
# length(critere1)
#
# critere2=which(10<table$degre)
# length(critere2)
# critere3=which(table$degre<12)
# length(critere3)
#
# critere4=which(4.5<table$acidite_totale)
# length(critere4)
#
# critere5=which(table$pH<3.3)
# length(critere5)
#
# A=setdiff(critere1,critere2)
# B=intersect(groupe3,critere5)
# C=setdiff(which(table$degre<100000000),critere1)
# D=union(A,B)
# sort(union(D,C))
# groupe11=intersect(critere1,critere2)
# groupe1=intersect(critere1,critere2)
# groupe12=setdiff(critere1,critere2)
# groupe13=setdiff(critere2,critere1)
# groupe14=setdiff(1:195, union(critere1,critere2))
# length(groupe1)
# length(setdiff(critere1,critere2))
# length(setdiff(critere2,critere1))
# length(union(critere2,critere1))
# length(table$pds_100baies)
# length(groupe1)
# groupe2=intersect(groupe1,critere3)
# length(groupe2)
# groupe3=intersect(groupe2,critere4)
# length(groupe3)
# groupe4=intersect(groupe3,critere5)
# length(groupe4)
# length(table$acidite_totale)
# table$condition=F
# table$condition[groupe4]=T
# table$condition[groupe12]="B"
# table$condition[groupe13]="C"
# table$condition[groupe14]="D"

#length(unique(table$code)) #192. 192/6=32.
names(table)
#length(which(table$forme=="ronde-ovoide"))
#fm=fitcv(table$spectre, table$arome_simpl, lwplsdalm, segm, print=T, ncomp=10, k=100, ncompdis=3)

# fm=fitcv(table$spectre, table$fertilite, plsr, segm, print=T, ncomp=20)
fm=cvfit(table$spectre, table$pH, plsr, segm, print=T, ncomp=20)


##Pour plsda
# z <- err(fm, ~ ncomp + rep)
#
# plotmse(z, nam = "errp", group="rep")
#
# fmextrait=lapply(fm,function (x) {x[x$ncomp==5,]})
#
# table(fmextrait$y$x1,fmextrait$fit$x1)

#plotsp(sp3$spectre)

##Pour plsr
z <- mse(fm, ~ ncomp + rep)
#str(z)


#plotmse(z, nam = "rmsep", group ="rep")
plotmse(z, nam = "r2", group ="rep")
plotmse(z, nam = "rmsep", group ="rep")


fm20=lapply(fm,function (x) {x[x$ncomp==8,]})
plot(fm20$y$y1,fm20$fit$y1)
hist(table$fertilite)





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



#pds 100 baies <200   degré 10 à 12    acidité totale > 4,5    pH < 3,3




