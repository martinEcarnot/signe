library(FactoMineR)
library(factoextra)
library(pls)
library(MASS)
#library(mixOmics)
library(signal)
library(plyr)
library(caret)
library(dplyr)
#library(nirs)
library("ellipse")
library(rgl)
library(ade4)
library("plotrix")



rm(list = ls())

source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')

brb4="C:/Users/avitvale/Documents/Test/globalmatrix"
load(file=brb4)
#dates=list("20180809T","20180724N","20170710P")

### FIXATION DES PARAMETRES UTILISES:
repet=4
p=2
n=11
m=1
ncmax=10
k=5


sp=globalmatrix
sp

iout=which(nchar(rownames(sp))>18)

sp <- sp[-iout,]
sp
###Pourquoi ???


## Filtrage des spectres aberrants
sp=sp[sp[,500]>0.6,]
sp=sp[sp[,1]<0.2,]
sp=sp[sp[,2000]<0.25,]

### Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
sp1=adj_asd(sp,c(602,1402))
## Reduction des variables (extremites bruitees)
sp2=sp1[,seq(51,ncol(sp1)-30,1)]
## SNV
sp3=t(scale(t(sp2)))
# ## Derivation Savitsky Golay
sp4=t(apply(sp3,1,sgolayfilt,p=p,n=n,m=m))

###Ajout du cepage au nom de la ligne.
##Filtre en fonction du cepage
titre=rownames(sp4)

ya=(substr(titre,11,13))== "015" |  (substr(titre,11,13))== "169" |  (substr(titre,11,13))== "685"
yb=(substr(titre,11,13))== "222" |  (substr(titre,11,13))== "509" |  (substr(titre,11,13))== "787"
yc=(substr(titre,11,13))== "471" |  (substr(titre,11,13))== "525" |  (substr(titre,11,13))== "747" |  (substr(titre,11,13))== "877"

##Separation en 3 matrices
spa=sp4[ya,]
spb=sp4[yb,]
spc=sp4[yc,]

##Rajoute le nom du cepage au nom de la ligne à la place de P/T/N
substr(rownames(spa),9,9)="C"
substr(rownames(spb),9,9)="G"
substr(rownames(spc),9,9)="S"
#substr(rownames(spa),9,9)="1"
#substr(rownames(spb),9,9)="2"
#substr(rownames(spc),9,9)="3"

##Recombine les 3 matrices pour reformer sp
sp4=rbind(spa,spb,spc)
rownames(sp4)


#### Chamanisme <- Ca change la couleur d'affichage. POur des raisons mystérieuses.

 data = read.table("cancers.txt",sep="\t",header=T)

 colnames(data)

 rownames(data) = data[,1] # On donne ? chaque ligne le nom de son cancer (Cancer 1, Cancer 2)

 x = data[1:7,-1] # On exclue la colonne 1 qui contenait les noms des cancers.


 acp3 = PCA(x, scale.unit=F, ncp=5, graph=T, axes=c(1,2))

####


## Creation de la matrice de classes
rownames(sp4)[1]
clas=as.factor(paste(as.factor(substr(rownames(sp4),1,4)),as.factor(substr(rownames(sp4),9,9)))) #cépages
clas=as.factor(substr(rownames(sp4),9,9))
# axes PLS
rplsda=caret::plsda(sp4, clas,ncomp=20)
axe1<-1
axe2<-2
axe3<-3
axeX <- rplsda$scores[,axe1] ; axeY <- rplsda$scores[,axe2] ; axeZ<- rplsda$scores[,axe3]

colo=as.numeric(clas)
colo

#acp = PCA(sp4, scale.unit=F, ncp=5, graph=F)
#axeX <- acp$ind$coord[,1] ; axeY <- acp$ind$coord[,2]

#Tracer le graphique

#plot3d(axeX,axeY,axeZ,pch=16, col=colo,
#     main=paste0("Représentation en f des cépages"),xlab=paste0("VL ",axe1),ylab=paste0("VL ",axe2));grid();



plot3d(axeX[clas==levels(clas)[1]],axeY[clas==levels(clas)[1]],axeZ[clas==levels(clas)[1]],

       col=which(levels(clas)==levels(clas)[1]),radius=0.2, type="p",xlab="Dim1",ylab="Dim2",zlab="Dim3")

for (i in levels(clas)) {
  points3d(axeX[clas==i],axeY[clas==i],axeZ[clas==i], col=which(levels(clas)==i),radius=0.2)
}



for (i in levels(clas)) {
  x = rplsda$scores[,axe1][clas==i] ; y = rplsda$scores[,axe2][clas==i] ; z = rplsda$scores[,axe3][clas==i]

  ellipse <- ellipse3d(cov(cbind(x,y,z)), centre=c(mean(x), mean(y), mean(z)), level = 0.95)

  plot3d(ellipse, col=which(levels(clas)==i), alpha = 0.2, add = TRUE, type="shade")
}











plot(axeX,axeY,pch=16, col=clas,
     main=paste0("Représentation en f des cépages"),xlab=paste0("VL ",axe1),ylab=paste0("VL ",axe2));grid();


legend(x="right", legend=unique(clas), col=unique(clas), pch=16,bg="white")

###Comment sait-on que la légende correspond ? Réponse : on sait pas. Hihi.
for (i in levels(clas)) {
  x_ell <- rplsda$scores[,axe1][clas==i] ; y_ell <- rplsda$scores[,axe2][clas==i] ; xy_ell <- data.frame(x_ell,y_ell)
  lines(ellipse(cov(xy_ell),centre=colMeans(xy_ell),level=0.95),type="l", lty=1, col=which(levels(clas)==i))
}








#text(axeX,axeY,rownames(x)) #Afficher un texte au dessus des points, par ex leur nom ou leur c?page...

#plot(acp$li[,2])























#plot(acp,choix="ind",habillage=1993)
#dimdesc(acp, axes = 1:2)


## Fonction pareto
#
# pareto = function(x, bar.col="cyan", line.col="red", pch=16, h=80, h.lty=3,main="",xlab="D?fauts",ylab="Fr?quence (%)", names.arg=c(), ylab2="Cumul",mar=c(5,4,3,4)) {
#
#   if (length(names.arg)>0) {names.arg=names.arg[order(x, decreasing = TRUE)]}
#
#   x = sort(x,decreasing=T); x = x*100/sum(x);
#
#   cumul = (cumsum(x)/sum(x))*100
#
#   simulation = barplot(x,col=bar.col) ; simulation
#
#   plot.new()
#
#   par(mar=mar)
#
#   barplot(x,col=bar.col,axes=F,ylim=c(0,100),main=main,xlab=xlab,ylab="",names.arg=names.arg)
#
#   #par(new=TRUE)
#
#   points(simulation,cumul,pch=pch,col=line.col,xlab="",ylab="",type="o")
#
#   abline(h=h,lty=h.lty) ; box()
#
#   axis(2) ; axis(4,c(0,20,40,60,80,100),col.axis=line.col,col=line.col)
#
#   mtext(ylab,side=2,line=2,cex=1.2) ; mtext(ylab2,side=4,col="red",line=2,cex=1.2)
#
#   result = c(x , cumul) ; result = matrix(result,nc=length(x), byrow=T)
#
#   if (length(names.arg)>0) {colnames(result) = names.arg }
#
#   rownames(result) = c("frequency","cumul")
#
#   return(result)}

####



#pareto(acp$eig[,2], h=95)
#acp
