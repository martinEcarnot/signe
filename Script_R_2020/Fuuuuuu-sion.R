

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
library(ggplot2)


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












### Pretraitements
## SNV
sp3=t(scale(t(sp)))
# ## Derivation Savitsky Golay
sp4=t(apply(sp3,1,sgolayfilt,p=p,n=n,m=m))



#### Chamanisme <- Ca change la couleur d'affichage. POur des raisons mystérieuses.

 data = read.table("cancers.txt",sep="\t",header=T)

 colnames(data)

 rownames(data) = data[,1] # On donne ? chaque ligne le nom de son cancer (Cancer 1, Cancer 2)

 x = data[1:7,-1] # On exclue la colonne 1 qui contenait les noms des cancers.


 acp3 = PCA(x, scale.unit=F, ncp=5, graph=T, axes=c(1,2))

####
rownames(sp4)[1]
tri=as.factor(substr(rownames(sp4),9,9))
levels(tri)
c=which(tri=="S")
sp5=sp4[c,]
sp5=sp4
length(sp5[,1])
# sp5=sp4
#iout=sample(5210,4000)
#sp5 <- sp5[-iout,]
length(sp5[1,])
## Creation de la matrice de classes
rownames(sp5)[1]
#clas=as.factor(paste(as.factor(substr(rownames(sp5),1,4)),as.factor(substr(rownames(sp5),15,16)))) #cépages
clas=as.factor(substr(rownames(sp5),11,13))

# axes PLS
rplsda=caret::plsda(sp5, clas,ncomp=20)
axe1<- 1
axe2<- 2
axe3<- 3
axeX <- rplsda$scores[,axe1] ; axeY <- rplsda$scores[,axe2] ; axeZ<- rplsda$scores[,axe3]

#acp = PCA(sp4, scale.unit=F, ncp=5, graph=F)
#axeX <- acp$ind$coord[,1] ; axeY <- acp$ind$coord[,2]

#Tracer le graphique

colo=as.factor(paste(as.factor(substr(rownames(sp5),1,4)),as.factor(substr(rownames(sp5),9,9)))) #cépages
colo=as.factor(substr(rownames(sp5),9,13))
#forme=as.factor(substr(rownames(sp5),18,18))
texte=as.factor(substr(rownames(sp5),9,9))

#2D

 plot(axeX,axeY,pch=16, col=coloclone[colo], #type="n",
      main=paste0("Représentation en f des cépages"),xlab=paste0("VL ",axe1),ylab=paste0("VL ",axe2));grid();

#text(axeX,axeY,texte, col=as.numeric(colo))

legend(x="right", legend=unique(colo), col=unique(coloclone), pch=15, bg="white")
#legend(x="left", legend=unique(forme), col=1, pch=unique(as.numeric(forme)), bg="white")

###Comment sait-on que la légende correspond ? Réponse : on sait pas. Hihi.
for (i in levels(colo)) {
  x_ell <- axeX[colo==i] ; y_ell <- axeY[colo==i] ; xy_ell <- data.frame(x_ell,y_ell)
  lines(ellipse(cov(xy_ell),centre=colMeans(xy_ell),level=0.95),type="l", lty=1, col=coloclone[which(levels(colo)==i)])
}



###3D



plot3d(axeX[colo==levels(colo)[1]],axeY[colo==levels(colo)[1]],axeZ[colo==levels(colo)[1]],

       col=which(levels(colo)==levels(colo)[1]),radius=0.2, type="p",xlab="Dim1",ylab="Dim2",zlab="Dim3")

for (i in levels(colo)) {
  points3d(axeX[colo==i],axeY[colo==i],axeZ[colo==i], col=which(levels(colo)==i),radius=0.2)
}



for (i in levels(colo)) {
  x = rplsda$scores[,axe1][colo==i] ; y = rplsda$scores[,axe2][colo==i] ; z = rplsda$scores[,axe3][colo==i]

  ellipse <- ellipse3d(cov(cbind(x,y,z)), centre=c(mean(x), mean(y), mean(z)), level = 0.95)

  plot3d(ellipse, col=which(levels(colo)==i), alpha = 0.2, add = TRUE, type="shade")
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


