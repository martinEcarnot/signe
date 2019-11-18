#package esquisse

library(FactoMineR)
library(factoextra)
library(pls)
library(MASS)
library(signal)
library(plyr)
library(caret)
library(dplyr)
library("ellipse")
library(rgl)
library(ade4)
library("plotrix")
library(ggplot2)
library(plotly)


rm(list = ls())

#source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')

brb4="C:/Users/avitvale/Documents/Test/globalmatrix" #Déplacer globalmatrix
load(file=brb4)
#dates=list("20180809T","20180724N","20170710P")

### FIXATION DES PARAMETRES UTILISES:
repet=4
p=2
n=11
m=1
ncmax=10
k=5


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

coloclone=c(rouge, rouge2, rouge3, bleu, bleu2, bleu3, vert, vert2, vert3, vert4)
colocepage=c(rouge,bleu,vert)


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
sp6=sp4

## FILTRE
date=as.factor(substr(rownames(sp6),5,8))
parc=as.factor(substr(rownames(sp6),18,18))
DATE=c("0703","0704","0710","0709","0702")
#DATE="0702"
sp3=sp6[which(date  %in%  DATE),]
parc=as.factor(substr(rownames(sp3),18,18))
PARC="G"
sp5=sp3[which(parc %in% PARC),]


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
colo=as.factor(substr(rownames(sp5),9,9))
annee=as.factor(substr(rownames(sp5),4,4))
#forme=as.factor(substr(rownames(sp5),18,18))
texte=as.factor(substr(rownames(sp5),9,9))

#2D

 plot(axeX,axeY,pch=16, col=coloclone[colo], type="n",
      main=paste0("Représentation en f des cépages"),xlab=paste0("VL ",axe1),ylab=paste0("VL ",axe2));grid();

text(axeX,axeY,annee, col=colocepage[colo])

text(axeX,axeY,levels(colo), col=as.numeric(annee))

#legend(x="right", legend=unique(colo), col=unique(coloclone), pch=15, bg="white")

 #legend(x="left", legend=unique(forme), col=1, pch=unique(as.numeric(forme)), bg="white")


# for (i in levels(colo)) {
#   x_ell <- axeX[colo==i] ; y_ell <- axeY[colo==i] ; xy_ell <- data.frame(x_ell,y_ell)
#   lines(ellipse(cov(xy_ell),centre=colMeans(xy_ell),level=0.95),type="l", lty=1, col=coloclone[which(levels(colo)==i)])
# }

for (i in levels(annee)) {
  x_ell <- axeX[annee==i] ; y_ell <- axeY[annee==i] ; xy_ell <- data.frame(x_ell,y_ell)
  lines(ellipse(cov(xy_ell),centre=colMeans(xy_ell),level=0.95),type="l", lty=1, col=which(levels(annee)==i))
}
legend(x="right", legend=unique(annee), col=c(1,2,3), pch=15, bg="white")

###3D



plot3d(axeX[colo==levels(colo)[1]],axeY[colo==levels(colo)[1]],axeZ[colo==levels(colo)[1]],

       col=colocepage[which(levels(colo)==levels(colo)[1])],radius=0.2, type="p",xlab="Dim1",ylab="Dim2",zlab="Dim3")

for (i in levels(colo)) {
  points3d(axeX[colo==i],axeY[colo==i],axeZ[colo==i], col=colocepage[which(levels(colo)==i)],radius=0.2)
}



for (i in levels(colo)) {
  x = rplsda$scores[,axe1][colo==i] ; y = rplsda$scores[,axe2][colo==i] ; z = rplsda$scores[,axe3][colo==i]

  ellipse <- ellipse3d(cov(cbind(x,y,z)), centre=c(mean(x), mean(y), mean(z)), level = 0.95)

  plot3d(ellipse, col=colocepage[which(levels(colo)==i)], alpha = 0.05, add = TRUE, type="shade")
}





###



sp6=data.frame(sp=sp5,axeX=rplsda$scores[,axe1],axeY=rplsda$scores[,axe2],jour=substr(rownames(sp5),4,8),annee=substr(rownames(sp5),1,4),cepage=substr(rownames(sp5),9,9),clone=substr(rownames(sp5),9,13),parcelle=substr(rownames(sp5),18,18))
ANNEE="2018"
sp7=sp6[which(sp6$annee  %in%  ANNEE),]
# Nuage de points avec estimation de la densité 2d
sp <- ggplot(sp7, aes(x=axeX, y=axeY,colour=clone)) +
  geom_point(size=2, alpha=1) +
  scale_color_manual(values = coloclone) #+
#  geom_density_2d(data = sp6, size=0.3, bins=3)
ggplotly(sp)
