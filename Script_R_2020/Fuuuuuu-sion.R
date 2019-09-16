library(FactoMineR)
library(factoextra) 
library(pls)
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
#source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha0.R')


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
### Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
sp1=adj_asd(sp,c(602,1402))
## Reduction des variables (extremites bruitees)
sp2=sp1[,seq(51,ncol(sp1)-30,1)]
## SNV
sp3=t(scale(t(sp2)))
# ## Derivation Savitsky Golay
sp4=t(apply(sp3,1,sgolayfilt,p=p,n=n,m=m))


#spok=cbind(spok1, )
#print(rownames(spok1)[1])
#spok.pca<-PCA(spok, graph=T, axes=c(3,4),quali.sup=1993)
#plot(spok.pca, habillage = 1993, axes=c(1,2))

###ICI déterminer les classes.

#spok1=cbind(sp4,as.numeric(substr(rownames(sp4),11,13))) #Pour les clones
#spok1=cbind(sp4,as.numeric(substr(rownames(sp4),1,8))) #Pour les dates
#spok1=cbind(sp4,as.numeric(substr(rownames(sp4),1,4))) #Pour l'année

## Creation de la matrice de classes
clas=as.factor(substr(rownames(sp4),1,8))
# axes PLS
rplsda=caret::plsda(sp4, clas,ncomp=20)
axeX <- rplsda$scores[,1] ; axeY <- rplsda$scores[,2] 

#acp = PCA(spok1, scale.unit=F, ncp=5, graph=F, quali.sup = 1993)

#axeX <- acp$ind$coord[,1] ; axeY <- acp$ind$coord[,2] 

#Tracer le graphique

plot(axeX,axeY,pch=16,xlab="COUCOU", col=spok1[,1993]+3);grid()



legend(x="topright", legend=unique(spok1[,1993]), col=unique(spok1[,1993]+3), pch=16,bg="white")

#text(axeX,axeY,rownames(x)) #Afficher un texte au dessus des points, par ex leur nom ou leur cépage...

#plot(acp$li[,2])

























#plot(acp,choix="ind",habillage=1993)
#dimdesc(acp, axes = 1:2)


## Fonction pareto
# 
# pareto = function(x, bar.col="cyan", line.col="red", pch=16, h=80, h.lty=3,main="",xlab="Défauts",ylab="Fréquence (%)", names.arg=c(), ylab2="Cumul",mar=c(5,4,3,4)) {
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
