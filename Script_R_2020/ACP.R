
# ACP 

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
source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha0.R')


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
# set.seed(1)

brb4="C:/Users/avitvale/Documents/Test/globalmatrix"
load(file=brb4)

## Filtrage des spectres aberrants
# globalmatrix=globalmatrix[globalmatrix[,500]>0.6,]
# globalmatrix=globalmatrix[globalmatrix[,1]<0.2,]
# globalmatrix=globalmatrix[globalmatrix[,2000]<0.25,]

# Data Filter
# Select dates
dates=list(
  # "20180619T"
  # ,"20180619P"
  # ,"20180619N"
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
 # "20180724T"
  "20180809T"
   ,"20180724N"
  # "20180816X"
  # ,"20180816E"
  # ,"20180710X"
  # ,"20180710E"
  # "20180816A"
  # ,"20180816E"
  # ,"20180731A"
  # ,"20180731E"
  
  
)
print(dates)
#Esca <- read.table("C:/Users/Noémie/Desktop/Test_poids_ok/Esca.csv",
#                   header=TRUE, sep=";",dec=".",row.names=1, check.names=FALSE,
#                   fileEncoding="latin1")

### FIXATION DES PARAMETRES UTILISES:
## Nombre de repetitions de la boucle de PLSDA:
repet=4
## Parametres du Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
p=2
n=11
m=1
## Nombre de VL max autorisees
ncmax=10
## Taille de l'echantillon de validation (1/v):
#v=3
## Nombre de groupes de CV
k=5
##Variable pour la moyenne des spectres
# z=6
# moy= 220
# print(globalmatrix)
# mean1_final=matrix(nrow = 1, ncol = ncol(globalmatrix))
# mean2_final=matrix(nrow = 220, ncol = ncol(globalmatrix))
# 
# mean1_final=colMeans(globalmatrix[(1:z),])

# for(j in 1:moy) {
#   mean2=colMeans(globalmatrix[((z+1):(z+6)),])
#   z=z+6
#   #print(mean2)
#   #print(z)
#   
#   mean2_final[j,]=mean2
# }

#sp=rbind(mean1_final,mean2_final)
sp=globalmatrix
### Pretraitements
## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
sp=adj_asd(sp,c(602,1402))
## Reduction des variables (extremites bruitees)
sp=sp[,seq(51,ncol(sp)-30,1)]
## SNV
sp=t(scale(t(sp)))
# ## Derivation Savitsky Golay
sp=t(apply(sp,1,sgolayfilt,p=p,n=n,m=m))

##Reunion des spectres et des valeurs de reference
#rownames(sp)=rownames(Esca)

#spok=cbind(sp,Esca)
#spok=cbind(sp,as.numeric(substr(rownames(sp),11,13)))  #pour les clones
spok=cbind(sp,as.numeric(substr(rownames(sp),1,8)))  #pour les jours
#spok=sp
print(spok[,1993])
str(sp)
# ## Filtrage des spectres aberrants
# spok=spok[spok[,500]>0.6,]
# spok=spok[spok[,1]<0.2,]



# globalmatrix2 <- read.table("C:/Users/Noémie/Desktop/SFE/Caracteristiques_agro/globalmatrix2.csv", 
#                             header=TRUE, sep=";",dec=".", row.names=1, check.names=FALSE, 
#                             fileEncoding="latin1")
# 
# globalmatrix2=globalmatrix2[globalmatrix2[,500]>0.6,]
# globalmatrix2=globalmatrix2[globalmatrix2[,1]<0.2,]
# globalmatrix2=globalmatrix2[globalmatrix2[,2000]<0.25,]

##AVEC PACKAGE FACTOEXTRA

#x=substr(rownames(spok),1,9)=="MPG 10/33"
#spok=spok[(x==FALSE),]

spok.pca <- PCA(spok, graph = TRUE)
print(spok.pca)
str(spok)
# Valeurs propres et graphique des valeurs propres
eig.val <- get_eigenvalue(spok.pca)
print (eig.val)
p <-fviz_eig(spok.pca, addlabels = TRUE) +
labs(title = "Graph des valeurs propres", x= "Valeurs propres" ,
y= "Pourcentage de variabilité expliquée (%)") +
 theme_minimal()
print (p)


## Graphique des individus

spok.pca<-PCA(spok, graph=T, axes=c(3,4),quali.sup=1993)
plot(spok.pca, habillage = 1993, axes=c(1,2))


##TESTS


acp3<-PCA(sp, scale.unit=F, ncp=5, graph=T, axes=c(1,2) )
axeX <- acp3$ind$coord[,1] ; axeY <- acp3$ind$coord[,2] 
plot(axeX,axeY,pch=16);grid()

###

d <- fviz_pca_ind (spok.pca, axes = c(1,2) , mean.point = FALSE, legend.title = "Modalités", axes.linetype = "solid",
              habillage= "none", col.ind=3, pointsize = 5) + 
             scale_color_brewer(palette="Set1")+
              labs(title = "ACP: graph des individus") +
                theme_minimal()
#Pb dans habillage
print(d)
str(spok[,1993])
## Graphique des individus
e <- fviz_pca_ind (spok.pca, axes = c(1,3) , mean.point = FALSE, legend.title = "Modalités", axes.linetype = "solid", label="none",
                     habillage =  spok$Degré, pointsize = 5) + 
  scale_color_brewer(palette="Set1")+
  labs(title = "ACP: graph des individus") +
  theme_minimal()
print(e)

## Graphique des individus
g <- fviz_pca_ind (spok.pca, axes = c(1,4) , mean.point = FALSE, legend.title = "Modalités", axes.linetype = "solid", label="none",
                   habillage =  spok$Degré, pointsize = 5) + 
  scale_color_brewer(palette="Set1")+
  labs(title = "ACP: graph des individus") +
  theme_minimal()
print(g)

####END#####





#rpca=prcomp(globalmatrix2) 
#scor=rpca$x 
# Representation graphique de l'ACP: 
# plot(scor, col=as.factor(substr(rownames(globalmatrix2),1,3))) 





#library(FactoMineR)

#globalmatrix2 <- read.table("C:/Users/Noémie/Desktop/SFE/Resultats/PTN1/globalmatrix.csv",
#                        header=TRUE, sep=";",dec=".", row.names=1, check.names=FALSE, fileEncoding="latin1")
#res <- PCA(globalmatrix2[,1:2072])
#plot(res, cex=0.8, invisible="quali", label="none", title="Graphe des individus")
#plot(res, choix= "ind", cex=0.8, invisible="quali", label="none", title="Graphe des individus", axes= 2:3)
#plot(res, cex=0.8, invisible="quali", title="Graphe des individus")



#brb="C:/Users/Noémie/Desktop/SFE/Resultats/PTN1/globalmatrix2"
#load(file=brb)
#rpca=PCA(globalmatrix,graph=FALSE)

#condition = as.factor(condition)
#df=data.frame(rpca=rpca$ind$coord,condition=condition)
#p=ggplot(data = globalmatrix2, aes(rpca.Dim.1,rpca.Dim.2, colour=dates)) + geom_point() + geom_text(aes(label=condition),hjust=0, vjust=0)
#plot(p)

