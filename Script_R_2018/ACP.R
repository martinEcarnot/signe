
# ACP

library(FactoMineR)
library(factoextra)

globalmatrix2 <- read.table("C:/Users/No?mie/Desktop/SFE/Resultats/PTN1/globalmatrix2.csv", header=TRUE, sep=";",dec=".", row.names=1, check.names=FALSE, fileEncoding="latin1")
globalmatrix2=globalmatrix2[globalmatrix2[,500]>0.6,]
globalmatrix2=globalmatrix2[globalmatrix2[,1]<0.2,]
globalmatrix2=globalmatrix2[globalmatrix2[,2000]<0.25,]
#matplot(t(globalmatrix2),pch='.')

##AVEC PACKAGE FACTOEXTRA

globalmatrix2.pca <- PCA(globalmatrix2 [, -2073], graph = FALSE)
print(globalmatrix2.pca)

## Valeurs propres et graphique des valeurs propres
#eig.val <- get_eigenvalue(globalmatrix2.pca)
#print (eig.val)
#p <-fviz_eig(globalmatrix2.pca, addlabels = TRUE) +
#labs(title = "Graph des valeurs propres", x= "Valeurs propres" , y= "Pourcentage de variabilit? expliqu?e (%)") +
#  theme_minimal()
#print (p)


## Graphique des individus
#d <- fviz_pca_ind (globalmatrix2.pca, axes = c(1,2) , axes.linetype = "solid",
#                habillage =  globalmatrix2$Conditions) +
#                 scale_color_brewer(palette="Set1")+
#                 labs(title = "ACP: graph des individus") +
#                 theme_minimal()

#print(d)

####END#####





#rpca=prcomp(globalmatrix2)
#scor=rpca$x
# Representation graphique de l'ACP:
# plot(scor, col=as.factor(substr(rownames(globalmatrix2),1,3)))





#library(FactoMineR)

#globalmatrix2 <- read.table("C:/Users/No?mie/Desktop/SFE/Resultats/PTN1/globalmatrix.csv",
#                        header=TRUE, sep=";",dec=".", row.names=1, check.names=FALSE, fileEncoding="latin1")
#res <- PCA(globalmatrix2[,1:2072])
#plot(res, cex=0.8, invisible="quali", label="none", title="Graphe des individus")
#plot(res, choix= "ind", cex=0.8, invisible="quali", label="none", title="Graphe des individus", axes= 2:3)
#plot(res, cex=0.8, invisible="quali", title="Graphe des individus")



#brb="C:/Users/No?mie/Desktop/SFE/Resultats/PTN1/globalmatrix2"
#load(file=brb)
#rpca=PCA(globalmatrix,graph=FALSE)

#condition = as.factor(condition)
#df=data.frame(rpca=rpca$ind$coord,condition=condition)
#p=ggplot(data = globalmatrix2, aes(rpca.Dim.1,rpca.Dim.2, colour=dates)) + geom_point() + geom_text(aes(label=condition),hjust=0, vjust=0)
#plot(p)

