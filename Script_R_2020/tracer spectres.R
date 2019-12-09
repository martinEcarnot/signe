library(MASS)
#library("mixOmics")
library(FactoMineR)
library(signal)
library(plyr)
library(caret)


rm(list = ls())

source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha.R')
source('C:/Users/avitvale/Documents/Script_R/SIGNE_maha0.R')


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
# set.seed(1)


#"C:/Users/No?mie/Desktop/SFE/Toutesdates/Globalmatrix/globalmatrix"
brb4="C:/Users/avitvale/Documents/Test/globalmatrix"
load(file=brb4)

# brb2="C:/Users/No?mie/Desktop/SFE/Toutesdates/Globalmatrix/globalmatrixB"
#  load(file=brb2)
#
# brb3="C:/Users/No?mie/Desktop/SFE/Toutesdates/Globalmatrix/globalmatrixG"
#  load(file=brb3)
#
# globalmatrix = rbind(globalmatrix1,globalmatrix2,globalmatrix3)
#
# brb4="C:\\Users\\No?mie\\Desktop\\SFE\\Toutesdates\\Globalmatrix\\"
#   save(globalmatrix, file=paste(brb4,"globalmatrix",sep=""))
#
# write.table(globalmatrix, file=paste(brb4,"globalmatrix.csv",sep=""),sep=";", quote=FALSE)



# # Data Filter
# # Select dates
 dates=list(

#   "20180823N"
#   ,"20180823P"
#   ,"20180823T"
#   # "20180619P"
#   # ,"20180627P"
#   #  "20180704P"
#   #"20180710P"
#   # ,"20180710E"
#   # "20180619T"
#   # "20180619N"
#   # ,"20180627T"
#   # ,"20180627N"
#   # ,"20180709T"
#   # ,"20180709N"
#
 )
#
#iok=substr(rownames(globalmatrix),1,8) %in% dates
#globalmatrix=globalmatrix[iok,]
#

 ### FIXATION DES PARAMETRES UTILISES:
 p=2
 n=11
 m=1

 sp=globalmatrix
 ## SNV
 sp3=t(scale(t(sp)))
 # ## Derivation Savitsky Golay
 sp4=t(apply(sp3,1,sgolayfilt,p=p,n=n,m=m))


matplot(t(globalmatrix),pch = ".",xlab = "Longueurs d'ondes (nm)", ylab = "Transflectance")
matplot(t(sp4),pch = ".",xlab = "Longueurs d'ondes (nm)", ylab = "Transflectance")
#
# ##Filtrer les feuilles par age
# titre=rownames(globalmatrix)
# #w=as.numeric(substr(titre,15,16))==1 | as.numeric(substr(titre,15,16))==2 | as.numeric(substr(titre,5,6))==3 | as.numeric(substr(titre,15,16))==7 | as.numeric(substr(titre,15,16))==8 | as.numeric(substr(titre,5,6))==9 | as.numeric(substr(titre,15,16))==13 | as.numeric(substr(titre,15,16))==14 | as.numeric(substr(titre,15,16))==15
# # ## retire les lignes correspondantes (a mettre en commentaire si pas de selection de feuilles ou cepages)
# #globalmatrix=globalmatrix[(w==TRUE),]
#
# ## FIXATION DES PARAMETRES UTILISES:
# # nombre de repetitions de la boucle de FDA:
# #repet= 50
# #Parametres du Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
# p=2
# n=11
# m=1
# # nombre de DV max autorisees
# ncmax=15
# #Taille de l'echantillon de validation (1/v):
# #v=3
# # Nombre de groupes de CV
# #k=6
#
# ## LDA ##
# sp=globalmatrix
#
# # creation de la matrice de classes
# class=as.factor(substr(rownames(sp),11,13))
# # variable qui mesure le nombre de classes
# c=length(levels(class))
#
# # ## Pretraitements
# #  # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
# # sp=adj_asd(sp,c(602,1402))
# #  # Reduction des variables (extremites bruitees)
# # sp=sp[,seq(1,ncol(sp)-30,1)]
# #  # SNV
# # sp=t(scale(t(sp)))
# #  # Derivation Savitsky Golay
# #  sp=t(apply(sp,1,sgolayfilt,p=p,n=n,m=m))
#
# # ##Selectionne les lignes de la matrice en fonction des modalites si pretraitements, remplacer "globalmatrix" par "sp"
# # titre=rownames(sp)
# # x1=substr(titre,9,9)== "N"
# # x2=substr(titre,9,9)== "P"
# # x3=substr(titre,9,9)== "T"
# #
# #
# # sp1=sp[(x1==TRUE),] #N
# # sp2=sp[(x2==TRUE),] #P
# # sp3=sp[(x3==TRUE),] #T
# #
#
#
# # #Plot
# # matplot(t(sp1),pch = ".",col = "red") #N
# # matplot(t(sp2), pch = ".",col = "blue", add = TRUE) #P
# # matplot(t(sp3),pch = "." ,col = "green",add = TRUE) #T
# # legend(x="topright", legend=c("N","P","T"), col=c("red","blue" ,"green"),
# #        lty =  c(1,1,1))

