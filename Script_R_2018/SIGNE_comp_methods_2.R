library(MASS)
library(mixOmics)
library(FactoMineR)
library(signal)
library(plyr)
library(caret)

rm(list = ls())

source('C:/Users/No?mie/Desktop/SFE/Script_R/adj_asd.R')
source('C:/Users/No?mie/Desktop/SFE/Script_R/SIGNE_load.R')
source('C:/Users/No?mie/Desktop/SFE/Script_R/SIGNE_maha.R')
source('C:/Users/No?mie/Desktop/SFE/Script_R/SIGNE_maha0.R')


# Choix de la fixation du tirage aleatoire (pour comparaison, rend les repetitions inutiles)
# set.seed(1)

brb="C:/Users/No?mie/Desktop/SFE/Resultats/PTN1/globalmatrix"
load(file=brb)

#globalmatrix[,500]>0.6
globalmatrix=globalmatrix[globalmatrix[,500]>0.6,]

# Data Filter
# Select dates
dates=list(
 "20180619P"
 ,"20180627P"
# ,"20180704P"
 ,"20180710P"
# ,"20180710E"
 ,"20180619T"
 ,"20180619N"
 ,"20180627T"
 ,"20180627N"
 ,"20180709T"
 ,"20180709N"
# "0529S"
# ,"0606S"
# ,"0612S"
# ,"0619S"
# ,"0626S"
# ,"0703S"
# ,"0710S"
# ,"0717S"
# ,"0724S"
# ,"0731S"
)
iok=substr(rownames(globalmatrix),1,9) %in% dates
globalmatrix=globalmatrix[iok,]

titre=rownames(globalmatrix)
## decommenter pour retirer les jeunes feuilles:
 #w=as.numeric(substr(titre,5,6))==4 | as.numeric(substr(titre,5,6))==5 | as.numeric(substr(titre,5,6))==6 | as.numeric(substr(titre,5,6))==10 | as.numeric(substr(titre,5,6))==11 | as.numeric(substr(titre,5,6))==12 | as.numeric(substr(titre,5,6))==16 | as.numeric(substr(titre,5,6))==17 | as.numeric(substr(titre,5,6))==18
# ## decommenter pour retirer les vieilles feuilles:
# # w=as.numeric(substr(titre,5,6))==1 | as.numeric(substr(titre,5,6))==2 | as.numeric(substr(titre,5,6))==3 | as.numeric(substr(titre,5,6))==7 | as.numeric(substr(titre,5,6))==8 | as.numeric(substr(titre,5,6))==9 | as.numeric(substr(titre,5,6))==13 | as.numeric(substr(titre,5,6))==14 | as.numeric(substr(titre,5,6))==15
# ## retire les lignes correspondantes (a mettre en commentaire si pas de selection de feuilles ou cepages)
# globalmatrix=globalmatrix[(w==TRUE),]
titre=rownames(globalmatrix)
## decommenter pour retirer la syrah:
# z=as.numeric(substr(titre,1,3))==471 | as.numeric(substr(titre,1,3))==525 | as.numeric(substr(titre,1,3))==747 | as.numeric(substr(titre,1,3))==877
## decommenter pour retirer le cabernet sauvignon:
# z=as.numeric(substr(titre,1,3))==015 | as.numeric(substr(titre,1,3))==169 | as.numeric(substr(titre,1,3))==685
## decommenter pour retirer le gamay:
# z=as.numeric(substr(titre,1,3))==787 | as.numeric(substr(titre,1,3))==509 | as.numeric(substr(titre,1,3))==222
## retire les lignes correspondantes (a mettre en commentaire si pas de selection de feuilles ou cepages)
# globalmatrix=globalmatrix

## FIXATION DES PARAMETRES UTILISES:
# nombre de repetitions de la boucle de FDA:
repet=50
#Parametres du Savitsky-Golay (p=degre du polynome, n= taille de la fenetre, m=ordre de derivation)
p=2
n=11
m=1
# nombre de DV max autorisees
ncmax=15
#Taille de l'echantillon de validation (1/v):
v=3
# Nombre de groupes de CV
k=6

## LDA ##
sp=globalmatrix

# creation de la matrice de classes
class=as.factor(substr(rownames(sp),11,13))
# variable qui mesure le nombre de classes
c=length(levels(class))

## Pretraitements
# # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651??me l.o.) et 1800 (1451))
sp=adj_asd(sp,c(651,1451))
# # Reduction des variables (extremites bruitees)
sp=sp[,seq(+1,ncol(sp)-30,1)]
# # SNV
sp=t(scale(t(sp)))
# # Derivation Savitsky Golay
sp=t(apply(sp,1,sgolayfilt,p=p,n=n,m=m))

## Boucle pour effectuer plusieurs FDA (reduire impact tirage aleatoire)

# initialisation vecteur de % de bons classments par DV
perok=vector(mode='numeric',length=ncmax)
# creation matrice de % de mauvais classements par clone
mc=matrix(nrow = ncmax,ncol = c)

## definition des matrices de resultat final
# creation matrice des perok finale
perok_final=matrix(nrow = repet, ncol = ncmax)
perok_finalm=matrix(nrow = repet, ncol = ncmax)
perok_finalm0=matrix(nrow = repet, ncol = ncmax)
# perok_final0=matrix(nrow = repet, ncol = ncmax)
perok_finalPLSDA=matrix(nrow = repet, ncol = ncmax)

# creation de la matrice des dv et perok maximaux
maxi_final=matrix(nrow= repet, ncol = 2)
# creation de la matrice de % de mauvais classements
mc_final=matrix(nrow= repet, ncol = length(levels(class)))
# creation d'un matrice cubique pour enregistrer les tables de contingence
t_final=array(dim=c(c,c,repet))
# noms des colonnes et lignes
colnames(t_final)=c(basename(levels(class)))
rownames(t_final)=c(basename(levels(class)))
colnames(maxi_final)= c("maxi.id","perok max")
colnames(mc_final)= c(basename(levels(class)))

for(j in 1:repet) {

  # creation des jeux d'apprentissage et validation
  flds <- createFolds(1:sum(iok), k = k)
  pred=as.data.frame(matrix(nrow = sum(iok), ncol = ncmax))  # As many sample as the whole set
  predm=as.data.frame(matrix(nrow = sum(iok), ncol = ncmax))  # As many sample as the whole set
  predm0=as.data.frame(matrix(nrow = sum(iok), ncol = ncmax))
  # pred0=as.data.frame(matrix(nrow = sum(iok), ncol = ncmax))
  predPLSDA=as.data.frame(matrix(nrow = sum(iok), ncol = ncmax))


  predpost=array(0,dim=c(sum(iok),c,ncmax))
  predmah=array(0,dim=c(sum(iok),c,ncmax))   # As many sample as the whole set


  # Boucle sur les groupes de CV
  for (i in 1:k) {
    id_val=sort(unlist(flds[i]))
    sp_val=sp[id_val,]
    class_val=class[id_val]
    sp_cal=sp[-id_val,]
    class_cal=class[-id_val]

    # PLSDA and application to have loadings and scores
    rplsda=caret::plsda(sp_cal, class_cal,ncomp=ncmax)
    sc_cal=rplsda$scores
    # sp_cal_c=scale(sp_cal,center=rplsda$Xmeans,scale = F)
    # sc_cal=sp_cal_c%*%rplsda$projection
    sp_val_c=scale(sp_val,center=rplsda$Xmeans,scale = F)
    sc_val=sp_val_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

    # rpca=prcomp(sp_cal)
    # spca_cal=rpca$x
    # spca_val=predict(rpca,sp_val)

    for (ii in 2:ncmax) {
      # Apprentissage
      rlda=lda(sc_cal[,1:ii], class_cal,tol=1.0e-5)
      # Validation
      pred[id_val,ii]=predict(rlda ,sc_val[,1:ii])$class
      predpost[id_val,,ii]=predict(rlda ,sc_val[,1:ii])$posterior
      # browser()
      predm[id_val,ii]=SIGNE_maha(rlda,sc_cal[,1:ii],class_cal,sc_val[,1:ii])$class
      predmah[id_val,,ii]=SIGNE_maha(rlda,sc_cal[,1:ii],class_cal,sc_val[,1:ii])$posterior

      predm0[id_val,ii]=SIGNE_maha0(sc_cal[,1:ii], class_cal, sc_val[,1:ii])$class

      # pred0[id_val,ii]=SIGNE_maha0(spca_cal[,1:ii], class_cal, spca_val[,1:ii])$class

      predPLSDA[id_val,ii]=predict(rplsda,sp_val,ncomp=ii)

      # p=predict(rlda ,sc_val[,1:ii])
      # predx[id_val,ii]=predict(rlda ,sc_val[,1:ii])$x
      # predp[id_val,ii]=predict(rlda ,sc_val[,1:ii])$posterior

    }
  }
    # Table de contingence


    ts=lapply(as.list(pred), class, FUN = table)
    tsm=lapply(as.list(predm), class, FUN = table)
    tsm0=lapply(as.list(predm0), class, FUN = table)
    # ts0=lapply(as.list(pred0), class, FUN = table)
    tsPLSDA=lapply(as.list(predPLSDA), class, FUN = table)

    # Matrice mauvais classements par clone
    #mc[i,]=(rowSums((as.matrix(t)))-diag(as.matrix(t)))/rowSums((as.matrix(t)))
    Sumpred=lapply(ts, FUN = rowSums)
    diags=lapply(ts, FUN = diag)
    diagsm=lapply(tsm, FUN = diag)
    diagsm0=lapply(tsm0, FUN = diag)
    # diags0=lapply(ts0, FUN = diag)
    diagsPLSDA=lapply(tsPLSDA, FUN = diag)

    # Pourcentage de bien classes
    perok=100*unlist(lapply(diags, FUN = sum))/length(class)
    perokm=100*unlist(lapply(diagsm, FUN = sum))/length(class)
    perokm0=100*unlist(lapply(diagsm0, FUN = sum))/length(class)
    # perok0=100*unlist(lapply(diags0, FUN = sum))/length(class)
    perokPLSDA=100*unlist(lapply(diagsPLSDA, FUN = sum))/length(class)

    # pred[id_val,]=predict(rplsda ,sp_val,ncomp = 1:ncmax)


  # # Table de contingence
   ts=lapply(as.list(pred), class, FUN = table)
  # # Matrice mauvais classements par clone
   #mc[i,]=(rowSums((as.matrix(t)))-diag(as.matrix(t)))/rowSums((as.matrix(t)))
   #Sumpred=lapply(ts, FUN = rowSums)
   #diags=lapply(ts, FUN = diag)
  # # Pourcentage de bien classes
   #perok=100*unlist(lapply(diags, FUN = sum))/length(class)

  maxi=max(perok)
  maxi.id=which.max(perok)


  ## Enregistrement des matrices de resultat final
  # remplissage de la matrice des perok finale
  perok_final[j,]=perok
  perok_finalm[j,]=perokm
  perok_finalm0[j,]=perokm0
  #perok_final0[j,]=perok0
  perok_finalPLSDA[j,]=perokPLSDA


  # remplissage de la dv max et de son % de bon classements globaux
  maxi_final[j,1]=maxi.id
  maxi_final[j,2]=maxi
  # remplissage de la matrice de mauvais classements par clone
   mc_final[j,]=mc[maxi.id,]
   #t_final[,,j]=ts
}

# affichage matrices de resultat final
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
# # Tracage des spectres
# # matplot(t(sp),pch='.')
# # Tracage de l'evolution des perok en fonction du nombre de DV utilisees
plot(colMeans(perok_finalm0))
points(colMeans(perok_finalm),col="red")
points(colMeans(perok_final),col="blue")
#points(colMeans(perok_final0),col="yellow")
points(colMeans(perok_finalPLSDA),col="green")
legend(ncmax*2/3,15,legend=c("Maha on PLSDA scores", "Maha on LDA+PLSDA scores","Predict on LDA+PLSDA scores","Predict PLSDA"),
       col=c("black","red","blue","green"), lty=1, cex=0.8)

# # # Tracage des moyennes de mauvais classements par clone
# # plot(colMeans(mc_final))
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
############ END ############

## Global model
rplsdag=caret::plsda(sp, class,ncomp=ncmax)
scg=rplsdag$scores
rldag=list()
for (i in 2:ncmax) {
rldag[[i]]=lda(scg[,1:i], class,tol=1.0e-5)
}

## PLot
 #plot(sc_cal[,1],sc_cal[,2], col=class_cal,pch=substr(rownames(sp_cal),3,3))
 #plot(sc_cal[,1],sc_cal[,2], pch=''
 #plot(sc_cal[,1],sc_cal[,2], col=dat_cal)
 #text(sc_val[,1],sc_val[,2], labels=substr(rownames(sp_val),11,12), col = "red")
 #plot(sc_cal[,1],sc_cal[,2], col=as.factor(v))

 #plot(scg[,1],scg[,2], col=class)
 #points(sc2[,1],sc2[,2],pch=3)
