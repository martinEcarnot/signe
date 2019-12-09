## COMPRESSION DES SPECTRES EN UNE MATRICE ##

# nettoyage des variables de l'espace de travail
rm(list = ls())

source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')
source('C:/Users/avitvale/Documents/Script_R/asd_read.R')
source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')


# #création de la matrice qui rassemble tous les dossiers:
globalmatrix=matrix(ncol=2072)

#print(globalmatrix)
setwd("C:/Users/avitvale/Documents/SIGNE")
for (a in 2017:2019){
  setwd("C:/Users/avitvale/Documents/SIGNE")
  setwd(as.character(a))
  for (n in 1:length(dir())){
    if ((dir()[n])=="Aude" | (dir()[n])=="Beaujolais" | (dir()[n])=="GDRautre" |(dir()[n])=="GDR96"){
      lieu=substr(dir()[n],1,1)
      if ((dir()[n])=="GDRautre"){
        lieu='g'
      }
      setwd(as.character(dir()[n]))
      for (m in 1:length(dir())){
        if(substr(dir()[m],9,9)=="P"){
          em=paste(dir()[m], sep="")
          w=SIGNE_load(em)
          rownames(w)=paste(dir()[m],rownames(w),lieu)
#          print(rownames(w))
          globalmatrix=rbind(globalmatrix, w)
        }
      }
      setwd("../")
    }
  }
}

setwd("C:/Users/avitvale/Documents/signeG")


#print(globalmatrix)
globalmatrix=globalmatrix[complete.cases(globalmatrix),]
str(globalmatrix)

#Filtrage des noms non standardisés
iout=which(nchar(rownames(globalmatrix))>18)


# rownames(globalmatrix[iout,])
# unique(rownames(globalmatrix[which(substr(rownames(globalmatrix),1,8)=="20170612"),]))



globalmatrix <- globalmatrix[-iout,]

## Filtrage des spectres aberrants
globalmatrix=globalmatrix[globalmatrix[,500]>0.6,]
globalmatrix=globalmatrix[globalmatrix[,1]<0.2,]
globalmatrix=globalmatrix[globalmatrix[,2000]<0.25,]


## Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
globalmatrix=adj_asd(globalmatrix,c(602,1402))

## Reduction des variables (extremites bruitees)
globalmatrix=globalmatrix[,seq(51,ncol(globalmatrix)-30,1)]



###Ajout du cepage au nom de la ligne.
##Filtre en fonction du cepage
titre=rownames(globalmatrix)

ya=(substr(titre,11,13))== "015" |  (substr(titre,11,13))== "169" |  (substr(titre,11,13))== "685"
yb=(substr(titre,11,13))== "222" |  (substr(titre,11,13))== "509" |  (substr(titre,11,13))== "787"
yc=(substr(titre,11,13))== "471" |  (substr(titre,11,13))== "525" |  (substr(titre,11,13))== "747" |  (substr(titre,11,13))== "877"

##Separation en 3 matrices
spa=globalmatrix[ya,]
spb=globalmatrix[yb,]
spc=globalmatrix[yc,]

##Rajoute le nom du cepage au nom de la ligne à la place de P/T/N
substr(rownames(spa),9,9)="C"
substr(rownames(spb),9,9)="G"
substr(rownames(spc),9,9)="S"
#substr(rownames(spa),9,9)="1"
#substr(rownames(spb),9,9)="2"
#substr(rownames(spc),9,9)="3"

##Recombine les 3 matrices pour reformer sp
globalmatrix=rbind(spa,spb,spc)


print(globalmatrix)
brb="C:/Users/avitvale/Documents/Test/"
save(globalmatrix, file=paste(brb,"globalmatrix",sep=""))
print(length(globalmatrix))
# write.table(globalmatrix, file=paste(brb,"globalmatrix.csv",sep=""),sep=";", quote=FALSE)

### END ###

