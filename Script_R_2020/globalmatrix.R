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



# ###Ajout du cepage au nom de la ligne. Retire aussi les cépages/clones qui sont pas les 3/10 notres.
# ##Filtre en fonction du cepage
# titre=rownames(globalmatrix)
#
# ya=(substr(titre,11,13))== "015"
# yb=(substr(titre,11,13))== "169"
# yc=(substr(titre,11,13))== "685"
#
# yd=(substr(titre,11,13))== "222"
# ye=(substr(titre,11,13))== "509"
# yf=(substr(titre,11,13))== "787"
#
# yg=(substr(titre,11,13))== "471"
# yh=(substr(titre,11,13))== "525"
# yi=(substr(titre,11,13))== "747"
# yj=(substr(titre,11,13))== "877"
#
# ##Separation en 3 matrices
# spa=globalmatrix[ya,]
# spb=globalmatrix[yb,]
# spc=globalmatrix[yc,]
# spd=globalmatrix[yd,]
# spe=globalmatrix[ye,]
# spf=globalmatrix[yf,]
# spg=globalmatrix[yg,]
# sph=globalmatrix[yh,]
# spi=globalmatrix[yi,]
# spj=globalmatrix[yj,]
#
# ##Rajoute le nom du cepage au nom de la ligne à la place de P/T/N
# substr(rownames(spa),9,9)="C"
# substr(rownames(spb),9,9)="C"
# substr(rownames(spc),9,9)="C"
#
# substr(rownames(spd),9,9)="G"
# substr(rownames(spe),9,9)="G"
# substr(rownames(spf),9,9)="G"
#
# substr(rownames(spg),9,9)="S"
# substr(rownames(sph),9,9)="S"
# substr(rownames(spi),9,9)="S"
# substr(rownames(spj),9,9)="S"
#
# ##Recombine les 3 matrices pour reformer sp
# globalmatrix=rbind(spa,spb,spc,spd,spe,spf,spg,sph,spi,spj)

for (i in 1:length(globalmatrix[,1])){
  if (substr(rownames(globalmatrix)[i],11,13)=="015" | substr(rownames(globalmatrix)[i],11,13)=="169" | substr(rownames(globalmatrix)[i],11,13)=="685"){
    substr(rownames(globalmatrix)[i],9,9)="C"
  }
  else if  (substr(rownames(globalmatrix)[i],11,13)=="222" | substr(rownames(globalmatrix)[i],11,13)=="509" | substr(rownames(globalmatrix)[i],11,13)=="787"){
    substr(rownames(globalmatrix)[i],9,9)="G"
  }
  else if  (substr(rownames(globalmatrix)[i],11,13)=="471" | substr(rownames(globalmatrix)[i],11,13)=="525" | substr(rownames(globalmatrix)[i],11,13)=="747" | substr(rownames(globalmatrix)[i],11,13)=="877"){
    substr(rownames(globalmatrix)[i],9,9)="S"
  }
  else {
    substr(rownames(globalmatrix)[i],9,9)="A"
  }
}

iout=which(substr(rownames(globalmatrix),9,9)=="A")
globalmatrix <- globalmatrix[-iout,]


print(globalmatrix)
brb="C:/Users/avitvale/Documents/Test/"
save(globalmatrix, file=paste(brb,"globalmatrix",sep=""))
print(length(globalmatrix))
# write.table(globalmatrix, file=paste(brb,"globalmatrix.csv",sep=""),sep=";", quote=FALSE)

### END ###

