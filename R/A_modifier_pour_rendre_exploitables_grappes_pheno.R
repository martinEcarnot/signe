library(rnirs)

#### Globalmatrix modifié et adapté à la situation.
## COMPRESSION DES SPECTRES EN UNE MATRICE ##

# nettoyage des variables de l'espace de travail
rm(list = ls())

source('C:/Users/avitvale/Documents/signeG/Script_R_2020/adj_asd.R')
source('C:/Users/avitvale/Documents/signeG/Script_R_2020/asd_read.R')
source('C:/Users/avitvale/Documents/signeG/Script_R_2020/SIGNE_load_modif.R')

# source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')       ##### N'arrive pas à trouver le chemin "normal (dans signeG), comprends pas pourquoi
# source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')
# source('C:/Users/avitvale/Documents/Script_R/asd_read.R')


# #création de la matrice qui rassemble tous les dossiers:
globalmatrix=matrix(ncol=2072)

setwd("C:/Users/avitvale/Documents/SIGNE/Prediction_pheno/Grappes")
  for (n in 1:length(dir())){
      lieu="B"
      for (m in 1:length(dir())){
          em=paste(dir()[m], sep="")
          w=SIGNE_load_modif(em)
          rownames(w)=paste(dir()[m],rownames(w),lieu) #On ajoute l'identifiant lieu au nom de chaque spectre
          globalmatrix=rbind(globalmatrix, w)
      }
}




globalmatrix=globalmatrix[complete.cases(globalmatrix),]
#str(globalmatrix)

rownames(globalmatrix)
paste((substr(rownames(globalmatrix),1,13)),(substr(rownames(globalmatrix),17,20)))

which(substr(rownames(globalmatrix),14,16)!="000")


# ## Filtrage des spectres aberrants        #Tous les spectres aberrants sont-ils filtrés par ce moyen ?
# globalmatrix2=globalmatrix[globalmatrix[,500]>0.6,] #Actuellement, 33 spectres jugés aberrants filtrés ici. (Dont 5 qui ne font pas partie de nos 10 clones étudiés..?)
# globalmatrix2=globalmatrix2[globalmatrix2[,1]<0.2,]
# globalmatrix2=globalmatrix2[globalmatrix2[,2000]<0.25,]


# ## Ajustement des sauts de detecteur (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
# globalmatrix=adj_asd(globalmatrix,c(602,1402))

## Reduction des variables (extremites bruitees)
globalmatrix=globalmatrix[,seq(51,ncol(globalmatrix)-30,1)]

## Filtrage des spectres aberrants
globalmatrix=globalmatrix[globalmatrix[,1]<0.6,]

#plotsp(globalmatrix[1:3593,])




globalmatrixconservatoiregrappes=globalmatrix



brb="C:/Users/avitvale/Documents/Test/"
save(globalmatrixconservatoiregrappes, file=paste(brb,"globalmatrixconservatoiregrappes",sep=""))
print(length(globalmatrixconservatoiregrappes))
# write.table(globalmatrix, file=paste(brb,"globalmatrix.csv",sep=""),sep=";", quote=FALSE)
setwd("C:/Users/avitvale/Documents/signeG")
### END ###


# gsub sert à remplacer un groupe de caractères par un autre dans une liste.
# gsub("A remplacer","Remplacement", Liste)



