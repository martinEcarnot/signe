## COMPRESSION DES SPECTRES EN UNE MATRICE ##

# nettoyage des variables de l'espace de travail
rm(list = ls())

source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')
source('C:/Users/avitvale/Documents/Script_R/asd_read.R')



# #création de la matrice qui rassemble tous les dossiers:
globalmatrix=matrix(ncol=2072)

#print(globalmatrix)
setwd("C:/Users/avitvale/Documents/SIGNE")
for (a in 2017:2019)
  {
  setwd("C:/Users/avitvale/Documents/SIGNE")
  setwd(as.character(a))
  for (n in 1:length(dir()))
    {

    if ((dir()[n])=="Aude" | (dir()[n])=="Beaujolais" | (dir()[n])=="GDRautre" |(dir()[n])=="GDR96")
    {
      setwd(as.character(dir()[n]))
      for (m in 1:length(dir()))
      {
        if(substr(dir()[m],9,9)=="P")
        {
          em=paste(dir()[m], sep="")
          w=SIGNE_load(em)
          rownames(w)=paste(dir()[m],rownames(w))
          globalmatrix=rbind(globalmatrix, w)
        }
      }
      setwd("../")
    }





    else if(substr(dir()[n],9,9)=="P")
      {
      em=paste(dir()[n], sep="")
      print(em)
      w=SIGNE_load(em)
      rownames(w)=paste(dir()[n],rownames(w))
      globalmatrix=rbind(globalmatrix, w)
    }
  }
}

setwd("C:/Users/avitvale/Documents/signeG")


#print(globalmatrix)
globalmatrix=globalmatrix[complete.cases(globalmatrix),]

print(globalmatrix)
brb="C:/Users/avitvale/Documents/Test/"
save(globalmatrix, file=paste(brb,"globalmatrix",sep=""))
print(length(globalmatrix))
# write.table(globalmatrix, file=paste(brb,"globalmatrix.csv",sep=""),sep=";", quote=FALSE)


### END ###
unique(nchar(rownames(globalmatrix)))
unique(rownames(globalmatrix))
nrow(globalmatrix)
5875/36
