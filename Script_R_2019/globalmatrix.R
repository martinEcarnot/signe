## COMPRESSION DES SPECTRES EN UNE MATRICE ##

# nettoyage des variables de l'espace de travail
rm(list = ls())

source('Script_R_2019/SIGNE_load.R')
#source('~/Documents/NICOLAS/Stage de fin annee au Grau du roi~/Documents/NICOLAS/Stage de fin annee au Grau du roi')

#creation de la matrice referencant les dossiers
dates=vector(mode = "logical", length = 30)

# dossiers ? assembler
# dates[1]="20180619N"
# dates[2]="20180627N"
# # dates[4]="20180619T"
# dates[5]="20180627P"
# # dates[6]="20180627T"
# # dates[7]="20180709N"
# # dates[8]="20180709T"
# dates[9]="20180710P"
# # dates[10]="20180816T"
# dates[11]="20180704P"
# # dates[12]="20180816N"
# # dates[13]="20190613P"
# # dates[14]="20190613T"
# # dates[15]="20190613N"
# # dates[10]="20180709N"
# # dates[15]="20180817T"
# dates[16]="20180619P"
# dates[17]="20180817N"
# dates[18]="20170524P"
# dates[19]="20170529P"
# dates[20]="20170606P"
# dates[21]="20170612P"
# dates[22]="20170619P"
dates[23]="20170626P"
dates[24]="20170703P"
dates[25]="20170710P"
dates[26]="20170717P"
dates[27]="20170724P"
dates[28]="20170731P"

dates=dates[dates!=F]
s=length(dates)
# #cr?ation de la matrice qui rassemble tous les dossiers:
globalmatrixN1=matrix(ncol=2072)

for (k in 1:s)
  {
  # em=paste("~/Documents/INRA/Projets/SIGNE/2018/GrauduRoi/Parcelles/",dates[k], sep="")
  em=paste("~/Documents/INRA/Projets/SIGNE/2018/GrauduRoi/Parcelles/",dates[k], sep="")
  w=SIGNE_load(em)
  rownames(w)=paste(dates[k],rownames(w))
  globalmatrixN1=rbind(globalmatrixN1, w)
}

globalmatrixN1=globalmatrixN1[complete.cases(globalmatrixN1),]


brb6="~/Documents/INRA/Projets/SIGNE/2018/GrauduRoi/PTN/"
save(globalmatrixN1, file=paste(brb6,"globalmatrixN1",sep=""))

write.table(globalmatrixN1, file=paste(brb6,"globalmatrixN1.csv",sep=""),sep=";", quote=FALSE)


### END ###

