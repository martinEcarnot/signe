## COMPRESSION DES SPECTRES EN UNE MATRICE ##

# nettoyage des variables de l'espace de travail
rm(list = ls())

source('Script_R_2019/SIGNE_load.R')

# SIGNE_load(d)



#source('C:/Users/No?mie/Desktop/SFE/Script_R/asd_read.R')

# fonction load
#load <- function (d) {
#  library(asdreader)
  #d=choose.dir()
#  l=Sys.glob(file.path(d, "*.asd"))
#  l=sort(l)

#  sp=matrix(nrow=length(l),ncol=2151)
#  spt=list()

#  for (i in 1:length(l)) {
    #sp1=get_spectra(l[i], type = "reflectance") # Fonction pour lire les fichiers .asd
#    sp1=asd_read(l[i])
#    sp1=sp1$spectrum/sp1$reference
#    sp[i,]=sp1
#  }


#  l1=basename(l)
#  l1=gsub(".asd","",l1)
#  l1=gsub("000","-",l1)
#  # Look for wavelength
#  md=get_metadata(l[i])
#  wl=seq(from=md$ch1_wavel, to=md$ch1_wavel+md$channels-1)

#  colnames(sp)=wl
#  row.names(sp)=l1

  # Create class file
#  clas=substr(l1,1,3)

#  uclas=unique(clas)
#  for (i in 1:length(uclas)) {
#    iok=which(clas==uclas[i])
#    clas[iok]=i
#  }

#  clas=data.frame(clas)
#  row.names(clas)=l1
#  colnames(clas)="clone"
#  sp=sp[,50:2121] # Remove noisy part of the spectra

#  return(sp)
#}


#nombre de dossiers de spectres ? assembler:
s=12
#creation de la matrice referencant les dossiers
dates=vector(mode = "logical", length = s)

# dossiers ? assembler
dates[1]="20180619N"
dates[2]="20180627N"
dates[3]="20180619P"
dates[4]="20180619T"
dates[5]="20180627P"
dates[6]="20180627T"
dates[7]="20180709N"
dates[8]="20180709T"
dates[9]="20180709P"
# dates[10]="20180709N"
# dates[11]="20180710P"
dates[10]="20180816T"
dates[11]="20180816P"
dates[12]="20180816N"
# dates[15]="20180817T"
# dates[16]="20180817P"
# dates[17]="20180817N"
# dates[18]="20170524P"
# dates[19]="20170529P"
# dates[20]="20170606P"
# dates[21]="20170612P"
# dates[22]="20170619P"
# dates[23]="20170626P"
# dates[24]="20170703P"
# dates[25]="20170710P"
# dates[26]="20170717P"
# dates[27]="20170724P"
# dates[28]="20170731P"



# #cr?ation de la matrice qui rassemble tous les dossiers:
globalmatrixN1=matrix(nrow=(s*1326), ncol=2072)

for (k in 1:s)
  {
  em=paste("~/Documents/NICOLAS/Stage de fin annee au Grau du roi/Pour Noemie/SIGNE_2018_LAFOUGE/Spectres_SIGNE/Parcelles/Clones/",dates[k], sep="")
  w=SIGNE_load(em)
  rownames(w)=paste(dates[k],rownames(w))
  globalmatrixN1=rbind(globalmatrixN1, w)
}

globalmatrixN1=globalmatrixN1[complete.cases(globalmatrixN1),]


strbrb6="~/Documents/NICOLAS/Stage de fin annee au Grau du roi/"
save(globalmatrixN1, file=paste(brb6,"globalmatrixN1",sep=""))

write.table(globalmatrixN1, file=paste(brb6,"globalmatrixN1.csv",sep=""),sep=";", quote=FALSE)


### END ###

