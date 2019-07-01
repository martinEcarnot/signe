## COMPRESSION DES SPECTRES EN UNE MATRICE ##

# nettoyage des variables de l'espace de travail
rm(list = ls())

user="ME"
if (user=="NL") {source('C:/Users/Noémie/Desktop/SFE/Script_R/SIGNE_load.R')} else {source("./Script_R_2018/SIGNE_load.R")}


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
s=9
#creation de la matrice referencant les dossiers
dates=vector(mode = "logical", length = s)

# dossiers ? assembler
dates[1]="20180619P"
dates[2]="20180627P"
#dates[3]="20180704P"
dates[3]="20180710P"
#dates[5]="20180710E"
dates[4]="20180619T"
dates[5]="20180619N"
dates[6]="20180627T"
dates[7]="20180627N"
dates[8]="20180709T"
dates[9]="20180709N"

# #cr?ation de la matrice qui rassemble tous les dossiers:
globalmatrix=matrix(nrow=(s*180), ncol=2072)

for (k in 1:s)
{
  if (user=="NL") {
    em=paste("C:\\Users\\Noémie\\Desktop\\SFE\\Spectres_SIGNE\\Fichiers_ASD\\PTN\\",dates[k], sep="")} else {
      em=paste("/home/ecarnot/Documents/INRA/Projets/SIGNE/2018/GrauduRoi/PTN/",dates[k], sep="")}

  w=SIGNE_load(em)
  rownames(w)=paste(dates[k],rownames(w))
  globalmatrix=rbind(globalmatrix, w)
}

globalmatrix=globalmatrix[complete.cases(globalmatrix),]

if (user=="NL") {brb="C:\\Users\\Noémie\\Desktop\\SFE\\Resultats\\PTN1\\"} else{ brb="/home/ecarnot/Documents/INRA/Projets/SIGNE/2018/GrauduRoi/PTN/"}

save(globalmatrix, file=paste(brb,"globalmatrix2",sep=""))

write.table(globalmatrix, file=paste(brb,"globalmatrix2.csv",sep=""),sep=";", quote=FALSE)


### END ###
