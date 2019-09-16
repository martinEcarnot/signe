## COMPRESSION DES SPECTRES EN UNE MATRICE ##

# nettoyage des variables de l'espace de travail
rm(list = ls())

source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')

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
s=7#creation de la matrice referencant les dossiers
dates=vector(mode = "logical", length = s)

# dossiers ? assembler
dates[1]="20180619P"
dates[2]="20180627P"
dates[3]="20180704P"
dates[4]="20180709P"
dates[5]="20180710P"
dates[6]="20180816P"
dates[7]="20180817P"
# dates[8]="20170524P"
# dates[9]="20170529P"
# dates[10]="20170606P"
# dates[11]="20170612P"
# dates[12]="20170619P"
# dates[13]="20170626P"
# dates[14]="20170703P"
# dates[15]="20170710P"
# dates[16]="20170717P"
# dates[17]="20170724P"
# dates[18]="20170731P"


# #cr?ation de la matrice qui rassemble tous les dossiers:
globalmatrix=matrix(ncol=2072)

#print(globalmatrix)
for (k in 1:s)
  {
  em=paste("C:/Users/avitvale/Documents/SIGNE/2018/GDR2018/",dates[k], sep="")
  w=SIGNE_load(em)
  rownames(w)=paste(dates[k],rownames(w))
  globalmatrix=rbind(globalmatrix, w)
}

#print(globalmatrix)
globalmatrix=globalmatrix[complete.cases(globalmatrix),]

print(globalmatrix)
brb="C:/Users/avitvale/Documents/Test/"
save(globalmatrix, file=paste(brb,"globalmatrix",sep=""))
print(globalmatrix)
# write.table(globalmatrix, file=paste(brb,"globalmatrix.csv",sep=""),sep=";", quote=FALSE)


### END ###
