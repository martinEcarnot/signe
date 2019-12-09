#' Function to import ASD spectra
#'
#' @param d character
#'
#' @return matrix
#' @author Martin Ecarnot
#' @export
#'
#'
source('Script_R_2020/asd_read.R')

SIGNE_loadMLO <- function (d) {
library(asdreader)
#d=choose.dir()
#d="C:\\Users\\Noemie\\Desktop\\SFE\\Pour_Nathan\\Documents_finaux\\test_spectres\\"
l=Sys.glob(file.path(d, "*.asd"))
l=sort(l)

sp=matrix(nrow=length(l),ncol=751)
spt=list()

for (i in 1:length(l)) {
 #sp1=get_spectra(l[i], type = "reflectance") # Fonction pour lire les fichiers .asd
  sp1=asd_read(l[i])
  sp1=sp1$spectrum/sp1$reference
  sp[i,]=sp1
}


l1=basename(l)
l1=gsub(".asd","",l1)
l1=gsub("000","-",l1)
# Look for wavelength
md=get_metadata(l[i])
wl=seq(from=md$ch1_wavel, to=md$ch1_wavel+md$channels-1)

colnames(sp)=wl
row.names(sp)=l1

# Create class file
 clas=substr(l1,1,3)

# uclas=unique(clas)
# for (i in 1:length(uclas)) {
#   iok=which(clas==uclas[i])
#   clas[iok]=i
# }

 clas=data.frame(clas)
 row.names(clas)=l1
 colnames(clas)="clone"
sp=sp[,50:751] # Remove noisy part of the spectra

return(sp)
# setwd(d)
# write.table(sp,file=paste(basename(d),"_sp.csv",sep=""),sep=';', quote=FALSE,col.names = NA)
# write.table(clas,file=paste(basename(d),"_clon.csv",sep=""),sep=';', quote=FALSE,col.names = NA)

}

