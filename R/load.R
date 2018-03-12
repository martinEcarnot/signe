#' Function to import ASD spectra
#'
#' @param d character
#'
#' @return matrix
#' @author Martin Ecarnot
#' @export
#' 
"15:38"

load <- function (d) {
library(asdreader)
l=Sys.glob(file.path(d, "*.asd"))

sp=matrix(nrow=length(l),ncol=2151)	
spt=list()

for (i in 1:length(l)) {
  # sp[i,]=get_spectra(l[i], type = "reflectance")
  # spt[i]=asd_read(l[i])  # Fonction pour lire les fichiers .asd
  sp1=asd_read(l[i])  # Fonction pour lire les fichiers .asd
  sp[i,]=t(sp1$spectrum/sp1$reference)

}
return(sp)
}
