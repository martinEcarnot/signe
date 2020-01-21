# Spectres SAADA 2019
# Essai de discri de genotypes

library(rnirs)
library(dplyr)
library(sampling)
source('Script_R_2020/SIGNE_maha0.R')
source("~/Documents/INRA/R/pre.R")


d='/home/ecarnot/Documents/INRA/melanges/2019'
sp=read.table(file.path(d,'saada2019_spectres.txt'), header=T, row.names=1, sep=';', dec='.')

# Enleve ce qui n'est pas de SAADA
sp=sp[-seq(1,29),]
sp=sp[grep("es",rownames(sp), invert=T),]

md=read.table(file.path(d,'SAADA prepa semis_code champ19.csv'),header=T,  sep=';')

# On duplique les infos de md sur sp
colnames(sp)=1e7/as.numeric(gsub("X","",colnames(sp)))
sp=data.frame(x=I(as.matrix(sp)))
class(sp$x)="matrix"
sp=cbind(code.champ=as.numeric(sub("\\-.*", "", rownames(sp))), plante=sub(".*-", "", rownames(sp)),sp)
sp$plante=as.factor(trimws(gsub("\\(1\\)","",sp$plante)))

sp=cbind(sp,left_join(sp,md,by ="code.champ")[,4:9])

mono=sp[grep("monogénotype",sp$culture),]




# Pré
# p=rbind(list('red',c(1300, 1, 2)),list('snv',''),list('sder',c(2,3,11)))
p=rbind(list('red',c(100,1,2)),list('snv',''),list('sder',c(2,3,41)))
mono$xp=pre(mono$x,p)
comp=15
rpca=pca(mono$xp,ncomp=comp)

z=rpca$explvarx
barplot(100 * z$pvar, names.arg = paste("comp", z$ncomp),ylab = "Pct. of variance explained")
plotxy(rpca$Tr[, c(1,2)],group=mono$rep)




# Discri
ngeno=nlevels(mono$Genotype1)
segm=segmFact(mono, var=list("Genotype1","plante"),nel=list(ngeno,rep(4,ngeno)), nrep=10 )
fm <- fitcv(mono$xp, mono$Genotype1,fun = plsdalm,  ncomp = 15,  segm = segm)

z <- err(fm, ~ ncomp)
plotmse(z, nam = "errp")
min(z$errp)

# Avec 2 geno
deux=levels(mono$Genotype1)[sample(65,2)]
ideux=mono$Genotype1 %in% deux
mono3=mono[itrois,]
mono3 <- droplevels(mono3)
ngeno=nlevels(mono3$Genotype1)
segm=segmFact(mono3, var=list("Genotype1","plante"),nel=list(ngeno,rep(4,ngeno)), nrep=10 )
fm <- fitcv(mono3$xp, mono3$Genotype1,fun = plsdalm,  ncomp = 15,  segm = segm)
z <- err(fm, ~ ncomp)
plotmse(z, nam = "errp")
i6=fm$y$ncomp==5
table(fm$y$x1[i6],fm$fit$x1[i6])
