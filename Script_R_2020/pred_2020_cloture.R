source('Script_R_2020/SIGNE_load.R')
library(rnirs)
library(caret)
library(MASS)
source('Script_R_2020/SIGNE_maha0.R')
source('Script_R_2020/sp2dfclo.R')


load("Script_R_2020/globalmatrix")
sp=globalmatrix
site=as.factor(substr(rownames(sp),18,18))
sp=sp[site=="G",]

# Chargement des données de l'été 2020
sp1=SIGNE_load("/home/ecarnot/Documents/INRA/Projets/SIGNE/2020/ASD_ete_2020/SIGNE_20200611/")
nsp1=nrow(sp1)
rownames(sp1)=gsub("524","525",rownames(sp1))
sp2=SIGNE_load("/home/ecarnot/Documents/INRA/Projets/SIGNE/2020/ASD_ete_2020/SIGNE_20200721/")
nsp2=nrow(sp2)
spn=rbind(sp1,sp2)

##Rajoute le nom du cepage au nom de la ligne
titre=rownames(spn)
yC= substr(titre,1,3)== "015" |  substr(titre,1,3)== "169" |  substr(titre,1,3)== "685"
yG= substr(titre,1,3)== "222" | substr(titre,1,3)== "509" |  substr(titre,1,3)== "787"
yS= substr(titre,1,3)== "471" |  substr(titre,1,3)== "525" |  substr(titre,1,3)== "747" |  substr(titre,1,3)== "877"
cep=character(nrow(spn))
cep[yC]="C"
cep[yG]="G"
cep[yS]="S"
rownames(spn)=paste0(c(rep("20200611",nsp1),rep("20200721",nsp2)),cep, " ",rownames(spn)," ","G")


sptot=rbind(sp,spn)
class=as.factor(substr(rownames(sptot),9,9)) #9,9
classclo=as.factor(substr(rownames(sptot),9,13)) #9,13
sptot=sp2dfclo(sptot,class,classclo)
colnames(sptot)[1:2]=c("cepage","clone")


# Pre
p=rbind(list('adj',c(552,1352)),list('red',c(1,1,2)),list('snv',''),list('sder',c(2,3,101))) # pour cepages
# p=rbind(list('adj',c(552,1352)),list('red',c(1,1,1)),list('snv',''))# pour C

# p=rbind(list('adj',c(552,1352)),list('red',c(1,1,2)),list('snv',''))
sptot$xp=pre(sptot$x,p)


sp=sptot[1:nrow(sp),]
spn=sptot[-(1:nrow(sp)),]

# spnr=spn[1:180,]
# spnr$xp[1:180,]=spnr$xp[1:180,]-t(replicate(nrow(sp1), colMeans(spn$xp[1:180,]))) +t(replicate(nrow(sp1), colMeans(sp$xp)))

# # PLSDA locale
# fm <- lwplsdalm(sp$xp, sp$cepage  , spn$xp, spn$cepage, ncomp = 25,k = nrow(sp),ncompdis = 25,h=1)
#
# fm <- rnirs::plsda(sp$xp, sp$cepage  , spn$xp, spn$cepage, ncomp = 25,da=dalm)
# err(fm, ~ ncomp)
#
# # CV en PLSDA locales
# nrep=3
# ncomp=25
# K=seq(30,nrow(sp), length.out = 5)
# h=c(0.5,1,2,5,10) #seq(0.5,10, length.out = 5)
# seg=segmcvkfold(nrow(sp),K=4,nrep=nrep)
# fmpls=fitcv(sp$xp, sp$cepage,fun=lwplsda,segm=seg, ncomp=ncomp, k=K,h=h, ncompdis=25)
# err(fmpls, ~ ncomp+rep)
#
# #PLSDA locale par cepage
# #Cabernet
spC=droplevels(sp[sp$cepage=="C",])
spnC=droplevels(spn[spn$cepage=="C",])
# fmC <- rnirs::plsda(spC$xp, spC$clone  , spnC$xp, spnC$clone, ncomp = 25,da=dalm)
# fmC <- lwplsdalm(spC$xp, spC$clone  , spnC$xp, spnC$clone, ncomp = 20,k = c(nrow(spC)),ncompdis = 20,h=c(5))
#
spG=droplevels(sp[sp$cepage=="G",])
spnG=droplevels(spn[spn$cepage=="G",])
# fmG <- rnirs::plsda(spG$xp, spG$clone  , spnG$xp, spnG$clone, ncomp = 25,da=dadis)
# fmG <- lwplsdalm(spG$xp, spG$clone  , spnG$xp, spnG$clone, ncomp = 25,k = c(nrow(spG)),ncompdis = 25,h=5)
#
spS=droplevels(sp[sp$cepage=="S",])
spnS=droplevels(spn[spn$cepage=="S",])
# fmS <- lwplsdalm(spS$xp, spS$clone  , spnS$xp, spnS$clone, ncomp = 25,k = c(nrow(spG)),ncompdis = 25,h=5)
#
# err(fmC, ~ ncomp)

# Reprise le l'ancien script PLSDA
spcal=rbind(spG,spnG[1:54,]) #spG[1:54,]  #
spval=spnG[55:108,] # spnG[1:54,]
classcal=spcal$clone
classval=spval$clone
ncmax=25

predmF=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))
distances=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))

rplsda=caret::plsda(spcal$xp, classcal,ncomp=ncmax)
sccal=rplsda$scores
spval_c=scale(spval$xp,center=rplsda$Xmeans,scale = F)
scval=spval_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

for (ii in 2:ncmax) {
  predmF[,ii]=SIGNE_maha0(sccal[,1:ii], classcal, scval[,1:ii])$class
  distances[,ii]=SIGNE_maha0(sccal[,1:ii], classcal, scval[,1:ii])$dist
}

tsm=lapply(as.list(predmF), classval, FUN = table)
diagsm=lapply(tsm, FUN = diag)
perokm =100*unlist(lapply(diagsm, FUN = sum))/length(classval)
print(perokm)

seuil=seq(1,1.8,0.05)
effseuil=matrix(nrow = length(seuil),ncol=2)
bestLV=which(perokm==max(perokm))[1]
d=distances[[bestLV]]
min1=t(apply(d,1,sort))
rapdista=min1[,2]/min1[,1]
for (i in 1:length(seuil)) {
iok=rapdista>seuil[i]
t=table(predmF[[bestLV]][iok],classval[iok])
effseuil[i,1]=sum(t)*100/nrow(d)
effseuil[i,2]=sum(diag(t))*100/sum(t)
}
plot(effseuil[,1],effseuil[,2])

#
# # ACP
# fm <- pca(spcal$xp, spval$xpr, ncomp = 10)
# Tr <- fm$Tr
# Tu <- fm$Tu
#
# T <- rbind(Tr, Tu)
# row.names(T) <- 1:nrow(T)
# group <- c(rep("Reference", nrow(Tr)), rep("Unknown", nrow(Tu)))
# plotxy(T[,c(3,4)], group = group, pch = 16, zeroes = TRUE)
#
#
#
#
# spn$xpr=spn$xp-t(replicate(nrow(spn), colMeans(spn$xp))) +t(replicate(nrow(spn), colMeans(sp$xp)))

