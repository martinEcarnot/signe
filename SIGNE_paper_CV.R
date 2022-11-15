# Avoir le dataframe sp avec toutes les meta donnée OK

datclone=substr(rownames(sp),1,13)
ndc=length(unique(datclone))
# On créé un facteur souche qui groupe les 6 spectres de chaque souche
numsp=as.numeric(substr(rownames(sp),15,16))
souche=cut(numsp, breaks = c(0,6,12,18),labels=c("s1","s2","s3"))
feuill=cut(numsp, breaks = c(0,3,6,9,12,15,19),labels=c("f1","f2","f3","f4","f5","f6"))
sp$feuill=paste(datclone,souche,feuill,sep=" ")

feuilles=c("f1","f2","f3","f4","f5","f6")
dat=data.frame(matrix(NA, nrow = ncmax-1, ncol = 0))
an=NULL
typ=NULL
LV=NULL
for (i in 1:length(feuilles)) {
  idcal= which(!grepl(feuilles[i], sp$feuill))

  spcal= droplevels(sp[idcal,])
  spval=droplevels(sp[-idcal,])
  predmF=as.data.frame(matrix(nrow = nrow(spval), ncol = ncmax))
  distances=as.data.frame(matrix(nrow = nrow(spval), ncol = ncmax))
  rapdista=as.data.frame(matrix(nrow = nrow(spval), ncol = ncmax))
  perokm =as.data.frame(matrix(nrow = 1, ncol = ncmax))
  lost_rap=as.data.frame(matrix(nrow = 1, ncol = ncmax))

  rplsda=caret::plsda(spcal$x, spcal$clo,ncomp=ncmax)
  sccal=rplsda$scores
  spval_c=scale(spval$x,center=rplsda$Xmeans,scale = F)
  scval=spval_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

  for (ii in 2:ncmax) {
    predmF[,ii]=SIGNE_maha0(sccal[,1:ii], spcal$clo, scval[,1:ii])$class
    distances[,ii]=SIGNE_maha0(sccal[,1:ii], spcal$clo, scval[,1:ii])$dist
  }
  tsm=lapply(as.list(predmF), spval$clo, FUN = table)
  diagsm=lapply(tsm, FUN = diag)
  peroktt =100*unlist(lapply(diagsm, FUN = sum))/nrow(spval)

  dat=cbind(dat,peroktt[-1])

  # an=c(an,rep(annees[i],ncmax-1),rep(annees[i],ncmax-1))
  typ=c(typ,rep("All smpls",ncmax-1),rep("Unambig. smpls",ncmax-1))
  LV=c(LV,rep(1:(ncmax-1),2))
}
