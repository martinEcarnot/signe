# From matrix cal matrix, cal classes and val matrix, compute Mahalanobis dist from each group

SIGNE_maha0 = function (scal,class_cal,sval) {
ng=length(unique(class_cal))
nm=unique(class_cal)

if (is.null(nrow(sval))){
  mdist=matrix(nrow = 1,ncol = ng)
} else {
  mdist=matrix(nrow = nrow(sval),ncol = ng)
}
for (i in 1:ng) {
  igroup=which(class_cal==nm[i])
  if (is.null(nrow(sval))){
    cm=cov(as.matrix(scal[igroup]))
    center=colMeans(as.matrix(scal[igroup]))  # rlda$means[i,]%*%rlda$scaling
  } else {
    cm=cov(scal[igroup,]) # Create covariance matrix     # xc=x-matrix(data=1, nrow=54)%*%colMeans(x), cm=t(xc)%*%xc
    center=colMeans(scal[igroup,])  # rlda$means[i,]%*%rlda$scaling
  }
  # mdist[,i]=mahalanobis(sval,center,cm)
  mdist[,i]=mahalanobis(sval,center,ginv(cm), inverted = TRUE)           #Pose toujours un pb pour VL=1
}

cl <- factor(nm[max.col(-mdist)],levels = nm)
pred=list(class=cl,dist=mdist)
return(pred)
}
