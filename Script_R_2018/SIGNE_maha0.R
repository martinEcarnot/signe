# From matrix cal matrix, cal classes and val matrix, compute Mahalanobis dist from each group 

SIGNE_maha0 = function (scal,class_cal,sval) {

ng=nlevels(class)
nm=levels(class)
mdist=matrix(nrow = nrow(sval),ncol = ng)

for (i in 1:ng) {
  igroup=which(class_cal==nm[i])
  cm=cov(scal[igroup,]) # Create covariance matrix     # xc=x-matrix(data=1, nrow=54)%*%colMeans(x), cm=t(xc)%*%xc
  center=colMeans(scal[igroup,])  # rlda$means[i,]%*%rlda$scaling
  # mdist[,i]=mahalanobis(sval,center,cm)
  mdist[,i]=mahalanobis(sval,center,ginv(cm), inverted = TRUE)
}

cl <- factor(nm[max.col(-mdist)],levels = nm)
pred=list(class=cl,dist=mdist)

return(pred)
}
