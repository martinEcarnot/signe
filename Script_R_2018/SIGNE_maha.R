# From score of LDA, compute Mahalanobis distances

SIGNE_maha = function (rlda,sc_cal,class_cal,sc_val) {

  # browser()
lsc_cal=predict(rlda,sc_cal)$x
lsc_val=predict(rlda ,sc_val)$x
ng=length(rlda$lev)
mdist=matrix(nrow = nrow(lsc_val),ncol = ng)

for (i in 1:ng) {
  igroup=which(class_cal==rlda$lev[i])
  # Create covariance matrix
  cm=cov(lsc_cal[igroup,]) # xc=x-matrix(data=1, nrow=54)%*%colMeans(x), cm=t(xc)%*%xc
  center=colMeans(lsc_cal[igroup,])  # rlda$means[i,]%*%rlda$scaling
  mdist[,i]=mahalanobis(lsc_val,center,cm)
}

nm <- rlda$lev
cl <- factor(nm[max.col(-mdist)],levels = nm)
pred=list(class=cl,posterior=mdist)

# perok=100*sum(diag(table(pred$class,class_val)))/length(class_val)
# print(perok)

return(pred)
}

# Plots
# cal=data.frame(calscor)
# val=data.frame(newsc)
# pl=ggplot(cal) + geom_point(aes(LD1, LD2, colour = class_cal), size = 2.5) +
#   geom_point(data=val,aes(LD1, LD2, colour = clas), size = 4, shape="+")
# # labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
# #      y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))
# plot(pl)
