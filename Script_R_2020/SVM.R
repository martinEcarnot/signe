sp5=sp4
#iout=sample(5210,5100)
#sp5 <- sp5[-iout,]
clas=as.factor(substr(rownames(sp5),9,9))
rplsda=caret::plsda(sp5, clas,ncomp=20)
axe1<- 1
axe2<- 2
axeX <- rplsda$scores[,axe1] ; axeY <- rplsda$scores[,axe2]

length(rplsda$scores[,1])


type=as.character(substr(rownames(sp5),9,9))
#type[which(type=="C")]="A"
#type[which(type=="S")]="A"
#type[which(type=="G")]="B"
df=cbind.data.frame(a=rplsda$scores[,1], b=rplsda$scores[,2], c=rplsda$scores[,3], d=rplsda$scores[,4], e=rplsda$scores[,5], f=rplsda$scores[,6], g=rplsda$scores[,7], h=rplsda$scores[,8], i=rplsda$scores[,9], j=rplsda$scores[,10], k=rplsda$scores[,11], l=rplsda$scores[,12], m=rplsda$scores[,13], n=rplsda$scores[,14], o=rplsda$scores[,15], p=rplsda$scores[,16], q=rplsda$scores[,17], r=rplsda$scores[,18], s=rplsda$scores[,19], t=rplsda$scores[,20], y=type)
rownames(df)=1:5210
summary(df)
#iout=sample(110,50)
#df <- df[-iout,]
df

plot(df$a,df$b,type="n")
text(df$a,df$b,rownames(df),col=c("blue","red","green")[df$y],cex=0.75)

#charger le packagee1071
library(e1071)


#svm.ovo <- SVM(x=iris[,1:4], y=iris[,5], class.type="one.versus.one", verbosity=0)

mpoly <-svm(y ~ a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t, data=df, class.type="one.versus.one", kernel="polynomial", scale=F, cost=100, coef0=1, degree=2)
print(mpoly)


#prédiction sur l’échantillon d’apprentissage, matrice de confusion
ypoly.df <-predict(mpoly,df)
T=table(df$y,ypoly.df)
M=as.matrix(T)
M
sum(diag(M))/sum(M)
M[1,1]/(M[1,1]+M[1,2]+M[1,3])
M[2,2]/(M[2,1]+M[2,2]+M[2,3])
M[3,3]/(M[3,1]+M[3,2]+M[3,3])

#représentation des points supports dans le plan
plot(df$a,df$b,col=c("blue","red","green")[df$y])
points(df$a[mpoly$index],df$b[mpoly$index],pch=5,cex=1.75,col=rgb(0,0,0))


#plotting
plot(df$a,df$b,col=c("blue","red","green")[ypoly.df])






#pour disposer toujours de la même séquence de données

set.seed(10)

#fonction générant le concept à modéliser
generate.data <-function(n){
  x2 <-runif(n)
  x1 <-runif(n)
  y <-factor(ifelse((x2>2*x1)|(x2>(2-2*x1)),1,2))
  return(data.frame(x1,x2,y))
}

#échantillon d’apprentissage -n = 30
dtrain <-generate.data(30)
#représentation graphique
plot(dtrain$x1,dtrain$x2,col=c("blue","red")[dtrain$y],xlim=c(0,1),ylim=c(0,1))
abline(0,2)
abline(2,-2)


#échantillon test -ntest= 5000
dtest <-generate.data(5000)
#représentation graphique
plot(dtest$x1,dtest$x2,col=c("blue","red")[dtest$y],xlim=c(0,1),ylim=c(0,1))
abline(0,2)
abline(2,-2)

#SVM linéaire
mlin <-svm(y ~ x1+x2, data=dtrain, kernel="linear",scale=F)
print(mlin)

#représentation des points supports dans le plan
plot(dtrain$x1,dtrain$x2,col=c("blue","red")[dtrain$y],xlim=c(0,1),ylim=c(0,1))
points(dtrain$x1[mlin$index],dtrain$x2[mlin$index],pch=5,cex=1.75,col=rgb(0,0,0))
abline(0,2)
abline(2,-2)


#constante 
beta.0 <--mlin$rho
#coefficients1et 
beta.1 <-sum(mlin$coefs*dtrain$x1[mlin$index])
beta.2 <-sum(mlin$coefs*dtrain$x2[mlin$index])
#tracé de la droite de séparation
plot(dtrain$x1,dtrain$x2,col=c("blue","red")[dtrain$y],xlim=c(0,1),ylim=c(0,1))
abline(-beta.0/beta.2,-beta.1/beta.2,col="green")

#prédiction sur l’échantillon d’apprentissage, matrice de confusion
ylin.train <-predict(mlin,dtrain)
table(dtrain$y,ylin.train)

#utilisation d’un noyau polynomial de degré 2
mpoly <-svm(y ~ x1+x2, data=dtrain, kernel="polynomial", scale=F, coef0=1, degree=2)
print(mpoly)

#prédictionsur l’échantillon test
ypoly.test <-predict(mpoly,dtest)
#matrice de confusion
mc.poly <-table(dtest$y,ypoly.test)
print(mc.poly)
#taux d’erreur en test
err.poly <-1-sum(diag(mc.poly))/sum(mc.poly)
print(err.poly)

#plotting
plot(dtest$x1,dtest$x2,col=c("blue","red")[ypoly.test],xlim=c(0,1),ylim=c(0,1))
abline(0,2,lwd=2)
abline(2,-2,lwd=2)

#nouveau paramétrage avec C = 1000
mpoly <-svm(y ~ x1+x2, data=dtrain, kernel="polynomial", scale=F, cost=1000, coef0=1, degree=2)
print(mpoly)
##Refaire tourner ce qu'il y au dessus avec ce nouveau paramètre


plot(dtrain$x1,dtrain$x2,col=c("blue","red")[dtrain$y],xlim=c(0,1),ylim=c(0,1))
points(dtrain$x1[mpoly$index],dtrain$x2[mpoly$index],pch=5,cex=1.75,col=rgb(0,0,0))
abline(0,2)
abline(2,-2)











#apprentissage
#standardisation automatique des données, cf. paramètre scale
m1 <-svm(class ~ ., data = dtrain)

#affichage
print(m1)

#prediction
y1 <-predict(m1,newdata=dtest)

#matrice de confusion
mc1 <-table(dtest$class,y1)
err1 <-1 -sum(diag(mc1))/sum(mc1)
print(err1)










## Not run:
# # train SVM from data in x and labels in y
# svm <- SVM(x, y, core="libsvm", kernel="linear", C=1)
#
# # train SVM using a dataset with both data and lables and a formula pointing to labels
# formula <- target ~ .
# svm <- SVM(formula, data, core="svmlight", kernel="rbf", gamma=1e3)
#
# # train a model with 2eSVM algorithm
# data(svm_breast_cancer_dataset)
# ds <- svm.breastcancer.dataset
# svm.2e <- SVM(x=ds[,-1], y=ds[,1], core="libsvm", kernel="linear", prep = "2e", C=10);
# # more at \url{http://r.gmum.net/samples/svm.2e.html}
#
# # train SVM on a multiclass data set
# data(iris)
# # with "one vs rest" strategy
# svm.ova <- SVM(Species ~ ., data=iris, class.type="one.versus.all", verbosity=0)
# # or with "one vs one" strategy
# svm.ovo <- SVM(x=iris[,1:4], y=iris[,5], class.type="one.versus.one", verbosity=0)
#
# # we can use svmlights sample weighting feature, suppose we have weights vector
# # with a weight for every sample in the traning data
# weighted.svm <- SVM(formula=y~., data=df, core="svmlight", kernel="rbf", C=1.0,
#                     gamma=0.5, example.weights=weights)
#
# # svmlight alows us to determine missing labels from a dataset
# # suppose we have a labels y with missing labels marked as zeros
# svm.transduction <- SVM(x, y, transductive.learning=TRUE, core="svmlight")
#
# # for more in-depth examples visit \url{http://r.gmum.net/getting_started.html}
# ## End(Not run)













