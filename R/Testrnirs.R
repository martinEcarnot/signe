# https://docs.google.com/document/d/1FwGD7RSJekT39E3wMwtO7nHWK7kGs-4-XMS8jupKh6w/edit

library(devtools)
install_github("mlesnoff/rnirs", dependencies = TRUE)



##### Exercice 1 ACP

library(rnirs)
library(gridExtra)

data(datforages)
names(datforages)
Xr <- datforages$Xr
yr <- datforages$yr
Xu <- datforages$Xu
yu <- datforages$yu


#rm(list = ls())
library(rnirs)

library(gridExtra)

data(datforages)
names(datforages)

#r pour référence, u pour unknown
Xr <- datforages$Xr
yr <- datforages$yr

Xu <- datforages$Xu
yu <- datforages$yu

headm(Xr)
headm(Xu)
table(yr)
table(yu)

n <- nrow(Xr)
m <- nrow(Xu)

zXr <- snv(Xr, center = FALSE)
zXu <- snv(Xu, center = FALSE)

################# PCA

ncomp <- 10
fm <- pca(zXr, ncomp = ncomp)

#fm <- pca(zXr, ncomp = ncomp, algo = pca.eigen)
#fm <- pca(zXr, ncomp = ncomp, algo = pca.nipals)
names(fm)

### Proportion of explained X-variance

z <- fm$explvarx
barplot(100 * z$pvar,
        names.arg = paste("comp", z$ncomp),
        ylab = "% Variance explained")

### Scores

headm(fm$Tr)

comp <- c(1, 2)
#comp <- c(2, 3)
z <-
  plotxy(fm$Tr[, comp])

plotxy(fm$Tr[, comp], label = TRUE)

plotxy(fm$Tr[, comp], group = yr)

plotxy(fm$Tr[, comp], group = yr, ellipse = TRUE)

p1 <- plotxy(fm$Tr[, c(1, 2)])
p2 <- plotxy(fm$Tr[, c(3, 4)])
grid.arrange(p1, p2, ncol = 2)

### Individual contributions (proportions)

headm(fm$contr.ind)
colSums(fm$contr.ind)

z <-fm$contr.ind ##Ici, on isole les individus qui participent le plus à la création de l'axe j
j <- 1           ## On dit >1/n, mais on peut les choisir différemment
#j <- 2
res <- z[z[, j] > 1 / n, j, drop = FALSE]
headm(res)
res <- z[order(z[, j], decreasing = TRUE), j, drop = FALSE]
headm(res)

### Loadings  #Les graphs suivants sont pas forcément très intéressants en spectro

headm(fm$P)

plotsp(t(fm$P), xlab = "wavelength (nm)", ylab = "Loadings")

plotsp(t(fm$P), text = TRUE, col = rainbow(ncomp)) #AUgment lisibilité

#plotsp1(t(fm$P), main = "Loadings") # "Interactif". Quand on fait "entrée", ca affiche un nouveau truc.

plotxy(fm$P, origin = c(0, 0))

### Variable coordonates

plotxy(fm$coord.var, origin = c(0, 0))

###  Correlation circle

comp <- c(1, 2)
#comp <- c(2, 3)
z <- fm$cor.circle[, comp]
plotxy(z, origin = c(0, 0), circle = TRUE)

################# PROJECTION OF Xu ON THE SCORE SPACE Tr

ncomp <- 10
fm <- pca(zXr, zXu, ncomp = ncomp)
headm(fm$Tr)
headm(fm$Tu)

res <- rbind(fm$Tr, fm$Tu)
group <- c(rep("Reference", n), rep("Unknown", m))
comp <- c(1, 2)
#comp <- c(2, 3)
plotxy(res[, comp], group = group, origin = c(0, 0))





##### Exercice 2

#On utilise SD et OD. OD distance entre point et point projeté dans l'espace utilisé,
#SD, distance entre point projeté et centre du nuage de cal.
#rm(list = ls())
library(rnirs)

library(gridExtra)

data(datcass)
names(datcass)

Xr <- datcass$Xr
Xu <- datcass$Xu

headm(Xr)
headm(Xu)

n <- nrow(Xr)
m <- nrow(Xu)

zXr <- detrend(Xr)
zXu <- detrend(Xu)

############# SD AND OD ON REFERENCE DATA SET

ncomp <- 15

fm <- pca(zXr, ncomp = ncomp)

res <- scordis(fm)        # <== SD
#res <- odis(fm, zXr)     # <== OD
names(res)

headm(res$dr)

d <- res$dr$d
summary(d)
# Extreme values
par(mfrow = c(1, 2))
plot(ecdf(d), xlab = "SD", pch = 1, ylab = "Fn(SD)", main = "")
abline(v = res$cut, col = "red", lty = 4)
hist(d, n = 100, col = "grey", xlab = "SD", main = "")
abline(v = res$cut, col = "red", lty = 4)
par(mfrow = c(1, 1))

z <- res$dr
z[z$dstand > 1, ]
res$cut

# For SD only:
# mean of the squared SD in the reference data set should be ~ a
d <- res$dr$d
mean(d^2)

############# SD-OD PLOT ON REFERENCE DATA SET
## Observations in the right-up side of the plot have "bad" leverages

ncomp <- 15
fm <- pca(zXr, ncomp = ncomp)

res.sd <- scordis(fm)
res.od <- odis(fm, zXr)

d.sd <- res.sd$dr$d
d.od <- res.od$dr$d

plot(d.sd, d.od, xlab = "SD", ylab = "OD")
abline(v = res.sd$cut, col = "red", lty = 4)
abline(h = res.od$cut, col = "red", lty = 4)

############# REFERENCE vs. TEST DATA SETS

ncomp <- 15
fm <- pca(zXr, zXu, ncomp = ncomp)

res.sd <- scordis(fm)
res.od <- odis(fm, zXr, zXu)

headm(res.sd$dr)
headm(res.sd$du)

headm(res.od$dr)
headm(res.od$du)

d.sdr <- res.sd$dr$d
d.odr <- res.od$dr$d
d.sdu <- res.sd$du$d
d.odu <- res.od$du$d

plot(d.sdr, d.odr, xlab = "SD", ylab = "OD",
     xlim = c(0, max(d.sdr, d.sdu)),
     ylim = c(0, max(d.odr, d.odu)))
points(d.sdu, d.odu, col = "purple", pch = 16)
abline(v = res.sd$cut, col = "red", lty = 4)
abline(h = res.od$cut, col = "red", lty = 4)





#####Exemple 3

#rm(list = ls())
library(rnirs)

data(datcass)
names(datcass)

Xr <- datcass$Xr
yr <- datcass$yr

Xu <- datcass$Xu
yu <- datcass$yu

headm(Xr)
headm(Xu)

n <- nrow(Xr)
m <- nrow(Xu)

zXr <- detrend(Xr, method = "lowess") #detrend : correction de la ligne de base. lowess, regression linéaire locale (fait des fenêtres, un peu comme Savitsky-Golay).
zXu <- detrend(Xu, method = "lowess")

plotsp(zXr)

#################### PLS model

ncomp <- 15
fm <- pls(zXr, yr, ncomp = ncomp)
names(fm)

### Proportion of explained X-variance

fm$explvarx

### Projection of Xu on the score space Tr ==> Scores Tu

fm <- pls(zXr, yr, zXu, ncomp = ncomp)
#fm <- pca(zXr, zXu, ncomp = ncomp) ##On constate ici que la fonction PLS a exactement les mêmes entrées/sorties que la fonction ACP.
headm(fm$Tr)
headm(fm$Tu)

res <- rbind(fm$Tr, fm$Tu)
group <- c(rep("Reference", n), rep("Unknown", m))
comp <- c(1, 2)
#comp <- c(2, 3)
plotxy(res[, comp], group = group, origin = c(0, 0))  ##Justifie ici l'utilisation de la locale, qui pourrait donner de meilleurs résultats qu'une globale, particulièrement ici

### b coefficients (model with ncomp components)

b <- bcoef(fm)
headm(b)

z <- t(b)
z <- z[, 2:ncol(z), drop = FALSE]
plotsp(z, xlab = "wavelength (nm)", ylab = "b-coefficient")

#################### PLSR model   ##PLSR, r pour regression. C'est ca qui permet de faire de la prédiciton quanti.

ncomp <-15
fm <- plsr(zXr, yr, zXu, yu, ncomp = ncomp)
#fm <- plsr(zXr, yr, zXu, ncomp = ncomp)
names(fm)
headm(fm$y)
headm(fm$fit)
headm(fm$r)
z <- mse(fm, ~ ncomp)
z
z[z$rmsep == min(z$rmsep), ]
plotmse(z)

### b coefficients (model with ncomp components)

b <- bcoef(fm)
headm(b)

bcoef(fm$fm)

### Cross-validations

segm <- segmcvkfold(n = n, K = 5, typ = "random", nrep = 20)
fm <- fitcv(
  zXr, yr,
  fun = plsr,
  ncomp = ncomp,
  segm = segm,
  print = TRUE
)
names(fm)
headm(fm$fit)
z <- mse(fm, ~ ncomp) ##"Moyenne sur l'ensemble des segments" ?
#z <- mse(fm, ~ ncomp + segm + rep)
headm(z)
z[z$rmsep == min(z$rmsep), ]
plotmse(z)

u <- selncomp.wold(z, alpha = .01)
u$res
u$sel
u$opt

segm <- segmcvkfold(n = n, K = 5, typ = "random", nrep = 20)
#segm <- segmcvmc(n = n, m = round(n / 5), nrep = 20)
fm <- fitcv(
  zXr, yr,
  fun = plsr,
  ncomp = ncomp,
  segm = segm,
  print = TRUE
)
headm(fm$fit)
z <- mse(fm, ~ ncomp + rep)
headm(z)
plotmse(z, group = "rep", col = "grey")






library(ggplot2)









#####Exercice 4 FDA (Analyse factorielle discriminante)

#rm(list = ls())
library(rnirs)

####################### fda

data(iris)

headm(iris)

Xr <- iris[, 1:4]
yr <- iris[, 5]

fm <- fda(Xr, yr)
names(fm)
# Xr-scores
headm(fm$Tr)
# Xr-class centers scores
fm$Tcenters
# Xr-loadings matrix
# = coefficients of linear discriminants
# = "LD" of function lda of MASS package
fm$P
# Explained variance by PCA of the class centers
# in transformed scale
fm$explvar

p <- plotxy(fm$Tr, group = yr)
p <- p + geom_point(
  data = data.frame(fm$Tc), aes(x = comp1, y = comp2),
  pch = 8, cex = 4, col = "blue"
)
p

# Object Tcenters is the projection of the class centers in the score space
fm <- fda(Xr, yr)
fm$Tcenters
centers <- centr(Xr, yr)$centers
centers
fm <- fda(Xr, yr, centers)
fm$Tu

###A partir d'ici, comparaison de différenets techniques de DA. (Pas important pour Benoît)

############################################################# dalm

data(iris)

X <- iris[, 1:4]
y <- iris[, 5]
N <- nrow(X)

m <- round(.25 * N)
n <- N - m
s <- sample(1:N, m)
Xr <- X[-s, ]
yr <- y[-s]
Xu <- X[s, ]
yu <- y[s]

fm <- dalm(Xr, yr, Xu, yu)
names(fm)
head(fm$y)
head(fm$fit)
head(fm$r)
head(fm$dummyfit)
fm$ni
err(fm)

############################################################# dadis ###AVec Mahalanobis, c'est la notre.

data(iris)

X <- iris[, 1:4]
y <- iris[, 5]
N <- nrow(X)

m <- round(.25 * N)
n <- N - m
s <- sample(1:N, m)
Xr <- X[-s, ]
yr <- y[-s]
Xu <- X[s, ]
yu <- y[s]

fm <- dadis(Xr, yr, Xu, yu)
names(fm)
head(fm$y)
head(fm$fit)
head(fm$r)
fm$ni
err(fm)

fm <- dadis(Xr, yr, Xu, yu, diss = "mahalanobis")
err(fm)

nclas <- length(unique(yr))
W <- matW(Xr, yr)$W * n / (n - nclas)
fm <- dadis(Xr, yr, Xu, yu, diss = "mahalanobis", sigma = W)
err(fm)

############################################################# daglm

data(iris)

X <- iris[, 1:4]
y <- iris[, 5]
N <- nrow(X)

m <- round(.25 * N) # Test
n <- N - m          # Training
s <- sample(1:N, m)
Xr <- X[-s, ]
yr <- y[-s]
Xu <- X[s, ]
yu <- y[s]

## Binomial model with logit link (logistic regression)
fm <- daglm(Xr, yr, Xu, yu)
names(fm)
head(fm$y)
head(fm$fit)
head(fm$r)
fm$ni
err(fm)

## Gaussian model with identity link (= usual linear model)
fm <- daglm(Xr, yr, Xu, yu, family = gaussian)
err(fm)

############################################################# daprob #"Les célèbres LDA."

data(iris)

X <- iris[, 1:4]
y <- iris[, 5]
N <- nrow(X)

m <- round(.25 * N) # Test
n <- N - m          # Training
s <- sample(1:N, m)
Xr <- X[-s, ]
yr <- y[-s]
Xu <- X[s, ]
yu <- y[s]

##### LDA (homogeneous covariances)

fm <- daprob(Xr, yr, Xu, yu, dens = dmnorm) #"dmnorm suppose que chaque classe a une distribution gaussienne multivariable"
names(fm)
head(fm$y)
head(fm$fit)
head(fm$r)
fm$ni
err(fm)

##### QDA (heterogeneous covariances)

fm <- daprob(Xr, yr, Xu, yu, dens = dmnorm, lda = FALSE)
err(fm)

##### Nonparametric DA

fm <- daprob(Xr, yr, Xu, yu, dens = dkern.gauss, h = .2)
err(fm)

############################################################# dasdod

data(iris)

X <- iris[, 1:4]
y <- iris[, 5]
N <- nrow(X)

m <- round(.25 * N) # Test
n <- N - m          # Training
s <- sample(1:N, m)
Xr <- X[-s, ]
yr <- y[-s]
Xu <- X[s, ]
yu <- y[s]

ncompcla <- 3
fm <- dasdod(Xr, yr, Xu, yu, ncompcla = ncompcla)
names(fm)
head(fm$y)
head(fm$fit)
head(fm$r)
head(fm$d)
head(fm$sd.stand)
head(fm$od.stand)
fm$ni
fm$ncompcla
fm$pvarcla
err(fm)
err(fm, ~ y1)

## With a preliminary PLS

fm <- dasdod(Xr, yr, Xu, yu, ncomp.pls = 3, ncompcla = 2)
fm$ncompcla
fm$pvarcla
err(fm)







##### Exercice 5

#rm(list = ls())
library(rnirs)

library(gridExtra)

data(datforages)

Xr <- datforages$Xr
yr <- datforages$yr

Xu <- datforages$Xu
yu <- datforages$yu

headm(Xr)
headm(Xu)

table(yr)
table(yu)

zXr <- savgol(snv(Xr), n = 21, p = 2, m = 2) ##J'ai utilisé n'imp quoi comme prétraitements, c'est ceux qui focntionnent le moins bien
zXu <- savgol(snv(Xu), n = 21, p = 2, m = 2)

######## PLS-DALM

ncomp <- 20
fm <- plsda(zXr, yr, zXu, yu, ncomp = ncomp)
# This is the same as:
# fm <- plsda(Xr, yr, Xu, yu, ncomp = ncomp,
#   da = dalm)
# By default in plsda, da = dalm

names(fm)
headm(fm$y)
headm(fm$fit)
headm(fm$r)

z <- err(fm, ~ ncomp)
plotmse(z, nam = "errp")
z[z$errp == min(z$errp), ]

## Same with plsdalm (faster)

ncomp <- 20
fm <- plsdalm(zXr, yr, zXu, yu, ncomp = ncomp)

z <- err(fm, ~ ncomp)
plotmse(z, nam = "errp")
z[z$errp == min(z$errp), ]

######## PLS-LDA

ncomp <- 20
fm <- plsda(zXr, yr, zXu, yu, ncomp,
            da = daprob)

z <- err(fm, ~ ncomp)
plotmse(z, nam = "errp")
z[z$errp == min(z$errp), ]

######## PLS-QDA

ncomp <- 20
fm <- plsda(zXr, yr, zXu, yu, ncomp,
            da = daprob, lda = FALSE)
names(fm)
headm(fm$y)
headm(fm$fit)
headm(fm$r)

z <- err(fm, ~ ncomp)
z
plotmse(z, nam = "errp")
z <- z[z$errp == min(z$errp), ]
z

######## PLS-KNNDA

ncomp <- 20
fm <- plsda(zXr, yr, zXu, yu, ncomp, da = knnwda,
            diss = "mahalanobis", k = 3) #Si on précise pas, h est infini, donc pas de plus proches voisins.

z <- err(fm, ~ ncomp + k)
plotmse(z, nam = "errp", group = "k")
z[z$errp == min(z$errp), ]





