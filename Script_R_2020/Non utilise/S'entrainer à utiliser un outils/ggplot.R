  sp6=data.frame(sp=sp5,axeX=rplsda$scores[,axe1],axeY=rplsda$scores[,axe2],jour=substr(rownames(sp5),4,8),annee=substr(rownames(sp5),1,4),cepage=substr(rownames(sp5),9,9),clone=substr(rownames(sp5),9,13),parcelle=substr(rownames(sp5),18,18))
sp6$axeX

# CAS 1
ggplot(sp6) +
  geom_point(aes(x = axeX, y = axeY, size = 1, color = cepage)) +
  scale_size(1)+scale_colour_discrete() # _ordinal
###Regarder ce que c'est que plotly, M est fan

library(plotly)

library(ggplot2)



# Changer la taille et la forme
g<-ggplot(sp6, aes(x=axeX, y=axeY, color=cepage)) +
  geom_point(size=0.5, alpha=1)+scale_color_discrete()

ggplotly(g)


unique(substr(rownames(sp5),9,13))




iout=sample(5210,5200)
sp1 <- sp[-iout,]
length(sp1[,1])
sp1
as.vector(sp1)





library(ggplot2)

spx=as.vector(sp1)
spp=data.frame(spx=as.vector(sp1),num=1:length(spx), colo=sp1)
spp

length(spx[1])

geom_point(size=5, alpha=0.5)

p <- ggplot(spp, aes(x = num, y = spx, colour=cepage)) + geom_point(size=0.5)
print(p)



p1 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet, group=Chick)) +
  geom_line() +
  ggtitle("Growth curve for individual chicks")
p1





aC= sp6$annee=="2017"
sp7 =sp6[(aC==TRUE),]



bleu=rgb(0.09,0.63,1)
bleu2=rgb(0,0.43,0.8)
bleu3=rgb(0.59,0.83,1)
vert=rgb(0.10,0.94,0.36)
vert2=rgb(0.12,0.75,0.10)
vert3=rgb(0.50,0.94,0.36)
vert4=rgb(0.50,0.75,0.36)
rouge=rgb(1,0.35,0.13)
rouge2=rgb(0.8,0.35,0.13)
rouge3=rgb(1,0.55,0.33)

coloclone=c(rouge, rouge2, rouge3, bleu, bleu2, bleu3, vert, vert2, vert3, vert4)

# Nuage de points avec estimation de la densité 2d
sp <- ggplot(sp6, aes(x=axeX, y=axeY,colour=clone)) +
  geom_point(size=0.5, alpha=1) +
  scale_color_manual(values = coloclone) +
  geom_density_2d(data = sp6, size=0.3, bins=3)
ggplotly(sp)


# Gradient de couleur
sp + stat_density_2d(aes(fill = ..level..), geom="polygon")
# Changer le gradient de couleur
sp + stat_density_2d(aes(fill = ..level..), geom="polygon")+
  scale_fill_gradient(low="blue", high="red")



# Une ellipse autour de tous les points
ggplot(faithful, aes(waiting, eruptions))+
  geom_point()+
  stat_ellipse()
# Ellipses par groupes
p <- ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3))+
  geom_point()
p + stat_ellipse()
# Changer le type d'ellipses:
# Valeurs possibles "t", "norm", "euclid"
p + stat_ellipse(type = "norm")









mypal

ggplot(sp6, aes(x=axeX, y=axeY, color=colo)) +
  geom_point(size=1) +
  geom_text(label=texte)
C=as.vector(c("red","blue","yellow","green"))
as.
rgb(0.1,0.1,0.2)

ggplot(sp6, aes(x=axeX, y=axeY, color=colo)) +
  geom_point(size=1) +
  geom_density2d(color = "red")# +
#  geom_smooth(method = "lm")


scale_color_brewer("Département", palette = "Set1")

# Changer le type de points en fonction des niveaux de cyl
ggplot(mtcars, aes(x=wt, y=mpg, shape=cyl)) +
  geom_point()
# Changer le type et la couleur
ggplot(mtcars, aes(x=wt, y=mpg, shape=cyl, color=cyl)) +
  geom_point()


library("RColorBrewer")
display.brewer.all()
colors <- brewer.pal(4, "BuPu")
colors



pal <- colorRampPalette(colors)
pal(5) # palette en 5 classes
pal(30) # palette en 30 classes




library("devtools")


mypal <- getpal("nicopal",10)
mypal




# 1 - CREATION DU FICHIER colors.RData
source_url("https://raw.githubusercontent.com/neocarto/R/master/colors/BuildColors.R")

blues<-buildpal("#8ECEE9","#146594","#DCF0F8","#0B2130","#368ECA",20)
greens<-buildpal("#ADD39D","#1B7D34","#E5F0DA","#142F1B","#3FA136",20)
reds<-buildpal("#EFA5A9","#991915","#FBDEE1","#1D0809","#DE4243",20)
purples<-buildpal("#FFB6F1","#87036E","#FFD9FA","#3F0036","#FF00C3",20)
oranges<-buildpal("#FFB888","#AA4500","#FFE1CC","#3D1900","#FF6700",20)
nicopal<-buildpal("#FCE37A","#FF0000","#FFEB97","#E50000","#FF6700",20)
nicopal2<-buildpal("#B6EFB6","#080E5B","#B6EFB6","#080E5B","#518996",20)

colors<-c(blues=list(blues),reds=list(reds),greens=list(greens),purples=list(purples),oranges=list(oranges),nicopal=list(nicopal),nicopal2=list(nicopal2))
save(colors, file = "colors.RData")






# 2 - FABRICATION DES PALETTES CARTOGRAPHIQUES
source_url("https://raw.githubusercontent.com/neocarto/R/master/colors/Palettes.R")
load("colors.RData")


# Simples
mypal <- getpal("nicopal",3)
displaypal(mypal)
mypal


# Doubles
mypal <- getpal("blues",4,"reds",4)
displaypal(mypal)
mypal


mypal <- getpal("blues",3,"reds",9,middle=T,alphaeffect=T)
displaypal(mypal)
mypal

??getpal


m <- ggplot(faithful, aes(x = eruptions, y = waiting)) +
  geom_point() +
  xlim(0.5, 6) +
  ylim(40, 110)
m + geom_density_2d()

m + stat_density_2d(aes(fill = stat(level)), geom = "polygon")

set.seed(4393)
dsmall <- diamonds[sample(nrow(diamonds), 1000), ]
d <- ggplot(dsmall, aes(x, y))
# If you map an aesthetic to a categorical variable, you will get a
# set of contours for each value of that variable
d + geom_density_2d(aes(colour = cut))

# Similarly, if you apply faceting to the plot, contours will be
# drawn for each facet, but the levels will calculated across all facets
d + stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
  facet_grid(. ~ cut) + scale_fill_viridis_c()
# To override this behavior (for instace, to better visualize the density
# within each facet), use stat(nlevel)
d + stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
  facet_grid(. ~ cut) + scale_fill_viridis_c()

# If we turn contouring off, we can use use geoms like tiles:
d + stat_density_2d(geom = "raster", aes(fill = stat(density)), contour = FALSE)
# Or points:
d + stat_density_2d(geom = "point", aes(size = stat(density)), n = 20, contour = FALSE)


