install.packages("cowplot")


library("cowplot")
library("gridGraphics")

df <- ToothGrowth
# Convertir la colonne dose en facteur
df$dose <- as.factor(df$dose)
head(df)


library(cowplot)
# Graphe par défaut
bp <- ggplot(df, aes(x=dose, y=len, color=dose)) +
  geom_boxplot() +
  theme(legend.position = "none")
bp
# Ajouter les grilles
bp + background_grid(major = "xy", minor = "none")


save_plot("mpg.pdf", bp,
          base_aspect_ratio = 1.3 #Laisser de la place pour la légende
)


# Nuage de points
sp <- ggplot(mpg, aes(x = cty, y = hwy, colour = factor(cyl)))+
  geom_point(size=2.5)
sp
# Bar plot
bp <- ggplot(diamonds, aes(clarity, fill = cut)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle=70, vjust=0.5))
bp

plot_grid(sp, bp, labels=c("A", "B"), ncol = 2, nrow = 1)

draw_plot(plot, x = 0, y = 0, width = 1, height = 1)





plot.iris <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
  geom_point() + facet_grid(. ~ Species) +
  stat_smooth(method = "lm") +
  background_grid(major = 'y', minor = "none") + # grilles horiz.
  panel_border() # Bordure autour de chaque panel

ggdraw() +
  draw_plot(plot.iris, 0, .5, 1, .5) +
  draw_plot(sp, 0, 0, .5, .5) +
  draw_plot(bp, .5, 0, .5, .5) +
  draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)



##### Facet


# Convertir la variable dose de type "numeric" au type "factor"
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
df <- ToothGrowth
head(df)

library(ggplot2)
bp <- ggplot(df, aes(x=dose, y=len, group=dose)) +
  geom_boxplot(aes(fill=dose))
bp

aff <- ggplot(Ldist, aes(x=Ldist[,6], y=(Ldist[,2]),colour=paste(Ldist[,3],Ldist[,5]),date=Ldist[,7],cepage=Ldist[,3],predit=Ldist[,4])) +
  geom_point(size=2, alpha=1) +
  scale_color_manual(values = syrah)
ggplotly(aff)

# Partitionnement vertical
bp + facet_grid(supp ~ .)
# Partitionnement horizontal
bp + facet_grid(. ~ supp)

# Facet avec deux variables: dose et supp.
# Les lignes sont "dose" et les colonnes sont "supp"
bp + facet_grid(dose ~ supp)
# Facet avec deux variables: inverser l'ordre des 2 variables
# Les lignes sont "supp" et les colonnes sont "dose"
bp + facet_grid(supp ~ dose)

bp + facet_grid(dose ~ supp, margins=TRUE)

bp + facet_grid(dose ~ supp, scales='free')

bp + facet_grid(dose ~ supp, labeller=label_both)


# Modifier le texte. Valeurs possibles pour le style de police:
#'plain', 'italic', 'bold', 'bold.italic'.
bp + facet_grid(dose ~ supp)+
  theme(strip.text.x = element_text(size=12, color="red",
                                    face="bold.italic"),
        strip.text.y = element_text(size=12, color="red",
                                    face="bold.italic"))
# Modifier l'apparence du rectangle autour
# des étiquettes de panneaux
bp + facet_grid(dose ~ supp)+
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=1.5, linetype="solid"))

bp + facet_wrap(~ dose)
bp + facet_wrap(~ dose, ncol=2)






















#set working directory
frames = 35
for(i in 1:frames){
  # creating a name for each plot file with leading zeros
  if (i < 10) {name = paste('000',i,'plot.png',sep='')}
  if (i < 100 && i >= 10) {name = paste('00',i,'plot.png', sep='')}
  if (i >= 100) {name = paste('0', i,'plot.png', sep='')}

  x = seq(0, i, 1)
  f.3 = dbinom(x, size = i, prob=.3)
  f.7 = dbinom(x, size = i, prob=.7)
  #saves the plot as a .png file in the working directory

  ggsave(name, plot =   ggplot(Ldist2, aes(x=numero, y=dista[,i],colour=paste(predclone[,i],succes[,i]),date=date,clone=clone,predit=predclone[,i])) +
           geom_point(size=2, alpha=1) +
           scale_color_manual(values = gamay) +
           theme(legend.position="none"))



  # ng(name)
  # ggplot(Ldist2, aes(x=numero, y=dista[,2],colour=paste(predclone[,2],succes[,2]),date=date,clone=clone,predit=predclone[,2])) +
  #   geom_point(size=2, alpha=1) +
  #   scale_color_manual(values = gamay) +
  #   theme(legend.position="none")
#  dev.off()
}




ggplot(Ldist2, aes(x=numero, y=dista[,2],colour=paste(predclone[,2],succes[,2]),date=date,clone=clone,predit=predclone[,2])) +
  geom_point(size=2, alpha=1) +
  scale_color_manual(values = gamay) +
  theme(legend.position="none")


plot(x, f.3, type='h', xlim = c(0,frames), ylim = c(0,.7), ylab ="probability",   main = paste("Binomial density with n = ", i), col = "red")
lines(x,f.7,type='h',col='blue')
text(45, .6, 'p = .3', col='red')
text(45, .6, 'p = .7', col='blue', pos=1)


# system("convert *.png -delay 3 -loop 0 binom.gif")
# system("convert -delay 80 *.png example_1.gif")







dir.create("examples")
setwd("examples")

# example 1: simple animated countdown from 10 to "GO!".
png(file="example%02d.png", width=200, height=200)
for (i in c(10:1, "G0!")){
  plot.new()
  text(.5, .5, i, cex = 6)
}
dev.off()

# convert the .png files to one .gif file using ImageMagick.
# The system() function executes the command as if it was done
# in the terminal. the -delay flag sets the time between showing
# the frames, i.e. the speed of the animation.
system("convert -delay 80 *.png example_1.gif")

# to not leave the directory with the single jpeg files
# I remove them.
file.remove(list.files(pattern=".png"))

