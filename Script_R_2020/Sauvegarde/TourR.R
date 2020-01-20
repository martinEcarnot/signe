library(tourr)
library(rnirs)
library(plotly)


idval=which(substr(rownames(sp),1,4)=="2019")
spval=sp[idval,]
# m=mstage(sp,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,ndc)))
# spval=getdata(sp,m)[[2]]
# idval=which(rownames(sp)  %in%  rownames(spval))

spval=sp[idval,]
spcal=sp[-idval,]
classval=class[idval]
classcal=class[-idval]
classvalclo=classclo[idval]
classcalclo=classclo[-idval]


predmF=as.data.frame(matrix(nrow = length(classval), ncol = ncmax))

rplsda=caret::plsda(spcal$x, classcal,ncomp=ncmax)
sccal=rplsda$scores
spval_c=scale(spval$x,center=rplsda$Xmeans,scale = F)
scval=spval_c%*%rplsda$projection

class(sccal)="matrix"
sc=rbind(sccal,scval)




# rplsda=caret::plsda(sp$x, class,ncomp=15)
# class(rplsda$scores)="matrix"
annee=substr(rownames(sp),1,4)
cepage=substr(rownames(sp),9,9)
cepageannee=paste(substr(rownames(sp),1,4),substr(rownames(sp),9,9))
Test=data.frame(sc[,1:15], cepageannee)


#CEPAGEANNEE=c("2017 C", "2018 C", "2019 C", "2017 G", "2018 G", "2019 G", "2017 S", "2018 S", "2019 S")
  CEPAGEANNEE=c("2017 G", "2018 G", "2019 G", "2017 S", "2018 S")

Test2=Test[which(Test$cepageannee %in% CEPAGEANNEE),]
Test2$cepageannee=droplevels(Test2$cepageannee)


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
colocepage=c(rouge, bleu, vert)
colocepageannee=c("#FF5B48", "#78C0FE", "#01FF7C", "#FA2700", "#6451DD", "#33CA22", "#C60803", "#190396", "#0A6C00")
#colocepageannee=c("#FF5B48", "#78C0FE", "#01FF7C", "#AC0F00", "#6451DD", "#33CA22")
colocepageannee2B=c("#78C0FE", "#6451DD", "#190396")
colocepageannee2=c("#78C0FE", "#01FF7C", "#6451DD", "#33CA22", "#190396")

colocepageannee2R=c("#FF5B48", "#FA2700", "#C60803")
colocepageannee2V=c("#01FF7C", "#33CA22", "#0A6C00")


#Violets #581845 #900C3F #C70039
#Bleus #132959 #234CA5 #457DBB
#Oranges #FF5733  #FFC30F #FFD940

#         Sombre      Moyen                   Clair
#Rouges : #C60803     #FA2700 / #AC0F00       #FF5B48
#Bleus : #190396      #6451DD                 #78C0FE
#Verts : #0A6C00      #33CA22                 #01FF7C

if(Sys.getenv("RSTUDIO") == "1" & # check if running in RStudio
   .Platform$OS.type == "unix") quartz() else X11()
animate_xy(Test2[, 1:8], col=colocepageannee2[Test2[, 16]])





# adapted from  https://github.com/rstudio/ggvis/blob/master/demo/tourr.r


mat <- rescale(as.matrix(flea[1:6]))
tour <- new_tour(mat, grand_tour(), NULL)

tour_dat <- function(step_size) {
  step <- tour(step_size)
  proj <- center(mat %*% step$proj)
  data.frame(x = proj[,1], y = proj[,2],
             species = flea$species)
}

proj_dat <- function(step_size) {
  step <- tour(step_size)
  data.frame(
    x = step$proj[,1], y = step$proj[,2], measure = colnames(mat)
  )
}

steps <- c(0, rep(1/15, 50))
stepz <- cumsum(steps)

# tidy version of tour data
tour_dats <- lapply(steps, tour_dat)
tour_datz <- Map(function(x, y) cbind(x, step = y), tour_dats, stepz)
tour_dat <- dplyr::bind_rows(tour_datz)

# tidy version of tour projection data
proj_dats <- lapply(steps, proj_dat)
proj_datz <- Map(function(x, y) cbind(x, step = y), proj_dats, stepz)
proj_dat <- dplyr::bind_rows(proj_datz)


ax <- list(
  title = "",
  range = c(-1, 1),
  zeroline = FALSE
)

# for nicely formatted slider labels
options(digits = 2)

proj_dat %>%
  plot_ly(x = ~x, y = ~y, frame = ~step, color = I("gray80")) %>%
  add_segments(xend = 0, yend = 0) %>%
  add_text(text = ~measure) %>%
  add_markers(color = ~species, data = tour_dat) %>%
  hide_legend() %>%
  layout(xaxis = ax, yaxis = ax) %>%
  animation_opts(33, redraw = FALSE)

