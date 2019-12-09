## COMPRESSION DES SPECTRES EN UNE MATRICE ##

# nettoyage des variables de l'espace de travail
rm(list = ls())

source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')       ##### N'arrive pas à trouver le chemin "normal (dans signeG), comprends pas pourquoi
source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')
source('C:/Users/avitvale/Documents/Script_R/asd_read.R')


# #création de la matrice qui rassemble tous les dossiers:
globalmatrix=matrix(ncol=2072)

setwd("C:/Users/avitvale/Documents/SIGNE")
for (a in 2017:2019){
  setwd("C:/Users/avitvale/Documents/SIGNE")                                     #Mettre ici le chemin d'accès au dossier contenant les spectres
  setwd(as.character(a)) #On balaie les trois dossiers année
  for (n in 1:length(dir())){
    if ((dir()[n])=="Aude" | (dir()[n])=="Beaujolais" | (dir()[n])=="GDRautre" |(dir()[n])=="GDR96"){
      lieu=substr(dir()[n],1,1) #On attribue au spectre un identifiant lieu selon sa parcelle d'origine : A, B, G ou g
      if ((dir()[n])=="GDRautre"){
        lieu='g'
      }
      setwd(as.character(dir()[n])) #On entre dans le sous-dossier.
      for (m in 1:length(dir())){
        if(substr(dir()[m],9,9)=="P"){ #on ne s'intéresse qu'au mesures effectuées en suivant la modalité P
          em=paste(dir()[m], sep="")
          w=SIGNE_load(em)
          rownames(w)=paste(dir()[m],rownames(w),lieu) #On ajoute l'identifiant lieu au nom de chaque spectre
          globalmatrix=rbind(globalmatrix, w)
        }
      }
      setwd("../")
    }
  }
}



globalmatrix=globalmatrix[complete.cases(globalmatrix),]
#str(globalmatrix)

#Filtrage des noms non standardisés
iout=which(nchar(rownames(globalmatrix))>18)
# rownames(globalmatrix[iout,])
# unique(rownames(globalmatrix[which(substr(rownames(globalmatrix),1,8)=="20170612"),]))
globalmatrix <- globalmatrix[-iout,]

## Filtrage des spectres aberrants        #Tous les spectres aberrants sont-ils filtrés par ce moyen ?
globalmatrix=globalmatrix[globalmatrix[,500]>0.6,] #Actuellement, 33 spectres jugés aberrants filtrés ici. (Dont 5 qui ne font pas partie de nos 10 clones étudiés..?)
globalmatrix=globalmatrix[globalmatrix[,1]<0.2,]
globalmatrix=globalmatrix[globalmatrix[,2000]<0.25,]


## Ajustement des sauts de detecteur (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
########globalmatrix=adj_asd(globalmatrix,c(602,1402))

## Reduction des variables (extremites bruitees)
globalmatrix=globalmatrix[,seq(51,ncol(globalmatrix)-30,1)]



# ###Ajout du cepage au nom de la ligne. Retire aussi les cépages/clones qui ne sont pas les 3/10 notres.
# ##Filtre en fonction du cepage
# titre=rownames(globalmatrix)
#
# ya=(substr(titre,11,13))== "015"
# yb=(substr(titre,11,13))== "169"
# yc=(substr(titre,11,13))== "685"
#
# yd=(substr(titre,11,13))== "222"
# ye=(substr(titre,11,13))== "509"
# yf=(substr(titre,11,13))== "787"
#
# yg=(substr(titre,11,13))== "471"
# yh=(substr(titre,11,13))== "525"
# yi=(substr(titre,11,13))== "747"
# yj=(substr(titre,11,13))== "877"
#
# ##Separation en 3 matrices
# spa=globalmatrix[ya,]
# spb=globalmatrix[yb,]
# spc=globalmatrix[yc,]
# spd=globalmatrix[yd,]
# spe=globalmatrix[ye,]
# spf=globalmatrix[yf,]
# spg=globalmatrix[yg,]
# sph=globalmatrix[yh,]
# spi=globalmatrix[yi,]
# spj=globalmatrix[yj,]
#
# ##Rajoute le nom du cepage au nom de la ligne à la place de P/T/N
# substr(rownames(spa),9,9)="C"
# substr(rownames(spb),9,9)="C"
# substr(rownames(spc),9,9)="C"
#
# substr(rownames(spd),9,9)="G"
# substr(rownames(spe),9,9)="G"
# substr(rownames(spf),9,9)="G"
#
# substr(rownames(spg),9,9)="S"
# substr(rownames(sph),9,9)="S"
# substr(rownames(spi),9,9)="S"
# substr(rownames(spj),9,9)="S"
#
# ##Recombine les 3 matrices pour reformer sp
# globalmatrix=rbind(spa,spb,spc,spd,spe,spf,spg,sph,spi,spj)


#Cette méthode est équivalente à ce qui précède, simplement les spectres dans la globalmatrix ne sont pas dans le même ordre.
for (i in 1:length(globalmatrix[,1])){
  if (substr(rownames(globalmatrix)[i],11,13)=="015" | substr(rownames(globalmatrix)[i],11,13)=="169" | substr(rownames(globalmatrix)[i],11,13)=="685"){
    substr(rownames(globalmatrix)[i],9,9)="C"
  }
  else if  (substr(rownames(globalmatrix)[i],11,13)=="222" | substr(rownames(globalmatrix)[i],11,13)=="509" | substr(rownames(globalmatrix)[i],11,13)=="787"){
    substr(rownames(globalmatrix)[i],9,9)="G"
  }
  else if  (substr(rownames(globalmatrix)[i],11,13)=="471" | substr(rownames(globalmatrix)[i],11,13)=="525" | substr(rownames(globalmatrix)[i],11,13)=="747" | substr(rownames(globalmatrix)[i],11,13)=="877"){
    substr(rownames(globalmatrix)[i],9,9)="S"
  }
  else {
    substr(rownames(globalmatrix)[i],9,9)="A"
  }
}

iout=which(substr(rownames(globalmatrix),9,9)=="A")
globalmatrix <- globalmatrix[-iout,]



#
# str(globalmatrix)
# a=matrix(ncol=1992)
# b=matrix(ncol=1992)
# c=matrix(ncol=1992)
# d=matrix(ncol=1992)
# e=matrix(ncol=1992)
# f=matrix(ncol=1992)
# g=matrix(ncol=1992)
# h=matrix(ncol=1992)
# i=matrix(ncol=1992)
# j=matrix(ncol=1992)
# k=matrix(ncol=1992)
# l=matrix(ncol=1992)
# m=matrix(ncol=1992)
# n=matrix(ncol=1992)
# o=matrix(ncol=1992)
# p=matrix(ncol=1992)
# q=matrix(ncol=1992)
# r=matrix(ncol=1992)
#
# colnames(a)=449:2440
# colnames(b)=449:2440
# colnames(c)=449:2440
# colnames(d)=449:2440
# colnames(e)=449:2440
# colnames(f)=449:2440
# colnames(g)=449:2440
# colnames(h)=449:2440
# colnames(i)=449:2440
# colnames(j)=449:2440
# colnames(k)=449:2440
# colnames(l)=449:2440
# colnames(m)=449:2440
# colnames(n)=449:2440
# colnames(o)=449:2440
# colnames(p)=449:2440
# colnames(q)=449:2440
# colnames(r)=449:2440
#
# # interm=as.data.frame(globalmatrix[which(substr(rownames(globalmatrix),1,13)==unique(substr(rownames(globalmatrix),1,13)[1])),])
# # if (length(interm[,1])==18){
# #   a=rbind(a,interm[1,])
# #   rownames(a)[length(a[,1])]=paste(substr(rownames(interm)[length(a[,1])],1,13),"-00 ", substr(rownames(interm)[length(a[,1])],18,18),1, sep="")
# #   # b=rbind(b,interm[2,])
# #   # c=rbind(c,interm[3,])
# #   # d=rbind(d,interm[4,])
# #   # e=rbind(e,interm[5,])
# #   # f=rbind(f,interm[6,])
# #   # g=rbind(g,interm[7,])
# #   # h=rbind(h,interm[8,])
# #   # i=rbind(i,interm[9,])
# #   # j=rbind(j,interm[10,])
# #   # k=rbind(k,interm[11,])
# #   # l=rbind(l,interm[12,])
# #   # m=rbind(m,interm[13,])
# #   # n=rbind(n,interm[14,])
# #   # o=rbind(o,interm[15,])
# #   # p=rbind(p,interm[16,])
# #   # q=rbind(q,interm[17,])
# #   # r=rbind(r,interm[18,])
# # }
#
#
#
# for (z in 1:length(unique(substr(rownames(globalmatrix),1,13)))) {
#   interm=globalmatrix[which(substr(rownames(globalmatrix),1,13)==unique(substr(rownames(globalmatrix),1,13))[z]),]
#
# #  truc1=as.data.frame(globalmatrix[which(substr(rownames(globalmatrix),1,13)==unique(substr(rownames(globalmatrix),1,13))[1]),])
# #  truc2=globalmatrix[which(substr(rownames(globalmatrix),1,13)==unique(substr(rownames(globalmatrix),1,13))[1]),]
#   if (length(interm[,1])==18){
#     print(i)
#     a=rbind(a,interm[1,])
#     K=as.character(paste(substr(rownames(interm)[1],1,13),"-01 ", substr(rownames(interm)[1],18,18), sep=""))
#     rownames(a)[length(a[,1])]=K
#
#     b=rbind(b,interm[2,])
#     K=as.character(paste(substr(rownames(interm)[2],1,13),"-01 ", substr(rownames(interm)[2],18,18), sep=""))
#     rownames(b)[length(b[,1])]=K
#
#     c=rbind(c,interm[3,])
#     K=as.character(paste(substr(rownames(interm)[3],1,13),"-01 ", substr(rownames(interm)[3],18,18), sep=""))
#     rownames(c)[length(c[,1])]=K
#
#     d=rbind(d,interm[4,])
#     K=as.character(paste(substr(rownames(interm)[4],1,13),"-01 ", substr(rownames(interm)[4],18,18), sep=""))
#     rownames(d)[length(d[,1])]=K
#
#     e=rbind(e,interm[5,])
#     K=as.character(paste(substr(rownames(interm)[5],1,13),"-01 ", substr(rownames(interm)[5],18,18), sep=""))
#     rownames(e)[length(e[,1])]=K
#
#     f=rbind(f,interm[6,])
#     K=as.character(paste(substr(rownames(interm)[6],1,13),"-01 ", substr(rownames(interm)[6],18,18), sep=""))
#     rownames(f)[length(f[,1])]=K
#
#     g=rbind(g,interm[7,])
#     K=as.character(paste(substr(rownames(interm)[7],1,13),"-07 ", substr(rownames(interm)[7],18,18), sep=""))
#     rownames(g)[length(g[,1])]=K
#
#     h=rbind(h,interm[8,])
#     K=as.character(paste(substr(rownames(interm)[8],1,13),"-07 ", substr(rownames(interm)[8],18,18), sep=""))
#     rownames(h)[length(h[,1])]=K
#
#     i=rbind(i,interm[9,])
#     K=as.character(paste(substr(rownames(interm)[9],1,13),"-07 ", substr(rownames(interm)[9],18,18), sep=""))
#     rownames(i)[length(i[,1])]=K
#
#     j=rbind(j,interm[10,])
#     K=as.character(paste(substr(rownames(interm)[10],1,13),"-07 ", substr(rownames(interm)[10],18,18), sep=""))
#     rownames(j)[length(j[,1])]=K
#
#     k=rbind(k,interm[11,])
#     K=as.character(paste(substr(rownames(interm)[11],1,13),"-07 ", substr(rownames(interm)[11],18,18), sep=""))
#     rownames(k)[length(k[,1])]=K
#
#     l=rbind(l,interm[12,])
#     K=as.character(paste(substr(rownames(interm)[12],1,13),"-07 ", substr(rownames(interm)[12],18,18), sep=""))
#     rownames(l)[length(l[,1])]=K
#
#     m=rbind(m,interm[13,])
#     K=as.character(paste(substr(rownames(interm)[13],1,13),"-13 ", substr(rownames(interm)[13],18,18), sep=""))
#     rownames(m)[length(m[,1])]=K
#
#     n=rbind(n,interm[14,])
#     K=as.character(paste(substr(rownames(interm)[14],1,13),"-13 ", substr(rownames(interm)[14],18,18), sep=""))
#     rownames(n)[length(n[,1])]=K
#
#     o=rbind(o,interm[15,])
#     K=as.character(paste(substr(rownames(interm)[15],1,13),"-13 ", substr(rownames(interm)[15],18,18), sep=""))
#     rownames(o)[length(o[,1])]=K
#
#     p=rbind(p,interm[16,])
#     K=as.character(paste(substr(rownames(interm)[16],1,13),"-13 ", substr(rownames(interm)[16],18,18), sep=""))
#     rownames(p)[length(p[,1])]=K
#
#     q=rbind(q,interm[17,])
#     K=as.character(paste(substr(rownames(interm)[17],1,13),"-13 ", substr(rownames(interm)[17],18,18), sep=""))
#     rownames(q)[length(q[,1])]=K
#
#     r=rbind(r,interm[18,])
#     K=as.character(paste(substr(rownames(interm)[18],1,13),"-13 ", substr(rownames(interm)[18],18,18), sep=""))
#     rownames(r)[length(r[,1])]=K
#   }
# }
#
#
#
#
#
# colnames(b)=(max(as.numeric(colnames(a)))+1):((max(as.numeric(colnames(a)))+1)+length(b[1,])-1)
# colnames(c)=(max(as.numeric(colnames(b)))+1):((max(as.numeric(colnames(b)))+1)+length(c[1,])-1)
# colnames(d)=(max(as.numeric(colnames(c)))+1):((max(as.numeric(colnames(c)))+1)+length(d[1,])-1)
# colnames(e)=(max(as.numeric(colnames(d)))+1):((max(as.numeric(colnames(d)))+1)+length(e[1,])-1)
# colnames(f)=(max(as.numeric(colnames(e)))+1):((max(as.numeric(colnames(e)))+1)+length(f[1,])-1)
#
# #colnames(g)=(max(as.numeric(colnames(f)))+1):((max(as.numeric(colnames(f)))+1)+length(g[1,])-1)
# colnames(h)=(max(as.numeric(colnames(g)))+1):((max(as.numeric(colnames(g)))+1)+length(h[1,])-1)
# colnames(i)=(max(as.numeric(colnames(h)))+1):((max(as.numeric(colnames(h)))+1)+length(i[1,])-1)
# colnames(j)=(max(as.numeric(colnames(i)))+1):((max(as.numeric(colnames(i)))+1)+length(j[1,])-1)
# colnames(k)=(max(as.numeric(colnames(j)))+1):((max(as.numeric(colnames(j)))+1)+length(k[1,])-1)
# colnames(l)=(max(as.numeric(colnames(k)))+1):((max(as.numeric(colnames(k)))+1)+length(l[1,])-1)
#
# #colnames(m)=(max(as.numeric(colnames(l)))+1):((max(as.numeric(colnames(l)))+1)+length(m[1,])-1)
# colnames(n)=(max(as.numeric(colnames(m)))+1):((max(as.numeric(colnames(m)))+1)+length(n[1,])-1)
# colnames(o)=(max(as.numeric(colnames(n)))+1):((max(as.numeric(colnames(n)))+1)+length(o[1,])-1)
# colnames(p)=(max(as.numeric(colnames(o)))+1):((max(as.numeric(colnames(o)))+1)+length(p[1,])-1)
# colnames(q)=(max(as.numeric(colnames(p)))+1):((max(as.numeric(colnames(p)))+1)+length(q[1,])-1)
# colnames(r)=(max(as.numeric(colnames(q)))+1):((max(as.numeric(colnames(q)))+1)+length(r[1,])-1)
#
#
# sptesta=cbind(a,b,c,d,e,f)
# sptestb=cbind(g,h,i,j,k,l)
# sptestc=cbind(m,n,o,p,q,r)
#
# sptest=rbind(sptesta,sptestb,sptestc)
#
#
#
# str(sptest)
# str(globalmatrix)

#
# interm=as.data.frame(globalmatrix[which(substr(rownames(globalmatrix),1,13)==unique(substr(rownames(globalmatrix),1,13)[1])   &   as.numeric(substr(rownames(globalmatrix),15,16)) %in% ((6*0)+1):(6*(0+1))     ),])
# if (length(interm[,1])==6){
#   #    I=interm[1,] #cbind(a=interm[1,], b=interm[2,], c=interm[3,], d=interm[4,], e=interm[5,], f=interm[6,])
#   #    for(k in 1:6){
#   #    rownames(interm)[1]=paste(substr(rownames(interm)[1],1,13),substr(rownames(interm)[1],18,18),j, sep="")
#   I=cbind(interm[1,], interm[2,], interm[3,], interm[4,], interm[5,], interm[6,])
#   spb=I
# }
#
#
# for (j in 1:2){
#   interm=as.data.frame(globalmatrix[which(substr(rownames(globalmatrix),1,13)==unique(substr(rownames(globalmatrix),1,13)[1])   &   as.numeric(substr(rownames(globalmatrix),15,16)) %in% ((6*j)+1):(6*(j+1))     ),])
#   if (length(interm[,1])==6){
# #    I=interm[1,] #cbind(a=interm[1,], b=interm[2,], c=interm[3,], d=interm[4,], e=interm[5,], f=interm[6,])
# #    for(k in 1:6){
# #    rownames(interm)[1]=paste(substr(rownames(interm)[1],1,13),substr(rownames(interm)[1],18,18),j, sep="")
#     I=cbind(interm[1,], interm[2,], interm[3,], interm[4,], interm[5,], interm[6,])
#     spb=rbind(spb,I)
#   }
# }




#
# for (i in 2:length(unique(substr(rownames(globalmatrix),1,13)))) {
#   for (j in 0:2){
#     interm=as.data.frame(globalmatrix[which(substr(rownames(globalmatrix),1,13)==unique(substr(rownames(globalmatrix),1,13)[i])   &   as.numeric(substr(rownames(globalmatrix),15,16)) %in% ((6*j)+1):(6*(j+1))     ),])
#     print(length(interm[,1]))
#     if (length(interm[,1])==6){           #6
#       I=cbind(interm[1,], interm[2,], interm[3,], interm[4,], interm[5,], interm[6,])
#       spb=rbind(spb,I)
#     }
#   }
# }
# str(spb) #Trop de lignes, non ?
# # L[which(L!=6)]
# # which(L!=6)
#
#
# interm=sp[which(paste(substr(rownames(sp),1,13),sp$souche)==unique(paste(substr(rownames(sp),1,13),sp$souche))[10]),]
#
# length(interm[,1])
#
#

# globalmatrix2=sptest
# globalmatrix2=globalmatrix2[complete.cases(globalmatrix2),]



str(globalmatrix)
# str(globalmatrix2)


brb="C:/Users/avitvale/Documents/Test/"
save(globalmatrix, file=paste(brb,"globalmatrix",sep=""))
print(length(globalmatrix))
# write.table(globalmatrix, file=paste(brb,"globalmatrix.csv",sep=""),sep=";", quote=FALSE)
setwd("C:/Users/avitvale/Documents/signeG")
### END ###

