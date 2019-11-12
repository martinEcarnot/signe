#NB : J'ai dû remplacer les "," par des "." dans le fichier d'origine.

data<-read.table("C:/Users/avitvale/Documents/Valentin Avit/Meteo_GDR.csv",sep=";",header=TRUE)

T0=10 #Zéro de végétation, caractéristique à chaque espèce. Pour la vigne, 10°C.


#Permet d'homogénéiser, afin que tous les nombres prennent 2 caractères (ce sera utile lorsqu'ils seront collés bout-à-bout pour former les rownames)
for (i in 1:length(data$JOUR)){
  if (data$JOUR[i] %in% 1:9){
    data$JOUR[i]=paste(0,data$JOUR[i], sep="")
  }
}

for (i in 1:length(data$MOIS)){
  if (data$MOIS[i] %in% 1:9){
    data$MOIS[i]=paste(0,data$MOIS[i], sep="")
  }
}






#On créé une data.frame avec une ligne par jour (=combinaison jour/mois/année unique)
D=unique(data[,2:4])

D2=as.data.frame(matrix(nrow=length(D[,1]), ncol=1), row.names=paste(D[,1],D[,2],D[,3], sep=""))
colnames(D2)="DJC"


#Remarque : les températures sont super élevées, quand même
D2[1,1]=0
compte=1
for (an in unique(data$AN)){
  dataan=data[which(data$AN==an),]
  for (mois in unique(dataan$MOIS)){
    datamois=dataan[which(dataan$MOIS==mois),]
    for (jour in unique(datamois$JOUR)){
      datajour=datamois[which(datamois$JOUR==jour),]
#      print(datajour)
      S=0
      compte=compte+1
      for (i in 1:length(datajour[,1])){
        s=max(0,as.numeric(as.character(datajour$T[i]))-10)
        # print(i)
        # print(datajour$T[i])
        # print(s)
        S=S+s
        # if (compte==30 & i ==20){ #pb à 68
        #   browser()
        # }
      }
      MJour=S/24
#      print(MJour)
#      print(compte)
      D2[compte,1]=as.numeric(D2[compte-1,1])+MJour
#      print(D2[compte,1])
    }
  }
}
mois=4
data[which(data$MOIS==mois),]





