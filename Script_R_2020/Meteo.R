#NB : J'ai dû remplacer les "," par des "." dans le fichier d'origine.
rm(list = ls())
#data<-read.table("C:/Users/avitvale/Documents/Valentin Avit/Meteo_GDR.csv",sep=";",header=TRUE)
data<-read.table("C:/Users/avitvale/Documents/Valentin Avit/2019-90/Donnees_quotidiennes_1.csv",sep=";",header=TRUE)
debut=4 #On commence au mois d'Avril


# #1
# data=data[which(as.numeric(data$MOIS)>=debut),] #On retire les mois précédent Avril


#2
data=data[which(as.numeric(substr(data$DATE,5,6))>=debut),] #On retire les mois précédent Avril

T0=10 #Zéro de végétation, caractéristique à chaque espèce. Pour la vigne, 10°C.


# #1
# #Permet d'homogénéiser, afin que tous les nombres prennent 2 caractères (ce sera utile lorsqu'ils seront collés bout-à-bout pour former les rownames)
# for (i in 1:length(data$JOUR)){
#   if (data$JOUR[i] %in% 1:9){
#     data$JOUR[i]=paste(0,data$JOUR[i], sep="")
#   }
# }
#
# for (i in 1:length(data$MOIS)){
#   if (data$MOIS[i] %in% 1:9){
#     data$MOIS[i]=paste(0,data$MOIS[i], sep="")
#   }
# }
#
# #On créé une data.frame avec une ligne par jour (=combinaison jour/mois/année unique)
# D=unique(data[,2:4])
# D2=as.data.frame(matrix(nrow=length(D[,1]), ncol=1), row.names=paste(D[,1],D[,2],D[,3], sep=""))


#2
D=data[which(data$POSTE=="34154001"),]
D2=as.data.frame(matrix(nrow=length(D[,1]), ncol=1), row.names=D$DATE)


colnames(D2)="DJC"

D2[1,1]=0
compte=1


# #1
# for (an in unique(data$AN)){
#   dataan=data[which(data$AN==an),]
#
#
#   for (mois in unique(dataan$MOIS)){
#     datamois=dataan[which(dataan$MOIS==mois),]
#     for (jour in unique(datamois$JOUR)){
#       datajour=datamois[which(datamois$JOUR==jour),]
#       S=0
#       compte=compte+1
#       for (i in 1:length(datajour[,1])){
#         s=max(0,as.numeric(as.character(datajour$T[i]))-10)
#         S=S+s
#       }
#       MJour=S/24
#       if (jour=="01" & mois=="04"){
#         D2[compte-1,1]=0
#       }
#       D2[compte,1]=as.numeric(D2[compte-1,1])+MJour
#     }
#   }
# }


#2
for (i in 2:length(D[,1])){
  MJour=max(0,((D$TN[i]+D$TX[i])/2)-10)
  D2[i,1]=as.numeric(D2[i-1,1])+MJour
  if (substr(D$DATE[i],5,8)=="0401"){
    D2[i,1]=0 #On réinitialise d'une année à l'autre
  }
}



#2017
#"0524G" "0529G" "0606G" "0612G" "0619G" "0626G" "0703G" "0710G" "0711g" "0717G" "0724G" "0731G"
L1=c("20170524", "20170529", "20170606", "20170612", "20170619", "20170626", "20170703", "20170710", "20170711", "20170717", "20170724", "20170731")

#2018 :
#"0619G" "0627G" "0704G" "0709G" "0710g" "0731A" "0724B" "0810B" "0816G" "0817g" "0823A"
L2=c("20180619", "20180627", "20180704", "20180709", "20180710", "20180731", "20180724", "20180810", "20180816", "20180817", "20180823")

#2019 :
#"0613G" "0617G" "0624G" "0628g" "0702G" "0703A" "0710G" "0718G" "0723A" "0726G" "0730G" "0822B"
L3=c("20190613", "20190617", "20190624", "20190628", "20190702", "20190703", "20190710", "20190718", "20190723", "20190726", "20190730", "20190822")

for (i in 1:length(L3)){
  print(D2[which(rownames(D2)==L3[i]),])
}


D2[which(rownames(D2)==L1[4]),]
D2[which(rownames(D2)==L2[2]),]
D2[which(rownames(D2)==L2[1]),]

# D2[which(rownames(D2)=="20170524"),]
# D2[which(rownames(D2)=="20170529"),]
# D2[which(rownames(D2)=="20170606"),]
# D2[which(rownames(D2)=="20170612"),]
# D2[which(rownames(D2)=="20170619"),]
# D2[which(rownames(D2)=="20170626"),]
# D2[which(rownames(D2)=="20170703"),]
# D2[which(rownames(D2)=="20170710"),]
# D2[which(rownames(D2)=="20170711"),]
# D2[which(rownames(D2)=="20170717"),]
# D2[which(rownames(D2)=="20170724"),]
# D2[which(rownames(D2)=="20170731"),]
#
# D2[which(rownames(D2)=="20180619"),]
# D2[which(rownames(D2)=="20180627"),]
# D2[which(rownames(D2)=="20180704"),]
# D2[which(rownames(D2)=="20180709"),]
# D2[which(rownames(D2)=="20180710"),]
# D2[which(rownames(D2)=="20180731"),]
# D2[which(rownames(D2)=="20180724"),]
# D2[which(rownames(D2)=="20180810"),]
# D2[which(rownames(D2)=="20180816"),]
# D2[which(rownames(D2)=="20180817"),]
# D2[which(rownames(D2)=="20180823"),]
#
# D2[which(rownames(D2)=="20190613"),]
# D2[which(rownames(D2)=="20190617"),]
# D2[which(rownames(D2)=="20190624"),]
# D2[which(rownames(D2)=="20190628"),]
# D2[which(rownames(D2)=="20190702"),]
# D2[which(rownames(D2)=="20190703"),]
# D2[which(rownames(D2)=="20190710"),]
# D2[which(rownames(D2)=="20190718"),]
# D2[which(rownames(D2)=="20190723"),]
# D2[which(rownames(D2)=="20190726"),]
# D2[which(rownames(D2)=="20190730"),]
# D2[which(rownames(D2)=="20190822"),]
