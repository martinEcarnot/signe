#### Globalmatrix modifié et adapté à la situation.
## COMPRESSION DES SPECTRES EN UNE MATRICE ##

# nettoyage des variables de l'espace de travail
rm(list = ls())

source('C:/Users/avitvale/Documents/signeG/Script_R_2020/adj_asd.R')
source('C:/Users/avitvale/Documents/signeG/Script_R_2020/asd_read.R')
source('C:/Users/avitvale/Documents/signeG/Script_R_2020/SIGNE_load_modif3.R')

# source('C:/Users/avitvale/Documents/Script_R/SIGNE_load.R')       ##### N'arrive pas à trouver le chemin "normal (dans signeG), comprends pas pourquoi
# source('C:/Users/avitvale/Documents/Script_R/adj_asd.R')
# source('C:/Users/avitvale/Documents/Script_R/asd_read.R')


# #création de la matrice qui rassemble tous les dossiers:
globalmatrix=matrix(ncol=2072)

setwd("C:/Users/avitvale/Documents/SIGNE/Prediction_pheno/Feuilles/2019")
lieu="B"
for (m in 1:length(dir())){
  em=paste(dir()[m], sep="")
  w=SIGNE_load_modif3(em)
  rownames(w)=paste(dir()[m],rownames(w),lieu) #On ajoute l'identifiant lieu au nom de chaque spectre
  globalmatrix=rbind(globalmatrix, w)
}



globalmatrix=globalmatrix[complete.cases(globalmatrix),]
#str(globalmatrix)

rownames(globalmatrix)


# ## Filtrage des spectres aberrants        #Tous les spectres aberrants sont-ils filtrés par ce moyen ?
# globalmatrix2=globalmatrix[globalmatrix[,500]>0.6,] #Actuellement, 33 spectres jugés aberrants filtrés ici. (Dont 5 qui ne font pas partie de nos 10 clones étudiés..?)
# globalmatrix2=globalmatrix2[globalmatrix2[,1]<0.2,]
# globalmatrix2=globalmatrix2[globalmatrix2[,2000]<0.25,]


## Ajustement des sauts de detecteur (Montpellier: sauts ?? 1000 (651 eme l.o.) et 1800 (1451))
globalmatrix=adj_asd(globalmatrix,c(602,1402))

## Reduction des variables (extremites bruitees)
globalmatrix=globalmatrix[,seq(51,ncol(globalmatrix)-30,1)]

globalmatrix=globalmatrix[complete.cases(globalmatrix),]


globalmatrix=globalmatrix[globalmatrix[,1]<0.3,]

globalmatrixconservatoirefeuilles2019=globalmatrix





brb="C:/Users/avitvale/Documents/Test/"
save(globalmatrixconservatoirefeuilles2019, file=paste(brb,"globalmatrixconservatoirefeuilles2019",sep=""))
print(length(Truc))
# write.table(globalmatrix, file=paste(brb,"globalmatrix.csv",sep=""),sep=";", quote=FALSE)
setwd("C:/Users/avitvale/Documents/signeG")
### END ###


# gsub sert à remplacer un groupe de caractères par un autre dans une liste.
# gsub("A remplacer","Remplacement", Liste)

# sp=globalmatrix
#
# trad=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Correspondance_code_conservatoire_gamay.csv")
# # setdiff(donnees1$clone, trad$clone)
#
# rownames(sp)=gsub("2153", "ARB 7 29", rownames(sp))
# sp=sp[-grep("787", rownames(sp)),]
#
#
#
#
#
# trad=read_csv2(file = "C:/Users/avitvale/Documents/Valentin Avit/Correspondance_code_conservatoire_gamay.csv")
#
# #rownames(sp)=gsub(" ", "", rownames(sp))
#
# A=gsub("-.*", "", rownames(sp))
# B=gsub(".*(N )|.*(T )", "", A)
# C=gsub(" ", "", B)
# C=gsub("0","",C)
#
#
# trad2=gsub("/", "", trad$clone)
# trad2=gsub(" ", "", trad2)
# trad2=gsub("0", "", trad2)
#
# setdiff(D, trad2)
# setdiff(trad2,C)
#
# D=gsub("CJB827","CJP827",C)
# D=gsub("CJB858","CJP858",D)
# D=gsub("CJP635","CJB635",D)
# D=gsub("CJPR2113","CJPR213",D)
# D=gsub("CJP635","CJB635",D)
# D=gsub("COG771","CPG771",D)
# D[1536]="GAL243"
# D[grep("SMHH",D)]="SMHH837"
#
# rownames(sp)[grep("2 43", rownames(sp))]
# unique(D[grep("BEC",D)])
# D[1536]
# setdiff(trad2,C)[grep("CJPR", setdiff(trad2,C))]
#
#
# length(rownames(globalmatrix))
#
# rownames(globalmatrix)[-grep("[0-9]{8}(N|T) ([0-9]+|[A-Z|a-z]+) ?[0-9]{0,3} ?[0-9]{0,3}-[0-9]{1,3} B", rownames(globalmatrix))]
# #20190730N MPG 10 33 -1 B. ont tous un pb d'espace avant le -.
#
# rownames(globalmatrix)[-grep("[0-9]{8}(N|T) ([0-9]+|[A-Z|a-z]+) [0-9]{0,3} ?[0-9]{0,3}-[0-9]{2} B", rownames(globalmatrix))]
#
#
#
#
# rownames(globalmatrix)=sub("--", "-0", rownames(globalmatrix))
#
# anomalies=grep("[0-9]{8}(N|T) ([0-9]+|[A-Z|a-z]+) ?[0-9]{0,3} ?[0-9]{0,3}-00[0-9]{2} B", rownames(globalmatrix)) #Zéros en trop dans le code d'identification individu
# rownames(globalmatrix)[anomalies]=sub("-00", "-", rownames(globalmatrix)[anomalies])
#
# anomalies=-grep("[0-9]{8}(N|T) ([0-9]+|[A-Z|a-z]+) ?[0-9]{0,3} ?[0-9]{0,3}-[0-9]{2,3} B", rownames(globalmatrix)) #Espace devant le "-"
# rownames(globalmatrix)[anomalies]=sub(" -", "-", rownames(globalmatrix)[anomalies])
#
# anomalies=grep("[0-9]{8}(N|T) ([0-9]+|[A-Z|a-z]+) ?[0-9]{0,3} ?[0-9]{0,3}-[0-9]{3} B", rownames(globalmatrix))
# rownames(globalmatrix)[anomalies]=sub("-0", "-", rownames(globalmatrix)[anomalies])
#
#
#
#
#
#
#
#
# indetermines=-grep("[0-9]{8}(N|T) ([0-9]+|[A-Z|a-z]+) ?[0-9]{0,3} ?[0-9]{1,3}-[0-9]{2} B", rownames(globalmatrix)) #La question à 9 millions : qui sont-ce ?
#
# rownames(globalmatrix)[-indetermines][-grep("[0-9]{8}(N|T) ([0-9]+|[A-Z|a-z]+) [0-9]{1,3} [0-9]{1,3}-[0-9]{2} B", rownames(globalmatrix)[-indetermines])]
#
# relous=grep("(2153)|(787)", rownames(globalmatrix)[-indetermines]) #Contient des chiffres là où tous les autres contiennent des lettres + autres anomalies.
#
# rownames(globalmatrix)[-indetermines][-relous][-grep("[0-9]{8}(N|T) ([A-Z|a-z]+) ?[0-9]{1,3} ?[0-9]{1,3}-[0-9]{2} B", rownames(globalmatrix)[-indetermines][-relous])]
#
#
# ### J'en suis là.
#
# sprintf("%05.2f",1)
# sprintf("%05d",1)
#
