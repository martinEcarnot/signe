# Add new spectra
# SIGNE 2017
#

### DECOMMENTER POUR LA COMPARAISON DU MODELE AVEC DE NOUVELLES DONNEES ###
#
# # On refait l'AFD en groupant les jeux d'apprentissage et de validation
# maxi.id_final=round(mean(maxi_final[,1]))
# rlda_final=lda(scor[,1:maxi.id_final], class)
#
# q='/run/user/1000/gvfs/smb-share:server=stocka2,share=agap-ble/Ble/ASD/SIGNE/2017/20170711_TEST'
# sp2=SIGNE_load(q)
# save(sp2,file='/home/ecarnot/Documents/INRA/Projets/SIGNE/2017/20170711_TEST')
load('/home/ecarnot/Documents/INRA/Projets/SIGNE/2017/20170711_TEST')
row.names(sp2)=gsub("C1","C1-015",row.names(sp2))
row.names(sp2)=gsub("C2","C2-015",row.names(sp2))
row.names(sp2)=gsub("C3","C3-169",row.names(sp2))
row.names(sp2)=gsub("C4","C4-169",row.names(sp2))
row.names(sp2)=gsub("C5","C5-685",row.names(sp2))
row.names(sp2)=gsub("C6","C6-685",row.names(sp2))
row.names(sp2)=gsub("G1","G1-222",row.names(sp2))
row.names(sp2)=gsub("G2","G2-509",row.names(sp2))
row.names(sp2)=gsub("G3","G3-509",row.names(sp2))
row.names(sp2)=gsub("G4","G4-787",row.names(sp2))
row.names(sp2)=gsub("S1","S1-471",row.names(sp2))
row.names(sp2)=gsub("S2","S2-471",row.names(sp2))
row.names(sp2)=gsub("S3","S3-525",row.names(sp2))
row.names(sp2)=gsub("S4","S4-747",row.names(sp2))
row.names(sp2)=gsub("S5","S5-747",row.names(sp2))
row.names(sp2)=gsub("S6","S6-877",row.names(sp2))
row.names(sp2)=gsub("S7","S7-877",row.names(sp2))
# creation de la matrice de classes
class2=as.factor(substr(rownames(sp2),4,6))

# #ajout des classes par cepage sur une deuxieme colonne de la matrice
# nom_cep=c("C","C","G","S","G","S","C","S","G","S")
# class2=cbind(class2,class2)
# class2[,2]=gsub("10","S",class2[,2])
# class2[,2]=gsub("9","G",class2[,2])
# class2[,2]=gsub("8","S",class2[,2])
# class2[,2]=gsub("7","C",class2[,2])
# class2[,2]=gsub("6","S",class2[,2])
# class2[,2]=gsub("5","G",class2[,2])
# class2[,2]=gsub("4","S",class2[,2])
# class2[,2]=gsub("3","G",class2[,2])
# class2[,2]=gsub("2","C",class2[,2])
# class2[,2]=gsub("1","C",class2[,2])
# #switch classement par clone ([,1]) ou par cepage ([,2])
# class2=as.factor(class2[,2])
#
## Pretraitements
p=2
n=11
m=2
# # Ajustement des sauts de detecteurs (Montpellier: sauts ?? 1000 (651??me l.o.) et 1800 (1451))
sp2=adj_asd(sp2,c(651,1451))
# # Reduction des variables (extremites bruitees)
sp2=sp2[,seq(100,2100,1)]
# # SNV
sp2=t(scale(t(sp2)))
iout=which(sp2[,1]>0)
# # Derivation Savitsky Golay
sp2=t(apply(sp2,1,sgolayfilt,p=p,n=n,m=m))
#
sp2=sp2[-iout,]
class2=class2[-iout]
# # Prediction des noueaux spectres

# sp_cal_c=scale(sp_cal,center=rplsda$Xmeans,scale = F)
# sc_cal=sp_cal_c%*%rplsda$projection 
sp2c=scale(sp2,center=rplsdag$Xmeans,scale = F)
sc2=sp2c%*%rplsdag$projection  


pred2=list()
for (i in 2:ncmax) {
  # pred2[[i]]=predict(rldag[[i]] ,sc2[,1:i], method="plug-in")
  pred2[[i]]=SIGNE_maha(rldag[[i]],scg[,1:i],class,sc2[,1:i])
}

pred2[[1]]=NULL
ts=lapply(lapply(pred2,"[[","class"), class2, FUN = table)

# Matrice mauvais classements par clone
# mc[i,]=(rowSums((as.matrix(t)))-diag(as.matrix(t)))/rowSums((as.matrix(t)))
# Sumpred=lapply(ts, FUN = rowSums)
diags=lapply(ts, FUN = diag)
# Poucentage de bien classes
perok=100*unlist(lapply(diags, FUN = sum))/length(class)


# scor_sp2=predict(rpca,sp2)
# pred_sp2=predict(rlda_final,scor_sp2[,1:maxi.id_final])
# print(pred_sp2$class)
# maxr=apply(pred_sp2$posterior, 1, max)
# iok=which(maxr>0.85)
# t2=table(pred_sp2$class[iok],class2[iok])

# print(t2)
# perok2=100*sum(diag(as.matrix(t2)))/length(class2[iok])
# # Matrice mauvais classements par clone
# mc2=(rowSums((as.matrix(t2)))-diag(as.matrix(t2)))/rowSums((as.matrix(t2)))
# print("Bon classements de la parcelle d'essai et nombre de DV utilisees:")
print(perok)
# print(maxi.id_final)
# print("moyenne de mauvais classements par clone:")
# print(mc2)
### FIN COMPARAISON DU MODELE AVEC DE NOUVELLES DONNEES ###

stop()

# Plots

# Scores après LDA
cal=data.frame(predict(rldag[[15]] ,scg[,1:15]))
calmeans=data.frame(predict(rldag[[15]] ,rldag[[15]]$means))
test = data.frame(pred2[[15]])

p1 <- ggplot(test) + geom_point(aes(x.LD1, x.LD2, colour = class2), size = 2.5) +
  geom_point(data=cal,aes(x.LD1, x.LD2, colour = class), size = 4, shape="+") +
  geom_point(data=calmeans,aes(x.LD1, x.LD2, colour = rownames(rldag[[15]]$means)), size = 15, shape="+")
  # labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
  #      y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))
plot(p1)

# p2 <- ggplot(dataset) + geom_point(aes(pca.PC1, pca.PC2, colour = species, shape = species), size = 2.5) #+
  # labs(x = paste("PC1 (", percent(prop.pca[1]), ")", sep=""),
  #      y = paste("PC2 (", percent(prop.pca[2]), ")", sep=""))

# grid.arrange(p1, p2)

# ACP
pcag=prcomp(sp)
pcag$scale=FALSE
scpca2=predict(pcag,sp2)
cal=data.frame(pcag$x)
mcal0=aggregate(pcag$x, list(class), mean)
mcal=data.frame(aggregate(pcag$x, list(class), mean)[,-1])
mval=data.frame(aggregate(scpca2, list(class2), mean)[,-1])

# calmeans=data.frame(predict(rldag[[15]] ,rldag[[15]]$means))
test = data.frame(scpca2)

pcx=1
pcy=2

p1 <- ggplot(test) + geom_point(aes(test[,pcx], test[,pcy], colour = class2), size = 2.5) +
  geom_point(data=cal,aes(cal[,pcx], cal[,pcy], colour = class), size = 4, shape="+") +
geom_point(data=mcal,aes(mcal[,pcx], mcal[,pcy], colour = mcal0[,1]), size = 18, shape="+")+
geom_point(data=mval,aes(mval[,pcx], mval[,pcy], colour = mcal0[,1]), size = 15, shape="o")


plot(p1)


# Scores après PLSDA
cal=data.frame(scg)
calmeans=data.frame(predict(rldag[[15]] ,rldag[[15]]$means))
test = data.frame(pred2[[15]])

p1 <- ggplot(test) + geom_point(aes(x.LD1, x.LD2, colour = class2), size = 2.5) +
  geom_point(data=cal,aes(x.LD1, x.LD2, colour = class), size = 4, shape="+") +
  geom_point(data=calmeans,aes(x.LD1, x.LD2, colour = rownames(rldag[[15]]$means)), size = 15, shape="+")
# labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
#      y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))
plot(p1)

