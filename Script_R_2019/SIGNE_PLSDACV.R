SIGNE_PLSDACV <- function(spcal,ncmax,repet,k)
{

  ## Definition des matrices de resultat final
  # Creation de la matrice des perok finale
  perok_final=matrix(nrow = repet, ncol = 1)
  perok_finalm0=matrix(nrow = repet, ncol = ncmax)

  ## Creation matrice de % de mauvais classements par clone
  mc=matrix(nrow = ncmax,ncol = c)

  ## Creation de la matrice des VL et perok maximaux
  maxi_final=matrix(nrow= repet, ncol = 2)

  ## Creation de la matrice de % de mauvais classements
  mc_final=matrix(nrow= repet, ncol = length(levels(spcal$y)))
  ## Creation d'un matrice cubique pour enregistrer les tables de contingence
  t_final=array(dim=c(c,c,repet))
  ## Noms des colonnes et des lignes
  colnames(t_final)=c(basename(levels(spcal$y)))
  rownames(t_final)=c(basename(levels(spcal$y)))
  colnames(maxi_final)= c("maxi.id","perok max")
  colnames(mc_final)= c(basename(levels(spcal$y)))

  ###s?paration validation calibration PLSDA###
  #set.seed(1) # fixe le tirage aleatoire
  for(j in 1:repet) {

    predm0=as.data.frame(matrix(nrow = nrow(spcal), ncol = ncmax))
    spcaldef=spcal # spcal deflaté du(des) groupe(s) de CV déjà validés
    ## Boucle CV
    for (i in 1:k) {
      ndc=length(unique(spcaldef$datclone))
      m=mstage(spcaldef,stage=list("cluster","cluster"), varnames=list("datclone","souche"),size=list(ndc,rep(1,30)))
      spvalCV=getdata(spcaldef,m)[[2]]

      idvalCV =which(rownames(spcal)  %in%  rownames(spvalCV))

      spcaldef=spcaldef[-(which(rownames(spcaldef)  %in%  rownames(spvalCV))),]

      # spvalCV=sp_cal[idvalCV,]       # matrice du jeu de validation
      # classvalCV=classcal[idvalCV]  #identifiants des classes du jeu de validation
      spcalCV=spcal[-idvalCV,]      #matrice du jeu de calibration compos?e de tout ce qui n'est pas en validation
      # classcalCV=classcal[-idvalCV] #identifiants des classes du jeu de calibration

      ## PLSDA and application to have loadings and scores
      rplsda=caret::plsda(spcalCV$x, spcalCV$y,ncomp=ncmax)
      sccalCV=rplsda$scores
      spvalCV_c=scale(spvalCV$x,center=rplsda$Xmeans,scale = F)
      scvalCV=spvalCV_c%*%rplsda$projection  # score_val=predict(rplsda,sc_val,type="scores") : ne marche pas

      for (ii in 2:ncmax) {

        ## Validation
        predm0[idvalCV,ii]=SIGNE_maha0(sccalCV[,1:ii], spcalCV$y, scvalCV[,1:ii])$class
      }
    }

    ## Table de contingence CV
    tsm0=lapply(as.list(predm0), spcal$y, FUN = table)

    ## Matrice mauvais classements par clone CV
    diagsm0=lapply(tsm0, FUN = diag)

    ## Pourcentage de bien classes CV
    perokm0 =100*unlist(lapply(diagsm0, FUN = sum))/nrow(spcal)

    ## Pourcentage de bien classes CV
    maxi=max(perokm0)
    maxi.id=which.max(perokm0)

    ### Enregistrement des matrices de resultat final
    ##Remplissage de la matrice des perok finale
    perok_finalm0[j,]=perokm0

    ## Remplissage de la VL max et de son % de bons classements globaux
    maxi_final[j,1]=maxi.id
    maxi_final[j,2]=maxi
    ## Remplissage de la matrice de mauvais classements par clone
    mc_final[j,]=mc[maxi.id,]
  }

  return(perok_finalm0)

}
