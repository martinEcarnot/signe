adj_asd <- function(Xi,iadj)
{
# Elimine les sauts dus au passage d'un detecteur  a un autre.
# On extrapole l'ajustement lineaire fait sur les [ws] points avant le saut, sur le point suivant
# La difference entre l'ancien et le nouveau "point suivant" est soustraite a tous les points apres le "point suivant"
#[651,1451]
ws=5
#browser()
for (i in 1:length(iadj)) {

	x= (iadj[i] - ws + 1):iadj[i]
	Y=Xi[, x]
	my=t(colMeans(t(Y)))

	sx=var(x) 
	mx=mean(x)
	b=cov(x,t(Y))/sx # Ajustement lineaire.
	b0=my-b*mx  # Ajustement lineaire.
	dif=Xi[,iadj[i]+1]-(b0+b*(iadj[i]+1)) 
	loremp=(iadj[i]+1):dim(Xi)[2]
	Xi[,loremp]=Xi[,loremp]-kronecker(matrix(1,1,dim(Xi)[2]-iadj[i]),t(dif))  # On utilise kronecker pour remplacer repmat
}
Xo=Xi
return(Xo)

}
