# Pour le calcul de SD et OD

# la fonction odis ne fonctionne pas avec les sorties de plsda (ou plsdalm)
# On utilise plsr en créant un tableau disjonctif avec la fonction dummy.
# Les od et sd sont calculé pour le nombre max de VL. Pour le faire varier, il faut refaire la plsr avec un nombre différent de composantes

rplsda=rnirs::plsr(spcal$x,dummy(classcal),spval$x,dummy(classval),ncmax)
od=odis(rpls,sp$x,spu$x)$du$d
sd=scordis(rpls)$du$d
