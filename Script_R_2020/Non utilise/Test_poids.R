Test_poids <- read.table("C:/Users/Noémie/Desktop/Test_poids.csv", 
                            header=TRUE, sep=";",dec=".", check.names=FALSE, 
                            fileEncoding="latin1")

Esca= Test_poids[rep(row.names(Test_poids), Test_poids$Poids),]


enr="C:\\Users\\Noémie\\Desktop\\Test_poids_ok\\"
write.table(Esca, file=paste(enr,"Esca.csv",sep=""),sep=";", quote=FALSE)