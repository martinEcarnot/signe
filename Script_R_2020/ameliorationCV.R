test = as.matrix(tsm0[["V9"]])
setwd("C:/Users/avitvale/Documents")
write.table(test, file = "test.txt", sep = "\t", col.names = F)

# ## creer les colonnes claires pour filtrer---------------------------
# library(tidyverse)
# data_names = as.data.frame(titre)
# data_names$clone = str_sub(data_names$titre, start = 11L, end = 13L)
# data_names$num_mesure = str_sub(data_names$titre, start = 15L, end = 16L)
# data_names = data_names %>% mutate(souche = ifelse(num_mesure == "01"|num_mesure == "02"|num_mesure == "03"|num_mesure == "04"|num_mesure == "05"|num_mesure == "06", 1, 
#                                                     ifelse(num_mesure == "07"|num_mesure == "08"|num_mesure == "09"|num_mesure == "10"|num_mesure == "11"|num_mesure == "12", 2, 3)))
# # test = as.data.frame(sp, row.names = rownames(sp))
# # test = cbind(data_names, test)
# 
# filtre_clone = as.factor(data_names$clone)
# filtre_souche = as.factor(data_names$souche)
# filtres = as.list(filtre_souche, filtre_clone)
# test1 = split(data_names, filtres)


                