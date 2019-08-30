
# Librairies -----------------------------------------------------------

library(data.table)
library(ggplot2)
# library(highcharter)
library(reshape2)
library(forecast)
library(gridExtra)
library(dplyr)
library(ggfortify)
library(plotly)
# library(cluster)
library(moments)
# library(NMOF)
library(tseries)
library(MASS)
# library(packagefinder)
# library(funtimes) #notament test pour savoir si tendance dans TS
library(gamlss) #ouf distributions
# library(skellam) #skellam distribution
# library(stabledist) #distribution stable
# library(GAS)
library(copula)
library(VineCopula) #select bivariate copulas
# library(Kendall)    #tau de kendall / man-kendall test
# library(energy)     #distance correlation
library(MTS) 
library(foreach)
library(doParallel)
library(purrr) #map(list,1) => tous les premiers passagers de chaque wagon
library(fBasics)
library(magrittr)
# library(GA) # Genetic algorithm
library(lmtest)
library(fExtremes)  # Hill estimator
library(cowplot) # plot_grid list ggplot
# library(rgl) # surface 3d
library(tidyverse)
# library(extRemes)
# library(ggpval) #https://rdrr.io/cran/ggpval/man/add_pval.html
library(combinat) # permet d'obtenir toutes les permutations possibles parmis deux vecteurs (permn)
library(fGarch)
# library(caTools) # Fonction trapz pour le calcul d'int?gral 
library(naniar) # replace value with NA
library(factoextra)
library(FactoMineR)

# #Données écriture-----------------------------------------
# path_ECB = "C:/Users/vdec/Desktop/M?moire/RStudio/Donn?es/ECB Yield"
# data_ECB = read.csv(paste(path_ECB,"data.csv",sep = "/"))
# nom_Toutes_Colonnes = colnames(data_ECB)
# data_ECB = data_ECB[,!unlist(lapply(1:ncol(data_ECB), function(x) length(unique(data_ECB[,x]))==1))]
# data_ECB = setDT(data_ECB)
# data_ECB_IF = data_ECB[1:1311712,]
# data_ECB = data_ECB[2647285:3958996,]
# data_ECB$TIME_PERIOD = as.Date(data_ECB$TIME_PERIOD, format = "%Y-%m-%d")
# data_ECB = data_ECB[,c(2,3,4)]
# data_ECB$"Maturite_an" = as.numeric(tstrsplit(
#   tstrsplit(data_ECB$DATA_TYPE_FM, split = c("_"))[[2]],split = "Y")[[1]])
# data_ECB[is.na(Maturite_an)]$Maturite_an=0
# data_ECB$'Maturite_mois' = as.numeric(tstrsplit(
#   tstrsplit(data_ECB$DATA_TYPE_FM, split = "Y")[[2]], split = "M")[[1]]
# )
# data_ECB[is.na(Maturite_mois)]$Maturite_mois=0
# val_intermediaire = as.numeric(tstrsplit(
#   tstrsplit(data_ECB$DATA_TYPE_FM, split = "M")[[1]], split = "_")[[2]])
# val_intermediaire[is.na(val_intermediaire)]=0
# data_ECB$Maturite_mois = data_ECB$Maturite_mois + val_intermediaire
# data_ECB = data_ECB[,-1]
# 
# data_ECB_annuel = data_ECB[data_ECB$Maturite_mois==0]
# data_Stats = dcast(data_ECB_annuel,  Maturite_an~TIME_PERIOD, value.var = "OBS_VALUE")[,-1]
# yield = data.table(Date = unique(data_ECB$TIME_PERIOD)[2:3664],yield = as.numeric(data_Stats[10,2:3664]),
#                    increment = as.numeric(data_Stats[10,2:3664])-as.numeric(data_Stats[10,1:3663]))
# yield_all_years = t(data_Stats)
# date_ECB = unique(data_ECB$TIME_PERIOD)
# 
# # \\SRV002-ms4\vdec$\Utilisateur\R\win-library\3.5
# 
# # CSV Yield all years ---------------------------------------------------------
# yield_all_years = cbind(yield_all_years, date_ECB)
# 
# write.csv(yield_all_years, file = "yield_all_years.csv")
# 
# # CSV Yield maturit? < 1 an ---------------------------------------------------
# yield_all_months = data_ECB[data_ECB$Maturite_an==0]
# yield_all_months = dcast(yield_all_months, Maturite_mois~TIME_PERIOD, value.var = "OBS_VALUE")[,-1] %>%
#   t %>% cbind(date_ECB) %>% as.data.table
# 
# write.csv(yield_all_months, file = "yield_all_months.csv")
# 
# # CSV All yield -----------------------------------------------------------------------
# all_yield = data_ECB %>% cbind(.$Maturite_an + .$Maturite_mois/12) %>% .[,-c("Maturite_an","Maturite_mois")] %>%
#   setnames(c("TIME_PERIOD","OBS_VALUE","Maturite"))
# write.csv(all_yield, file = "all_yield.csv")
# 
# # all_yield %>% ggplot(aes(TIME_PERIOD,OBS_VALUE)) + geom_line(aes(col = Maturite)) + 
# #   theme(legend.position = "none")


# Chargement data -------------------------------------------------------------

# Donn?es ECB
path_data = "C:/Users/vdec/Desktop/Memoire/02 - Redaction/09 - R memoire/Donnees/"

yield_all_years = read.csv(paste0(path_data,"yield_all_years.csv"))

yield_all_years = yield_all_years[,-32] ; yield_all_years[,1] = as.Date(yield_all_years[,1],format = "%Y-%m-%d")
yield_all_years = cbind.data.frame(yield_all_years[,2:31],yield_all_years[,1])
colnames(yield_all_years) = c("1_an",paste0(2:30,"_ans"),"date")
# yield_all_years %>% melt(id.vars = "date") %>% ggplot(aes(date,value)) + geom_line(aes(col = variable)) +
#   theme(legend.position = 'none')

yield_all_months = read.csv(paste0(path_data,"yield_all_months.csv"))

yield_all_months = yield_all_months[,-c(1,11)] 
yield_all_months = cbind.data.frame(yield_all_months, yield_all_years$date)
colnames(yield_all_months) = c(paste0(3:11,"_mois"),"date")
# yield_all_months %>% melt(id.vars = "date") %>% ggplot(aes(date,value)) + geom_line(aes(col = variable)) +
#   theme(legend.position = "none")

# Donn?es Bank of England
england_data = read.csv(paste0(path_data,"conc_data.csv"), sep = ";")
colnames(england_data) = c("date",seq(0.5,25,by = .5))

england_EIOPA = england_data[!is.na(england_data$`25`) & !is.na(england_data$`0.5`),] ; 
england_EIOPA$date %<>% as.character %>% as.Date("%d/%m/%Y")
england_EIOPA = england_EIOPA[england_EIOPA$date >= as.Date("29/01/1998","%d/%m/%Y"),]

# Données MSCI 
MSCI = read.csv(paste0(path_data,"data_MSCI.csv"), sep = ";")

# Donnees EIOPA ----------------------------------------------------------

EIOPA = read.csv(paste0(path_data,"EIOPA.csv"), sep = ";", dec = ",")
EIOPA[,2:ncol(EIOPA)] = EIOPA[,2:ncol(EIOPA)]*100
