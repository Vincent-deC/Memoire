
# Données ----------------------------------------------------------------------
source("C:/Users/vdec/Desktop/Memoire/02 - Redaction/09 - R memoire/Donnees/chargementDonnees.R")

# Attention changement de fichier ici
source("E:/Memoire/02 - Redaction/10 - Boulot/Fonctions_modeles_test.R")

# Ecart résidus - incréments -------------


(incr-res)^2 %>% mean
mean(abs(incr))

# KS test residus MA(1) ---------------

fit_ma_moment = MA1_estimate(incr)

gamlssML(fit_ma_moment$residus, family = "NO") %>% {ks.test(fit_ma_moment$residus, "pnorm", .$mu,.$sigma)}

gamlssML(fit_ma_moment$residus, family = "JSUo") %>% {ks.test(fit_ma_moment$residus, "pJSUo", .$mu, .$sigma, .$nu, .$tau)}

# Moments résidus et théorique JSUo -----------------

paramJSUo = gamlssML(fit_ma_moment$residus, family = "JSUo")

muJSUo = paramJSUo$mu ; sigmaJSUo = paramJSUo$sigma ; nuJSUo = paramJSUo$nu ; tauJSUo = paramJSUo$tau
w = exp(1/tauJSUo^2)

meanJSUo = muJSUo - sigmaJSUo*(w^.5)*sinh(nuJSUo/tauJSUo) 

varJSUo = sigmaJSUo^2 * (w-1)*(w*cosh(2*nuJSUo/tauJSUo)+1)/2

skewnessJSUo = (-sigmaJSUo^3 * w^.5 * (w-1)^2 * (w*(w+2)*sinh(3*nuJSUo/tauJSUo) + 
                                                   3*sinh(nuJSUo/tauJSUo))/4)/(varJSUo^1.5)

kurtosisJSUo = (sigmaJSUo^4 * (w-1)^2 * (w^2 * (w^4 + 2*w^3 + 3*w^2 -3)*cosh(4*nuJSUo/tauJSUo) + 4*w^2 * (w+2) *
                                           cosh(2*nuJSUo/tauJSUo) + 3*(2*w+1)) / 8)/varJSUo^2 - 3 


meanJSUo ; varJSUo ; skewnessJSUo ; kurtosisJSUo
mean(fit_ma_moment$residus) ; var(fit_ma_moment$residus) ; skewness(fit_ma_moment$residus) 
kurtosis(fit_ma_moment$residus)-3

gamlssML(fit_ma_moment$residus, family = "NO")$mu


# Valeurs statistiques projection In Sample -----------

proj_IS_naive_gaussien = projection_single_modele(m = 10^4, modele = "Naive Gaussien")
proj_IS_naive_JSUo = projection_single_modele(m = 10^4, modele = "Naive JSUo")

proj_IS_naive_gaussien$incr %>% mean
proj_IS_naive_gaussien$incr %>% var
proj_IS_naive_gaussien$incr %>% skewness
proj_IS_naive_gaussien$incr %>% kurtosis

proj_IS_naive_JSUo$incr %>% mean
proj_IS_naive_JSUo$incr %>% var
proj_IS_naive_JSUo$incr %>% skewness
proj_IS_naive_JSUo$incr %>% kurtosis

mean(incr) ; var(incr) ; skewness(incr) ; kurtosis(incr)

# Quantile historique projection In Sample -----------

quantile(incr,probs = c(0.025,0.0025,0.975,0.9975))



# 4 - Approche générale de la structure par terme : modèle de Diebold-Li et extensions ---------------------------------------------

# 4.1 - De Nelson-Siegel à Diebold-Li : contexte -----------------------------------------------------------------------------------

# Première estimation des betas ---------------------------------------------------------------------------------------------------

load_2_5_10 = loadings_DNS(t = c(2,5,10)) %>% t 
load_2_5_10 %<>% inv

beta = load_2_5_10 %*% (cbind(yield_all_years$`2_ans`,yield_all_years$`5_ans`,yield_all_years$`10_ans`) %>% t)

# 4.2 - Altération du modèle de Diebold-Li suite au retour d’expérience des précédentes modélisations ------------------------------

# Tests de stationnarités sur les beta et résidus ----------------------------------------------------------------------------------

beta_optim %>% {lapply(1:3,function(i) adf.test(.[,i])$p.value)} %>% unlist
beta_optim %>% {lapply(1:3,function(i) pp.test(.[,i])$p.value)} %>% unlist
beta_optim %>% {lapply(1:3,function(i) kpss.test(.[,i])$p.value)} %>% unlist

residus_NS %>% {adf.test(.)$p.value} ; residus_NS %>% {pp.test(.)$p.value} ; residus_NS %>% {kpss.test(.)$p.value} 

incr_beta = (beta_optim[2:3664,] - beta_optim[1:3663,])
incr_beta %>% {lapply(1:3,function(i) adf.test(.[,i])$p.value)} %>% unlist
incr_beta %>% {lapply(1:3,function(i) pp.test(.[,i])$p.value)} %>% unlist
incr_beta %>% {lapply(1:3,function(i) kpss.test(.[,i])$p.value)} %>% unlist

incr_residus = residus_NS[2:3664] - residus_NS[1:3663]
incr_residus %>% {adf.test(.)$p.value} ; incr_residus %>% {pp.test(.)$p.value} ; incr_residus %>% {kpss.test(.)$p.value} 


# Ajustement ARMA sur les incréments beta
arma1 = auto.arima(incr_beta[,1])
arma2 = auto.arima(incr_beta[,2])
arma3 = auto.arima(incr_beta[,3])
armaRes = auto.arima(incr_residus)

arima(incr_beta[,1], order = c(0,0,0))$aic ; arma1$aic
arima(incr_beta[,2], order = c(0,0,0))$aic ; arma2$aic
arima(incr_beta[,3], order = c(0,0,0))$aic ; arma3$aic
arima(incr_residus, order = c(0,0,0))$aic ; armaRes$aic

# AIC JSUo
incr_beta[,1] %>% gamlssML(family = "JSUo") %>% {dJSUo(incr_beta[,1], mu = .$mu, sigma = .$sigma, nu = .$nu, tau = .$tau)} %>% log %>% 
{2*4-2*sum(.)} 


incr_beta[,2] %>% gamlssML(family = "JSUo") %>% {dJSUo(incr_beta[,2], mu = .$mu, sigma = .$sigma, nu = .$nu, tau = .$tau)} %>% log %>% 
{2*4-2*sum(.)} 


incr_beta[,3] %>% gamlssML(family = "JSUo") %>% {dJSUo(incr_beta[,3], mu = .$mu, sigma = .$sigma, nu = .$nu, tau = .$tau)} %>% log %>% 
{2*4-2*sum(.)} 


incr_residus %>% gamlssML(family = "JSUo") %>% {dJSUo(incr_residus, mu = .$mu, sigma = .$sigma, nu = .$nu, tau = .$tau)} %>% log %>% 
{2*4-2*sum(.)} 

# Test ajustement ARMA fenêtre roulante ----------------------------------------------------------------------------------------------
pas = 1000
iteration = seq(1,(3663-pas),by=1)

ma = paste0("ma",1:10)
ar = paste0("ar",1:10)

nbre_cores = detectCores(logical = T) - 1
registerDoParallel(cores = nbre_cores)

modeles = foreach(i = 1:(3663-pas), .packages = "forecast") %dopar% {
  j = which(iteration==i)
  incr = incr_beta[i:(pas+i),] ; incr_res = incr_residus[i:(pas+i)]
  armaPas1 = auto.arima(incr[,1])
  armaPas2 = auto.arima(incr[,2])
  armaPas3 = auto.arima(incr[,3])
  armaPasRes = auto.arima(incr_res)
  
  list(c(armaPas1$coef[!is.na(armaPas1$coef["intercept"])] , armaPas1$coef[!is.na(armaPas1$coef[ma])], 
         armaPas1$coef[!is.na(armaPas1$coef[ar])]), 
    c(armaPas2$coef[!is.na(armaPas2$coef["intercept"])] , armaPas2$coef[!is.na(armaPas2$coef[ma])], 
      armaPas2$coef[!is.na(armaPas2$coef[ar])]),
    c(armaPas3$coef[!is.na(armaPas3$coef["intercept"])] , armaPas3$coef[!is.na(armaPas3$coef[ma])], 
      armaPas3$coef[!is.na(armaPas3$coef[ar])]),
    c(armaPasRes$coef[!is.na(armaPasRes$coef["intercept"])] , armaPasRes$coef[!is.na(armaPasRes$coef[ma])], 
      armaPasRes$coef[!is.na(armaPasRes$coef[ar])]))
}

modelAR1 = NULL ; modelMA1 = NULL 
modelAR2 = NULL ; modelMA2 = NULL 
modelAR3 = NULL ; modelMA3 = NULL 
modelARres = NULL ; modelMAres = NULL 
for(i in 1:2663){
  print(i)
  modelAR1 = c(modelAR1,sum(!is.na(modeles[[i]][[1]][ar])))
  modelMA1 = c(modelMA1,sum(!is.na(modeles[[i]][[1]][ma])))
  
  modelAR2 = c(modelAR2,sum(!is.na(modeles[[i]][[2]][ar])))
  modelMA2 = c(modelMA2,sum(!is.na(modeles[[i]][[2]][ma])))
  
  modelAR3 = c(modelAR3,sum(!is.na(modeles[[i]][[3]][ar])))
  modelMA3 = c(modelMA3,sum(!is.na(modeles[[i]][[3]][ma])))
  
  modelARres = c(modelARres,sum(!is.na(modeles[[i]][[4]][ar])))
  modelMAres = c(modelMAres,sum(!is.na(modeles[[i]][[4]][ma])))
  
}

modelRep1 = data.table(AR = modelAR1, MA = modelMA1, Date = unique(yield_all_years$date)[1001:3663])
modelRep2 = data.table(AR = modelAR2, MA = modelMA2, Date = unique(yield_all_years$date)[1001:3663])
modelRep3 = data.table(AR = modelAR3, MA = modelMA3, Date = unique(yield_all_years$date)[1001:3663])
modelRepRes = data.table(AR = modelARres, MA = modelMAres, Date = unique(yield_all_years$date)[1001:3663])


figure2.2.2 = ggplot(data = modelRepRes, aes(x = AR, y = MA)) +
  geom_count(alpha = .4, show.legend = F, color = "darkblue") 
figure2.2.2 = figure2.2.2 + geom_text(data = ggplot_build(figure2.2.2)$data[[1]], 
                                      aes(x, y, label = round(prop,2)), color = "darkred",size = 15, hjust=0.5, vjust=-1) + 
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30))
figure2.2.2 + ylim(c(0,6)) + scale_size_continuous(range = c(5, 20))


# Valeurs d'AIC dans les différentes configurations ----------------------------------------------------

# Beta1 : MA(1) et AR(2)
# Beta2 : MA(2) MA(1) AR(1) AR(2)
# Beta3 : MA(1) AR(1) AR(2)
# Residus : AR(1) AR(2) AR(3)


pas = 1000
iteration = seq(1,(3663-pas),by=1)

AIC_1 = foreach(i = 1:(3663-pas), .packages = "forecast", .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_beta[i:(pas+i),1]
  
  ma1 = arima(incr, order = c(0,0,1))
  ar2 = arima(incr, order = c(2,0,0))
  normal = arima(incr, order = c(0,0,0))
  
  c(ma1$aic, ar2$aic,normal$aic)
}

AIC_1 = as.data.table(AIC_1) ; colnames(AIC_1) = c("MA1","AR2","Normal") ; 
AIC_1$Date = unique(yield_all_years$date)[1001:3663]


AIC_2 = foreach(i = 1:(3663-pas), .packages = "forecast", .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_beta[i:(pas+i),2]
  
  ma1 = arima(incr, order = c(0,0,1))
  ma2 = arima(incr, order = c(0,0,2))
  ar1 = arima(incr, order = c(1,0,0))
  ar2 = arima(incr, order = c(2,0,0))
  normal = arima(incr, order = c(0,0,0))
  
  c(ma1$aic,ma2$aic, ar1$aic, ar2$aic,normal$aic)
}

AIC_2 = as.data.table(AIC_2) ; colnames(AIC_2) = c("MA1","MA2","AR1","AR2","Normal") ; 
AIC_2$Date = unique(yield_all_years$date)[1001:3663]

AIC_3 = foreach(i = 1:(3663-pas), .packages = "forecast", .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_beta[i:(pas+i),3]
  
  ma1 = arima(incr, order = c(0,0,1))
  ar1 = arima(incr, order = c(1,0,0))
  ar2 = arima(incr, order = c(2,0,0))
  normal = arima(incr, order = c(0,0,0))
  
  c(ma1$aic, ar1$aic, ar2$aic,normal$aic)
}

AIC_3 = as.data.table(AIC_3) ; colnames(AIC_3) = c("MA1","AR1","AR2","Normal") ; 
AIC_3$Date = unique(yield_all_years$date)[1001:3663]


AIC_res = foreach(i = 1:(3663-pas), .packages = "forecast", .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_residus[i:(pas+i)]
  
  ar1 = arima(incr, order = c(1,0,0))
  ar2 = arima(incr, order = c(2,0,0))
  ar3 = arima(incr, order = c(3,0,0))
  normal = arima(incr, order = c(0,0,0))
  
  c(ar1$aic, ar2$aic, ar3$aic,normal$aic)
}

AIC_res = as.data.table(AIC_res) ; colnames(AIC_res) = c("AR1","AR2","AR3","Normal") ; 
AIC_res$Date = unique(yield_all_years$date)[1001:3663]


AIC_1 %>% melt(id.vars = "Date") %>% ggplot(aes(Date,value)) + geom_line(aes( col = variable))
AIC_2 %>% melt(id.vars = "Date") %>% ggplot(aes(Date,value)) + geom_line(aes( col = variable))
AIC_3 %>% melt(id.vars = "Date") %>% ggplot(aes(Date,value)) + geom_line(aes( col = variable))
AIC_res %>% melt(id.vars = "Date") %>% ggplot(aes(Date,value)) + geom_line(aes( col = variable))

# Independance des beta -------------------------------------------------------------------------------------------

independance = foreach(i = 1:(3663-pas), .packages = "forecast", .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_beta[i:(pas+i),] ; incr_res = incr_residus[i:(pas+i)]
  
  ind1 = Box.test(incr[,1], lag = 1)$p.value
  ind2 = Box.test(incr[,2], lag = 1)$p.value
  ind3 = Box.test(incr[,3], lag = 1)$p.value
  indRes = Box.test(incr_res, lag = 1)$p.value
  
  c(ind1, ind2, ind3, indRes)
}

independance %>% cbind(1:nrow(.)) %>% as.data.table %>% setnames(c("B1","B2","B3","Res","Index")) %>% melt(id.vars = "Index") %>%
  ggplot(aes(Index,value)) + geom_line(aes(col = variable))









# 4.2 - Dépendance DNS --------------------------------------------------------------------------------------------------------------

# Exemple IS -------------------------------------------------

aleas = foreach(i=1:3, .combine = cbind) %do% {
  incr_beta[,i] %>% gamlssML(family = "JSUo") %>% list("mu" = .$mu, "sigma" = .$sigma, "nu" = .$nu, "tau" = .$tau) %>%
    {alea_JSUo(par_JSUo = ., indice = incr_beta[,i])}
}

aleas = cbind(aleas, (incr_residus %>% gamlssML(family = "JSUo") %>% list("mu" = .$mu, "sigma" = .$sigma, "nu" = .$nu, "tau" = .$tau) %>%
{alea_JSUo(par_JSUo = ., indice = incr_residus)}))

aleas %>% cor %>% chol













# 6 - Etude et réplication des chocs de taux Solvabilité 2 -------------------------------------------------
# 6.1 - Méthode retenue par l'EIOPA -------------------------------------------------------------------------


data_england_all = england_EIOPA[,-1]

acp_england = prcomp(data_england_all, scale = T, center = T)


comp = acp_england$x[,1:4] 
changes_england = (data_england_all[261:(nrow(data_england_all)),] - data_england_all[1:(nrow(data_england_all)-260),])/
  data_england_all[1:(nrow(data_england_all)-260),]

comp_stand = apply(comp,2, function(x) (x-mean(x))/sd(x))
# comp_stand = comp_stand[1:(nrow(data_england_all)-260),]
comp_stand = comp_stand[261:(nrow(data_england_all)),]

index = seq(2,50,by = 2)
choc_each_england = foreach(i=1:25, .combine = rbind) %do% {
  
  fit_10 = lm(changes_england[,index[i]] ~ 0 + comp_stand[,1] + comp_stand[,2] + comp_stand[,3] + comp_stand[,4])
  
  f1 = (fit_10$coefficients[1]*(comp_stand[,1])) %>% quantile(probs = c(0.0025,0.9975))
  f2 = (fit_10$coefficients[2]*(comp_stand[,2])) %>% quantile(probs = c(0.0025,0.9975))
  f3 = (fit_10$coefficients[3]*(comp_stand[,3])) %>% quantile(probs = c(0.0025,0.9975))
  f4 = (fit_10$coefficients[4]*(comp_stand[,4])) %>% quantile(probs = c(0.0025,0.9975))
  
  c(round(f1[1] + f2[1] + f3[1] + f4[1],2)*100,round(f1[2] + f2[2] + f3[2] + f4[2],2)*100)
}

# 6.2 - Réplication chocs EIOPA donnes BCE et modeles --------------------------------------------------------

data_BCE = yield_all_years[1:2200,-31]

acp_bce = prcomp(data_BCE, scale = T, center = T)

(acp_bce$sdev^2/sum(acp_bce$sdev^2)) %>% .[1:4] %>% sum
  
comp_bce = acp_bce$x[,1:4]


changes_bce = (data_BCE[261:(nrow(data_BCE)),] - data_BCE[1:(nrow(data_BCE)-260),])/
  data_BCE[1:(nrow(data_BCE)-260),]

comp_stand_bce = apply(comp_bce,2, function(x) (x-mean(x))/sd(x))
comp_stand_bce = comp_stand_bce[1:(nrow(data_BCE)-260),]
# comp_stand_bce = comp_stand_bce[261:(nrow(data_BCE)),]


choc_each_bce = foreach(i=1:30, .combine = rbind) %do% {
  
  fit_10 = lm(changes_bce[,i] ~ 0 + comp_stand_bce[,1] + comp_stand_bce[,2] + comp_stand_bce[,3] + comp_stand_bce[,4])
  
  f1 = (fit_10$coefficients[1]*(comp_stand_bce[,1])) %>% quantile(probs = c(0.0025,0.9975))
  f2 = (fit_10$coefficients[2]*(comp_stand_bce[,2])) %>% quantile(probs = c(0.0025,0.9975))
  f3 = (fit_10$coefficients[3]*(comp_stand_bce[,3])) %>% quantile(probs = c(0.0025,0.9975))
  f4 = (fit_10$coefficients[4]*(comp_stand_bce[,4])) %>% quantile(probs = c(0.0025,0.9975))
  
  c(round(f1[1] + f2[1] + f3[1] + f4[1],2)*100,round(f1[2] + f2[2] + f3[2] + f4[2],2)*100)
}



















# 6.2 - Chocs modeles -------------------------------------------

# Fin 2015 : 30/12/2015
# Fin 2016 : 30/12/2016
# Fin 2017 : 29/12/2017
# Fin 2018 : 28/12/2018

yc_BCE_31_12_2015 = yield_all_years[yield_all_years$date == as.Date("30/12/2015",format = "%d/%m/%Y"),-31] %>% as.numeric
yc_BCE_31_12_2016 = yield_all_years[yield_all_years$date == as.Date("30/12/2016",format = "%d/%m/%Y"),-31] %>% as.numeric
yc_BCE_31_12_2017 = yield_all_years[yield_all_years$date == as.Date("29/12/2017",format = "%d/%m/%Y"),-31] %>% as.numeric
yc_BCE_31_12_2018 = yield_all_years[yield_all_years$date == as.Date("28/12/2018",format = "%d/%m/%Y"),-31] %>% as.numeric


y = EIOPA$spot_2015[1:20]
par_eiopa_2015 = estim_par_EIOPA(yc = EIOPA$spot_2015[1:20])

y = EIOPA$spot_2016[1:20]
par_eiopa_2016 = estim_par_EIOPA(yc = EIOPA$spot_2016[1:20])

y = EIOPA$spot_2017[1:20]
par_eiopa_2017 = estim_par_EIOPA(yc = EIOPA$spot_2017[1:20])

y = EIOPA$spot_2018[1:20]
par_eiopa_2018 = estim_par_EIOPA(yc = EIOPA$spot_2018[1:20])

choc_2015_original = foreach(i=1:20, .combine = rbind) %do% {
  param_DNS = calibration_DNS_original(x = yield_all_years[yield_all_years$date<=as.Date("30/12/2015",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_original(param_EIOPA = par_eiopa_2015, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}


choc_2015_copule = foreach(i=1:20, .combine = rbind) %do% {
  param_DNS = calibration_DNS_cor(x = yield_all_years[yield_all_years$date<=as.Date("30/12/2015",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_cor(param_EIOPA = par_eiopa_2015, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}

cbind(1:20, EIOPA$spot_2015[1:20], EIOPA$up_2015[1:20], EIOPA$down_2015[1:20], choc_2015_original[,1], choc_2015_original[,2], 
      EIOPA$spot_2016[1:20],choc_2015_copule[,1],choc_2015_copule[,2]) %>%
  as.data.table %>% setnames(c("Maturite","EIOPA_2015","EIOPA_UP","EIOPA_DOWN","DNS_orig_DOWN","DNS_orig_UP","EIOPA_2016",
                               "DNS_cor_DOWN","DNS_cor_UP")) %>%
  ggplot(aes(Maturite,EIOPA_2015)) + geom_line() + geom_line(aes(y=EIOPA_UP, col = "Choc EIOPA")) + 
  geom_line(aes(y = EIOPA_DOWN, col = "Choc EIOPA")) + geom_line(aes(y = DNS_orig_DOWN, col = "Quantile DNS Original")) + 
  geom_line(aes(y = DNS_orig_UP, col = "Quantile DNS Original")) + geom_line(aes(y = EIOPA_2016, col = "EIOPA 2016")) +
  geom_line(aes(y = DNS_cor_DOWN, col = "Quantile DNS Corrélation")) + 
  geom_line(aes(y = DNS_cor_UP, col = "Quantile DNS Corrélation"))


s_up = c(0.61,0.53,0.49,0.46,0.45,0.41,0.37,0.34,0.32,0.30,0.3,0.3,0.3,0.29,0.28,0.28,0.27,0.26,0.26,0.25)
s_down = c(0.58,0.51,0.44,0.4,0.4,0.38,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.47,0.48,0.49,0.49,0.5)

b_up = c(2.14,1.86,1.72,1.61,1.58,1.44,1.3,1.19,1.12,1.05,1.05,1.05,1.05,1.02,0.98,0.98,0.95,0.91,0.91,0.88)
b_down = c(1.16,0.99,0.83,0.74,0.71,0.67,0.63,0.62,0.61,0.61,0.6,0.6,0.59,0.58,0.57,0.56,0.55,0.54,0.52,0.5)


EIOPA_new_stress = matrix(rep(0,20*8),ncol = 8) %>% as.data.table %>% 
  setnames(c("down_2015","up_2015","down_2016","up_2016","down_2017","up_2017","down_2018","up_2018"))

index = rep(seq(2,11,by =3),each = 2)
for(i in seq(1,8,by = 2)){
  EIOPA_new_stress[,i] = (lapply(1:20, function(j) EIOPA[j,index[i]]*(1-s_down[j])-b_down[j]) %>% unlist) 
}

for(i in seq(2,8,by = 2)){
  EIOPA_new_stress[,i] = (lapply(1:20, function(j) EIOPA[j,index[i]]*(1+s_up[j])+b_up[j]) %>% unlist) 
}





















# Annexes -----------------------------------------------------------------------------------

# Annexe B OoS ------------------------------------------------------------------------------

# 1 à 10 ans ---------------------------------------------------------------------------------

DNS_copule = matrix(rep(0,7*30), ncol = 7)
annee = 0
for(i in (annee+1):(annee+30)){
  DNS_copule[i-annee,] = c(i,(projection_modeles_stats_OoS(x = yield_all_years[,i], m = 10^3, modeles = "DNS copule", taux = i, pas = 10) %>% unlist)) 
  print(i)
}
DNS_copule %>% cbind.data.frame(rep("DNS copule",nrow(DNS_copule)))


# 11 à 20 ans -------------------------------------------------------------------------------- 

DNS_copule = matrix(rep(0,7*10), ncol = 7)
for(i in 11:20){
  DNS_copule[i-10,] = c(i,(projection_modeles_stats_OoS(x = yield_all_years[,i], m = 10^3, modeles = "DNS copule", taux = i, pas = 10) %>% unlist)) 
  print(i)
}
DNS_copule %>% cbind.data.frame(rep("DNS copule",nrow(DNS_copule)))




# 21 à 30 ans -------------------------------------------------------------------------------- 


hyb = matrix(rep(0,7*10), ncol = 7)
annee = 24
for(i in (annee+1):(annee+10)){
  hyb[i-annee,] = c(i,(projection_modeles_stats_OoS(x = yield_all_years[,i], m = 10^3, modeles = "DNS copule", taux = i, pas = 10) %>% unlist)) 
  print(i)
}
hyb %>% cbind.data.frame(rep("DNS copule",nrow(hyb)))













