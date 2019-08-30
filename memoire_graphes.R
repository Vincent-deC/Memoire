
# Donnees ---------------------------------------------------------------------------------------------
source("C:/Users/vdec/Desktop/Memoire/02 - Redaction/09 - R memoire/Donnees/chargementDonnees.R")
source("C:/Users/vdec/Desktop/Memoire/02 - Redaction/09 - R memoire/Donnees/Fonctions_modeles_graphes.R")
path_save = "C:/Users/vdec/Desktop/Memoire/02 - Redaction/09 - R memoire/FigureFinal/"

yield = cbind.data.frame("date" = yield_all_years$date[2:3664], "yield" = yield_all_years$`10_ans`[2:3664],
                         increment = yield_all_years$`10_ans`[2:3664]-yield_all_years$`10_ans`[1:3663]) %>% as.data.table


# -----------------------------------------------------------------------------------------------------
# 1 - Préambule ---------------------------------------------------------------------------------------

# 1.4 - Courbes de taux et pr?sentation des données utilisées -----------------------------------------

# Figure 1.4.1 - Forme accessible ? la forme paramétrique NS ----------
y = function(){
  a = seq(-6,12,length.out = 7) ; m = seq(0,10,length.out = 10^4)
  ts = matrix(rep(0,length(a)*(length(m)-1)), ncol = length(a))
  ts = rbind(rep(0,7),ts)
  for(i in 1:length(a)){
    ts[2:(length(m)),i] = 1- (1-a[i])*( (1-exp(-m[2:10^4])))/m[2:10^4] - a[i]*exp(-m[2:10^4])
  }
  return(ts)
}

data = as.data.table(y()) ; data$m =seq(0,10,length.out = 10^4)

figure1.4.1 = ggplot(data, aes(x=m)) + geom_line(aes(y=V1, col = "s")) + geom_line(aes(y=V2, col = "d")) +
  geom_line(aes(y=V3, col = "f")) + geom_line(aes(y=V4, col = "g")) + 
  geom_line(aes(y=V5, col = "q")) + geom_line(aes(y=V6, col = "h")) +
  geom_line(aes(y=V7, col = "v")) + xlab("") + ylab("") +
  theme(axis.text = element_blank(), legend.position = "none") 

# ggsave(paste0(path_save,"Figure1_4_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure1.4.1)




# Figure 1.4.3 - Plot des taux sur les maturités ------------------------------------------------------

figure1.4.3 = yield_all_years %>% melt(id.vars = "date") %>% ggplot(aes(date,value)) + geom_line(aes(col = variable)) + 
  theme_cowplot() +
  theme(legend.position = "none",axis.title = element_text(size = 30),
        axis.text = element_text(size = 40)) + ylab("Taux") + xlab("Date")

# ggsave(paste0(path_save,"Figure1_4_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure1.4.3)


# 1.5 - Critères de validation et de séléction de modèle ----------------------------------------------

# Figures 1.5.1&2&3 : Illustration critère IS --------------------------------------

proj_IS = projection_modeles(m=10^3, modeles = c("JSUo","naive gaussien", "DNS gaussien"))
gg_In = gg_IS(proj_IS)

figure1.5.1 = gg_In$g_IC + theme_cowplot() + 
  theme(legend.position = "none", axis.text = element_text(size = 40),
        axis.title = element_text(size = 40))  + labs(title = "")

figure1.5.2 = gg_In$g_log_inf + labs(title = "") + ylab("Log fonction répartition") + xlab("Log valeur absolue support") + theme_cowplot() +
  theme(legend.position = "none", axis.text = element_text(size = 45),
        axis.title = element_text(size = 45))

figure1.5.3 = gg_In$g_log_sup + labs(title = "") + ylab("Log fonction de survie") + xlab("Log support") + theme_cowplot() +
  theme(legend.position = "none", axis.text = element_text(size = 45),
        axis.title = element_text(size = 45))

# ggsave(paste0(path_save,"Figure1_5_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure1.5.1)
# ggsave(paste0(path_save,"Figure1_5_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure1.5.2)
# ggsave(paste0(path_save,"Figure1_5_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure1.5.3)


# Figures 1.5.4&5&6&7 : Illustration critère OoS --------------------------------------------

proj_OoS = proj_OoS_par(m=10^3,pas = 10, modeles = c("Gaussien"))

figure1.5.4 = (proj_OoS$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")

figure1.5.5 = (proj_OoS$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")


log_log_Gaussien = gg_all_log_log_plot(pas = 10, modele = "Gaussien", m = 10^3)

figure1.5.6 = log_log_Gaussien$g_inf + xlim(c(-6,0)) + ylim(c(-6,0)) + xlab("") + ylab("") + labs(title = "") + 
  theme(axis.text = element_text(size = 40))
figure1.5.7 = log_log_Gaussien$g_sup + xlim(c(-6,0)) + ylim(c(-6,0)) + xlab("") + ylab("") + labs(title = "") + 
  theme(axis.text = element_text(size = 40))

# ggsave(paste0(path_save,"Figure1_5_4.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure1.5.4)
# ggsave(paste0(path_save,"Figure1_5_5.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure1.5.5)
# ggsave(paste0(path_save,"Figure1_5_6.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure1.5.6)
# ggsave(paste0(path_save,"Figure1_5_7.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure1.5.7)





# -----------------------------------------------------------------------------------------------------
# 2 - APPROCHE FONDEE SUR L'ETUDE EMPIRIQUE DES VARIATIONS JOURNALIERES DU TAUx 10 ANS ----------------

# 2.1 - Caractéristiques statistiques des variations et étude de la série temporelle ------------------

#Figure 2.1.1 - Historique taux 10 ans  -------------------------------------

gIncr = ggplot(yield, aes(date,increment)) + geom_line(color = "darkred")  +
  xlab("") + ylab("Incréments") + theme(axis.text = element_text(size = 45), plot.margin = unit(c(0,2,0,0),"mm"),
                                        axis.text.x = element_blank(), axis.title = element_text(size = 0))


gYield = ggplot(yield, aes(date,yield)) + geom_line(color = "darkblue") +
  xlab("") + ylab("Taux") + theme(axis.text = element_text(size = 45), plot.margin = unit(c(3,2,5,8),"mm"),
                                  axis.title = element_text(size = 0 ))

figure2.1.1 = arrangeGrob(gIncr,gYield, ncol = 1)

# ggsave(paste0(path_save,"Figure2_1_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.1.1)




# Figure 2.1.2 - Autocorrélogramme ----------------------------------

figure2.1.2 = autoplot(acf(yield$increment, type  ="correlation", trace = F)) + 
  theme(axis.title = element_text(size = 45), axis.text = element_text(size = 45)) + ylab("Corrélation")

# ggsave(paste0(path_save,"Figure2_1_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.1.2)



  

# 2.2 - Modélisation de la série temporelle par un processus ARMA ---------------------------

#Figure 2.2.1 - Increments et residus ------------------------------

residus = arima(yield$increment, order = c(0,0,1))$residuals

figure2.2.1 = ggplot(yield, aes(date,increment)) + geom_line(color = "darkblue", alpha = .5) + 
  geom_line(aes(y = residus, col = "darkred"), alpha = .5) +
  xlab("") + ylab("") + theme(axis.text = element_text(size = 40), plot.margin = unit(c(-3.5,2,5,5.2),"mm"),
                                  axis.title = element_text(size = 0 ), legend.position = "none")

# ggsave(paste0(path_save,"Figure2_2_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.2.1)



#Figure 2.2.2 - Proportion de modèles ARMA fittés -------------------

pas = 1000
iteration = seq(1,(3663-pas),by=1)

ma = paste0("ma",1:10)
ar = paste0("ar",1:10)

nbre_cores = detectCores(logical = T) - 1
registerDoParallel(cores = nbre_cores)

model = foreach(i = 1:(3663-pas), .packages = "forecast") %dopar% {
    j = which(iteration==i)
    incr = yield$increment[i:(pas+i)]
    armaPas = auto.arima(incr)
    
    c(armaPas$coef[!is.na(armaPas$coef["intercept"])] , armaPas$coef[!is.na(armaPas$coef[ma])], 
                       armaPas$coef[!is.na(armaPas$coef[ar])])
}

modelAR = NULL ; modelMA = NULL 
for(i in 1:2663){
  print(i)
  modelAR = c(modelAR,sum(!is.na(model[[i]][ar])))
  modelMA = c(modelMA,sum(!is.na(model[[i]][ma])))
  
}

modelRep = data.table(AR = modelAR, MA = modelMA, Date = unique(yield$date)[1001:3663])

figure2.2.2 = ggplot(data = modelRep, aes(x = AR, y = MA)) +
  geom_count(alpha = .4, show.legend = F, color = "darkblue") 
figure2.2.2 = figure2.2.2 + geom_text(data = ggplot_build(figure2.2.2)$data[[1]], 
                  aes(x, y, label = round(prop,2)), color = "darkred",size = 15, hjust=0.5, vjust=-1) + 
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45))
figure2.2.2 = figure2.2.2 + ylim(c(0,6)) + scale_size_continuous(range = c(5, 20))

# ggsave(paste0(path_save,"Figure2_2_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.2.2)


#Figure 2.2.3 - AIC pour MA(1), AR(1) et AR(2) ------------------------------


pas = 1000
iteration = seq(1,(3663-pas),by=1)

AIC = foreach(i = 1:(3663-pas), .packages = "forecast", .combine = rbind) %dopar% {

  j = which(iteration==i)
  incr = yield$increment[i:(pas+i)]
  
  ma1 = arima(incr, order = c(0,0,1))
  ar1 = arima(incr, order = c(1,0,0))
  ar2 = arima(incr, order = c(2,0,0))
  
  c(ma1$aic, ar1$aic, ar2$aic)
}

AIC = as.data.table(AIC) ; colnames(AIC) = c("MA1","AR1","AR2") ; 
AIC$Simulation = 1:nrow(AIC)

figure2.2.3 = ggplot(AIC, aes(x=Simulation)) + geom_line(aes(y=MA1, col = "MA(1)"), size = 2) + geom_line(aes(y=AR1, col = "AR(1)"), size = 2) +
  geom_line(aes(y=AR2, col = "AR(2)"), size = 2) + ylab("AIC") + 
  scale_color_manual(name = "Modèles :", values = gg_color_hue(3)) +
  theme(legend.title = element_text(size = 45), legend.text = element_text(size = 45),
        axis.title = element_text(size = 45), axis.text = element_text(size = 45)) + guides(colour = guide_legend(override.aes = list(size=5)))

# ggsave(paste0(path_save,"Figure2_2_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.2.3)



#Calcul - Valeurs des paramètres MA(1), AR(1) et AR(2) -------------------------------------------------------


pas = 1000
iteration = seq(1,(3663-pas),by=1)

param = foreach(i = 1:(3663-pas), .packages = "forecast", .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = yield$increment[i:(pas+i)]
  
  ma1 = arima(incr, order = c(0,0,1))$coef[1]
  ar1 = arima(incr, order = c(1,0,0))$coef[1]
  ar2 = arima(incr, order = c(2,0,0))$coef[1:2]
  
  c(ma1, ar1, ar2)
}

param = as.data.table(param) ; colnames(param) = c("MA1","AR1","AR2_1", "AR2_2") ; 
param$Index = 1:length(iteration)

param$MA1 %>% mean ; param$MA1 %>% min ; param$MA1 %>% max
param$AR1 %>% mean ; param$AR1 %>% min ; param$AR1 %>% max
param$AR2_1 %>% mean ; param$AR2_1 %>% min ; param$AR2_1 %>% max
param$AR2_2 %>% mean ; param$AR2_2 %>% min ; param$AR2_2 %>% max

# figure = ggplot(param, aes(x=Index)) + geom_line(aes(y=MA1, col = "MA(1)"), size = 2) + geom_line(aes(y=AR1, col = "AR(1)"), size = 2) +
#   geom_line(aes(y=AR2_1, col = "AR(2)_1"), size = 2) + geom_line(aes(y=AR2_2, col = "AR(2)_2"), size = 2) + ylab("Valeur de paramètre") + 
#   scale_color_manual(name = "Paramètres :", values = gg_color_hue(4)) +
#   theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30),
#         axis.title = element_text(size = 30), axis.text = element_text(size = 30))

#Figure 2.2.4 - Indépendance des résidus ------------------------------

pas = 1000
iteration = seq(1,(3663-pas),by=1)

independance = foreach(i = 1:(3663-pas), .packages = "forecast", .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = yield$increment[i:(pas+i)]
  
  ma1 = arima(incr, order = c(0,0,1))$residuals
  ar1 = arima(incr, order = c(1,0,0))$residuals
  ar2 = arima(incr, order = c(2,0,0))$residuals
  
  c(Box.test(ma1, lag = 5)$p.value, Box.test(ar1, lag = 5)$p.value, Box.test(ar2, lag = 5)$p.value)
}

independance = as.data.table(independance) ; colnames(independance) = c("MA1","AR1","AR2") ; 
independance$Simulation = 1:nrow(independance)

figure2.2.4 = ggplot(independance, aes(x=Simulation)) + geom_line(aes(y=MA1, col = "MA(1)"), size = 2) + 
  geom_line(aes(y=AR1, col = "AR(1)"), size = 2) +
  geom_line(aes(y=AR2, col = "AR(2)"), size = 2) + ylab("P-value") + 
  scale_color_manual(name = "Modèles :", values = gg_color_hue(3)) +
  theme(legend.title = element_text(size = 45), legend.text = element_text(size = 45),
        axis.title = element_text(size = 45), axis.text = element_text(size = 45)) + ylim(c(0,1)) + 
  guides(colour = guide_legend(override.aes = list(size=5)))

# ggsave(paste0(path_save,"Figure2_2_4.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.2.4)

#Figure 2.2.5 - Densité KDE résidus MA(1) -------------------------
arma = auto.arima(yield$increment)

densites = data.table(Bruit = arma$residuals)

figure2.2.5 = ggplot(densites, aes(x=as.numeric(Bruit))) + geom_density(fill = "darkred", alpha = .2) + ylab("Densité") +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30),
        axis.title = element_text(size = 30)) + 
  scale_color_manual(name = "Densités : ", values = c(gg_color_hue(2)))  + xlab("")

ggsave(paste0(path_save,"Figure2_2_5.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.2.5)

#Figure 2.2.7 - Ajustement d'une loi normale et JSUo sur densité des résidus -------------------------
arma = auto.arima(yield$increment)

fitNormal = fitdistr(arma$residuals, densfun = "normal")
fitJSUo = gamlssML(arma$residuals, family = "JSUo") 

dNormal = dnorm(arma$residuals, mean = fitNormal$estimate[1], sd = fitNormal$estimate[2])
dJSUo = dJSUo(x=arma$residuals, mu = fitJSUo$mu, sigma = fitJSUo$sigma, 
              nu = fitJSUo$nu, tau = fitJSUo$tau)

densites = data.table(Bruit = arma$residuals, Normal = dNormal, JSUo = dJSUo)

figure2.2.7 = ggplot(densites, aes(x=as.numeric(Bruit))) + geom_density(fill = "darkred", alpha = .2) + 
  geom_line(aes(y=Normal, col = "Normal"), size = 2) +
  geom_line(aes(y=JSUo, col = "JSUo"), size = 2) + ylab("Densité") +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30),
        axis.title = element_text(size = 30)) + 
  scale_color_manual(name = "Densités : ", values = c(gg_color_hue(2)))  + xlab("")

ggsave(paste0(path_save,"Figure2_2_7.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.2.7)


#Figure 2.2. 8&9 - Queues de distribution de l'ajustement  ------------------

arma = auto.arima(yield$increment)

fitNormal = fitdistr(arma$residuals, densfun = "normal")
fitJSUo = gamlssML(arma$residuals, family = "JSUo") 

dNormal = dnorm(arma$residuals, mean = fitNormal$estimate[1], sd = fitNormal$estimate[2])
dJSUo = dJSUo(x=arma$residuals, mu = fitJSUo$mu, sigma = fitJSUo$sigma, 
              nu = fitJSUo$nu, tau = fitJSUo$tau)

densites = data.table(Bruit = arma$residuals, Normal = dNormal, JSUo = dJSUo)

figure2.2.8.9 = ggplot(densites, aes(x=as.numeric(Bruit))) + geom_density(fill = "darkred", alpha = .2) + 
  geom_line(aes(y=Normal, col = "Normal"), size = 1) +
  geom_line(aes(y=JSUo, col = "JSUo"), size = 1) + ylab("") + theme_cowplot() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 30)) + 
  scale_color_manual(name = "Densités : ", values = c(gg_color_hue(2)))  + xlab("") 
ggplotly(figure2.2.8.9)

# 2.3 - Enseignements généraux de cette première modélisation ----------------------------------------
#Figure 2.3. 1&2&3&4 - Projection IS approche naïve -----------------------

proj_IS_naive_gaussien = projection_single_modele(m = 10^4, modele = "Naive Gaussien")
proj_IS_naive_JSUo = projection_single_modele(m = 10^4, modele = "Naive JSUo")

figure2.3.1 = gg_faisceau(traj = proj_IS_naive_gaussien$yield, m_faisceau = 10^2, double_intervalle = T, probs = c(0.025,0.975,0.0025,0.9975)) + 
  ylim(c(-7,10)) + theme(axis.text = element_text(size = 45))

figure2.3.2 = gg_faisceau(traj = proj_IS_naive_JSUo$yield, m_faisceau = 10^2, double_intervalle = T, probs = c(0.025,0.975,0.0025,0.9975)) + 
  ylim(c(-7,10)) + theme(axis.text = element_text(size = 45))

y_list = list("MA(1) gaussien" = proj_IS_naive_gaussien$incr, "MA(1) JSUo" = proj_IS_naive_JSUo$incr)
gg_log = gg_log_multi_comp(x = incr, y_list = y_list)

log_incr = log_log_dt(incr)
x_inf_025 = (round(log_incr$queue_inf$y_inf,2) == round(log10(0.025),2)) %>% {log_incr$queue_inf$x_inf[which(.)]} %>% mean
x_inf_0025 = (round(log_incr$queue_inf$y_inf,1) == round(log10(0.0025),1)) %>% {log_incr$queue_inf$x_inf[which(.)]} %>% mean

x_sup_975 = (round(log_incr$queue_sup$y_sup,2) == round(log10(0.025),2)) %>% {log_incr$queue_sup$x_sup[which(.)]} %>% mean
x_sup_9975 = (round(log_incr$queue_sup$y_sup,1) == round(log10(0.0025),1)) %>% {log_incr$queue_sup$x_sup[which(.)]} %>% mean


figure2.3.3 = gg_log$g_inf + theme_cowplot() + theme(legend.position = "none", axis.text = element_text(size = 45)) + xlab("") + ylab("") + 
  labs(tilte = "0") + geom_segment(aes(x = min(log_incr$queue_inf$x_inf), xend = x_inf_025, y = log10(0.025), yend = log10(0.025)), 
                                   linetype = "dashed") + 
  geom_segment(aes(x =x_inf_025, xend = x_inf_025, y = -7, yend = log10(0.025)), linetype = "dashed") + 
  geom_segment(aes(x = min(log_incr$queue_inf$x_inf), xend = x_inf_0025, y = log10(0.0025), yend = log10(0.0025)), linetype = "dashed") + 
  geom_segment(aes(x =x_inf_0025, xend = x_inf_0025, y = -7, yend = log10(0.0025)), linetype = "dashed") + labs(title = "")

figure2.3.4 = gg_log$g_sup + theme_cowplot() + theme(legend.position = "none", axis.text = element_text(size = 45)) + xlab("") + ylab("") + 
  labs(tilte = "0")  + geom_segment(aes(x = min(log_incr$queue_sup$x_sup), xend = x_sup_975, y = log10(0.025), yend = log10(0.025)), 
                                    linetype = "dashed") + 
  geom_segment(aes(x =x_sup_975, xend = x_sup_975, y = -7, yend = log10(0.025)), linetype = "dashed") + 
  geom_segment(aes(x = min(log_incr$queue_sup$x_sup), xend = x_sup_9975, y = log10(0.0025), yend = log10(0.0025)), linetype = "dashed") + 
  geom_segment(aes(x =x_sup_9975, xend = x_sup_9975, y = -7, yend = log10(0.0025)), linetype = "dashed") + labs(title = "")

# ggsave(paste0(path_save,"Figure2_3_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.1)
# ggsave(paste0(path_save,"Figure2_3_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.2)
# 
# ggsave(paste0(path_save,"Figure2_3_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.3)
# ggsave(paste0(path_save,"Figure2_3_4.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.4)




#Figure 2.3. 5&6&7&8&9&10&11&12 - Résultats OoS approche naïve - Erreurs ---------------------------------------------------------

proj_OoS_naive_gaussien = proj_OoS_par(m = 10^4, pas = 10, modeles = "Naive Gaussien")
proj_OoS_naive_JSUo = proj_OoS_par(m = 10^4, pas = 10, modeles = "Naive JSUo")
proj_OoS_Gaussien = proj_OoS_par(m = 10^4, pas = 10, modeles = "Gaussien")
proj_OoS_JSUo = proj_OoS_par(m = 10^4, pas = 10, modeles = "JSUo")

# Naive gaussien
figure2.3.5 = (proj_OoS_naive_gaussien$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") +
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("") + geom_hline(yintercept = 0.05, linetype = "dashed", col = "darkred")

# Naive JSUo
figure2.3.6 = (proj_OoS_naive_JSUo$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") +
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("") + geom_hline(yintercept = 0.05, linetype = "dashed", col = "darkred")

# Gaussien
figure2.3.7 = (proj_OoS_Gaussien$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") +
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("") + geom_hline(yintercept = 0.05, linetype = "dashed", col = "darkred")

# JSUo
figure2.3.8 = (proj_OoS_JSUo$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") +
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("") + geom_hline(yintercept = 0.05, linetype = "dashed", col = "darkred")

# ggsave(paste0(path_save,"Figure2_3_5.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.5)
# ggsave(paste0(path_save,"Figure2_3_6.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.6)
# ggsave(paste0(path_save,"Figure2_3_7.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.7)
# ggsave(paste0(path_save,"Figure2_3_8.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.8)

# Naive gaussien
figure2.3.9 = (proj_OoS_naive_gaussien$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") +
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("") 

# Naive JSUo
figure2.3.10 = (proj_OoS_naive_JSUo$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") +
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("") 

# Gaussien
figure2.3.11 = (proj_OoS_Gaussien$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") +
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")

# JSUo
figure2.3.12 = (proj_OoS_JSUo$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") +
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")

# ggsave(paste0(path_save,"Figure2_3_9.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.9)
# ggsave(paste0(path_save,"Figure2_3_10.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.10)
# ggsave(paste0(path_save,"Figure2_3_11.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.11)
# ggsave(paste0(path_save,"Figure2_3_12.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.12)

proj_OoS_Bootstrap = proj_OoS_par(m = 10^4, pas = 10, modeles = "Bootstrap")

proj_OoS_naive_gaussien$Nbre_erreur_95 %>% apply(2,mean)/1000
proj_OoS_naive_JSUo$Nbre_erreur_95 %>% apply(2,mean)/1000
proj_OoS_Gaussien$Nbre_erreur_95 %>% apply(2,mean)/1000
proj_OoS_JSUo$Nbre_erreur_95 %>% apply(2,mean)/1000

proj_OoS_naive_gaussien$Nbre_erreur_99_5 %>% apply(2,mean)/1000
proj_OoS_naive_JSUo$Nbre_erreur_99_5 %>% apply(2,mean)/1000
proj_OoS_Gaussien$Nbre_erreur_99_5 %>% apply(2,mean)/1000
proj_OoS_JSUo$Nbre_erreur_99_5 %>% apply(2,mean)/1000

proj_OoS_Bootstrap$Nbre_erreur_95 %>% apply(2,mean)/1000
proj_OoS_Bootstrap$Nbre_erreur_99_5 %>% apply(2,mean)/1000


# Figure 2.3. 13&14&15&16&17&18&19&20&21&22 - Faisceaux de queues de distribution --------------------------------------------

queues_Gaussien = gg_all_log_log_plot(m = 10^3, modele = "Gaussien", pas = 10)
queues_JSUo = gg_all_log_log_plot(m = 10^3, modele = "JSUo", pas = 10)
queues_naive_Gaussien = gg_all_log_log_plot(m = 10^3, modele = "Naive Gaussien", pas = 10)
queues_naive_JSUo = gg_all_log_log_plot(m = 10^3, modele = "Naive JSUo", pas = 10)
queues_empirique = gg_all_log_log_plot(m = 10^3, modele = "Empirique", pas = 10)

figure2.3.13 = queues_Gaussien$g_inf + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
  theme(axis.text = element_text(size = 45))
# figure2.3.14 = queues_Gaussien$g_sup + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
#   theme(axis.text = element_text(size = 40))

figure2.3.15 = queues_JSUo$g_inf + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
  theme(axis.text = element_text(size = 45))
# figure2.3.16 = queues_JSUo$g_sup + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
#   theme(axis.text = element_text(size = 40))

figure2.3.17 = queues_naive_Gaussien$g_inf + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
  theme(axis.text = element_text(size = 45))
# figure2.3.18 = queues_naive_Gaussien$g_sup + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
#   theme(axis.text = element_text(size = 40))

figure2.3.19 = queues_naive_JSUo$g_inf + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
  theme(axis.text = element_text(size = 45))
# figure2.3.20 = queues_naive_JSUo$g_sup + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
#   theme(axis.text = element_text(size = 40))

figure2.3.21 = queues_empirique$g_inf + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
  theme(axis.text = element_text(size = 45))
# figure2.3.22 = queues_empirique$g_sup + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
#   theme(axis.text = element_text(size = 40))

ggsave(paste0(path_save,"Figure2_3_13.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.13)
# ggsave(paste0(path_save,"Figure2_3_14.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.14)
ggsave(paste0(path_save,"Figure2_3_15.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.15)
# ggsave(paste0(path_save,"Figure2_3_16.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.16)
ggsave(paste0(path_save,"Figure2_3_17.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.17)
# ggsave(paste0(path_save,"Figure2_3_18.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.18)
ggsave(paste0(path_save,"Figure2_3_19.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.19)
# ggsave(paste0(path_save,"Figure2_3_20.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.20)
ggsave(paste0(path_save,"Figure2_3_21.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.21)
# ggsave(paste0(path_save,"Figure2_3_22.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure2.3.22)


# ---------------------------------------------------------------------------------------------------
# 3 - Approche distributionnelle des variations du taux 10 ans --------------------------------------

# 3.1 - Ajustement d'une loi hybride sur les variations ---------------------------------------------

# Figure 3.1. 1&2 - Zones de recherche des droites de Pareto ----------------------------------------------
dt_log_incr = log_log_dt(yield$increment)

dt_log_incr$queue_sup = dt_log_incr$queue_sup[-nrow(dt_log_incr$queue_sup),]

figure3.1.1 = dt_log_incr$queue_inf %>% ggplot(aes(x_inf,y_inf)) + geom_line() + geom_point(size = 3, col =  "darkblue", alpha = .5) + 
  xlim(c(-6,0)) + theme_cowplot() + 
  theme(axis.text = element_text(size = 45)) +
  annotate("rect", xmax=-.8, xmin=-2, ymin=-Inf, ymax=Inf, alpha=0.2, fill="darkred") + xlab("") + ylab("") + labs(title = "")
figure3.1.2 = dt_log_incr$queue_sup %>% ggplot(aes(x_sup,y_sup)) + geom_line() + geom_point(size = 3, col =  "darkblue", alpha = .5) + 
  xlim(c(-6,0)) + theme_cowplot() + 
  theme(axis.text = element_text(size = 45)) +
  annotate("rect", xmax=-.8, xmin=-2, ymin=-Inf, ymax=Inf, alpha=0.2, fill="darkred") + xlab("") + ylab("") + labs(title = "")

# ggsave(paste0(path_save,"Figure3_1_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.1)
# ggsave(paste0(path_save,"Figure3_1_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.2)


# Figure 3.1. 3&4 - Droites de Pareto et points de coupure ----------------------------------------------

intervalle_recherche_inf = dt_log_incr$queue_inf$x_inf[(dt_log_incr$queue_inf$y_inf<= -.8 & 
                                                          dt_log_incr$queue_inf$y_inf >= -2)]
intervalle_recherche_sup = dt_log_incr$queue_sup$x_sup[(dt_log_incr$queue_sup$y_sup<= -.8 & 
                                                          dt_log_incr$queue_sup$y_sup >= -2)]


R2_reg_inf = unlist(lapply(1:length(intervalle_recherche_inf), function(i) 
  dt_log_incr$queue_inf %>% .[.$x_inf >= intervalle_recherche_inf[i], ] %>% {lm(.$y_inf ~ .$x_inf)} %>%
    summary %>% .$adj.r.squared)) %>% {list("R2" = max(.), "Seuil" = intervalle_recherche_inf[which.max(.)])}

R2_reg_sup = unlist(lapply(1:length(intervalle_recherche_sup), function(i) 
  dt_log_incr$queue_sup %>% .[.$x_sup >= intervalle_recherche_sup[i], ] %>% {lm(.$y_sup ~ .$x_sup)} %>%
    summary %>% .$adj.r.squared)) %>% {list("R2" = max(.), "Seuil" = intervalle_recherche_sup[which.max(.)])}

x_inf_seuil = -1.091891 ; 
x_sup_seuil = -1.195444 ; 
reg_inf = dt_log_incr$queue_inf %>% .[.$x_inf >= x_inf_seuil, ] %>% {lm(.$y_inf ~ .$x_inf)}
reg_sup = dt_log_incr$queue_sup %>% .[.$x_sup >= x_sup_seuil, ] %>% {lm(.$y_sup ~ .$x_sup)}  # 0.9822
slope_sup = reg_sup$coefficients[2] ; int_sup =  reg_sup$coefficients[1]
slope_inf = reg_inf$coefficients[2] ; int_inf = reg_inf$coefficients[1]

figure3.1.3 = figure3.1.1 + geom_abline(slope = slope_inf, intercept = int_inf, col = "darkblue", size = 1.5) + 
  geom_vline(xintercept = x_inf_seuil, col = "darkred", size = 1.5) + theme_cowplot() + 
  theme(axis.text = element_text(size = 45)) 
figure3.1.4 = figure3.1.2 + geom_abline(slope = slope_sup, intercept = int_sup, col = "darkblue", size = 1.5) + 
  geom_vline(xintercept = x_sup_seuil, col = "darkred", size = 1.5) + theme_cowplot() + 
  theme(axis.text = element_text(size = 45)) 

# ggsave(paste0(path_save,"Figure3_1_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.3)
# ggsave(paste0(path_save,"Figure3_1_4.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.4)

# Figure 3.1. 5&6&7 - Fonction de répartition de la loi hybride -------------------------------------------


x_inf_lien =  dt_log_incr$queue_inf$x_inf[dt_log_incr$queue_inf$x_inf >= x_inf_seuil][which.min(
  abs((slope_inf * (dt_log_incr$queue_inf %>% {.[.$x_inf >= x_inf_seuil, ]$x_inf}) + int_inf) -
        dt_log_incr$queue_inf$y_inf[dt_log_incr$queue_inf$x_inf >= x_inf_seuil]))]

x_sup_lien =  dt_log_incr$queue_sup$x_sup[dt_log_incr$queue_sup$x_sup >= x_sup_seuil][which.min(
  abs((slope_sup * (dt_log_incr$queue_sup %>% {.[.$x_sup >= x_sup_seuil, ]$x_sup}) + int_sup) -
        dt_log_incr$queue_sup$y_sup[dt_log_incr$queue_sup$x_sup >= x_sup_seuil]))]


fit_JSUo =  gamlssML(yield$increment, family = "JSUo") %>% {list("mu" = .$mu, "sigma" = .$sigma, "nu" = .$nu,
                                                      "tau" = .$tau)}

q = seq(-.5,.5,length.out = 10^4)

ecdf_hybride = pHybrid(q = q, x_inf_lien_fn = x_inf_lien, x_sup_lien_fn = x_sup_lien,
                        int_inf_fn = int_inf, int_sup_fn = int_sup, slope_inf_fn = slope_inf, slope_sup_fn = slope_sup,par_JSUo_fn = fit_JSUo)

ecdf_JSUo = pJSUo(q = q,mu = fit_JSUo$mu, sigma = fit_JSUo$sigma, nu = fit_JSUo$nu, tau = fit_JSUo$tau)

figure3.1.5 = cbind(q,ecdf_hybride,ecdf_JSUo) %>% as.data.table %>% setnames(c("Support","Proba","JSUo")) %>%
  ggplot(aes(Support,Proba)) + geom_line() + geom_line(aes(y = JSUo), col = "darkred") +
  geom_vline(xintercept = -exp(x_inf_lien*log(10))) + geom_vline(xintercept = exp(x_sup_lien*log(10))) + 
  ylab("Probabilité cumulée") + theme(axis.text = element_text(size = 40), axis.title = element_text(size = 30))

figure3.1.6 = cbind(q,ecdf_hybride,ecdf_JSUo) %>% as.data.table %>% setnames(c("Support","Proba","JSUo")) %>%
  ggplot(aes(Support,Proba)) + geom_line(size = 2) + geom_line(aes(y = JSUo), col = "darkred", size = 2) +
  geom_vline(xintercept = -exp(x_inf_lien*log(10)),size = 2) + geom_vline(xintercept = exp(x_sup_lien*log(10)), size = 2) + 
  ylab("Probabilité cumulée") + theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30)) + 
  xlim(c(0.070,0.072)) + ylim(c(0.965,0.97)) + xlab('') + ylab("")

figure3.1.7 = cbind(q,ecdf_hybride,ecdf_JSUo) %>% as.data.table %>% setnames(c("Support","Proba","JSUo")) %>%
  ggplot(aes(Support,Proba)) + geom_line(size = 2) + geom_line(aes(y = JSUo), col = "darkred", size = 2) +
  geom_vline(xintercept = -exp(x_inf_lien*log(10)),size = 2) + geom_vline(xintercept = exp(x_sup_lien*log(10)), size = 2) + 
  ylab("Probabilité cumulée") + theme(axis.text = element_text(size = 45), axis.title = element_text(size = 30)) + 
  xlim(c(-0.1,-0.09855)) + ylim(c(7.63*10^-3,8.25*10^-3)) + xlab('') + ylab("")


# ggsave(paste0(path_save,"Figure3_1_5.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.5)
# ggsave(paste0(path_save,"Figure3_1_6.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.6)
# ggsave(paste0(path_save,"Figure3_1_7.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.7)




# Figure 3.1.8 - Projection IS loi hybride comparée à JSUo -------------------------------------------------------

proj_IS_hybride = projection_single_modele(m = 10^4, modele = "Hybride")
proj_IS_JSUo = projection_single_modele(m = 10^4, modele = "Gaussien")

y_list = list("Hybride" = proj_IS_hybride$yield, "JSUo" = proj_IS_JSUo$yield)

IC_95 = gg_multi_IC(y = yield_all_years$`10_ans`, y_list = y_list, date = yield_all_years$date, 
                          probs = c(0.025,0.975))
IC_99_5 = gg_multi_IC(y = yield_all_years$`10_ans`, y_list = y_list, date = yield_all_years$date, 
                      probs = c(0.0025,0.9975))


figure3.1.8 = cbind(IC_95$IC,IC_99_5$IC) %>% setnames(c("l1_hyb","l1_j","u1_hyb","u1_j","l2_hyb","l2_j","u2_hyb","u2_j")) %>%
  ggplot(aes(x = yield_all_years$date)) + geom_line(aes(y=l1_hyb, col = "Hybride")) + geom_line(aes(y=l1_j, col = "Gaussien")) +
  geom_line(aes(y=u1_hyb, col = "Hybride")) + geom_line(aes(y=u1_j, col = "Gaussien")) + geom_line(aes(y=l2_hyb, col = "Hybride")) +
  geom_line(aes(y=l2_j, col = "Gaussien")) + geom_line(aes(y=u2_hyb, col = "Hybride")) + geom_line(aes(y=u2_j, col = "Gaussien")) + 
  geom_line(aes(y = yield_all_years$`10_ans`)) + xlab("Date") + ylab("Taux") + 
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 35), legend.text = element_text(size = 40),
        legend.title = element_blank()) + guides(colour = guide_legend(override.aes = list(size=10)))

# ggsave(paste0(path_save,"Figure3_1_8.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.8)



# Figure 3.1. 9&10&11&12 - Projection OoS loi hybride comparée à normale -------------------------------------------------------

a1 = Sys.time() # 1H
proj_OoS_hybride = proj_OoS_par(m = 10^4, pas = 10, modeles = "Hybride")
a2 = Sys.time() - a1

proj_OoS_normale = proj_OoS_par(m = 10^4, pas = 10, modeles = "Gaussien")

figure3.1.9 = (proj_OoS_hybride$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")

figure3.1.10 = (proj_OoS_normale$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")

figure3.1.11 = (proj_OoS_hybride$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")

figure3.1.12 = (proj_OoS_normale$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")

# ggsave(paste0(path_save,"Figure3_1_9.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.9)
# ggsave(paste0(path_save,"Figure3_1_10.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.10)
#   ggsave(paste0(path_save,"Figure3_1_11.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.11)
#   ggsave(paste0(path_save,"Figure3_1_12.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.12)

proj_OoS_hybride$Nbre_erreur_95 %>% apply(2,mean)/1000
proj_OoS_JSUo$Nbre_erreur_95 %>% apply(2,mean)/1000

proj_OoS_hybride$Nbre_erreur_99_5 %>% apply(2,mean)/1000
proj_OoS_JSUo$Nbre_erreur_99_5 %>% apply(2,mean)/1000

# Figure 3.1. 13&14&15 - Faisceaux de queues de distribution --------------------------------------------

# queues_JSUo = gg_all_log_log_plot(m = 10^3, modele = "JSUo", pas = 10) # Figure2.3.15
queues_Hybride = gg_all_log_log_plot(m = 10^3, modele = "Hybride", pas = 10)
# queues_empirique = gg_all_log_log_plot(m = 10^3, modele = "Empirique", pas = 10) # Figure2.3.21


figure3.1.15 = queues_Hybride$g_inf + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + 
  theme(axis.text = element_text(size = 40))

ggsave(paste0(path_save,"Figure3_1_15.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure3.1.15)





# 3.3 - Du taux 10 ans aux autres maturités --------------------------------------------------

# Figure 3.3.1 - Structure de corrélation ---------------------------------------------------


maturite = 1:30

f = list(family = "Courier New, monospace", size = 18, color = "#7f7f7f")
xlab = list(title = "Maturité")
ylab = list(title = "Maturité")
zlab = list(title = "Corrélation")

data = as.matrix(yield_all_years[,-31]) %>% cor

plot_ly(x = maturite, y = maturite ,z = data) %>% 
  add_surface(contours = list( z = list(show=TRUE, usecolormap=TRUE, highlightcolor="#ff0000", project=list(z=TRUE)))) %>%
  layout(scene = list(camera=list(eye = list(x=1.87, y=0.88, z=-0.64)), xaxis = xlab,
      yaxis = ylab, zaxis = zlab), title = "Structure par terme - Corrélation")
# 
# data = as.matrix(yield_all_years[,-31]) %>% {.[2:3664,] - .[1:3663,]} %>% cor
# 
# plot_ly(x = maturite, y = maturite ,z = data) %>% 
#   add_surface(contours = list( z = list(show=TRUE, usecolormap=TRUE, highlightcolor="#ff0000", project=list(z=TRUE)))) %>%
#   layout(scene = list(camera=list(eye = list(x=1.87, y=0.88, z=-0.64)), xaxis = xlab,
#                       yaxis = ylab, zaxis = zlab), title = "Incréments - Corrélation")

# -----------------------------------------------------------------------------------------------------
# 4 - Approche g?n?rale de la structure par terme : mod?le de Diebold-Li et extensions ----------------
# 4.1 - De Nelson Siegel ? Diebold-Li : contexte


# 4.1 - De Nelson-Siegel ? Diebold-Li : contexte ------------------------------------------------------------------


# Figure 4.1. 1&2&3 - Level twist et butterfly ---------------------------------------------

x = seq(0,5,length.out = 10^4)

y = 1-exp(-x*1.5/2)
y_bis = 1-exp(-x^2*1.5)
y_level = y + 0.1
y_twist = 1 - exp(-x)
y_but = y_twist + sin(x*1.4)/7

figure4.4.1 = cbind(x,y,y_level) %>% as.data.table %>% setnames(c("x","y1","y2")) %>% ggplot(aes(x,y1)) + 
  geom_line() + geom_line(aes(y=y2)) + theme_cowplot() + theme(axis.text = element_blank()) + xlab("") + ylab("")

figure4.4.2 = cbind(x,y_bis,y_twist) %>% as.data.table %>% setnames(c("x","y1","y2")) %>% ggplot(aes(x,y1)) + 
  geom_line() + geom_line(aes(y=y2)) + theme_cowplot() + theme(axis.text = element_blank()) + xlab("") + ylab("")

figure4.4.3 = cbind(x,y_twist,y_but) %>% as.data.table %>% setnames(c("x","y1","y2")) %>% ggplot(aes(x,y1)) + 
  geom_line() + geom_line(aes(y=y2)) + theme_cowplot() + theme(axis.text = element_blank()) + xlab("") + ylab("")

# ggsave(paste0(path_save,"Figure4_4_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.1)
# ggsave(paste0(path_save,"Figure4_4_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.2)
# ggsave(paste0(path_save,"Figure4_4_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.3)



# Figure 4.1.4 - Loadings ------------------------------------------

loadings = function(t,lambda){
  load_b1 = rep(1,max(length(t),length(lambda),length(t)*length(lambda)))
  load_b2 = (1-exp(-lambda*t))/(lambda*t)
  load_b3 = (1-exp(-lambda*t))/(lambda*t) -exp(-lambda*t)
  return(matrix(c(load_b1,load_b2,load_b3),ncol = 3))
}

t = seq(0,30,length.out = 10^4) ; lambda = 0.598

figure4.1.4 = cbind(t, loadings(t, lambda)) %>% as.data.table %>% setnames(c('Maturite',"B1","B2","B3")) %>%
  ggplot(aes(Maturite,B1)) + geom_line(col = c("black")) + geom_line(aes(y=B2), col = "darkblue") + 
  geom_line(aes(y=B3), col = "darkred")   + xlab("Maturité") + theme(axis.text = element_text(size = 40),
                                                                     axis.title = element_text(size = 40)) + ylab("Loading")

# ggsave(paste0(path_save,"Figure4_1_4.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.1.4)

# Figure 4.1. 5&6&7 - Betas empiriques vs Beta estim?s -----------------------------------------------

figure_4.1.5 = cbind.data.frame(yield_all_years$date, beta_estim[,1],beta_optim[,1]) %>% as.data.table %>% 
  setnames(c("Date","Beta_estim","Beta_optim")) %>%
  ggplot(aes(Date,Beta_estim)) + geom_line(aes(col = "Beta_estim")) + geom_line(aes(y = Beta_optim,col = "Beta_optim")) + theme_cowplot() + 
  theme(axis.text = element_text(size = 45)) + xlab("") + ylab("")

figure_4.1.6 = cbind.data.frame(yield_all_years$date, beta_estim[,2],beta_optim[,2]) %>% as.data.table %>% 
  setnames(c("Date","Beta_estim","Beta_optim")) %>%
  ggplot(aes(Date,Beta_estim)) + geom_line(aes(col = "Beta_estim")) + geom_line(aes(y = Beta_optim,col = "Beta_optim")) + theme_cowplot() + 
  theme(legend.position = "none", axis.text = element_text(size = 45)) + xlab("") + ylab("")

figure_4.1.7 = cbind.data.frame(yield_all_years$date, beta_estim[,3],beta_optim[,3]) %>% as.data.table %>% 
  setnames(c("Date","Beta_estim","Beta_optim")) %>%
  ggplot(aes(Date,Beta_estim)) + geom_line(aes(col = "Beta_estim")) + geom_line(aes(y = Beta_optim,col = "Beta_optim")) + theme_cowplot() + 
  theme(legend.position = "none", axis.text = element_text(size = 45)) + xlab("") + ylab("")

# ggsave(paste0(path_save,"Figure4_1_5.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure_4.1.5)
# ggsave(paste0(path_save,"Figure4_1_6.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure_4.1.6)
# ggsave(paste0(path_save,"Figure4_1_7.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure_4.1.7)


# Figure 4.1.bis - Surface des résidus NS --------------------------------------------



f = list(family = "Courier New, monospace", size = 18, color = "#7f7f7f")
xlab = list(title = "Année")
ylab = list(title = "Maturité")
zlab = list(title = "Résidus")

dates = yield_all_years$date
maturite = 1:30

data = residus_NS %>% unlist %>% matrix(ncol = 30) %>% t

figure4_1_bis = plot_ly(x = dates, y = maturite ,z = as.matrix(data)) %>% add_surface(contours = list(z = list(show=TRUE, usecolormap=TRUE,
    highlightcolor="#ff0000", project=list(z=TRUE)))) %>% layout(scene = list(camera=list(eye = list(x=1.87, y=0.88, z=-0.64))
      , xaxis = xlab, yaxis = ylab, zaxis = zlab), title = "")
figure4_1_bis



# Figure 4.1. 8&9 - Yield 10 ans reconstruit plus r?sidus ------------------------------------------

y_th = unlist(lapply(1:3664, function(i) NS_YC(t = 10, par = beta_optim[i,], lambda = 0.598)))

figure4.1.8 = cbind.data.frame(yield_all_years$date, yield_all_years$`10_ans`,y_th) %>% as.data.table %>% 
  setnames(c("Date","Y_emp","Y_th")) %>%
  ggplot(aes(Date,Y_emp)) + geom_line(aes(col = "Y_emp")) + geom_line(aes(y = Y_th,col = "Y_th")) + theme_cowplot() + 
  theme(legend.position = "none", axis.text = element_text(size = 45)) + xlab("") + ylab("")

figure4.1.9 = cbind.data.frame(yield_all_years$date, residus_NS) %>% as.data.table %>% 
  setnames(c("Date","residus_NS")) %>%
  ggplot(aes(Date,residus_NS)) + geom_line() + theme_cowplot() + 
  theme(legend.position = "none", axis.text = element_text(size = 45)) + xlab("") + ylab("")

# ggsave(paste0(path_save,"Figure4_1_8.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.1.8)
# ggsave(paste0(path_save,"Figure4_1_9.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.1.9)




# 4.2 - Altération du modèle de Diebold-Li ------------------------------------------------------------------------------

  # Figure 4.2. 1&2&3&4 - Proportions de modèles fittés betas et résidus --------------------------------------------------
  
  incr_beta = (beta_optim[2:3664,] - beta_optim[1:3663,])
  incr_residus = residus_NS[2:3664] - residus_NS[1:3663]
  
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
  
  # Beta1
  figure4.2.1 = ggplot(data = modelRep1, aes(x = AR, y = MA)) + theme_cowplot() +
    geom_count(alpha = .4, show.legend = F, color = "darkblue") 
  figure4.2.1 = figure4.2.1 + geom_text(data = ggplot_build(figure4.2.1)$data[[1]], 
                                        aes(x, y, label = round(prop,2)), color = "darkred",size = 15, hjust=0.5, vjust=-1) + 
    theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45))
  figure4.2.1 = figure4.2.1 + ylim(c(0,6)) + scale_size_continuous(range = c(5, 20))
  
  # Beta2
  figure4.2.2 = ggplot(data = modelRep2, aes(x = AR, y = MA)) + theme_cowplot() +
    geom_count(alpha = .4, show.legend = F, color = "darkblue") 
  figure4.2.2 = figure4.2.2 + geom_text(data = ggplot_build(figure4.2.2)$data[[1]], 
                                        aes(x, y, label = round(prop,2)), color = "darkred",size = 15, hjust=0.5, vjust=-1) + 
    theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45))
  figure4.2.2 = figure4.2.2 + ylim(c(0,6)) + scale_size_continuous(range = c(5, 20))
  
  # Beta3
  figure4.2.3 = ggplot(data = modelRep3, aes(x = AR, y = MA)) + theme_cowplot() +
    geom_count(alpha = .4, show.legend = F, color = "darkblue") 
  figure4.2.3 = figure4.2.3 + geom_text(data = ggplot_build(figure4.2.3)$data[[1]], 
                                        aes(x, y, label = round(prop,2)), color = "darkred",size = 15, hjust=0.5, vjust=-1) + 
    theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45))
  figure4.2.3 = figure4.2.3 + ylim(c(0,6)) + scale_size_continuous(range = c(5, 20))
  
  # Residus
  figure4.2.4 = ggplot(data = modelRepRes, aes(x = AR, y = MA)) + theme_cowplot() +
    geom_count(alpha = .4, show.legend = F, color = "darkblue") 
  figure4.2.4 = figure4.2.4 + geom_text(data = ggplot_build(figure4.2.4)$data[[1]], 
                                        aes(x, y, label = round(prop,2)), color = "darkred",size = 15, hjust=0.5, vjust=-1) + 
    theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45))
  figure4.2.4 = figure4.2.4 + ylim(c(0,6)) + scale_size_continuous(range = c(5, 20))
  
  # ggsave(paste0(path_save,"Figure4_2_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.2.1)
  # ggsave(paste0(path_save,"Figure4_2_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.2.2)
  # ggsave(paste0(path_save,"Figure4_2_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.2.3)
  # ggsave(paste0(path_save,"Figure4_2_4.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.2.4)
  
  
# Figure 4.2. 5&6&7&8 - Valeurs d'AIC dans les différentes configurations ----------------------------------------------------

pas = 1000
iteration = seq(1,(3663-pas),by=1)

AIC_1 = foreach(i = 1:(3663-pas), .packages = c("forecast", "gamlss"), .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_beta[i:(pas+i),1]
  
  ma1 = arima(incr, order = c(0,0,1))
  ar2 = arima(incr, order = c(2,0,0))
  normal = arima(incr, order = c(0,0,0))
  JSUo = incr%>% gamlssML(family = "JSUo") %>% {dJSUo(incr, mu = .$mu, sigma = .$sigma, nu = .$nu, tau = .$tau)} %>% log %>% 
  {2*4-2*sum(.)} 
  
  c(ma1$aic, ar2$aic,normal$aic, JSUo)
}

AIC_1 = as.data.table(AIC_1) ; colnames(AIC_1) = c("MA1","AR2","Normal","JSUo") ; 
AIC_1$Date = unique(yield_all_years$date)[1001:3663]


AIC_2 = foreach(i = 1:(3663-pas), .packages = c("forecast", "gamlss"), .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_beta[i:(pas+i),2]
  
  ma1 = arima(incr, order = c(0,0,1))
  ma2 = arima(incr, order = c(0,0,2))
  ar1 = arima(incr, order = c(1,0,0))
  ar2 = arima(incr, order = c(2,0,0))
  normal = arima(incr, order = c(0,0,0))
  JSUo = incr%>% gamlssML(family = "JSUo") %>% {dJSUo(incr, mu = .$mu, sigma = .$sigma, nu = .$nu, tau = .$tau)} %>% log %>% 
  {2*4-2*sum(.)} 
  
  c(ma1$aic,ma2$aic, ar1$aic, ar2$aic,normal$aic, JSUo)
}

AIC_2 = as.data.table(AIC_2) ; colnames(AIC_2) = c("MA1","MA2","AR1","AR2","Normal","JSUo") ; 
AIC_2$Date = unique(yield_all_years$date)[1001:3663]

AIC_3 = foreach(i = 1:(3663-pas), .packages = c("forecast", "gamlss"), .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_beta[i:(pas+i),3]
  
  ma1 = arima(incr, order = c(0,0,1))
  ar1 = arima(incr, order = c(1,0,0))
  ar2 = arima(incr, order = c(2,0,0))
  normal = arima(incr, order = c(0,0,0))
  JSUo = incr%>% gamlssML(family = "JSUo") %>% {dJSUo(incr, mu = .$mu, sigma = .$sigma, nu = .$nu, tau = .$tau)} %>% log %>% 
  {2*4-2*sum(.)} 
  
  c(ma1$aic, ar1$aic, ar2$aic,normal$aic, JSUo)
}

AIC_3 = as.data.table(AIC_3) ; colnames(AIC_3) = c("MA1","AR1","AR2","Normal","JSUo") ; 
AIC_3$Date = unique(yield_all_years$date)[1001:3663]


AIC_res = foreach(i = 1:(3663-pas), .packages = c("forecast", "gamlss"), .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_residus[i:(pas+i)]
  
  ar1 = arima(incr, order = c(1,0,0))
  ar2 = arima(incr, order = c(2,0,0))
  ar3 = arima(incr, order = c(3,0,0))
  normal = arima(incr, order = c(0,0,0))
  JSUo = incr%>% gamlssML(family = "JSUo") %>% {dJSUo(incr, mu = .$mu, sigma = .$sigma, nu = .$nu, tau = .$tau)} %>% log %>% 
  {2*4-2*sum(.)} 
  
  c(ar1$aic, ar2$aic, ar3$aic,normal$aic, JSUo)
}

AIC_res = as.data.table(AIC_res) ; colnames(AIC_res) = c("AR1","AR2","AR3","Normal","JSUo") ; 
AIC_res$Date = unique(yield_all_years$date)[1001:3663]

AIC_1[,-"Date"] %>% apply(2,mean)
AIC_1[,-"Date"] %>% apply(1,which.min) %>% {list("MA1" = sum(.==1), "AR2" = sum(.==2), "Normal" = sum(.==3), "JSUo" = sum(.==4))}

AIC_2[,-"Date"] %>% apply(2,mean)
AIC_2[,-"Date"] %>% apply(1,which.min) %>% {list("MA1" = sum(.==1), "MA2" = sum(.==2), "AR1" = sum(.==3), "AR2" = sum(.==4),
                                                 "Normal" = sum(.==5),"JSUo" = sum(.==6))}

AIC_3[,-"Date"] %>% apply(2,mean)
AIC_3[,-"Date"] %>% apply(1,which.min) %>% {list("MA1" = sum(.==1), "AR1" = sum(.==2), "AR2" = sum(.==3), "Normal" = sum(.==4),
                                                 "JSUo" = sum(.==5))}

AIC_res[,-"Date"] %>% apply(2,mean)
AIC_res[,-"Date"] %>% apply(1,which.min) %>% {list("AR1" = sum(.==1), "AR2" = sum(.==2), "AR3" = sum(.==3), "Normal" = sum(.==4),
                                                   "JSUo" = sum(.==5))}

# Figure 4.2.9 - Independance des beta et résidus -------------------------------------------------------------------------------------------

independance = foreach(i = 1:(3663-pas), .packages = "forecast", .combine = rbind) %dopar% {
  
  j = which(iteration==i)
  incr = incr_beta[i:(pas+i),] ; incr_res = incr_residus[i:(pas+i)]
  
  ind1 = Box.test(incr[,1], lag = 1)$p.value
  ind2 = Box.test(incr[,2], lag = 1)$p.value
  ind3 = Box.test(incr[,3], lag = 1)$p.value
  indRes = Box.test(incr_res, lag = 1)$p.value
  
  c(ind1, ind2, ind3, indRes)
}

figure4.2.9 = independance %>% cbind(1:nrow(.)) %>% as.data.table %>% setnames(c("B1","B2","B3","Residus","Simulations")) %>% 
  melt(id.vars = "Simulations") %>%
  ggplot(aes(Simulations,value)) + geom_line(aes(col = variable)) + theme_cowplot() + ylab("p-value") + theme(axis.text = element_text(size = 40),
                                                                                                        axis.title = element_text(size = 40),
                                                                                                        legend.title = element_text(size = 40), 
                                                                                                        legend.text = element_text(size = 40)) +
  scale_color_manual(name = "", values = gg_color_hue(4))  + guides(colour = guide_legend(override.aes = list(size=10)))

# ggsave(paste0(path_save,"Figure4_2_9.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.2.9)









# 4.3 - Ajout d'une structure de d?pendance sur les facteurs et les r?sidus -----------------------------------------

# Figure 4.3.1 - Correlations des facteurs ------------------------------------------------------

data_cross = beta_optim
data_cross_plot = beta_optim[500:3664,] ; data_cross_plot = as.data.table(data_cross_plot)

cor = c("Beta_1_2", "Beta_2_3", "Beta_1_3")
data_cross_plot = cbind(data_cross_plot,matrix(rep(0,3*nrow(data_cross_plot)),ncol=3))
colnames(data_cross_plot) = c(c("Beta1","Beta2","Beta3"),cor)

for(i in 1:(3664-500)){
  cor_mat = cor(data_cross[i:(i+99),])
  data_cross_plot[i,4] = cor_mat[1,2]
  data_cross_plot[i,5] = cor_mat[1,3]
  data_cross_plot[i,6] = cor_mat[2,3]
  
}


data_cross_plot$"Date" = yield_all_years$date[500:3664]
data_cross_plot_yield = cbind.data.frame(beta_optim,yield_all_years$date) %>% setnames(c(1:3,'Date')) %>%
  melt(id.vars = "Date", variable.name = "series")
data_cross_plot_c = data_cross_plot[,-c(1,2,3)] %>% melt(id.vars = "Date", variable.name = "series") 

data_cross_plot_cor = data.frame("Date" = rep(yield_all_years$date[1:500],each = 3),"series" = rep(NA,1500),"value" = rep(NA,1500)) %>% 
  rbind(data_cross_plot_c)



figureCor = ggplot(data_cross_plot_yield, aes(Date,value)) + geom_line(aes(color = series)) +
  xlab("") + ylab("") + ggtitle("")  + ylab("") + 
  geom_vline(xintercept = yield_all_years$date[1255+100], col = "darkgreen", size = 1) + 
  geom_vline(xintercept = yield_all_years$date[756+100], col = "darkgreen", size = 1) +
  annotate("rect", xmax=yield_all_years$date[1255+100], xmin=yield_all_years$date[756+100],
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="darkgreen") + theme_cowplot() +
  theme(plot.margin = unit(c(-10,2,10,4.1),"mm"),legend.position = "none", axis.text = element_text(size = 40)) + 
  scale_color_manual(name = "", values = c("black",gg_color_hue(2)))  

figureYield = ggplot(data_cross_plot_cor, aes(Date,value)) + geom_line(aes(color = series)) +
  xlab("") + ylab("") + ggtitle("") + ylab("") + theme_cowplot() +
  geom_vline(xintercept = yield_all_years$date[1255+100], col = "darkgreen", size = 1) + 
  theme(plot.margin = unit(c(4,2,-5,4),"mm"), axis.text.x = element_blank(), legend.position = "none", axis.text = element_text(size = 40))  + 
  scale_color_manual(name = "", values = c("black",gg_color_hue(2))) 

figure4.3.1 = arrangeGrob(figureYield,figureCor, ncol = 1)

# ggsave(paste0(path_save,"Figure4_3_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.3.1)




# 4.4 - R?sultats IS / Out of Sample de ces mod?lisations ------------------------------------------------------

# Figure 4.4.1 - Resultats IS 3 mod?les DNS ----------------------------------------------------------
proj_IS_DNS_original = projection_single_modele(m = 10^4, modele = "DNS original")
proj_IS_DNS_JSUo = projection_single_modele(m = 10^4, modele = "DNS JSUo")
proj_IS_DNS_Copule = projection_single_modele(m = 10^4, modele = "DNS copule")


y_list = list("Original" = proj_IS_DNS_original$yield, "JSUo" = proj_IS_DNS_JSUo$yield, "Copule" = proj_IS_DNS_Copule$yield)

IC_95 = gg_multi_IC(y = yield_all_years$`10_ans`, y_list = y_list, date = yield_all_years$date, 
                    probs = c(0.025,0.975))
IC_99_5 = gg_multi_IC(y = yield_all_years$`10_ans`, y_list = y_list, date = yield_all_years$date, 
                      probs = c(0.0025,0.9975))

figure4.4.1 = cbind(IC_95$IC,IC_99_5$IC) %>% 
  setnames(c("l1_original","l1_j","l1_cop","u1_original","u1_j","u1_cop","l2_original","l2_j","l2_cop","u2_original","u2_j","u2_cop")) %>%
  ggplot(aes(x = yield_all_years$date)) + geom_line(aes(y=l1_original, col = "Original")) + geom_line(aes(y=l1_j, col = "JSUo")) +
  geom_line(aes(y=u1_original, col = "Original")) + geom_line(aes(y=u1_j, col = "JSUo")) + geom_line(aes(y=u2_original, col = "Original")) +
  geom_line(aes(y=l2_j, col = "JSUo")) + geom_line(aes(y=l2_original, col = "Original")) + geom_line(aes(y=u2_j, col = "JSUo")) + 
  geom_line(aes(y=l1_cop, col = "Copule")) + geom_line(aes(y=u1_cop, col = "Copule")) + geom_line(aes(y=l2_cop, col = "Copule")) +
  geom_line(aes(y=u2_cop, col = "Copule")) +
  geom_line(aes(y = yield_all_years$`10_ans`)) + xlab("Date") + ylab("Taux") + theme_cowplot() + 
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40),
        legend.title = element_blank()) + guides(colour = guide_legend(override.aes = list(size=5)))

# ggsave(paste0(path_save,"Figure4_4_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.1)

# Figure 4.4. 2&3 - Queues de distribution IS -----------------------------------------------------------------------

y_list = list("Original" = proj_IS_DNS_original$incr, "JSUo" = proj_IS_DNS_JSUo$incr, "Copule" = proj_IS_DNS_Copule$incr)

gg_log = gg_log_multi_comp(x = incr, y_list = y_list)


log_incr = log_log_dt(incr)
x_inf_025 = (round(log_incr$queue_inf$y_inf,2) == round(log10(0.025),2)) %>% {log_incr$queue_inf$x_inf[which(.)]} %>% mean
x_inf_0025 = (round(log_incr$queue_inf$y_inf,1) == round(log10(0.0025),1)) %>% {log_incr$queue_inf$x_inf[which(.)]} %>% mean

x_sup_975 = (round(log_incr$queue_sup$y_sup,2) == round(log10(0.025),2)) %>% {log_incr$queue_sup$x_sup[which(.)]} %>% mean
x_sup_9975 = (round(log_incr$queue_sup$y_sup,1) == round(log10(0.0025),1)) %>% {log_incr$queue_sup$x_sup[which(.)]} %>% mean


figure4.4.2 = gg_log$g_inf + theme_cowplot() + theme(legend.position = "none", axis.text = element_text(size = 40)) + xlab("") + ylab("") + 
  labs(tilte = "0") + geom_segment(aes(x = min(log_incr$queue_inf$x_inf), xend = x_inf_025, y = log10(0.025), yend = log10(0.025)), 
                                   linetype = "dashed") + 
  geom_segment(aes(x =x_inf_025, xend = x_inf_025, y = -7, yend = log10(0.025)), linetype = "dashed") + 
  geom_segment(aes(x = min(log_incr$queue_inf$x_inf), xend = x_inf_0025, y = log10(0.0025), yend = log10(0.0025)), linetype = "dashed") + 
  geom_segment(aes(x =x_inf_0025, xend = x_inf_0025, y = -7, yend = log10(0.0025)), linetype = "dashed") + labs(title = "")

figure4.4.3 = gg_log$g_sup + theme_cowplot() + theme(legend.position = "none", axis.text = element_text(size = 40)) + xlab("") + ylab("") + 
  labs(tilte = "0")  + geom_segment(aes(x = min(log_incr$queue_sup$x_sup), xend = x_sup_975, y = log10(0.025), yend = log10(0.025)), 
                                    linetype = "dashed") + 
  geom_segment(aes(x =x_sup_975, xend = x_sup_975, y = -7, yend = log10(0.025)), linetype = "dashed") + 
  geom_segment(aes(x = min(log_incr$queue_sup$x_sup), xend = x_sup_9975, y = log10(0.0025), yend = log10(0.0025)), linetype = "dashed") + 
  geom_segment(aes(x =x_sup_9975, xend = x_sup_9975, y = -7, yend = log10(0.0025)), linetype = "dashed") + labs(title = "")

# ggsave(paste0(path_save,"Figure4_4_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.2)
# ggsave(paste0(path_save,"Figure4_4_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.3)




# Figure 4.4. 4&5&6&7&8&9 - Resultats OoS 3 mod?les DNS ------------------------------------------------------------

proj_OoS_DNS_original = proj_OoS_par(m = 10^4, pas = 10, modeles = "DNS original")
proj_OoS_DNS_JSUo = proj_OoS_par(m = 10^3, pas = 10, modeles = "DNS JSUo")
proj_OoS_DNS_copule = proj_OoS_par(m = 10^4, pas = 10, modeles = "DNS copule")

# DNS original
figure4.4.4 = (proj_OoS_DNS_original$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() +
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("") + geom_hline(yintercept = 0.05, linetype = "dashed", col = "darkred")

figure4.4.5 = (proj_OoS_DNS_original$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() + 
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")


# DNS JSUo
figure4.4.6 = (proj_OoS_DNS_JSUo$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() + 
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("") + geom_hline(yintercept = 0.05, linetype = "dashed", col = "darkred")

figure4.4.7 = (proj_OoS_DNS_JSUo$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() + 
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")

# DNS copule
figure4.4.8 = (proj_OoS_DNS_copule$Nbre_erreur_95/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() + 
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("") + geom_hline(yintercept = 0.05, linetype = "dashed", col = "darkred")

figure4.4.9 = (proj_OoS_DNS_copule$Nbre_erreur_99_5/1000) %>% cbind(1:nrow(.)) %>% as.data.table %>% 
  melt(id.vars = "V2",  variable.name = "variable") %>% ggplot(aes(V2,value)) + geom_line() + xlab("") + theme_cowplot() + 
  theme(axis.text = element_text(size = 45), axis.title = element_text(size = 45),legend.position = "none") + geom_point(size = 2, alpha = .5) +
  ylim(c(0,1)) + ylab("")


# ggsave(paste0(path_save,"Figure4_4_4.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.4)
# ggsave(paste0(path_save,"Figure4_4_5.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.5)
# ggsave(paste0(path_save,"Figure4_4_6.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.6)
# ggsave(paste0(path_save,"Figure4_4_7.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.7)
# ggsave(paste0(path_save,"Figure4_4_8.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.8)
# ggsave(paste0(path_save,"Figure4_4_9.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.9)

proj_OoS_DNS_original$Nbre_erreur_95 %>% apply(2,mean)/1000
proj_OoS_DNS_JSUo$Nbre_erreur_95 %>% apply(2,mean)/1000
proj_OoS_DNS_copule$Nbre_erreur_95 %>% apply(2,mean)/1000

proj_OoS_DNS_original$Nbre_erreur_99_5 %>% apply(2,mean)/1000
proj_OoS_DNS_JSUo$Nbre_erreur_99_5 %>% apply(2,mean)/1000
proj_OoS_DNS_copule$Nbre_erreur_99_5 %>% apply(2,mean)/1000



# Figure 4.4. 10&11&12 - Faisceaux queues de distribution ------------------------------------------------------------

queues_DNS_orignal = gg_all_log_log_plot(m = 10^3, modele = "DNS original", pas = 10)
queues_DNS_JSUo = gg_all_log_log_plot(m = 10^3, modele = "DNS JSUo", pas = 10)
queues_DNS_copule = gg_all_log_log_plot(m = 10^3, modele = "DNS copule", pas = 10)

figure4.4.10 = queues_DNS_orignal$g_inf + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + theme_cowplot() + 
  theme(axis.text = element_text(size = 40), legend.position = "none")

figure4.4.11 = queues_DNS_JSUo$g_inf + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + theme_cowplot() +  
  theme(axis.text = element_text(size = 40), legend.position = "none")

figure4.4.12 = queues_DNS_copule$g_inf + xlim(c(-7,2)) + ylim(c(-7,0)) + labs(title = '') + xlab('') + ylab('') + theme_cowplot() +  
  theme(axis.text = element_text(size = 40), legend.position = "none")

# ggsave(paste0(path_save,"Figure4_4_10.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.10)
# ggsave(paste0(path_save,"Figure4_4_11.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.11)
# ggsave(paste0(path_save,"Figure4_4_12.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure4.4.12)


# -----------------------------------------------------------------------------------------------------

# 5 - Discussions sur la partie mod?lisation ----------------------------------------------------------

# 5.1 - Br?ve revue des mod?les -----------------------------------------------------------------------

# 5.1. 1&2&3 - Synth?se des r?sultats IS ------------------------------------------------------------------

couleur = c("darkred",gg_color_hue(6)[c(1:4,6)])

proj_IS_tous = projection_modeles(m = 10^4, modeles = c("Gaussien", "naive gaussien", "Hybride", "DNS original", "DNS JSUo","DNS copule"))


figure5.1.1 = proj_IS_tous$g_log_inf + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.1.2 = proj_IS_tous$g_IC +  theme_cowplot() +
  theme(plot.title = element_text(size = 45), axis.text = element_text(size = 45), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.1.3 = proj_IS_tous$g_log_sup +  theme_cowplot() +
  theme(plot.title = element_text(size = 45), axis.text = element_text(size = 45), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

# ggsave(paste0(path_save,"Figure5_1_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.1.1)
# ggsave(paste0(path_save,"Figure5_1_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.1.2)
# ggsave(paste0(path_save,"Figure5_1_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.1.3)

# 5.2 - Comparaisons et impacts des diff?rentes mod?lisations ------------------------------------------------------------



# Figure5.2. 1&2&3&4 - Projection 2009 ----------------------------------------------------------------------------------

# Noir / Bleu / Orange / vert / Rouge / Rouge foncé
couleur = c("darkred",gg_color_hue(6)[c(1:4,6)])

proj_OoS_2009 = projection_modeles(start_calib = 1, n_calib = 1000, IS = F, n_proj = 1000, m =10^3, trace = T,
                                   modeles = c("Gaussien", "naive gaussien", "Hybride", "DNS original", "DNS JSUo","DNS copule"))

figure5.2.1 = proj_OoS_2009$g_log_inf + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.2 = proj_OoS_2009$g_IC +theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.3 = proj_OoS_2009$g_log_sup + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.4 = proj_OoS_2009$g_fen +  theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

# ggsave(paste0(path_save,"Figure5_2_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.1)
# ggsave(paste0(path_save,"Figure5_2_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.2)
# ggsave(paste0(path_save,"Figure5_2_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.3)
# ggsave(paste0(path_save,"Figure5_2_4.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.4)


# Figure5.2. 5&6&7&8 - Projection 2013 --------------------------------------------------------------------------------

proj_OoS_2013 = projection_modeles(start_calib = 1000, n_calib = 1000, IS = F, n_proj = 1000, m =10^3, trace = T,
                                   modeles = c("Gaussien", "naive gaussien", "Hybride", "DNS original", "DNS JSUo","DNS copule"))

figure5.2.5 = proj_OoS_2013$g_log_inf + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.6 = proj_OoS_2013$g_IC + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.7 = proj_OoS_2013$g_log_sup + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.8 = proj_OoS_2013$g_fen + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

# ggsave(paste0(path_save,"Figure5_2_5.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.5)
# ggsave(paste0(path_save,"Figure5_2_6.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.6)
# ggsave(paste0(path_save,"Figure5_2_7.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.7)
# ggsave(paste0(path_save,"Figure5_2_8.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.8)

# Figure5.2. 9&10&11&12 - Projection 2013 - 8 ans -------------------------------------------------------------------
proj_OoS_2013_max = projection_modeles(start_calib = 1000, n_calib = 1000, IS = F, n_proj = 1663, m =10^3, trace = T,
                                       modeles = c("Gaussien", "naive gaussien", "Hybride", "DNS original", "DNS JSUo","DNS copule"))

figure5.2.9 = proj_OoS_2013_max$g_log_inf + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.10 = proj_OoS_2013_max$g_IC + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.11 = proj_OoS_2013_max$g_log_sup + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.12 = proj_OoS_2013_max$g_fen + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)


# ggsave(paste0(path_save,"Figure5_2_9.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.9)
# ggsave(paste0(path_save,"Figure5_2_10.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.10)
# ggsave(paste0(path_save,"Figure5_2_11.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.11)
# ggsave(paste0(path_save,"Figure5_2_12.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.12)

# Figure5.2. 13&14&15&16 - Projection 2007 - 4 ans -----------------------------------------------------------------------
proj_OoS_2007_4 = projection_modeles(start_calib = 1, n_calib = 500, IS = F, n_proj = 1000, m =10^3, trace = T,
                                     modeles = c("Gaussien", "naive gaussien", "Hybride", "DNS original", "DNS JSUo","DNS copule"))

figure5.2.13 = proj_OoS_2007_4$g_log_inf + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.14 = proj_OoS_2007_4$g_IC + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.15 = proj_OoS_2007_4$g_log_sup + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

figure5.2.16 = proj_OoS_2007_4$g_fen + theme_cowplot() +
  theme(plot.title = element_text(size = 40), axis.text = element_text(size = 40), legend.position = "none") + 
  xlab("") + ylab("") + labs(title = "") + scale_color_manual(values = couleur)

# ggsave(paste0(path_save,"Figure5_2_13.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.13)
# ggsave(paste0(path_save,"Figure5_2_14.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.14)
# ggsave(paste0(path_save,"Figure5_2_15.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.15)
# ggsave(paste0(path_save,"Figure5_2_16.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.2.16)


# Figurer MSCI -----------------------------------

MSCI$Dernier = as.numeric(MSCI$Dernier)
MSCI$Date = as.Date(MSCI$Date, "%d/%m/%Y")

incr_log = 
  (MSCI$Dernier[2:length(MSCI$Dernier)]/MSCI$Dernier[1:(length(MSCI$Dernier)-1)]) %>%
  log

data = cbind.data.frame(MSCI$Date[2:nrow(MSCI)],incr_log) %>% 
  setnames(c("Date","Rdt")) 
ggplot(data,aes(Date,Rdt)) + geom_line()  


dt_log_exp = function(indice){
  incr = (indice[2:(length(indice))] / indice[1:(length(indice)-1)]) %>% log
  ecdf_incr = ecdf(incr)
  incr_pos = incr[incr>0] %>% sort ; incr_neg = incr[incr<0] %>% sort 
  
  res = list("Queue_inf" = data.table(
    "x_inf" = -incr_neg %>% log10,
    "y_inf" = ecdf_incr(incr_neg) %>% log10
  ),
  "Queue_sup" = data.table(
    "x_sup" = incr_pos %>% log10,
    "y_sup" = (1-ecdf_incr(incr_pos)) %>% log10
  ))
  
  res$Queue_sup = res$Queue_sup[-nrow(res$Queue_sup),]
  return(res)
}

a = dt_log_exp(MSCI$Dernier)
ggplot(a$Queue_inf, aes(x_inf,y_inf)) + geom_line()
ggplot(a$Queue_sup, aes(x_sup,y_sup)) + geom_line()


# Recherche  meilleure r??gression
recherche = seq(-.8,-2,length.out = 10^4)

reg_inf = lapply(1:length(recherche), function(i) 
  a$Queue_inf[a$Queue_inf$y_inf<recherche[i],] %>%
  {lm(.$y_inf ~ .$x_inf)} %>% 
  {list("R2" = summary(.)$adj.r.squared,
        "Coefficients" = list("slope" = .$coefficients[2],
                              "intercept" = .$coefficients[1]))})

reg_sup = lapply(1:length(recherche), function(i) 
  a$Queue_sup[a$Queue_sup$y_sup<recherche[i],] %>%
  {lm(.$y_sup ~ .$x_sup)} %>% 
  {list("R2" = summary(.)$adj.r.squared,
        "Coefficients" = list("slope" = .$coefficients[2],
                              "intercept" = .$coefficients[1]))})

reg_inf_bis= reg_inf
reg_inf = reg_inf_bis[which.max(map(reg_inf,1))-50]
reg_sup = reg_sup[which.max(map(reg_sup,1))]

# Recherche du meilleur point de rencontre :

int_inf_0 = a$Queue_inf$x_inf[which.min(abs((reg_inf[[1]]$Coefficients$intercept +
                                               a$Queue_inf$x_inf*
                                               reg_inf[[1]]$Coefficients$slope) -
                                              a$Queue_inf$y_inf))]

int_sup_0 = a$Queue_sup$x_sup[which.min(abs((reg_sup[[1]]$Coefficients$intercept + 
                                               a$Queue_sup$x_sup*
                                               reg_sup[[1]]$Coefficients$slope) - 
                                              a$Queue_sup$y_sup))]

y_0_inf = a$Queue_inf$y_inf[which.min(abs((reg_inf[[1]]$Coefficients$intercept +
                                             a$Queue_inf$x_inf*
                                             reg_inf[[1]]$Coefficients$slope) -
                                            a$Queue_inf$y_inf))]


y_0_sup = a$Queue_sup$y_sup[which.min(abs((reg_sup[[1]]$Coefficients$intercept + 
                                             a$Queue_sup$x_sup*
                                             reg_sup[[1]]$Coefficients$slope) - 
                                            a$Queue_sup$y_sup))]

fit_JSUo = fitDist(incr_log, type = 'realline')


pHybrid = function(q, x_inf = int_inf_0,x_sup = int_sup_0,coef_inf = reg_inf[[1]],
                   coef_sup = reg_sup[[1]],fit = fit_JSUo){
  
  
  x_0_inf = - exp(-log(10)*coef_inf$Coefficients$intercept/coef_inf$Coefficients$slope)
  x_0_sup =  exp(-log(10)*coef_sup$Coefficients$intercept/coef_sup$Coefficients$slope)
  
  x_inf_delog = - exp(log(10)*x_inf)
  x_sup_delog = exp(log(10)*x_sup)
  r = rep(0,length(q))
  for(i in 1:length(q)){
    if(q[i] <= x_inf_delog){
      r[i] = (q[i]/x_0_inf)^coef_inf$Coefficients$slope
    } else if(q[i] >= x_sup_delog){
      r[i] = 1 - (q[i]/x_0_sup)^coef_sup$Coefficients$slope
    } else {
      r[i] = pSHASH(q[i], mu = fit$mu, sigma = fit$sigma, nu = fit$nu,
                    tau = fit$tau)
    }
  }
  return(r)
}


qHybrid = function(p, y_inf = y_0_inf,y_sup = y_0_sup,coef_inf = reg_inf[[1]],
                   coef_sup = reg_sup[[1]],fit = fit_JSUo){
  
  x_0_inf = - exp(-log(10)*coef_inf$Coefficients$intercept/coef_inf$Coefficients$slope)
  x_0_sup =  exp(-log(10)*coef_sup$Coefficients$intercept/coef_sup$Coefficients$slope)
  
  y_inf_delog = exp(log(10)*y_inf)
  y_sup_delog = 1-exp(log(10)*y_sup)
  r = rep(0,length(p))
  for(i in 1:length(p)){
    if(p[i] <= y_inf_delog){
      r[i] = (p[i]^(1/coef_inf$Coefficients$slope))*x_0_inf
    } else if(p[i] >= y_sup_delog){
      r[i] = x_0_sup*((1-p[i])^(1/coef_sup$Coefficients$slope))
    } else {
      r[i] = qSHASH(p[i], mu = fit$mu, sigma = fit$sigma, nu = fit$nu,
                    tau = fit$tau)
    }
  }
  return(r)
}

rHybrid = function(n, y_inf = y_0_inf,y_sup = y_0_sup,coef_inf = reg_inf[[1]],
                   coef_sup = reg_sup[[1]],fit = fit_JSUo){
  
  r_unif = runif(n)
  x_0_inf = - exp(-log(10)*coef_inf$Coefficients$intercept/coef_inf$Coefficients$slope)
  x_0_sup =  exp(-log(10)*coef_sup$Coefficients$intercept/coef_sup$Coefficients$slope)
  
  y_inf_delog = exp(log(10)*y_inf)
  y_sup_delog = 1-exp(log(10)*y_sup)
  r = rep(0,n)
  for(i in 1:n){
    if(r_unif[i] <= y_inf_delog){
      r[i] = (r_unif[i]^(1/coef_inf$Coefficients$slope))*x_0_inf
    } else if(r_unif[i] >= y_sup_delog){
      r[i] = x_0_sup*((1-r_unif[i])^(1/coef_sup$Coefficients$slope))
    } else {
      r[i] = qSHASH(r_unif[i], mu = fit$mu, sigma = fit$sigma, nu = fit$nu,
                    tau = fit$tau)
    }
  }
  return(r)
}


x_inf_delog = - exp(log(10)*int_inf_0)
x_sup_delog = exp(log(10)*int_sup_0)

p = seq(0,1,length.out = 10^4)
q = qHybrid(p)
plot(q,p,type = 'l')
hist(q,breaks = 500, xlim = c(-.1,.1))
hist(incr_log, breaks = 500, xlim = c(-.1,.1))


dens_hyb = cbind(p,q) %>% data.table %>% setnames(c("", "Support")) %>%
  ggplot(aes(Support)) + geom_density() + 
  geom_vline(xintercept = x_inf_delog, col = "darkred") + geom_vline(xintercept = x_sup_delog,
                                                                     col = "darkred") + ylab("") + xlab("") +
  xlim(c(-.05,.05))

dens_incr = cbind(incr_log) %>% data.table %>% setnames(c("Support")) %>%
  ggplot(aes(Support)) + geom_density() + 
  geom_vline(xintercept = x_inf_delog, col = "darkred") + geom_vline(xintercept = x_sup_delog,
                                                                     col = "darkred") +   ylab("") + xlab("") +
  xlim(c(-.05,.05)) 

dens_norm = cbind(p,qnorm(p,mean = mean(incr_log), sd = sd(incr_log))) %>% data.table %>% setnames(c("","Support")) %>%
  ggplot(aes(Support)) + geom_density() + 
  geom_vline(xintercept = x_inf_delog, col = "darkred") + geom_vline(xintercept = x_sup_delog,
                                                                     col = "darkred") +
  xlim(c(-.05,.05)) + ylab("") + xlab("")

figure5.3.1 = dens_incr
figure5.3.2 = dens_hyb
figure5.3.3 = dens_norm

ggsave(paste0(path_save,"Figure5_3_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.3.1)
ggsave(paste0(path_save,"Figure5_3_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.3.2)
ggsave(paste0(path_save,"Figure5_3_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure5.3.3)



# -----------------------------------------------------------------------------------------------------
# 6 - Etude et réplication des chocs de taux Solvabilité 2 --------------------------------------------

# 6.1 - La méthode retenue par l'EIOPA ----

# Figure 6.1.1&2&3 - PCA Bank of England ------------------------------

pca_england = PCA(england_EIOPA[,-1], graph = F)

figure6.1.1 = fviz_eig(pca_england, addlabels = TRUE, ylim = c(0, 50)) + labs(title = "") + ylab("Pourcentage de variance expliqu?e") + 
  xlab("Dimension") + theme_cowplot() + theme(axis.text = element_text(size = 40), axis.title = element_text(size = 35))


# ggsave(paste0(path_save,"Figure6_6_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6.1.1)

pca_england = prcomp(england_EIOPA[,-1], scale = T, center = T)

figure6.1.2 = england_EIOPA[,c(sample(2:51,50),1)] %>% melt(id.vars = "date") %>% ggplot(aes(date,value)) + geom_line(aes(col = variable)) + 
  theme_cowplot() +
  theme(legend.position = "none",axis.title = element_text(size = 40),
        axis.text = element_text(size = 40)) + ylab("Taux") + xlab("Date")

# ggsave(paste0(path_save,"Figure6_6_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6.1.2)


figure6.1.3 = autoplot(pca_england, data = england_EIOPA, colour = "date", loadings = T,
                    loadings.label = T, loadings.label.size = 4, loadings.colour = "black", frame = F, loadings.label.colour = "darkred",
         loadings.label.hjust = 2, loadings.label.vjust = 0) +
  theme_cowplot() + theme(axis.title = element_text(size = 40), legend.position = "none", axis.text = element_text(size = 40))

# ggsave(paste0(path_save,"Figure6_6_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6.1.3)

# Figure 6.1. 4&5&6 - Composante PCA BCE ----------------------------------------------------------------

pca_BCE = prcomp(yield_all_years[,-31], scale = T, center = T)

comp = pca_BCE$x[,1:3] 

figure6_1_4 = cbind.data.frame(yield_all_years$date,-comp[,1]) %>% as.data.table %>% setnames(c("index","value")) %>% 
  ggplot(aes(index,value)) + geom_line() + xlab("") + ylab("")

figure6_1_5 = cbind.data.frame(yield_all_years$date,-comp[,2]) %>% as.data.table %>% setnames(c("index","value")) %>% 
  ggplot(aes(index,value)) + geom_line() + xlab("") + ylab("")

figure6_1_6 = cbind.data.frame(yield_all_years$date,-comp[,3]) %>% as.data.table %>% setnames(c("index","value")) %>% 
  ggplot(aes(index,value)) + geom_line() + xlab("") + ylab("")

ggsave(paste0(path_save,"figure6_1_4.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_1_4)
ggsave(paste0(path_save,"figure6_1_5.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_1_5)
ggsave(paste0(path_save,"figure6_1_6.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_1_6)


# 6.2 - Application de la méthode et production de chocs ------------------------------------------------

# Figure 6.2.1 - PCA pour les taux BCE ----------------------------------------------------------
pca_BCE = prcomp(yield_all_years[,-31], scale = T, center = T)


figure6.2.1 = autoplot(pca_BCE, data = yield_all_years, colour = "date", loadings = T,
                       loadings.label = T, loadings.label.size = 4, loadings.colour = "black", frame = F, loadings.label.colour = "darkred",
                       loadings.label.hjust = 2, loadings.label.vjust = 0) +
  theme_cowplot() + theme(axis.title = element_text(size = 40), legend.position = "none", axis.text = element_text(size = 40))

# ggsave(paste0(path_save,"Figure6_2_1.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6.2.1)



# Figure 6.2.2&3 - Problématique des taux nuls ----------------------------------------------------------

data_BCE = yield_all_years[,-31]

changes_bce = (data_BCE[261:(nrow(data_BCE)),] - data_BCE[1:(nrow(data_BCE)-260),])/
  data_BCE[1:(nrow(data_BCE)-260),]

# 1 an et 10 ans

figure6_2_2 = cbind.data.frame(yield_all_years$date[261:nrow(data_BCE)],changes_bce[,1]) %>% as.data.table %>% setnames(c("Date","Variations")) %>%
  ggplot(aes(Date,Variations)) + geom_line() + 
  theme(axis.title = element_text(size = 40), legend.position = "none", axis.text = element_text(size = 45))
  

figure6_2_3 = cbind.data.frame(yield_all_years$date[261:nrow(data_BCE)],changes_bce[,10]) %>% as.data.table %>% setnames(c("Date","Variations")) %>%
  ggplot(aes(Date,Variations)) + geom_line() + 
  theme(axis.title = element_text(size = 40), legend.position = "none", axis.text = element_text(size = 45))

# ggsave(paste0(path_save,"Figure6_2_2.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_2)
# ggsave(paste0(path_save,"Figure6_2_3.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_3)
  
  





















































# Figure 6.4 - Sensibilité des maturités aux composantes --------------------------------------



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
  
  fit_10$coefficients
}

choc_each_england %>% cbind(1:25) %>% as.data.table %>% setnames(c(paste0("PC",1:4),"Maturité")) %>% melt(id.vars = "Maturité") %>%
  ggplot(aes(Maturité,value)) + geom_line(aes(col = variable))




# Figure 6.2 5&6&7&8 - Courbes EIOPA et BCE -----------------------------------------------------------


yc_BCE_31_12_2015 = yield_all_years[yield_all_years$date == as.Date("30/12/2015",format = "%d/%m/%Y"),-31] %>% as.numeric
yc_BCE_31_12_2016 = yield_all_years[yield_all_years$date == as.Date("30/12/2016",format = "%d/%m/%Y"),-31] %>% as.numeric
yc_BCE_31_12_2017 = yield_all_years[yield_all_years$date == as.Date("29/12/2017",format = "%d/%m/%Y"),-31] %>% as.numeric
yc_BCE_31_12_2018 = yield_all_years[yield_all_years$date == as.Date("28/12/2018",format = "%d/%m/%Y"),-31] %>% as.numeric

figure6_2_5 = cbind(1:30,EIOPA$spot_2015[1:30],yc_BCE_31_12_2015) %>% as.data.table %>% setnames(c("index","EIOPA","BCE")) %>%
  melt(id.vars = "index") %>% ggplot(aes(index,value)) + geom_line(aes(col = variable)) + geom_vline(xintercept = 20) + 
  theme(legend.position = "none", plot.title = element_text(size = 45), axis.text = element_text(size = 45)) + 
  xlab("") + ylab("") + labs(title = "") 

figure6_2_6 = cbind(1:30,EIOPA$spot_2016[1:30],yc_BCE_31_12_2016) %>% as.data.table %>% setnames(c("index","EIOPA","BCE")) %>%
  melt(id.vars = "index") %>% ggplot(aes(index,value)) + geom_line(aes(col = variable)) + geom_vline(xintercept = 20) + 
  theme(legend.position = "none", plot.title = element_text(size = 45), axis.text = element_text(size = 45)) + 
  xlab("") + ylab("") + labs(title = "") 

figure6_2_7 = cbind(1:30,EIOPA$spot_2017[1:30],yc_BCE_31_12_2017) %>% as.data.table %>% setnames(c("index","EIOPA","BCE")) %>%
  melt(id.vars = "index") %>% ggplot(aes(index,value)) + geom_line(aes(col = variable)) + geom_vline(xintercept = 20) + 
  theme(legend.position = "none", plot.title = element_text(size = 45), axis.text = element_text(size = 45)) + 
  xlab("") + ylab("") + labs(title = "") 

figure6_2_8 = cbind(1:30,EIOPA$spot_2018[1:30],yc_BCE_31_12_2018) %>% as.data.table %>% setnames(c("index","EIOPA","BCE")) %>%
  melt(id.vars = "index") %>% ggplot(aes(index,value)) + geom_line(aes(col = variable)) + geom_vline(xintercept = 20) + 
  theme(legend.position = "none", plot.title = element_text(size = 45), axis.text = element_text(size = 45)) + 
  xlab("") + ylab("") + labs(title = "") 

# ggsave(paste0(path_save,"figure6_2_5.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_5)
# ggsave(paste0(path_save,"figure6_2_6.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_6)
# ggsave(paste0(path_save,"figure6_2_7.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_7)
# ggsave(paste0(path_save,"figure6_2_8.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_8)


# Figure 6.2 9&10&11&12 - Courbes EIOPA et BCE -----------------------------------------------------------

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


choc_2015_original = foreach(i=1:20, .combine = rbind, .export = "purrr") %do% {
  param_DNS = calibration_DNS_original(x = yield_all_years[yield_all_years$date<=as.Date("30/12/2015",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_original(param_EIOPA = par_eiopa_2015, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}

choc_2015_copule = foreach(i=1:20, .combine = rbind, .packages = c("purrr","gamlss")) %dopar% {
  param_DNS = calibration_DNS_cor(x = yield_all_years[yield_all_years$date<=as.Date("30/12/2015",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_cor(param_EIOPA = par_eiopa_2015, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}
# choc_2015_copule[18,] = choc_2015_copule[17,]+((choc_2015_copule[20,]-choc_2015_copule[17,])/3)
# choc_2015_copule[19,] = choc_2015_copule[17,]+2*((choc_2015_copule[20,]-choc_2015_copule[17,])/3)

choc_2016_original = foreach(i=1:20, .combine = rbind, .export = "purrr") %do% {
  param_DNS = calibration_DNS_original(x = yield_all_years[yield_all_years$date<=as.Date("30/12/2016",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_original(param_EIOPA = par_eiopa_2016, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}

choc_2016_copule = foreach(i=1:20, .combine = rbind, .packages = c("purrr","gamlss")) %dopar% {
  param_DNS = calibration_DNS_cor(x = yield_all_years[yield_all_years$date<=as.Date("30/12/2016",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_cor(param_EIOPA = par_eiopa_2016, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}
# choc_2016_copule[18,] = choc_2016_copule[17,]+((choc_2016_copule[20,]-choc_2016_copule[17,])/3)
# choc_2016_copule[19,] = choc_2016_copule[17,]+2*((choc_2016_copule[20,]-choc_2016_copule[17,])/3)

choc_2017_original = foreach(i=1:20, .combine = rbind, .export = "purrr") %do% {
  param_DNS = calibration_DNS_original(x = yield_all_years[yield_all_years$date<=as.Date("29/12/2017",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_original(param_EIOPA = par_eiopa_2017, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}

choc_2017_copule = foreach(i=1:20, .combine = rbind, .packages = c("purrr","gamlss")) %dopar% {
  param_DNS = calibration_DNS_cor(x = yield_all_years[yield_all_years$date<=as.Date("29/12/2017",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_cor(param_EIOPA = par_eiopa_2017, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}
# choc_2017_copule[18,] = choc_2017_copule[17,]+((choc_2017_copule[20,]-choc_2017_copule[17,])/3)
# choc_2017_copule[19,] = choc_2017_copule[17,]+2*((choc_2017_copule[20,]-choc_2017_copule[17,])/3)

choc_2018_original = foreach(i=1:20, .combine = rbind, .export = "purrr") %do% {
  param_DNS = calibration_DNS_original(x = yield_all_years[yield_all_years$date<=as.Date("28/12/2018",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_original(param_EIOPA = par_eiopa_2018, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}

choc_2018_copule = foreach(i=1:20, .combine = rbind, .packages = c("purrr","gamlss")) %dopar% {
  param_DNS = calibration_DNS_cor(x = yield_all_years[yield_all_years$date<=as.Date("28/12/2018",format = "%d/%m/%Y"),i], taux = i)
  choc_S2_DNS_cor(param_EIOPA = par_eiopa_2018, param_DNS = param_DNS,m = 10^4,n = 260,taux = i)
}
# choc_2018_copule[18,] = choc_2018_copule[17,]+((choc_2018_copule[20,]-choc_2018_copule[17,])/3)
# choc_2018_copule[19,] = choc_2018_copule[17,]+2*((choc_2018_copule[20,]-choc_2018_copule[17,])/3)

figure6_2_9 = cbind(1:20, EIOPA$spot_2015[1:20], EIOPA$up_2015[1:20], EIOPA$down_2015[1:20], choc_2015_original[,1], choc_2015_original[,2], 
      EIOPA$spot_2016[1:20],choc_2015_copule[,1],choc_2015_copule[,2]) %>%
  as.data.table %>% setnames(c("Maturite","EIOPA_2015","EIOPA_UP","EIOPA_DOWN","DNS_orig_DOWN","DNS_orig_UP","EIOPA_2016",
                               "DNS_cor_DOWN","DNS_cor_UP")) %>%
  ggplot(aes(Maturite,EIOPA_2015)) + geom_line(size = 2) + geom_line(aes(y=EIOPA_UP, col = "Choc EIOPA"), size = 2) + 
  geom_line(aes(y = EIOPA_DOWN, col = "Choc EIOPA"), size = 2) + geom_line(aes(y = DNS_orig_DOWN, col = "Quantile DNS Original"), size = 2) + 
  geom_line(aes(y = DNS_orig_UP, col = "Quantile DNS Original"), size = 2) + geom_line(aes(y = EIOPA_2016, col = "EIOPA 2016"), size = 2) +
  geom_line(aes(y = DNS_cor_DOWN, col = "Quantile DNS Corrélation"), size = 2) + 
  geom_line(aes(y = DNS_cor_UP, col = "Quantile DNS Corrélation"), size = 2) + 
  theme(legend.position = "none", plot.title = element_text(size = 40), axis.text = element_text(size = 45), axis.title = element_text(size =40)) + 
  xlab("Maturite") + ylab("Taux") + labs(title = "") 

figure6_2_10 = cbind(1:20, EIOPA$spot_2016[1:20], EIOPA$up_2016[1:20], EIOPA$down_2016[1:20], choc_2016_original[,1], choc_2016_original[,2], 
                    EIOPA$spot_2017[1:20],choc_2016_copule[,1],choc_2016_copule[,2]) %>%
  as.data.table %>% setnames(c("Maturite","EIOPA_2015","EIOPA_UP","EIOPA_DOWN","DNS_orig_DOWN","DNS_orig_UP","EIOPA_2016",
                               "DNS_cor_DOWN","DNS_cor_UP")) %>%
  ggplot(aes(Maturite,EIOPA_2015)) + geom_line(size = 2) + geom_line(aes(y=EIOPA_UP, col = "Choc EIOPA"), size = 2) + 
  geom_line(aes(y = EIOPA_DOWN, col = "Choc EIOPA")) + geom_line(aes(y = DNS_orig_DOWN, col = "Quantile DNS Original"), size = 2) + 
  geom_line(aes(y = DNS_orig_UP, col = "Quantile DNS Original"), size = 2) + geom_line(aes(y = EIOPA_2016, col = "EIOPA 2017"), size = 2) +
  geom_line(aes(y = DNS_cor_DOWN, col = "Quantile DNS Corrélation"), size = 2) + 
  geom_line(aes(y = DNS_cor_UP, col = "Quantile DNS Corrélation"), size = 2) + 
  theme(legend.position = "none", plot.title = element_text(size = 40), axis.text = element_text(size = 45), axis.title = element_text(size =40)) + 
  xlab("Maturite") + ylab("Taux") + labs(title = "") 

figure6_2_11 = cbind(1:20, EIOPA$spot_2017[1:20], EIOPA$up_2017[1:20], EIOPA$down_2017[1:20], choc_2017_original[,1], choc_2017_original[,2], 
                     EIOPA$spot_2018[1:20],choc_2017_copule[,1],choc_2017_copule[,2]) %>%
  as.data.table %>% setnames(c("Maturite","EIOPA_2015","EIOPA_UP","EIOPA_DOWN","DNS_orig_DOWN","DNS_orig_UP","EIOPA_2016",
                               "DNS_cor_DOWN","DNS_cor_UP")) %>%
  ggplot(aes(Maturite,EIOPA_2015)) + geom_line(size = 2) + geom_line(aes(y=EIOPA_UP, col = "Choc EIOPA"), size = 2) + 
  geom_line(aes(y = EIOPA_DOWN, col = "Choc EIOPA"), size = 2) + geom_line(aes(y = DNS_orig_DOWN, col = "Quantile DNS Original"), size = 2) + 
  geom_line(aes(y = DNS_orig_UP, col = "Quantile DNS Original"), size = 2) + geom_line(aes(y = EIOPA_2016, col = "EIOPA 2018"), size = 2) +
  geom_line(aes(y = DNS_cor_DOWN, col = "Quantile DNS Corrélation"), size = 2) + 
  geom_line(aes(y = DNS_cor_UP, col = "Quantile DNS Corrélation"), size = 2) + 
  theme(legend.position = "none", plot.title = element_text(size = 40), axis.text = element_text(size = 45), axis.title = element_text(size =40)) + 
  xlab("Maturite") + ylab("Taux") + labs(title = "") 

figure6_2_12 = cbind(1:20, EIOPA$spot_2018[1:20], EIOPA$up_2018[1:20], EIOPA$down_2018[1:20], choc_2018_original[,1], choc_2018_original[,2], 
                     rep(NA,20),choc_2018_copule[,1],choc_2018_copule[,2]) %>%
  as.data.table %>% setnames(c("Maturite","EIOPA_2018","EIOPA_UP","EIOPA_DOWN","DNS_orig_DOWN","DNS_orig_UP","EIOPA_2016",
                               "DNS_cor_DOWN","DNS_cor_UP")) %>%
  ggplot(aes(Maturite,EIOPA_2018)) + geom_line(size = 2) + geom_line(aes(y=EIOPA_UP, col = "Choc EIOPA"), size = 2) + 
  geom_line(aes(y = EIOPA_DOWN, col = "Choc EIOPA"), size = 2) + geom_line(aes(y = DNS_orig_DOWN, col = "Quantile DNS Original"), size = 2) + 
  geom_line(aes(y = DNS_orig_UP, col = "Quantile DNS Original"), size = 2) + geom_line(aes(y = EIOPA_2016, col = "EIOPA 2018"), size = 2) +
  geom_line(aes(y = DNS_cor_DOWN, col = "Quantile DNS Corrélation"), size = 2) + 
  geom_line(aes(y = DNS_cor_UP, col = "Quantile DNS Corrélation"), size = 2) + 
  theme(legend.position = "none", plot.title = element_text(size = 40), axis.text = element_text(size = 45), axis.title = element_text(size =40)) + 
  xlab("Maturite") + ylab("Taux") + labs(title = "") 


# ggsave(paste0(path_save,"figure6_2_9.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_9)
# ggsave(paste0(path_save,"figure6_2_10.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_10)
# ggsave(paste0(path_save,"figure6_2_11.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_11)
# ggsave(paste0(path_save,"figure6_2_12.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_12)




# Figure 6.2 13&14&15&16 - Nouveaux chocs EIOPA -----------------------------------------------------------------------



s_up = c(0.61,0.53,0.49,0.46,0.45,0.41,0.37,0.34,0.32,0.30,0.3,0.3,0.3,0.29,0.28,0.28,0.27,0.26,0.26,0.25)
s_down = c(0.58,0.51,0.44,0.4,0.4,0.38,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.47,0.48,0.49,0.49,0.5)

b_up = c(2.14,1.86,1.72,1.61,1.58,1.44,1.3,1.19,1.12,1.05,1.05,1.05,1.05,1.02,0.98,0.98,0.95,0.91,0.91,0.88)
b_down = c(1.16,0.99,0.83,0.74,0.71,0.67,0.63,0.62,0.61,0.61,0.6,0.6,0.59,0.58,0.57,0.56,0.55,0.54,0.52,0.5)


EIOPA_new_stress = matrix(rep(0,20*8),ncol = 8) %>% as.data.table %>% 
  setnames(c("down_2015","up_2015","down_2016","up_2016","down_2017","up_2017","down_2018","up_2018"))

index = rep(seq(2,11,by =3),each = 2)
for(i in seq(1,8,by = 2)){
  EIOPA_new_stress[,i] = (lapply(1:20, function(j) ifelse(EIOPA[j,index[i]]<=0,EIOPA[j,index[i]]-b_down[j],
                                                          EIOPA[j,index[i]]*(1-s_down[j])-b_down[j])) %>% unlist) 
}

for(i in seq(2,8,by = 2)){
  EIOPA_new_stress[,i] = (lapply(1:20, function(j) EIOPA[j,index[i]]*(1+s_up[j])+b_up[j]) %>% unlist) 
}



figure6_2_13 = cbind(1:20, EIOPA$spot_2015[1:20], EIOPA_new_stress$up_2015, EIOPA_new_stress$down_2015, choc_2015_original[,1], choc_2015_original[,2], 
                    EIOPA$spot_2016[1:20],choc_2015_copule[,1],choc_2015_copule[,2]) %>%
  as.data.table %>% setnames(c("Maturite","EIOPA_2015","EIOPA_UP","EIOPA_DOWN","DNS_orig_DOWN","DNS_orig_UP","EIOPA_2016",
                               "DNS_cor_DOWN","DNS_cor_UP")) %>%
  ggplot(aes(Maturite,EIOPA_2015)) + geom_line(size = 2) + geom_line(aes(y=EIOPA_UP, col = "Choc EIOPA"), size = 2) + 
  geom_line(aes(y = EIOPA_DOWN, col = "Choc EIOPA"), size = 2) + geom_line(aes(y = DNS_orig_DOWN, col = "Quantile DNS Original"), size = 2) + 
  geom_line(aes(y = DNS_orig_UP, col = "Quantile DNS Original"), size = 2) + geom_line(aes(y = EIOPA_2016, col = "EIOPA 2016"), size = 2) +
  geom_line(aes(y = DNS_cor_DOWN, col = "Quantile DNS Corrélation"), size = 2) + 
  geom_line(aes(y = DNS_cor_UP, col = "Quantile DNS Corrélation"), size = 2) + 
  theme(legend.position = "none", plot.title = element_text(size = 40), axis.text = element_text(size = 45), axis.title = element_text(size =40)) + 
  xlab("") + ylab("") + labs(title = "") 

figure6_2_14 = cbind(1:20, EIOPA$spot_2016[1:20],EIOPA_new_stress$up_2016, EIOPA_new_stress$down_2016, choc_2016_original[,1], choc_2016_original[,2], 
                     EIOPA$spot_2017[1:20],choc_2016_copule[,1],choc_2016_copule[,2]) %>%
  as.data.table %>% setnames(c("Maturite","EIOPA_2015","EIOPA_UP","EIOPA_DOWN","DNS_orig_DOWN","DNS_orig_UP","EIOPA_2016",
                               "DNS_cor_DOWN","DNS_cor_UP")) %>%
  ggplot(aes(Maturite,EIOPA_2015)) + geom_line(size = 2) + geom_line(aes(y=EIOPA_UP, col = "Choc EIOPA"), size = 2) + 
  geom_line(aes(y = EIOPA_DOWN, col = "Choc EIOPA"), size = 2) + geom_line(aes(y = DNS_orig_DOWN, col = "Quantile DNS Original"), size = 2) + 
  geom_line(aes(y = DNS_orig_UP, col = "Quantile DNS Original"), size = 2) + geom_line(aes(y = EIOPA_2016, col = "EIOPA 2017"), size = 2) +
  geom_line(aes(y = DNS_cor_DOWN, col = "Quantile DNS Corrélation"), size = 2) + 
  geom_line(aes(y = DNS_cor_UP, col = "Quantile DNS Corrélation"), size = 2) + 
  theme(legend.position = "none", plot.title = element_text(size = 40), axis.text = element_text(size = 45), axis.title = element_text(size =40)) + 
  xlab("") + ylab("Taux") + labs(title = "") 

figure6_2_15 = cbind(1:20, EIOPA$spot_2017[1:20], EIOPA_new_stress$up_2017, EIOPA_new_stress$down_2017, choc_2017_original[,1], choc_2017_original[,2], 
                     EIOPA$spot_2018[1:20],choc_2017_copule[,1],choc_2017_copule[,2]) %>%
  as.data.table %>% setnames(c("Maturite","EIOPA_2015","EIOPA_UP","EIOPA_DOWN","DNS_orig_DOWN","DNS_orig_UP","EIOPA_2016",
                               "DNS_cor_DOWN","DNS_cor_UP")) %>%
  ggplot(aes(Maturite,EIOPA_2015)) + geom_line(size = 2) + geom_line(aes(y=EIOPA_UP, col = "Choc EIOPA"), size = 2) + 
  geom_line(aes(y = EIOPA_DOWN, col = "Choc EIOPA"), size = 2) + geom_line(aes(y = DNS_orig_DOWN, col = "Quantile DNS Original"), size = 2) + 
  geom_line(aes(y = DNS_orig_UP, col = "Quantile DNS Original"), size = 2) + geom_line(aes(y = EIOPA_2016, col = "EIOPA 2018"), size = 2) +
  geom_line(aes(y = DNS_cor_DOWN, col = "Quantile DNS Corrélation"), size = 2) + 
  geom_line(aes(y = DNS_cor_UP, col = "Quantile DNS Corrélation"), size = 2) + 
  theme(legend.position = "none", plot.title = element_text(size = 40), axis.text = element_text(size = 45), axis.title = element_text(size =40)) + 
  xlab("") + ylab("") + labs(title = "") 

figure6_2_16 = cbind(1:20, EIOPA$spot_2018[1:20], EIOPA_new_stress$up_2018, EIOPA_new_stress$down_2018, choc_2018_original[,1], choc_2018_original[,2], 
                     rep(NA,20),choc_2018_copule[,1],choc_2018_copule[,2]) %>%
  as.data.table %>% setnames(c("Maturite","EIOPA_2018","EIOPA_UP","EIOPA_DOWN","DNS_orig_DOWN","DNS_orig_UP","EIOPA_2016",
                               "DNS_cor_DOWN","DNS_cor_UP")) %>%
  ggplot(aes(Maturite,EIOPA_2018)) + geom_line(size = 2) + geom_line(aes(y=EIOPA_UP, col = "Choc EIOPA"), size = 2) + 
  geom_line(aes(y = EIOPA_DOWN, col = "Choc EIOPA"), size = 2) + geom_line(aes(y = DNS_orig_DOWN, col = "Quantile DNS Original"), size = 2) + 
  geom_line(aes(y = DNS_orig_UP, col = "Quantile DNS Original"), size = 2) + geom_line(aes(y = EIOPA_2016, col = "EIOPA 2018"), size = 2) +
  geom_line(aes(y = DNS_cor_DOWN, col = "Quantile DNS Corrélation"), size = 2) + 
  geom_line(aes(y = DNS_cor_UP, col = "Quantile DNS Corrélation"), size = 2) + 
  theme(legend.position = "none", plot.title = element_text(size = 45), axis.text = element_text(size = 45), axis.title = element_text(size =45)) + 
  xlab("") + ylab("") + labs(title = "") 


# ggsave(paste0(path_save,"figure6_2_13.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_13)
# ggsave(paste0(path_save,"figure6_2_14.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_14)
# ggsave(paste0(path_save,"figure6_2_15.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_15)
# ggsave(paste0(path_save,"figure6_2_16.pdf"), height = unit(c(15.49), "cm"), width = unit(c(25.4), "cm"), plot = figure6_2_16)

















