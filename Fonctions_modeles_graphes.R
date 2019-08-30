
# Fonctions génériques --------------------------------------------------------------------

mise_en_forme_stats = function(stats){
  
  stat_95 = map(stats,1) %>% unlist
  stat_99_5 = map(stats,2) %>% unlist
  moy = map(stats,3) %>% unlist
  var = map(stats,4) %>% unlist
  skew = map(stats,5) %>% unlist
  kurt = map(stats,6) %>% unlist
  maturite = 1:length(stat_95)
  
  res = cbind(maturite, stat_95, stat_99_5, moy, var, skew, kurt)
  
  return(res)
}


# Calibration MA(1)
MA1_estimate = function(incr){
  mu_estimate = mean(incr)
  rho_1 = (sum( (incr[1:(length(incr)-1)]-mu_estimate)*(incr[2:length(incr)]-mu_estimate)))/((length(incr)-1)*
                                                                                               var(incr))
  psi_estimate = (1-sqrt(1-4*rho_1^2) )/(2*rho_1)
  
  residus = c(incr[1] - mu_estimate,rep(0,(length(incr)-1)))
  for(i in 2:length(incr)){
    residus[i] = incr[i] - mu_estimate - psi_estimate*residus[i-1]
  }
  return(list("mu" = mu_estimate, "psi" = psi_estimate, "residus" = residus))
}

# G?n?rations de trajectoires de MA(1) avec bruit de loi JSUo
traj_MA1_JSUo = function(n, par_residus, mu_MA1, psi, y_0){

  bruit = rJSUo(n, mu = par_residus[1], sigma = par_residus[2], nu = par_residus[3], tau = par_residus[4])
  
  incr_MA = unlist(lapply(2:n, function(i) mu_MA1 + bruit[i] + psi*bruit[i-1]))
  
  traj = cumsum(c(y_0,incr_MA))
  return(traj)
}

simu_DL = function(t = 10,beta_et_res, lambda =.598){
  traj = beta_et_res[[1]] + beta_et_res[[2]] * (1-exp(-lambda*t))/(lambda*t) + beta_et_res[[3]] * 
    ( (1-exp(-lambda*t))/(lambda*t) - exp(-lambda*t)) + beta_et_res[[4]]
  return(traj)
}

# G?n?rations de trajectoires de AR(1) avec bruit de loi normale (sans moyenne)
# traj_MA1_JSUo = function(n, par_residus, phi, y_0){
#   
#   bruit = rnorm(n, mu = par_residus[1], sigma = par_residus[2])
#   
#   incr_AR = unlist(lapply(2:n, function(i)  bruit[i] + psi*bruit[i-1]))
#   
#   traj = cumsum(c(y_0,incr_MA))
#   return(traj)
# }

log_log_dt = function(incr){
  incr_pos = incr[incr>0] %>% sort ; incr_neg = - incr[incr<0] %>% sort ; emp_cdf_x = ecdf(incr)
  
  fn_prob_x = list("queue_inf" = data.table(x_inf = incr_neg %>% log10,
                                            y_inf = emp_cdf_x(-incr_neg) %>% log10),
                   "queue_sup" = data.table(x_sup = incr_pos %>% log10,
                                            y_sup = (1-emp_cdf_x(incr_pos)) %>% log10  ))
  return(fn_prob_x)
}

# Fonctions DNS ---------------------------------------------------------------------------

loadings_DNS = function(t, lambda = 0.598){
  load_b1 = rep(1,length(t))
  load_b2 = (1-exp(-lambda*t))/(lambda*t)
  load_b3 = (1-exp(-lambda*t))/(lambda*t) -exp(-lambda*t)
  return(matrix(c(load_b1,load_b2,load_b3), byrow = T, ncol = length(t)))
}

reg_BETA_per_obs_y = function(theta){
  t = 1:30 ; lambda = 0.598
  beta_1 = theta[1] ; beta_2 = theta[2] ; beta_3 = theta[3]
  
  y_th = beta_1 + beta_2*((1-exp(-lambda*t))/(lambda*t)) + beta_3*((1-exp(-lambda*t))/(lambda*t) - 
                                                                            exp(-lambda*t))
  return( sum((y_th-y)^2))
}

NS_YC = function(t,par,lambda = 0.598){
  beta_1 = par[1] ; beta_2 = par[2] ; beta_3 = par[3]
  y_th = beta_1 + beta_2*((1-exp(-lambda*t))/(lambda*t)) + beta_3*((1-exp(-lambda*t))/(lambda*t) - 
                                                                     exp(-lambda*t))
  return(y_th)
}

alea_JSUo = function(par_JSUo, indice){
  indice = indice[2:length(indice)] - indice[1:(length(indice)-1)]
  epsilon = par_JSUo$nu + par_JSUo$tau*asinh( (indice - par_JSUo$mu)/par_JSUo$tau)
  return(epsilon)
}

cor_epsilon = function(matrice_cor, n){
  L = chol(matrice_cor) ; n_var = nrow(matrice_cor)
  L_t = t(L)
  # Y = matrix(replicate(n, t(L) %*% matrix(rnorm(n_var),nrow = n_var)), ncol = n_var, byrow = T)
  Y = matrix(replicate(n, L_t %*% matrix(rnorm(n_var),nrow = n_var)), ncol = n_var, byrow = T)
  return(Y)
}

simu_JSUo = function(par, epsilon){
  mu = par$mu; sigma = par$sigma; nu = par$nu; tau = par$tau
  rJSUo = mu + sigma*sinh((epsilon-nu)/tau)
  return(rJSUo)
}

# Fonctions loi Hybride ------------------------------------------------------------------

# Recherche des param?tres
param_Hybrid = function(incr){
  
  dt_log_incr = log_log_dt(incr)
  dt_log_incr$queue_sup = dt_log_incr$queue_sup[-nrow(dt_log_incr$queue_sup),]
  # Meilleures r?gressions sur les queues de distribution
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
  
  # R?gressions et r?cup?ration des param?tres
  x_inf_seuil = R2_reg_inf$Seuil ; x_sup_seuil = R2_reg_sup$Seuil
  reg_inf = dt_log_incr$queue_inf %>% .[.$x_inf >= x_inf_seuil, ] %>% {lm(.$y_inf ~ .$x_inf)}
  reg_sup = dt_log_incr$queue_sup %>% .[.$x_sup >= x_sup_seuil, ] %>% {lm(.$y_sup ~ .$x_sup)}  
  
  slope_inf = reg_inf$coefficients[2] ; int_inf = reg_inf$coefficients[1]
  slope_sup = reg_sup$coefficients[2] ; int_sup =  reg_sup$coefficients[1]
  
  # Recherche des meilleurs seuils pour changement de loi
  x_inf_lien =  dt_log_incr$queue_inf$x_inf[dt_log_incr$queue_inf$x_inf >= x_inf_seuil][which.min(
    abs((slope_inf * (dt_log_incr$queue_inf %>% {.[.$x_inf >= x_inf_seuil, ]$x_inf}) + int_inf) -
          dt_log_incr$queue_inf$y_inf[dt_log_incr$queue_inf$x_inf >= x_inf_seuil]))]
  
  x_sup_lien =  dt_log_incr$queue_sup$x_sup[dt_log_incr$queue_sup$x_sup >= x_sup_seuil][which.min(
    abs((slope_sup * (dt_log_incr$queue_sup %>% {.[.$x_sup >= x_sup_seuil, ]$x_sup}) + int_sup) -
          dt_log_incr$queue_sup$y_sup[dt_log_incr$queue_sup$x_sup >= x_sup_seuil]))]
  
  # Proba associ?es aux meilleures seuils
  y_inf_lien =  dt_log_incr$queue_inf$y_inf[dt_log_incr$queue_inf$x_inf >= x_inf_seuil][which.min(
    abs((slope_inf * (dt_log_incr$queue_inf %>% {.[.$x_inf >= x_inf_seuil, ]$x_inf}) + int_inf) -
          dt_log_incr$queue_inf$y_inf[dt_log_incr$queue_inf$x_inf >= x_inf_seuil]))]
  
  y_sup_lien =  dt_log_incr$queue_sup$y_sup[dt_log_incr$queue_sup$x_sup >= x_sup_seuil][which.min(
    abs((slope_sup * (dt_log_incr$queue_sup %>% {.[.$x_sup >= x_sup_seuil, ]$x_sup}) + int_sup) -
          dt_log_incr$queue_sup$y_sup[dt_log_incr$queue_sup$x_sup >= x_sup_seuil]))]
  
  # Ajustement loi JSUos
  par_JSUo = gamlssML(incr, family = "JSUo") %>% {list("mu" = .$mu, "sigma" = .$sigma, "nu" = .$nu,
                                                       "tau" = .$tau)}
  
  param = list("x_inf_lien" = x_inf_lien, "x_sup_lien" = x_sup_lien, "slope_inf" = slope_inf,
               "slope_sup" = slope_sup, "int_inf" = int_inf, "int_sup" = int_sup, 
               "par_JSUo" = par_JSUo, "y_inf_lien" = y_inf_lien, "y_sup_lien" = y_sup_lien)
  return(param)
}

# Fonction pHybrid, qHybrid et rHybrid
pHybrid = function(q,x_inf_lien_fn, x_sup_lien_fn, slope_inf_fn, slope_sup_fn, int_inf_fn, 
                   int_sup_fn, par_JSUo_fn){
  
  x_0_inf = - exp(-int_inf_fn*log(10)/slope_inf_fn)
  x_0_sup = exp(-int_sup_fn*log(10)/slope_sup_fn)
  
  x_inf_lien_fn = -exp(x_inf_lien_fn*log(10))
  x_sup_lien_fn = exp(x_sup_lien_fn*log(10))
  
  p = rep(0,length(q))
  for(i in 1:length(q)){
    if(q[i] <=  x_inf_lien_fn){
      p[i] = ((q[i]/x_0_inf)^slope_inf_fn) 
    } else if(q[i] >= x_sup_lien_fn){
      p[i] = (1-(q[i]/x_0_sup)^slope_sup_fn)
    } else {
      p[i] =(pJSUo(q[i], mu=par_JSUo_fn$mu, sigma = par_JSUo_fn$sigma, nu = par_JSUo_fn$nu, tau = par_JSUo_fn$tau))
    }
  }  
  return(p)
}

qHybrid = function(p,y_inf_lien_fn, y_sup_lien_fn, slope_inf_fn, slope_sup_fn, int_inf_fn, int_sup_fn, 
                   par_JSUo_fn){
  
  x_0_inf_fn = - exp(-int_inf_fn*log(10)/slope_inf_fn)
  x_0_sup_fn = exp(-int_sup_fn*log(10)/slope_sup_fn)
  
  x_inf_lien_fn = exp(y_inf_lien_fn*log(10))
  y_sup_lien_fn = 1-exp(y_sup_lien_fn*log(10))
  
  x = rep(0,length(p))
  for(i in 1 :length(p)){
    if(p[i] <= y_inf_lien_fn){
      x[i] = x_0_inf_fn*(p[i]^(1/slope_inf_fn))
    } else if(p[i] >= y_sup_lien_fn){
      x[i] = x_0_sup_fn*((1-p[i])^(1/slope_sup_fn))
    } else{
      x[i] = qJSUo(p[i],mu=par_JSUo_fn$mu, sigma = par_JSUo_fn$sigma, nu = par_JSUo_fn$nu, tau = par_JSUo_fn$tau)
    }
  }
  return(x)
}

rHybrid = function(n,y_inf_lien_fn, y_sup_lien_fn, slope_inf_fn, slope_sup_fn, int_inf_fn, int_sup_fn, 
                   par_JSUo_fn){
  u = runif(n)
  r_hyb = rep(0,n)
  for(i in 1:n){
    r_hyb[i] = qHybrid(u[i], y_inf_lien_fn = y_inf_lien_fn, y_sup_lien_fn = y_sup_lien_fn,
                       slope_inf_fn = slope_inf_fn, slope_sup_fn = slope_sup_fn,
                       int_inf_fn = int_inf_fn, int_sup_fn = int_sup_fn, par_JSUo_fn = par_JSUo_fn)
  }
  return(r_hyb)
}

# Optimisation des séries betas issues du modèle DNS --------------------------------------

# Premi?re estimation :

loadings = loadings_DNS(t = c(2,5,10)) %>% inv
yield_estim = matrix(c(yield_all_years$`2_ans`, yield_all_years$`5_ans`, yield_all_years$`10_ans`), ncol = 3)
beta_estim = yield_estim %*% loadings 

# Estimation par optimisation (MCO) :

nbre_cores = detectCores(logical = T) - 1
registerDoParallel(cores = nbre_cores)

beta_optim = foreach(i=1:3664, .combine = rbind) %dopar% {
  y = as.numeric(yield_all_years[i,-31])
  optim(par = beta_estim[i,], fn = reg_BETA_per_obs_y)$par
}

# R?sidus 10 ans
residus_NS = foreach(i=1:3664, .combine = c) %dopar% {
  y = yield_all_years[i,10]
  y - NS_YC(t=10, par = beta_optim[i,])
}


# Fonctions graphiques ---------------------------------------------------------------------

# Couleurs ggplots
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# Régression lin?aire sur la queue de la distribution
incr = yield_all_years$`10_ans`[2:3664] - yield_all_years$`10_ans`[1:3663] 

# Diagrammes log - log : 
gg_log_multi_comp = function(x, y_list){
  incr_pos = x[x>0] ; incr_neg = x[x<0] ; emp_cdf_x = ecdf(x)
  
  fn_prob_x = list("queue_inf" = data.table(x_inf = -incr_neg %>% log10,
                                            y_inf = emp_cdf_x(incr_neg) %>% log10),
                   "queue_sup" = data.table(x_sup = incr_pos %>% log10,
                                            y_sup = (1-emp_cdf_x(incr_pos)) %>% log10  ))
  
  # fn_prob_x$queue_sup = fn_prob_x$queue_sup[- (fn_prob_x$queue_sup$y_sup %>% is.infinite %>% which),]
  
  y_list = lapply(y_list, function(x) ecdf(x)) %>% lapply(function(x) 
    list("queue_inf" = data.table(x_inf = -incr_neg %>% log10,
                                  y_inf = x(incr_neg) %>% log10),
         "queue_sup" = data.table(x_sup = incr_pos %>% log10,
                                  y_sup = (1-x(incr_pos)) %>% log10)))
  
  # Log - log inf :
  y_list_inf = mapply(function(u,v)  paste0("geom_line(aes(y =", 
                                            paste0("c(",paste0(u$queue_inf$y_inf, collapse = ","),")"),
                                            " ,col =",paste0("\"",
                                                             v,"\""),"), size = 3)"), 
                      u = y_list, v = names(y_list), SIMPLIFY = F) %>% 
    unlist %>% paste(collapse = "+")
  
  g_inf = quote(ggplot(fn_prob_x$queue_inf, aes(x_inf,y_inf)) + geom_line() + 
                  geom_point(col = "darkred", size = 3) + xlab("Log support") +
                  ylab("Log fonction répartition") + labs(title = "Diagramme log - log queue inférieure"))
  
  g_inf %<>% deparse %>% paste(collapse = " ") %>% paste("+", y_list_inf) %>% {parse(text = .)} %>% eval
  
  # Log - log sup :
  y_list_sup = mapply(function(u,v)  paste0("geom_line(aes(y =", 
                                            paste0("c(",paste0(u$queue_sup$y_sup, collapse = ","),")"),
                                            " ,col =",paste0("\"",
                                                             v,"\""),"), size = 2)"), 
                      u = y_list, v = names(y_list), SIMPLIFY = F) %>% 
    unlist %>% paste(collapse = "+")
  
  g_sup = quote(ggplot(fn_prob_x$queue_sup, aes(x_sup,y_sup)) + geom_line() + 
                  geom_point(col = "darkred", size = 3) + xlab("Log support") +
                  ylab("Log fonction survie") + labs(title = "Diagramme log - log queue supérieure"))
  
  g_sup %<>% deparse %>% paste(collapse = " ") %>% paste("+", y_list_sup) %>% {parse(text = .)} %>% eval
  
  return(list("g_inf" = g_inf, "g_sup" = g_sup))
}

# Plot des intervalles de confiance
gg_multi_IC = function(y, y_list, date, probs = c(0.025,0.975)){
  
  lower = mapply(function(u,v) unlist(lapply(1:ncol(u), function(i) quantile(u[,i], probs = probs[1]))), 
                 u = y_list, v = names(y_list)) 
  
  colnames(lower) = paste0("lower_",colnames(lower))
 
  upper = mapply(function(u,v) unlist(lapply(1:ncol(u), function(i) quantile(u[,i], probs = probs[2]))), 
                 u = y_list, v = names(y_list))
  colnames(upper) = paste0("upper_",colnames(upper))
  
  IC = cbind.data.frame(lower, upper) 
  geom_IC = mapply(function(u,v)  paste0("geom_line(aes(y =", paste0("c(",paste0(IC[,u], collapse = ","),")"),
                                                       " ,col =",paste0("\"",
                                                                        v,"\""),"), size = 3)"), 
                                 u = 1:ncol(IC), v = map(strsplit(colnames(IC),"_"),2), SIMPLIFY = F) %>% 
    unlist %>% paste(collapse = "+")

  dt = cbind.data.frame(y,date) %>% setnames(c("yield","date")) 
  
  g_IC = quote(ggplot(dt,aes(date,yield)) + geom_line() + xlab("Date") + ylab("Taux") + 
                 labs(title = "Intervalle(s) de confiance")) %>% 
    deparse %>% paste(collapse = " ") %>% paste("+", geom_IC)%>% {parse(text = .)} %>% eval
  
  res = list("g_IC" = g_IC, "IC" = IC)
  return(res)
}

# Plot fenêtre de calibration et de projection

gg_illust_fenetres = function(x_total, date_total, calib, proj, start){
  
  yield_calib = c(rep(NA,length(1:(start-1))),x_total[calib]) %>% c(rep(NA,length((length(.)+1):length(x_total))))
  date_calib =c(rep(NA,length(1:(start-1))),date_total[calib]) %>% c(rep(NA,length((length(.)+1):length(x_total))))
  
  yield_proj = c(rep(NA,length(1:(start+length(calib)-1))),x_total[proj]) %>% 
    c(rep(NA,length((length(.)+1):length(x_total))))
  date_proj = c(rep(NA,length(1:(start+length(calib)-1))),date_total[proj]) %>% 
    c(rep(NA,length((length(.)+1):length(x_total))))
  
  dt = data.table(yield = x_total, yield_calib = yield_calib, yield_proj = yield_proj ,date = date_total, 
                  date_calib = date_calib, date_proj = date_proj)
  
  g_fen = ggplot(dt, aes(date,yield)) + geom_line() + ylab("Taux") + xlab("Date") +
    geom_vline(xintercept = date_total[start], col = "darkgreen", size = 1) + 
    geom_vline(xintercept = date_total[start + length(calib)-1], col = "darkgreen", size = 1) +
    annotate("rect", xmax=date_total[start + length(calib)-1], xmin=date_total[start],
             ymin=-Inf, ymax=Inf, alpha=0.2, fill="darkgreen") +
    geom_vline(xintercept = date_total[start+ length(calib)], col = "darkred", size = 1) + 
    geom_vline(xintercept = date_total[start + length(calib) + length(proj)-1], col = "darkred", size = 1) +
    annotate("rect", xmax=date_total[start + length(calib) + length(proj)-1], 
             xmin=date_total[start+ length(calib)],
             ymin=-Inf, ymax=Inf, alpha=0.2, fill="darkred") + 
    labs(title = "Fenêtres de calibration et de projection")
}

# Plot des régressions sur les queues de distribution

gg_regression = function(yield, pas = 1000, n_reg = 50, n_elim = 10, saut = 500){
  
  index_1 = seq(1,(length(yield)-pas),by = saut)
  index_2 = (seq(1,(length(yield) - pas),by = saut)+pas)
  index = lapply(1:length(index_1), function(j) index_1[j]:index_2[j])
    
  incr_pas = lapply(index, function(i) yield[i+1] - yield[i])
  
  # Sous sections queues de distribution
  reg = lapply(incr_pas, function(y) list("queue_inf" = 
                                            data.table(x_inf = - y[y<0] %>% sort %>%
                                                         .[(length(.)-n_reg-n_elim):(length(.)-n_elim)] %>%
                                                                      log10,
                                                              y_inf = ecdf(y)(y[y<0] %>% sort(decreasing = T) %>%
                                                            .[(length(.)-n_reg-n_elim):(length(.)-n_elim)]) %>%
                                                                    log10),
                                          "queue_sup" = data.table(x_sup =  y[y>0] %>% sort %>%
                                                                .[(length(.)-n_reg-n_elim):(length(.)-n_elim)] %>%
                                                                     log10,
                                                                   y_inf = (1-ecdf(y)(y[y>0] %>% sort %>%
                                                            .[(length(.)-n_reg-n_elim):(length(.)-n_elim)]))%>%
                                                                     log10
                                                                   )))
  
  # R?gressions sur les queues
  g_reg_inf = lapply(reg, function(dt) dt$queue_inf %>% {lm(.$y_inf ~ .$x_inf)} %>% 
                       {list("intercept" = .$coef[1], "pente" = .$coef[2])} %>%
                       {paste0("geom_abline( intercept = ",.$intercept,", slope = ",.$pente,")")}) %>%
    paste0(collapse = "+")
  
  g_reg_sup = lapply(reg, function(dt) dt$queue_inf %>% {lm(.$y_inf ~ .$x_inf)} %>% 
  {list("intercept" = .$coef[1], "pente" = .$coef[2])} %>%
  {paste0("geom_abline( intercept = ",.$intercept,", slope = ",.$pente,")")}) %>%
    paste0(collapse = "+")
  
  g_reg_neg = quote(ggplot(reg[[1]]$queue_inf,aes(x_inf,y_inf)) + geom_line(col = "white") + xlab("Date") + 
                      ylab("Taux") + 
            labs(title = "Intervalle(s) de confiance")) %>% deparse %>% paste(collapse = " ") %>% 
    paste("+", g_reg_inf)%>% {parse(text = .)} %>% eval
  
  # Plot des log - log sur tout l'historique
  
  incr = yield[2:length(yield)] - yield[1:(length(yield)-1)]
  incr_pos = incr[incr>0] ; incr_neg = incr[incr<0] ; emp_cdf_x = ecdf(incr)
  
  fn_prob_x = list("queue_inf" = data.table(x_inf = -incr_neg %>% log10,
                                            y_inf = emp_cdf_x(incr_neg) %>% log10),
                   "queue_sup" = data.table(x_sup = incr_pos %>% log10,
                                            y_sup = (1-emp_cdf_x(incr_pos)) %>% log10  ))
  
  g_reg_neg = quote(ggplot(fn_prob_x$queue_inf,aes(x_inf,y_inf)) + geom_line() + xlab("Date") + 
                      ylab("Taux") + geom_point(col = "darkred", alpha = .2) + xlab("Log support") + 
                      ylab("Log fonction répartition") + xlim(-12,0) + 
    labs(title = "Régressions sur les queues de distribution inférieures")) %>% deparse %>% 
  paste(collapse = " ") %>% paste("+", g_reg_inf)%>% {parse(text = .)} %>% eval
  
  g_reg_pos = quote(ggplot(fn_prob_x$queue_sup,aes(x_sup,y_sup)) + geom_line() + xlab("Date") + 
                      ylab("Taux") + geom_point(col = "darkred", alpha = .2) +
                      xlab("Log support") + ylab("Log fonction répartition") + xlim(-12,0) +
                      labs(title = "Régressions sur les queues de distribution supérieures")) %>% deparse %>% 
    paste(collapse = " ") %>% paste("+", g_reg_sup)%>% {parse(text = .)} %>% eval
  
  
  g_queue = list("g_inf" = g_reg_neg, "g_sup" = g_reg_pos)
  return(g_queue)
}


gg_all_log_log_plot = function(start_calib = 1, m = nbre_simu, IS = F, n_calib = 1000, n_proj = 1000,
                               pas = 100, yield = yield_all_years$`10_ans`, date = yield_all_years$date,
                               modele =c("Gaussien" , "JSUo" , "Bootstrap", "Naive Gaussien", "Naive JSUo", 
                                         "DNS original", "DNS Gaussien" , "DNS JSUo",
                                         "GARCH norm", "GARCH std","GARCH sstd","DNS copule",
                                         "Empirique","Hybride")){
  
  index_1 = seq(start_calib,(length(yield)-n_proj-n_calib),by = pas)
  
  res = lapply(1:length(index_1), function(i) projection_single_modele(start_calib = index_1[i], 
                                                                       n_calib = n_calib, m = m, IS = F, 
                                                                       n_proj = n_proj, modele = modele)$incr %>%
                 as.vector)
  res = lapply(res, function(x) log_log_dt(incr = x))
  for(i in 1:length(res)){
    res[[i]]$queue_inf$x_inf = round(res[[i]]$queue_inf$x_inf,3)
    res[[i]]$queue_sup$x_sup = round(res[[i]]$queue_sup$x_sup,3)
  }
  
  sup_unique_inf = unlist(lapply(res, function(x) x$queue_inf$x_inf)) %>% unique %>% sort
  sup_unique_sup = unlist(lapply(res, function(x) x$queue_sup$x_sup)) %>% unique %>% sort
  
  mat_inf = sup_unique_inf %>% cbind(matrix(rep(NA,length(index_1)*length(sup_unique_inf)),ncol = length(index_1)))
  mat_sup = sup_unique_sup %>% cbind(matrix(rep(NA,length(index_1)*length(sup_unique_sup)),ncol = length(index_1)))
  
  for(i in 1:length(index_1)){
    res[[i]]$queue_inf = res[[i]]$queue_inf[!duplicated(res[[i]]$queue_inf$x_inf),]
    res[[i]]$queue_sup = res[[i]]$queue_sup[!duplicated(res[[i]]$queue_sup$x_sup),]
  }
  
  for(i in 2:(length(index_1)+1)){
    mat_inf[,i][mat_inf[,1] %in% res[[i-1]]$queue_inf$x_inf] =
      res[[i-1]]$queue_inf$y_inf %>% round(3)
    mat_sup[,i][mat_sup[,1] %in% res[[i-1]]$queue_sup$x_sup] =
      res[[i-1]]$queue_sup$y_sup %>% round(3)
  }
  
  g_inf = mat_inf %>% as.data.table %>% setnames(c("Sup", paste0("1_",2:ncol(mat_inf)))) %>% melt(id.vars = "Sup") %>%
    ggplot(aes(Sup,value)) + geom_line(aes(col = variable)) + theme(legend.position = "none")
  
  g_sup = mat_sup %>% as.data.table %>% setnames(c("Sup", paste0("1_",2:ncol(mat_inf)))) %>% melt(id.vars = "Sup") %>%
    ggplot(aes(Sup,value)) + geom_line(aes(col = variable)) + theme(legend.position = "none")
  
  g_inf = g_inf + theme(plot.title = element_text(size = 20)) +
    labs(title = paste0("Inf - ",modele)) + xlim(c(-15,2.5)) + ylim(c(-12.5,0))

  g_sup = g_inf  + theme(plot.title = element_text(size = 20)) +
    labs(title = paste0("Sup - ",modele)) + xlim(c(-15,2.5)) + ylim(c(-12.5,0))
  
  return(list("g_inf" = g_inf, "g_sup" = g_sup))
}

gg_qq_plot = function(x = incr, modeles = list("Gaussien" , "JSUo" , "Bootstrap", "Naive Gaussien", 
                                                         "Naive JSUo", "DNS Gaussien" , "DNS JSUo" ),
                      probs = seq(0,1, by = .01)){
  
  modeles_incr = lapply(modeles, function(x) projection_single_modele(modele = x)$incr)
  
  quant_th = modeles_incr %>% lapply(function(x) quantile(x, probs = probs)) %>% unlist %>%
    matrix(ncol = length(modeles))
  
  quantile_incr = quantile(incr, probs = probs)
  
  dt_quant = cbind(quantile_incr,quant_th) %>% as.data.table %>% setnames(c("Incr",unlist(modeles))) 
  plot = dt_quant %>% melt(id.vars = "Incr") %>% ggplot(aes(value,Incr)) + geom_line(aes(col = variable), 
                                                                                     alpha = .2) +
    geom_point(aes(col = variable)) + geom_abline(slope = 1, intercept = 0) + xlab("Quantiles théoriques") +
    ylab("Quantiles empiriques")
  
  return(list("qq_plot" = plot, "dt_qq" = dt_quant))
}




gg_qq_plot_fen = function(x = incr, modeles = list("Gaussien" , "JSUo" , "Bootstrap", "Naive Gaussien", 
                                                   "Naive JSUo", "DNS Gaussien" , "DNS JSUo" ),
                          probs = seq(0,1, by = .0005), fen = 1000, pas = 1000){
  
  index_1 = seq(1,(length(incr)-fen), by = pas) ; index_2 = index_1 + fen
  calib = lapply(1:length(index_1), function(j) index_1[j]:index_2[j])
  
  res = rep(list(rep(NA,(length(modeles)+1))),length(index_1))
  for(i in 1:length(res)){
    print(i)
    res[[i]] = gg_qq_plot(x = x[calib], modeles = modeles, probs = probs)$qq_plot
  }
  
  return(res)
  
}


moments_evol = function(incr, fen = 1000, pas = 10){
  index_1 = seq(1,(length(incr)-fen), by = pas) ; index_2 = index_1 + fen
  calib = lapply(1:length(index_1), function(j) index_1[j]:index_2[j])
  
  moy = unlist(lapply(calib, function(x) mean(incr[x])))
  sd = unlist(lapply(calib, function(x) var(incr[x])))
  skew = unlist(lapply(calib, function(x) skewness(incr[x])))
  kurt = unlist(lapply(calib, function(x) kurtosis(incr[x])))
  
  g_moy = cbind(1:length(index_1), moy) %>% as.data.table %>% setnames(c("Index","Moyenne")) %>%
    ggplot(aes(Index,Moyenne)) + geom_line() + labs(title = "Moyenne") + ylim(c(-0.005,0.005))
  g_sd = cbind(1:length(index_1), sd) %>% as.data.table %>% setnames(c("Index","Ecart_type")) %>%
    ggplot(aes(Index,Ecart_type)) + geom_line() + labs(title = "Ecart_type") + ylim(c(0,0.003))
  g_skew = cbind(1:length(index_1), skew) %>% as.data.table %>% setnames(c("Index","Skewness")) %>%
    ggplot(aes(Index,Skewness)) + geom_line() + labs(title = "Skewness") + ylim(c(-0.8,0.8))
  g_kurt = cbind(1:length(index_1), kurt) %>% as.data.table %>% setnames(c("Index","Kurtosis")) %>%
    ggplot(aes(Index,Kurtosis)) + geom_line() + labs(title = "Excès kurtosis") + ylim(c(-1,6))
  
  plots = list("g_moy" = g_moy, "g_sd" = g_sd, "g_skew" = g_skew, "g_kurt" = g_kurt)
  return(plots)
}

gg_faisceau = function(traj, probs = c(0.025,0.975), double_intervalle = F , date = yield_all_years$date, IS = T, yield = yield_all_years$`10_ans`, m_faisceau){
  
  traj = as.matrix(traj)
  m = nrow(traj)
  
  lower = unlist(lapply(1:ncol(traj), function(i) quantile(traj[,i], prob = probs[1])))
  upper = unlist(lapply(1:ncol(traj), function(i) quantile(traj[,i], prob = probs[2])))
  
  traj_sample = traj[sample(1:m,m_faisceau),] %>% as.data.table
  trial = paste0("Traj ",1:m_faisceau) ; colnames(traj_sample) = as.character(date) ; traj_sample$trial = trial
  
  traj_melt = melt(traj_sample, id.vars = "trial") 
  traj_melt$variable = as.Date(traj_melt$variable, format = "%Y-%m-%d")
  
  traj_melt$lower = rep(lower,each = m_faisceau)
  traj_melt$upper = rep(upper,each = m_faisceau)
  
  
  gg_faisceau = ggplot(traj_melt, aes(x=variable,y=value, group = trial)) + geom_line(alpha = .5, size = .2)  + 
    geom_line(aes(y=lower), color = "darkred", size = 1) +
    geom_line(aes(y = upper), color = "darkred",size = 1) +
    xlab("") + ylab("") + theme(axis.text = element_text(size = 30))
  
  if(double_intervalle & length(probs) == 4){
    lower_bis = unlist(lapply(1:ncol(traj), function(i) quantile(traj[,i], prob = probs[3]))) %>% rep(each = m_faisceau)
    upper_bis = unlist(lapply(1:ncol(traj), function(i) quantile(traj[,i], prob = probs[4]))) %>% rep(each = m_faisceau)
    
    gg_faisceau = gg_faisceau + geom_line(aes(y = lower_bis), col = "darkred") + geom_line(aes(y = upper_bis), col = "darkred")
  }
  
  if(IS){
    gg_faisceau = gg_faisceau + geom_line(aes(y=rep(yield,each = m_faisceau)), color ="darkblue",size = 1)
  }
  
  return(gg_faisceau)
}

gg_IS = function(proj_IS){
  
  g_log_inf = proj_IS$g_log_inf + theme(plot.title = element_text(size = 30), axis.text = element_text(size = 30),
                                        legend.text = element_text(size = 20), legend.title = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=5)))
  g_IC = proj_IS$g_IC + theme(plot.title = element_text(size = 30), axis.text = element_text(size = 30),
                              legend.text = element_text(size = 20), legend.title = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=5)))
  g_log_sup = proj_IS$g_log_sup + theme(plot.title = element_text(size = 30), axis.text = element_text(size =30),
                                        legend.text = element_text(size = 20), legend.title = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  res = list("g_log_inf" = g_log_inf, "g_IC" = g_IC, "g_log_sup" = g_log_sup)
  return(res)
}

gg_OoS = function(proj_OoS){
  
  # Erreur
  a_MW_erreur = proj_OoS$data_erreur %>% as.data.frame  %>%
    melt(id.vars = "index", variable.name = "Modeles")  %>%
    ggplot(aes(index,value)) + geom_line(aes(col = Modeles), size = 1.2) + xlab("Simulations") +
    ylab("Erreurs") + labs(title = "Erreurs sur 1663 simulations de 1000 points")
  a_MW_erreur = a_MW_erreur +  theme(plot.title = element_text(size = 25), axis.text = element_text(size = 20),
                               legend.text = element_text(size = 20), legend.title = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  # Ecart moyen
  a_MW_sub_ecart = proj_OoS$data_ecart_moyen %>% as.data.frame  %>%
    melt(id.vars = "index", variable.name = "Modeles")  %>%
    ggplot(aes(index,value)) + geom_line(aes(col = Modeles), size = 1) + xlab("Simulations") +
    ylab("Ecart moyen (%)") + labs(title = "Ecarts moyens de sortie d'intervalle")
  a_MW_sub_ecart = a_MW_sub_ecart + theme(plot.title = element_text(size = 25), axis.text = element_text(size = 20),
                                          legend.text = element_text(size = 20), legend.title = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  # IC terminauxs
  a_MW_sub_IC = proj_OoS$data_largeur %>% as.data.frame 
  
  if("GARCH_min" %in% colnames(a_MW_sub_IC)){
    a_MW_sub_IC[a_MW_sub_IC$GARCH_min< -7,5] = NA
    a_MW_sub_IC[a_MW_sub_IC$GARCH_max> 10,6] = NA
  }
  
  names_largeur = colnames(a_MW_sub_IC)[-length(colnames(a_MW_sub_IC))]
  
  g_largeur_geom = mapply(function(u,v)  paste0("geom_line(aes(y =", paste0("c(",paste0(a_MW_sub_IC[,u], 
                                                                                        collapse = ","),")"),
                                                " ,col =",paste0("\"",
                                                                 v,"\""),"))"), 
                          u = 1:(ncol(a_MW_sub_IC)-1), v = map(strsplit(names_largeur,split = "_"),1), SIMPLIFY = F) %>% 
    unlist %>% paste(collapse = "+")
  
  g_largeur = quote(ggplot(a_MW_sub_IC,aes(1:nrow(a_MW_sub_IC),yield_last)) + geom_line() + xlab("Simulations") + 
                      ylab("Taux (%)") + labs(title = "Largeur intervalles confiance fin de simulation")) %>% 
    deparse %>% paste(collapse = " ") %>% paste("+", g_largeur_geom)%>% {parse(text = .)} %>% eval
  
  g_largeur = g_largeur + theme(plot.title = element_text(size = 25), axis.text = element_text(size = 20),
                                legend.text = element_text(size = 20), legend.title = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  res = list("gg_erreur" = a_MW_erreur, "gg_ecart" = a_MW_sub_ecart, "gg_IC" = g_largeur)
  return(res)
}

# Projections globales -------------------------------------------------------------------------------

proj_densite = function(x = yield, m = nbre_simu, IS = T, n = NULL, densfun = c("norm","JSUo")){
 
  n = ifelse(IS,(length(x)-1),n)
  x_incr = x[2:length(x)] - x[1:(length(x)-1)]
  
  if(densfun == "norm"){
    mean = mean(x_incr) ; sd = sd(x_incr)
    incr_th = matrix(replicate(m,rnorm(n, mean = mean, sd = sd)), byrow = T, nrow = m)
  } else if(densfun == "JSUo"){
    fit_JSUo = gamlssML(x_incr, family = "JSUo") %>% {list("mu" = .$mu, "sigma" = .$sigma, "tau" = .$tau,
                                                           "nu" = .$nu)}
    incr_th = matrix(replicate(m,rJSUo(n, mu = fit_JSUo$mu, sigma = fit_JSUo$sigma, tau = fit_JSUo$tau,
                                         nu = fit_JSUo$nu)), byrow = T, nrow = m)}
  
  if(IS == T) {
    yield_th = cbind(rep(x[1],m),incr_th) %>% apply(1,cumsum) %>% t
  } else{
    yield_th = cbind(rep(x[length(x)],m),incr_th) %>% apply(1,cumsum) %>% t}
  
  return(list("incr" = as.vector(incr_th), "yield" = yield_th))
}


proj_bootstrap = function(x = yield, m = nbre_simu, IS = T, n = NULL){
  
  n = ifelse(IS,(length(x)-1),n)
  x_incr = x[2:length(x)] - x[1:(length(x)-1)]
  
  incr_th = matrix(replicate(m,sample(x_incr, size = n, replace = T)) ,byrow = T, nrow = m)

  if(IS == T) {
    yield_th= cbind(rep(x[1],m),incr_th) %>% apply(1,cumsum) %>% t}
  
  else{
    yield_th = cbind(rep(x[length(x)],m),incr_th) %>% apply(1,cumsum) %>% t}
  
  return(list("incr" = as.vector(incr_th), "yield" = yield_th))
}


proj_naive = function(x = yield, m = nbre_simu, IS = T, n = NULL, densfun = "norm"){
  
  n = ifelse(IS,(length(x)-1),n)
  x_incr = x[2:length(x)] - x[1:(length(x)-1)]
  
  par_x = MA1_estimate(x_incr)
  
  if(densfun == "norm"){
    mean = mean(par_x$residus) ; sd = sd(par_x$residus)
    res_th = matrix(replicate(m,rnorm(n+1, mean = mean, sd = sd)), byrow = T, nrow = m)
  } else if(densfun == "JSUo"){
    fit_JSUo = gamlssML(par_x$residus, family = "JSUo")
    res_th = matrix(replicate(m,rJSUo(n+1, mu = fit_JSUo$mu, sigma = fit_JSUo$sigma, tau = fit_JSUo$tau,
                                       nu = fit_JSUo$nu)), byrow = T, nrow = m)
    }
  
  incr_th = apply(res_th, 1, function(y) par_x$mu + y[2:length(y)] + par_x$psi*y[1:(length(y)-1)]) %>% t
  
  if(IS == T) {
    yield_th = cbind(rep(x[1],m),incr_th) %>% apply(1,cumsum) %>% t}
  
  else{
    yield_th = cbind(rep(x[length(x)],m),incr_th) %>% apply(1,cumsum) %>% t}
  
  
  return(list("incr" = as.vector(incr_th), "yield" = yield_th))
}


proj_DNS_remanie = function(x = yield, start_calib = NULL, m = nbre_simu, IS = T, n = NULL, 
                            densfun = c("norm","JSUo")){
  # Les betas otpimis?s sont beta_optim
  # Les résidus du modèle NS sont residus_NS
  
  # On récupère ici seulement les loadings pour le taux 10 ans
  loadings = as.numeric(loadings_DNS(t = 10))
  
  # Nombre de points de projection et fenêtre de calibration
  n = ifelse(IS,(length(x)-1),n)
  if(IS) calib = (1:length(x)) else calib = (start_calib:(start_calib+length(x)-1))
  
  # On r?cup?re les betas incr?ments sur la p?riode de calibration ainsi que les r?sidus de NS
  x_incr = x[2:length(x)] - x[1:(length(x)-1)]
  beta = beta_optim[calib,] ; beta_incr = beta[2:nrow(beta),] - beta[1:(nrow(beta)-1),]
  residus_10_ans = residus_NS[calib]
  
  #MA Betas
  MA_beta = lapply(1:3, function(i) MA1_estimate(beta_incr[,i])) ; names(MA_beta) = c("B1","B2","B3")
  
  # Ici on projette les r?sidus des mod?les MA(1) sur les b?tas et les r?sidus de NS
  if(densfun == "norm"){
    par = lapply(MA_beta, function(x) list( 'mean' = mean(x$residus), "sd" = sd(x$residus)))
    res_beta_th = lapply(par, function(x) matrix(replicate(m,rnorm(n+1, mean = x$mean, sd = x$sd)), byrow = T, 
                                                 nrow = m))
    fit_residus = list("mean" = mean(residus_10_ans), "sd" = sd(residus_10_ans))
    residus_th = matrix(replicate(m,rnorm(n, mean = fit_residus$mean, sd = fit_residus$sd)),byrow = T, nrow = m) 
    
  } else if(densfun == "JSUo"){
    par = lapply(MA_beta, function(x) gamlssML(x$residus, family = "JSUo"))
    res_beta_th = lapply(par, function(x) matrix(replicate(m,rJSUo(n+1, mu = x$mu, sigma = x$sigma, 
                                                              tau = x$tau, nu = x$nu)), 
                                            byrow = T, nrow = m))
    fit_residus = gamlssML(residus_NS,family = "JSUo")
    residus_th = matrix(replicate(m,rJSUo(n, mu = fit_residus$mu, sigma = fit_residus$sigma,
                                          tau = fit_residus$tau, nu = fit_residus$nu)),byrow = T, nrow = m) 
  }
  
  # On calcule les projections des incr?ments selon des mod?les MA(1) et on rajoute une colonne de r?sidus nuls
  # vu qu'on part du dernier point de la p?riode de calibration
  residus_th = cbind(rep(0,m),residus_th)
  incr_beta_th = lapply(1:3, function(i) apply(res_beta_th[[i]], 1 , function(y) MA_beta[[i]]$mu + 
                                                 y[2:length(y)] + MA_beta[[i]]$psi*y[1:(length(y)-1)]) %>% t)
  
  # Ici on calcule la projection des betas
  if(IS == T){
    beta_th = lapply(1:3, function(i) cbind(rep(beta[1,i],m),incr_beta_th[[i]]) %>% apply(1,cumsum) %>% t)
  } else{
    beta_th = lapply(1:3, function(i) cbind(rep(beta[nrow(beta),i],m),incr_beta_th[[i]]) %>% apply(1,cumsum) %>% t)
  }
  
  # Yield projet?
  yield_th = beta_th[[1]]*loadings[1] + beta_th[[2]]*loadings[2] + beta_th[[3]]*loadings[3] + residus_th
  incr_th = yield_th[,2:ncol(yield_th)] - yield_th[,1:(ncol(yield_th)-1)]
  
  return(list("incr" = as.vector(incr_th), "yield" = yield_th))
}

proj_DNS = function(x = yield, start_calib = NULL, m = nbre_simu, IS = T, n = NULL){
  # Les betas otpimis?s sont beta_optim
  # Les r?sidus du mod?le NS sont residus_NS
  
  # On r?cup?re ici seulement les loadings pour le taux 10 ans
  loadings = as.numeric(loadings_DNS(t = 10))
  
  # Nombre de points de projection et fen?tre de calibration
  n = ifelse(IS,(length(x)-1),n)
  if(IS) calib = (1:length(x)) else calib = (start_calib:(start_calib+length(x)-1))
  
  # On r?cup?re les betas sur la p?riode de calibration ainsi que les r?sidus de NS
  beta = beta_optim[calib,] 
  residus_10_ans = residus_NS[calib]
  
  #AR Betas
  AR_beta = lapply(1:3, function(i) lm(beta[2:nrow(beta),i]~beta[1:(nrow(beta)-1),i]) %>%
                     {list("intercept" = .$coef[1], "phi" = .$coef[2], "sd" = sd(.$residuals))}) 
  
  # Ici on projette les r?sidus des mod?les MA(1) sur les b?tas et les r?sidus de NS
  
  res_beta_th = lapply(AR_beta, function(x) matrix(replicate(m,rnorm(n, mean = 0, sd = x$sd)), byrow = T, 
                                                 nrow = m))
  fit_residus = list("mean" = mean(residus_10_ans), "sd" = sd(residus_10_ans))
  residus_th = matrix(replicate(m,rnorm(n, mean = fit_residus$mean, sd = fit_residus$sd)),byrow = T, nrow = m) 
    
  
  # On calcule les projections selon des mod?les AR(1) et on rajoute une colonne de r?sidus nuls
  # vu qu'on part du dernier point de la p?riode de calibration
  residus_th = cbind(rep(0,m),residus_th)
  beta_th = list("B1" = matrix(rep(0,m*(n+1)), nrow = m), "B2" = matrix(rep(0,m*(n+1)), nrow = m),
                 "B3" = matrix(rep(0,m*(n+1)), nrow = m))
  
  if(IS == T){
    for(i in 1:3){
      beta_th[[i]][,1] = rep(beta[1,i],m)
    }
  } else{
    for(i in 1:3){
      beta_th[[i]][,1] = rep(beta[nrow(beta),i],m)
    }
  }
  
  for(i in 1:3){
    for(j in 2:(n+1)){
      beta_th[[i]][,j] = AR_beta[[i]]$intercept + AR_beta[[i]]$phi*beta_th[[i]][,j-1] + res_beta_th[[i]][,j-1]
    }
  }
  
  # Ici on calcule la projection des betas
  # Yield projet?
  yield_th = beta_th[[1]]*loadings[1] + beta_th[[2]]*loadings[2] + beta_th[[3]]*loadings[3] + residus_th
  incr_th = yield_th[,2:ncol(yield_th)] - yield_th[,1:(ncol(yield_th)-1)]
  
  return(list("incr" = as.vector(incr_th), "yield" = yield_th))
}

proj_GARCH = function(x = yield, m = nbre_simu, IS = T, n = NULL, densfun = c("norm","std","sstd")){
  
  n = ifelse(IS,(length(x)-1),n)
  x_incr = x[2:length(x)] - x[1:(length(x)-1)]
  
  fit_garch = garchFit(formula = ~ garch(1,1), data = x_incr, cond.dist = "sstd",trace = F)
  alpha =  fit_garch@fit$par[3] ; beta = fit_garch@fit$par[4]
 
   # N?cessite alpha + beta < 1
  if(alpha + beta >= 1){
    beta = beta - ((alpha+beta)-1 + .000001) 
  }
  
  # alpha + beta
  # beta - ((alpha+beta)-1 + .000001) + alpha
  
  garch_spec = garchSpec(model = list(mu = fit_garch@fit$par[1],
                                      omega = fit_garch@fit$par[2],
                                      alpha = alpha,
                                      beta = beta,
                                      skew = fit_garch@fit$par[5],
                                      shape = fit_garch@fit$par[6]), cond.dist = densfun)
  
  
  garch_sim  = matrix(replicate(m,garchSim(spec = garch_spec, n = n)),byrow = T, nrow = m)
  
  if(IS == T) {
    yield_th = cbind(rep(x[1],m),garch_sim) %>% apply(1,cumsum) %>% t
  } else{
    yield_th = cbind(rep(x[length(x)],m),garch_sim) %>% apply(1,cumsum) %>% t}
  
  return(list("incr" = as.vector(garch_sim), "yield" = yield_th))
}

proj_DNS_copule = function(x = yield, start_calib = NULL, m = nbre_simu, IS = T, n = NULL){
  # Les betas otpimis?s sont beta_optim
  # Les r?sidus du mod?le NS sont residus_NS
  
  # On r?cup?re ici seulement les loadings pour le taux 10 ans
  loadings = as.numeric(loadings_DNS(t = 10))
  
  # Nombre de points de projection et fen?tre de calibration
  n = ifelse(IS,(length(x)-1),n)
  if(IS) calib = (1:length(x)) else calib = (start_calib:(start_calib+length(x)-1))
  
  # On r?cup?re les betas incr?ments sur la p?riode de calibration ainsi que les r?sidus de NS
  beta = beta_optim[calib,] ; beta_incr = beta[2:nrow(beta),] - beta[1:(nrow(beta)-1),]
  residus_10_ans = residus_NS[calib]
  
  # Ajustement de la copule et des marginales Beta incr?ments 
  beta_incr_pobs = beta_incr %>% {cbind(pobs(.[,1]),pobs(.[,2]),pobs(.[,3]))}
  cop_t = tCopula(dim = 3)
  fit_cop_t = fitCopula(cop_t,beta_incr_pobs, method  = "mpl")
  cop_fittee = tCopula(param = fit_cop_t@estimate[1], dim = 3, df = fit_cop_t@estimate[2])
  
  B1_incr_JSUo = gamlssML(beta_incr[,1], family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                    "tau" = .$tau,"nu" = .$nu)}
  B2_incr_JSUo = gamlssML(beta_incr[,2], family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                    "tau" = .$tau,"nu" = .$nu)}
  B3_incr_JSUo = gamlssML(beta_incr[,3], family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                    "tau" = .$tau,"nu" = .$nu)}
  
  mvdc_object = mvdc(cop_fittee, margins = c("JSUo","JSUo","JSUo"), paramMargins = list(B1_incr_JSUo,
                                                                                        B2_incr_JSUo,
                                                                                        B3_incr_JSUo))
  rBeta_incr = lapply(1:m, function(i) rMvdc(n,mvdc_object))
  
  
  
  # Ajustement des r?siuds NS
  res_NS_incr = residus_10_ans[2:length(residus_10_ans)] - residus_10_ans[1:(length(residus_10_ans)-1)]
  fit_res_incr = gamlssML(res_NS_incr,family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                 "tau" = .$tau,"nu" = .$nu)}
  
  if(IS == T){
    beta_th = lapply(1:m, function(i) rMvdc(n,mvdc_object) %>% {rbind(beta[1,],.)} %>% apply(2,cumsum))
    res_th_NS = lapply(1:m, function(i) rJSUo(n, mu = fit_res_incr$mu, sigma = fit_res_incr$sigma, 
                                              tau = fit_res_incr$tau,
                                              nu = fit_res_incr$nu) %>% {c(residus_10_ans[1],.)} %>%
                         cumsum)
    
  } else{
    beta_th = lapply(1:m, function(i) rMvdc(n,mvdc_object) %>% {rbind(beta[nrow(beta),],.)} %>% apply(2,cumsum))
    res_th_NS = lapply(1:m, function(i) rJSUo(n, mu = fit_res_incr$mu, sigma = fit_res_incr$sigma, 
                                              tau = fit_res_incr$tau,
                                              nu = fit_res_incr$nu) %>% 
                                              {c(residus_10_ans[length(residus_10_ans)],.)} %>%
                         cumsum)
  }
  
  yield_th = matrix(unlist(lapply(1:m, function(i) loadings[1]*beta_th[[i]][,1] + loadings[2]*beta_th[[i]][,2] +
                                    loadings[3]*beta_th[[i]][,3] + res_th_NS[[i]])), byrow = F, ncol = m) %>% t
  incr_th = yield_th[,2:ncol(yield_th)] - yield_th[,1:(ncol(yield_th)-1)] 
  
  return(list("incr" = as.vector(incr_th), "yield" = yield_th))
}


proj_DNS_copule_cor = function(x = yield, start_calib = NULL, m = nbre_simu, IS = T, n = NULL){
  # Les betas otpimis?s sont beta_optim
  # Les r?sidus du mod?le NS sont residus_NS
  
  # On r?cup?re ici seulement les loadings pour le taux 10 ans
  loadings = as.numeric(loadings_DNS(t = 10))
  
  # Nombre de points de projection et fen?tre de calibration
  n = ifelse(IS,(length(x)-1),n)
  if(IS) calib = (1:length(x)) else calib = (start_calib:(start_calib+length(x)-1))
  
  # On r?cup?re les betas incr?ments sur la p?riode de calibration ainsi que les r?sidus de NS
  beta = beta_optim[calib,] ; beta_incr = beta[2:nrow(beta),] - beta[1:(nrow(beta)-1),]
  residus_10_ans = residus_NS[calib] 
  resi_incr = residus_10_ans[2:length(residus_10_ans)] - residus_10_ans[1:(length(residus_10_ans)-1)]
  
  # Ajustement de la copule et des marginales Beta incr?ments 
  
  B1_incr_JSUo = gamlssML(beta_incr[,1], family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                    "tau" = .$tau,"nu" = .$nu)}
  B2_incr_JSUo = gamlssML(beta_incr[,2], family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                    "tau" = .$tau,"nu" = .$nu)}
  B3_incr_JSUo = gamlssML(beta_incr[,3], family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                    "tau" = .$tau,"nu" = .$nu)}
  fit_resi_incr = gamlssML(resi_incr,family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                 "tau" = .$tau,"nu" = .$nu)}
  
  
  a_B1 = alea_JSUo(par_JSUo = B1_incr_JSUo, indice = beta_incr[,1])
  a_B2 = alea_JSUo(par_JSUo = B2_incr_JSUo, indice = beta_incr[,2])
  a_B3 = alea_JSUo(par_JSUo = B3_incr_JSUo, indice = beta_incr[,3])
  a_R = alea_JSUo(par_JSUo = fit_resi_incr, indice = resi_incr)
  
  cor_mat = cbind(a_B1,a_B2,a_B3,a_R) %>% cor
  
  loadings = loadings_DNS(t = 10)
  
  # IS ou non IS
  if(IS == T){
    res_th = lapply(1:m, function(i) cor_epsilon(matrice_cor = cor_mat, n=n) %>%
    {list("B1" = c(beta[1,1],simu_JSUo(par = B1_incr_JSUo, epsilon = .[,1])) %>% cumsum,
          "B2" = c(beta[1,2],simu_JSUo(par = B2_incr_JSUo, epsilon = .[,2])) %>% cumsum,
          "B3" = c(beta[1,3],simu_JSUo(par = B3_incr_JSUo, epsilon = .[,3])) %>% cumsum,
          "R" = c(residus_NS[1],simu_JSUo(par = fit_resi_incr, epsilon = .[,4])) %>% cumsum)})
    
  } else{
    res_th = lapply(1:m, function(i) cor_epsilon(matrice_cor = cor_mat, n=n) %>%
    {list("B1" = c(last(beta[,1]),simu_JSUo(par = B1_incr_JSUo, epsilon = .[,1])) %>% cumsum,
          "B2" = c(last(beta[,2]),simu_JSUo(par = B2_incr_JSUo, epsilon = .[,2])) %>% cumsum,
          "B3" = c(last(beta[,3]),simu_JSUo(par = B3_incr_JSUo, epsilon = .[,3])) %>% cumsum,
          "R" = c(last(residus_NS),simu_JSUo(par = fit_resi_incr, epsilon = .[,4])) %>% cumsum)})
  }
  
  yield_th = matrix(unlist(lapply(1:m, function(i) res_th[[i]][[1]] + loadings[2]*res_th[[i]][[2]] + 
                                 loadings[3]*res_th[[i]][[3]] + res_th[[i]][[4]])), byrow = T, nrow = m)
  incr_th = yield_th[,2:ncol(yield_th)] - yield_th[,1:(ncol(yield_th)-1)] 
  
  return(list("incr" = as.vector(incr_th), "yield" = yield_th))
}


proj_Hybrid = function(x = yield, m = nbre_simu, IS = T, n = NULL){
  
  n = ifelse(IS,(length(x)-1),n)
  x_incr = x[2:length(x)] - x[1:(length(x)-1)]
  
  param_hybrid = param_Hybrid(incr = x_incr)

  incr_th = matrix(replicate(m,
                             rHybrid(n, y_inf_lien_fn = param_hybrid$y_inf_lien, y_sup_lien_fn = param_hybrid$y_sup_lien,
                     slope_inf_fn = param_hybrid$slope_inf, slope_sup_fn = param_hybrid$slope_sup,
                     int_inf_fn = param_hybrid$int_inf, int_sup_fn = param_hybrid$int_sup,
                     par_JSUo_fn = param_hybrid$par_JSUo)), 
                   byrow = T, nrow = m)
  
  if(IS == T) {
    yield_th = cbind(rep(x[1],m),incr_th) %>% apply(1,cumsum) %>% t
  } else{
    yield_th = cbind(rep(x[length(x)],m),incr_th) %>% apply(1,cumsum) %>% t}
  
  return(list("incr" = as.vector(incr_th), "yield" = yield_th))
}


# Fonctions génériques de projection -----------------------------------------------------------------------

# x = yield_all_years$`10_ans` ; date = yield_all_years$date; start_calib = 1000
# n_calib = 2000; m = nbre_simu; IS = F; n_proj = 600; trace = T

projection_modeles = function(x = yield_all_years$`10_ans` , date = yield_all_years$date, start_calib = NULL, 
                              n_calib = NULL, m = nbre_simu, IS = T, n_proj = NULL, trace = T,
                              modeles = c("Gaussien","JSUo","Bootstrap","naive gaussien",
                              "naive JSUo", "DNS original", "DNS gaussien", "DNS JSUo",
                                          "GARCH sstd", "DNS copule","Hybride")){
  
  if(IS) periode_calib = (1:length(x)) else periode_calib = (start_calib:(start_calib+n_calib-1))
  if(IS) periode_proj = (1:length(x)) else periode_proj = ((start_calib+n_calib):(start_calib+n_calib+n_proj-1))
  
  x_calib = x[periode_calib]; x_proj = x[periode_proj]
  date_calib = date[periode_calib]; date_proj = date[periode_proj]
  
  y_list = list()
  
  print(paste0("IS = ",IS))
  if(IS){
    if("Gaussien" %in% modeles){
      norm = proj_densite(densfun = "norm", x = x_calib, m = m)
      ifelse(trace,print("Mod?le Gaussien"),NULL)
      y_list[["Gaussien"]] = norm
      
    }  
    
    if("JSUo" %in% modeles){ 
      JSUo = proj_densite(densfun = "JSUo", x = x_calib, m = m)
      ifelse(trace,print("Mod?le JSUo"),NULL)
      y_list[["JSUo"]] = JSUo
    
    }  
    
    if("Bootstrap" %in% modeles){
      boot = proj_bootstrap(x = x_calib, m = m)
      ifelse(trace,print("Mod?le Bootstrap"),NULL)
      y_list[["Bootstrap"]] = boot
    
    }  
    
    if("naive gaussien" %in% modeles){
      naive_gaussien = proj_naive(densfun = "norm", x = x_calib, m = m)
      ifelse(trace,print("Naive Gaussien"),NULL)
      y_list[["naive gaussien"]] = naive_gaussien
    
    }  
    
    if("naive JSUo" %in% modeles){
      naive_JSUo = proj_naive(densfun = "JSUo", x = x_calib, m = m)
      ifelse(trace,print("Naive JSUo"),NULL)
      y_list[["naive JSUo"]] = naive_JSUo
    
    }  
    
    if("DNS original" %in% modeles){
      DNS_original = proj_DNS(x = x_calib, m=m)
      ifelse(trace,print("DNS original"),NULL)
      y_list[["DNS original"]] = DNS_original
    
    }  
    
    if("DNS gaussien" %in% modeles){
      DNS_gaussien = proj_DNS_remanie(densfun = "norm", x = x_calib, m = m)
      ifelse(trace,print("DNS Gaussien"),NULL)
      y_list[["DNS gaussien"]] = DNS_gaussien
    
    }  
    
    if("DNS JSUo" %in% modeles){
      DNS_JSUo = proj_DNS_remanie(densfun = "JSUo", x = x_calib, m = m)
      ifelse(trace,print("DNS JSUo"),NULL)
      y_list[["DNS JSUo"]] = DNS_JSUo
    
    }  
    
    if("GARCH sstd" %in% modeles){
      GARCH = proj_GARCH(densfun = "sstd", x = x_calib, m = m)
      ifelse(trace,print("GARCH"),NULL)
      y_list[["GARCH sstd"]] = GARCH
    
    }  
    
    if("DNS copule" %in% modeles){
      DNS_copule = proj_DNS_copule_cor(x = x_calib, m=m)
      ifelse(trace,print("DNS Copule"),NULL)
      y_list[["DNS copule"]] = DNS_copule
    
    }  
    
    if("Hybride" %in% modeles){
      loi_hybride = proj_Hybrid(x = x_calib, m=m)
      ifelse(trace,print("Loi hybride"),NULL)
      y_list[["Hybride"]] = loi_hybride
    }
    
  } else{
    if("Gaussien" %in% modeles){
      norm = proj_densite(x = x_calib, IS = F, n = n_proj, densfun = "norm", m = m)
      ifelse(trace,print("Mod?le Gaussien"),NULL)
      y_list[["Gaussien"]] = norm
    
    } 
    
    if("JSUo" %in% modeles){
      JSUo = proj_densite(x = x_calib, IS = F, n = n_proj, densfun = "JSUo", m = m)
      ifelse(trace,print("Mod?le JSUo"),NULL)
      y_list[["JSUo"]] = JSUo
    
    } 
    
    if("Bootstrap" %in% modeles){
      boot = proj_bootstrap(x = x_calib, IS = F, n = n_proj, m = m)
      ifelse(trace,print("Mod?le Bootstrap"),NULL)
      y_list[["Bootstrap"]] = boot
    
    } 
    
    if("naive gaussien" %in% modeles){
      naive_gaussien = proj_naive(densfun = "norm", x = x_calib, IS = F, n = n_proj, m = m)
      ifelse(trace,print("Naive Gaussien"),NULL)
      y_list[["naive gaussien"]] = naive_gaussien
    
    } 
    
    if("naive JSUo" %in% modeles){
      naive_JSUo = proj_naive(densfun = "JSUo", x = x_calib, IS = F, n = n_proj, m = m)
      ifelse(trace,print("Naive JSUo"),NULL)
      y_list[["naive JSUo"]] = naive_JSUo
    
    } 
    
    if("DNS original" %in% modeles){
      DNS_original = proj_DNS(x = x_calib, IS = F, n = n_proj,start_calib = start_calib, m = m)
      ifelse(trace,print("DNS original"),NULL)
      y_list[["DNS original"]] = DNS_original
    
    } 
    
    if("DNS Gaussien" %in% modeles){
      DNS_gaussien = proj_DNS_remanie(densfun = "norm", x = x_calib, IS = F, n = n_proj,
                                       start_calib = start_calib, m = m)
      ifelse(trace,print("DNS Gaussien"),NULL)
      y_list[["DNS Gaussien"]] = DNS_gaussien
    
    } 
    
    if("DNS JSUo" %in% modeles){
      DNS_JSUo = proj_DNS_remanie(densfun = "JSUo", x = x_calib, IS = F, n = n_proj,
                                   start_calib = start_calib, m = m)
      ifelse(trace,print("DNS JSUo"),NULL)
      y_list[["DNS JSUo"]] = DNS_JSUo
    
    } 
    
    if("GARCH sstd" %in% modeles){
      GARCH = proj_GARCH(x = x_calib, IS = F, n = n_proj, densfun = "sstd", m = m)
      ifelse(trace,print("GARCH"),NULL)
      y_list[["GARCH sstd"]] = GARCH
    
    } 
    
    if("DNS copule" %in% modeles){
      DNS_copule = proj_DNS_copule_cor(x = x_calib, IS = F, n = n_proj, start_calib = start_calib, m = m)
      ifelse(trace,print("DNS Copule"),NULL)
      y_list[["DNS copule"]] = DNS_copule
    
    } 
    
    if("Hybride" %in% modeles){
      loi_hybride = proj_Hybrid(x = x_calib, IS = F, n = n_proj, m = m)
      ifelse(trace,print("Loi hybride"),NULL)
      y_list[["Hybride"]] = loi_hybride
    }
  }
  
  y_l_yield = lapply(y_list, function(x) x$yield)
  
  y_l_incr = lapply(y_list, function(x) x$incr)
  
  if(IS) g_fen = NULL else{
    g_fen = gg_illust_fenetres(x_total = x, date_total = date, calib = periode_calib, proj = periode_proj,
                               start = start_calib)
  }
  
  print("Intervalles de confiance")
  if(IS) y = x_proj else y = c(last(x_calib),x_proj)
  if(IS) date = date_proj else date = c(last(date_calib),date_proj)
  g_IC = gg_multi_IC(y = y, y_list = y_l_yield, date = date)$g_IC
  
  print("Diagrammes Log - Log")
  g_log = gg_log_multi_comp(x = x_proj[2:length(x_proj)] - x_proj[1:(length(x_proj)-1)], y_list = y_l_incr)
  
  
  
  plots = list("g_IC" = g_IC, "g_log_inf" = g_log$g_inf, "g_log_sup" = g_log$g_sup, "g_fen" = g_fen)
  # plots = list("g_IC" = g_IC, "g_log_inf" = g_log$g_inf, "g_log_sup" = g_log$g_sup)
  return(plots)
}  


projection_single_modele = function(x = yield_all_years$`10_ans` , date = yield_all_years$date, start_calib = NULL, 
                                    n_calib = NULL, m = nbre_simu, IS = T, n_proj = NULL, 
                                    modele = c("Gaussien" , "JSUo" , "Bootstrap", "Naive Gaussien", "Naive JSUo", 
                                               "DNS original", "DNS Gaussien" , "DNS JSUo","GARCH norm",
                                               "GARCH std", "GARCH sstd","DNS copule","Empirique",
                                               "Hybride"),
                                    trace = F, taux){
  
  ifelse(trace,print(paste0("Modele : " ,modele)),NA)
  
  if(IS) periode_calib = (1:length(x)) else periode_calib = (start_calib:(start_calib+n_calib-1))
  if(IS) periode_proj = (1:length(x)) else periode_proj = ((start_calib+n_calib):(start_calib+n_calib+n_proj-1))
  
  x_calib = x[periode_calib]; x_proj = x[periode_proj]
  date_calib = date[periode_calib]; date_proj = date[periode_proj]
  
  ifelse(trace,print(paste0("IS = ",IS)),NA)
  if(IS==T){
    if(modele == "Gaussien"){
      proj = proj_densite(densfun = "norm", x = x_calib, m = m)
      print("norm")
    } else if (modele == "JSUo"){
      proj = proj_densite(densfun = "JSUo", x = x_calib, m = m)
      print("JSUo")
    } else if(modele == "Bootstrap"){
      proj = proj_bootstrap(x = x_calib, m = m)
      print("Bootstrap")
    } else if(modele == "Naive Gaussien"){
      proj = proj_naive(densfun = "norm", x = x_calib, m = m)
      print("naive gaussien")
    } else if(modele == "Naive JSUo"){
      proj = proj_naive(densfun = "JSUo", x = x_calib, m = m)
      print("naive JSUo")
    } else if(modele == "DNS original"){
      proj = proj_DNS(x = x_calib, m=m, taux = taux)
      print("DNS original")
    } else if (modele == "DNS Gaussien"){
      proj = proj_DNS_remanie(densfun = "norm", x = x_calib, m = m, taux = taux)
      print("DNS gaussien")
    } else if (modele == "DNS JSUo"){
      proj = proj_DNS_remanie(densfun = "JSUo", x = x_calib, m = m, taux = taux)
      print("DNS JSUo")
    } else if (modele == "GARCH norm"){
      proj = proj_GARCH(densfun = "norm", x = x_calib, m = m)
      print("GARCH norm")
    } else if (modele == "GARCH sstd"){
      proj = proj_GARCH(densfun = "sstd", x = x_calib, m = m)
      print("GARCH sstd")
    } else if (modele == "GARCH std"){
      proj = proj_GARCH(densfun = "std", x = x_calib, m = m)
      print("GARCH std")
    } else if (modele == "DNS copule"){
      proj = proj_DNS_copule_cor(x = x_calib, m = m, taux = taux)
      print("DNS copule")
    } else if (modele == "Empirique"){
      proj = list("yield" = x, "incr" = x[2:3664]-x[1:3663])
      print("empirique")
    } else if(modele == "Hybride"){
      proj = proj_Hybrid(x = x_calib, m = m)
      print("hybride")
    }
    
  } else{
    if(modele == "Gaussien"){
      proj = proj_densite(x = x_calib, IS = F, n = n_proj, densfun = "norm", m = m)
    } else if (modele == "JSUo"){
      proj = proj_densite(x = x_calib, IS = F, n = n_proj, densfun = "JSUo", m = m)
    } else if(modele == "Bootstrap"){
      proj = proj_bootstrap(x = x_calib, IS = F, n = n_proj, m = m)
    } else if(modele == "Naive Gaussien"){
      proj = proj_naive(densfun = "norm", x = x_calib, IS = F, n = n_proj, m = m)
    } else if(modele == "Naive JSUo"){
      proj = proj_naive(densfun = "JSUo", x = x_calib, IS = F, n = n_proj, m = m)
    } else if(modele == "DNS original"){
      proj = proj_DNS(x = x_calib, IS = F, n = n_proj,start_calib = start_calib, m = m, taux = taux)
    } else if (modele == "DNS Gaussien"){
      proj = proj_DNS_remanie(densfun = "norm", x = x_calib, IS = F, n = n_proj, 
                              start_calib = start_calib, m = m, taux = taux)
    } else if (modele == "DNS JSUo"){
      proj = proj_DNS_remanie(densfun = "JSUo", x = x_calib, IS = F, n = n_proj, 
                              start_calib = start_calib, m = m, taux = taux)
    } else if (modele == "GARCH norm"){
      proj = proj_GARCH(densfun = "norm", x = x_calib, m = m, IS = F, n = n_proj)
    } else if (modele == "GARCH std"){
      proj = proj_GARCH(densfun = "std", x = x_calib, m = m, IS = F, n = n_proj)
    } else if (modele == "GARCH sstd"){
      proj = proj_GARCH(densfun = "sstd", x = x_calib, m = m, IS = F, n = n_proj)
    } else if (modele == "DNS copule"){
      proj = proj_DNS_copule_cor(x = x_calib, IS = F, n = n_proj,start_calib = start_calib, m = m, taux = taux)
    } else if (modele == "Empirique"){
      proj = list("yield" = x[periode_proj], "incr" = x[(min(periode_proj)+1):max(periode_proj)]-
                    x[min(periode_proj):(max(periode_proj)-1)])
    } else if(modele == "Hybride"){
      proj = proj_Hybrid(x = x_calib, IS = F, n = n_proj, m = m)
    }
    
  }
  res = list("incr" = proj$incr, "yield" = proj$yield)
  return(res)
}  

  

res_OoS = function(x = yield_all_years$`10_ans` , date = yield_all_years$date, start_calib = NULL, 
                       n_calib = NULL, m = nbre_simu, n_proj = NULL, 
                   modeles = c("Gaussien" , "JSUo" , "Bootstrap", "Naive Gaussien", "Naive JSUo", "DNS original", "DNS Gaussien" , 
                              "DNS JSUo","GARCH norm","GARCH std", "GARCH sstd","DNS copule","Empirique", "Hybride")){
  
  
  # Base calib et proj :
  periode_calib = (start_calib:(start_calib+n_calib-1))
  periode_proj = ((start_calib+n_calib):(start_calib+n_calib+n_proj-1))
  
  x_calib = x[periode_calib]; x_proj = x[periode_proj]
  date_calib = date[periode_calib]; date_proj = date[periode_proj]
  
  y_list = list()
  
  # Mod?les :
  if("Gaussien" %in% modeles){
    norm = proj_densite(x = x_calib, IS = F, n = n_proj, densfun = "norm", m = m)$yield
    y_list[["Gaussien"]] = norm
  }
  
  if("JSUo" %in% modeles){
    JSUo = proj_densite(x = x_calib, IS = F, n = n_proj, densfun = "JSUo", m = m)$yield
    y_list[["JSUo"]] = JSUo
  }
  
  if("Bootstrap" %in% modeles){
    boot = proj_bootstrap(x = x_calib, IS = F, n = n_proj, m = m)$yield
    y_list[["Bootstrap"]] = boot
  }
  
  if("Naive Gaussien" %in% modeles){
    naive_gaussien = proj_naive(densfun = "norm", x = x_calib, IS = F, n = n_proj, m = m)$yield
    y_list[["Naive Gaussien"]] = naive_gaussien
  }
  
  if("Naive JSUo" %in% modeles){
    naive_JSUo = proj_naive(densfun = "JSUo", x = x_calib, IS = F, n = n_proj, m = m)$yield
    y_list[["Naive JSUo"]] = naive_JSUo
  }
  
  if("DNS original" %in% modeles){
    DNS_original = proj_DNS(x = x_calib, IS = F, n = n_proj,start_calib = start_calib, m = m)$yield
    y_list[["DNS original"]] = DNS_original
  }
  
  if("DNS Gaussien" %in% modeles){
    DNS_gaussien = proj_DNS_remanie(densfun = "norm", x = x_calib, IS = F, n = n_proj,
                                    start_calib = start_calib, m = m)$yield
    y_list[["DNS Gaussien"]] = DNS_gaussien
  }
  
  if("DNS JSUo" %in% modeles){
    DNS_JSUo = proj_DNS_remanie(densfun = "JSUo", x = x_calib, IS = F, n = n_proj,
                                start_calib = start_calib, m = m)$yield
    y_list[["DNS JSUo"]] = DNS_JSUo
  }
  
  if("GARCH norm" %in% modeles){
    GARCH_norm = proj_GARCH(x = x_calib, IS = F, n = n_proj, densfun = "norm", m = m)$yield
    y_list[["GARCH norm"]] = GARCH_norm
  }
  
  if("GARCH std" %in% modeles){
    GARCH_std = proj_GARCH(x = x_calib, IS = F, n = n_proj, densfun = "std", m = m)$yield
    y_list[["GARCH std"]] = GARCH_std
  }
  
  if("GARCH sstd" %in% modeles){
    GARCH_sstd = proj_GARCH(x = x_calib, IS = F, n = n_proj, densfun = "sstd", m = m)$yield
    y_list[["GARCH sstd"]] = GARCH_sstd
  }
  
  if("DNS copule" %in% modeles){
    DNS_copule = proj_DNS_copule_cor(x = x_calib, IS = F, n = n_proj,
                                     start_calib = start_calib, m = m)$yield
    y_list[["DNS copule"]] = DNS_copule
  }
  
  if("Hybride" %in% modeles){
    Hybride = proj_Hybrid(x = x_calib, IS = F, n = n_proj, m = m)$yield
    y_list[["Hybride"]] = Hybride
  }
  
  # Nombre d'erreurs, ?cart moyen, largeur intervalle terminal, points sous 2
  IC_95 = gg_multi_IC(y = x_proj, y_list = y_list, date = date_proj, probs = c(0.025,0.975))$IC %>% .[2:(n_proj+1),]
  IC_99_5 = gg_multi_IC(y = x_proj, y_list = y_list, date = date_proj, probs = c(0.0025,0.9975))$IC %>% .[2:(n_proj+1),]
  
  res_95 = lapply(1:length(y_list), function(i) c(sum((x_proj<IC_95[,i] | x_proj>IC_95[,i+length(y_list)])),
                                               (sum(abs((x_proj[x_proj<IC_95[,i]]-IC_95[x_proj<IC_95[,i],i]))) + 
                                                  sum(abs(x_proj[x_proj>IC_95[,i+length(y_list)]]-
                                                            IC_95[x_proj>IC_95[,i+length(y_list)],i+length(y_list)])))/
                                                 sum((x_proj<IC_95[,i] | x_proj>IC_95[,i+length(y_list)])),
                                               sum(IC_95[,i]< -2),
                                               c(last(IC_95[,i]),last(IC_95[,i+length(y_list)]))))
  
  res_99_5 = lapply(1:length(y_list), function(i) c(sum((x_proj<IC_99_5[,i] | x_proj>IC_99_5[,i+length(y_list)])),
                                               (sum(abs((x_proj[x_proj<IC_99_5[,i]]-IC_99_5[x_proj<IC_99_5[,i],i]))) + 
                                                  sum(abs(x_proj[x_proj>IC_99_5[,i+length(y_list)]]-
                                                            IC_99_5[x_proj>IC_99_5[,i+length(y_list)],i+length(y_list)])))/
                                                 sum((x_proj<IC_99_5[,i] | x_proj>IC_99_5[,i+length(y_list)])),
                                               sum(IC_99_5[,i]< -2),
                                               c(last(IC_99_5[,i]),last(IC_99_5[,i+length(y_list)]))))
  res_95 = c(res_95,last(x_proj))
  res_99_5 = c(res_99_5,last(x_proj))
  names(res_95) = c(names(y_list),"y_last") ; names(res_99_5) = c(names(y_list),"y_last")
  return(list(res_95,res_99_5))
}



projection_modeles_stats_IS = function(x , m ,  trace = T, taux ,
                                       modeles = c("Gaussien","JSUo","Bootstrap","naive gaussien",
                                                   "naive JSUo", "DNS original", "DNS gaussien", "DNS JSUo",
                                                   "GARCH sstd", "DNS copule","Hybride")){
  
  
  x_calib = x
  
  y_list = list()
  
  if("Gaussien" %in% modeles){
    norm = proj_densite(densfun = "norm", x = x_calib, m = m)
    ifelse(trace,print("Mod?le Gaussien"),NULL)
    y_list[["Gaussien"]] = norm
    
  }  
  
  if("JSUo" %in% modeles){ 
    JSUo = proj_densite(densfun = "JSUo", x = x_calib, m = m)
    ifelse(trace,print("Mod?le JSUo"),NULL)
    y_list[["JSUo"]] = JSUo
    
  }  
  
  if("Bootstrap" %in% modeles){
    boot = proj_bootstrap(x = x_calib, m = m)
    ifelse(trace,print("Mod?le Bootstrap"),NULL)
    y_list[["Bootstrap"]] = boot
    
  }  
  
  if("naive gaussien" %in% modeles){
    naive_gaussien = proj_naive(densfun = "norm", x = x_calib, m = m)
    ifelse(trace,print("Naive Gaussien"),NULL)
    y_list[["naive gaussien"]] = naive_gaussien
    
  }  
  
  if("naive JSUo" %in% modeles){
    naive_JSUo = proj_naive(densfun = "JSUo", x = x_calib, m = m)
    ifelse(trace,print("Naive JSUo"),NULL)
    y_list[["naive JSUo"]] = naive_JSUo
    
  }  
  
  if("DNS original" %in% modeles){
    DNS_original = proj_DNS(x = x_calib, m=m, taux = taux)
    ifelse(trace,print("DNS original"),NULL)
    y_list[["DNS original"]] = DNS_original
    
  }  
  
  if("DNS gaussien" %in% modeles){
    DNS_gaussien = proj_DNS_remanie(densfun = "norm", x = x_calib, m = m, taux = taux)
    ifelse(trace,print("DNS Gaussien"),NULL)
    y_list[["DNS gaussien"]] = DNS_gaussien
    
  }  
  
  if("DNS JSUo" %in% modeles){
    DNS_JSUo = proj_DNS_remanie(densfun = "JSUo", x = x_calib, m = m, taux = taux)
    ifelse(trace,print("DNS JSUo"),NULL)
    y_list[["DNS JSUo"]] = DNS_JSUo
    
  }  
  
  if("GARCH sstd" %in% modeles){
    GARCH = proj_GARCH(densfun = "sstd", x = x_calib, m = m)
    ifelse(trace,print("GARCH"),NULL)
    y_list[["GARCH sstd"]] = GARCH
    
  }  
  
  if("DNS copule" %in% modeles){
    DNS_copule = proj_DNS_copule_cor(x = x_calib, m=m, taux = taux)
    ifelse(trace,print("DNS Copule"),NULL)
    y_list[["DNS copule"]] = DNS_copule
    
  }  
  
  if("Hybride" %in% modeles){
    loi_hybride = proj_Hybrid(x = x_calib, m=m)
    ifelse(trace,print("Loi hybride"),NULL)
    y_list[["Hybride"]] = loi_hybride
  }
  
  
  
  
  y_l_yield = lapply(y_list, function(x) x$yield)
  
  y_l_incr= lapply(y_list, function(x) x$incr)
  
  moy = lapply(y_l_incr, function(x) mean(x))
  var = lapply(y_l_incr, function(x) var(x))
  skew = lapply(y_l_incr, function(x) skewness(x))
  kurt = lapply(y_l_incr, function(x) kurtosis(x))
  
  stat_95 = lapply(y_l_yield, function(x) 
    lapply(1:ncol(x), function(i) quantile(x[,i], probs = c(0.025,0.975)))) %>% unlist %>% matrix(byrow = T, ncol = 2) %>%
    {sum( (x < as.numeric(.[,1])) | (x > as.numeric(.[,2])) )/length(x)}
  
  stat_99_5 = lapply(y_l_yield, function(x) 
    lapply(1:ncol(x), function(i) quantile(x[,i], probs = c(0.0025,0.9975)))) %>% unlist %>% matrix(byrow = T, ncol = 2) %>%
    {sum( (x < as.numeric(.[,1])) | (x > as.numeric(.[,2])) )/length(x)}
  
  res = list("stat_95" = stat_95, "stat_99_5" = stat_99_5, "mean" = moy, "var" = var, "skew" = skew, "kurt" = kurt)
  return(res)
}  

projection_modeles_stats_OoS = function(x , m ,  trace = T , n_proj = 1000, n_calib = 1000, pas = 10,
                                        modeles = c("Gaussien","JSUo","Bootstrap","naive gaussien",
                                                    "naive JSUo", "DNS original", "DNS gaussien", "DNS JSUo",
                                                    "GARCH sstd", "DNS copule","Hybride"), taux){
  print(modeles)
  index = seq(1,(length(x)-(n_proj+n_calib)), by = pas)
  
  nbre_cores = detectCores(logical = T) - 1
  registerDoParallel(cores = nbre_cores)
  
  stats = foreach(i=1:length(index), .packages = c("gamlss","magrittr","purrr","data.table","ggplot2","fGarch",
                                                   "copula"), 
                  .export = c("res_OoS","proj_densite","proj_bootstrap","proj_naive","proj_DNS",
                              "gg_multi_IC","projection_single_modele","MA1_estimate","loadings_DNS","beta_optim",
                              "residus_NS","proj_GARCH","proj_DNS_copule_cor","simu_JSUo","proj_Hybrid",
                              "cor_epsilon","alea_JSUo","param_Hybrid","rHybrid","qHybrid","log_log_dt","yield_all_years")) %dopar% {
                                
                                sim = projection_single_modele(x=x,date=yield_all_years$date,start_calib = index[i],n_calib = n_calib, 
                                                               n_proj = n_proj, 
                                                               m= m, modele = modeles,IS =F, taux = taux)
                                moy = mean(sim$incr)
                                var = var(sim$incr)
                                skew = skewness(sim$incr)
                                kurt = kurtosis(sim$incr)
                                
                                yield_proj = x[(index[i]+n_calib):(index[i] +n_calib+n_proj-1)]
                                
                                sim$yield = sim$yield[,-1]
                                stat_95 = lapply(1:ncol(sim$yield), function(i) quantile(sim$yield[,i], probs = c(0.025,0.975))) %>% unlist %>% 
                                  matrix(byrow = T, ncol = 2) %>% 
                                  {sum( (yield_proj < as.numeric(.[,1])) | (yield_proj > as.numeric(.[,2])) )/length(yield_proj)}
                                
                                stat_99_5 = lapply(1:ncol(sim$yield), function(i) quantile(sim$yield[,i], probs = c(0.0025,0.9975))) %>% unlist %>% 
                                  matrix(byrow = T, ncol = 2) %>% 
                                  {sum( (yield_proj < as.numeric(.[,1])) | (yield_proj > as.numeric(.[,2])) )/length(yield_proj)}
                                
                                list("stat_95" = stat_95, "stat_99_5" = stat_99_5, "mean" = moy, "var" = var, "skew" = skew, "kurt" = kurt)
                                
                              }
  
  stat_95 = map(stats,1) %>% unlist %>% mean
  stat_99_5 = map(stats,2) %>% unlist %>% mean
  moy = map(stats,3) %>% unlist %>% mean
  var = map(stats,4) %>% unlist %>% mean
  skew = map(stats,5) %>% unlist %>% mean
  kurt = map(stats,6) %>% unlist %>% mean
  
  res = list("stat_95" = stat_95, "stat_99_5" = stat_99_5, "moy" = moy, "var" = var, "skew" = skew, "kurt" = kurt)
  
  return(res)
}  


# x = yield_all_years$`10_ans` ; date = yield_all_years$date;
# n_calib = 1000; m = 100; n_proj = 1000; pas = 1000 ; modeles = "Gaussien"
# proj_OoS_par()

proj_OoS_par = function(x = yield_all_years$`10_ans` , date = yield_all_years$date, 
                        n_calib = 1000, m = nbre_simu, n_proj = 1000, pas = 100, 
                        modeles = c("Gaussien" , "JSUo" , "Bootstrap", "Naive Gaussien", "Naive JSUo", "DNS original", "DNS Gaussien" , 
                                    "DNS JSUo","GARCH norm","GARCH std", "GARCH sstd","DNS copule","Empirique", "Hybride")){
  
  index = seq(1,(length(x)-(n_proj+n_calib)), by = pas)
  
  # a1 = Sys.time() #35 sec pour 17 simu
  nbre_cores = detectCores(logical = T) - 1
  registerDoParallel(cores = nbre_cores)
  
  res = foreach(i=1:length(index), .packages = c("gamlss","magrittr","purrr","data.table","ggplot2","fGarch",
                                                 "copula"), 
                .export = c("res_OoS","proj_densite","proj_bootstrap","proj_naive","proj_DNS","proj_DNS_remanie",
                            "gg_multi_IC","projection_modeles","MA1_estimate","loadings_DNS","beta_optim",
                            "residus_NS","proj_GARCH","proj_DNS_copule_cor","simu_JSUo","proj_Hybrid",
                            "cor_epsilon","alea_JSUo","param_Hybrid","rHybrid","qHybrid","log_log_dt")) %dopar% {
    sim = res_OoS(x=x,date=date,start_calib = index[i],n_calib = n_calib, n_proj = n_proj, m= m, modeles = modeles)
    sim
  }
  
  label_modeles = names(res[[1]][[1]])[1:(length(modeles))]
  
  # Stats par modèle
  #95%
  erreur= foreach(i=1:(length(modeles))) %do% {
    erreur_point = unlist(lapply(1:length(index),function(j) unlist(map(map(res[[j]],i),1)))) %>% matrix(byrow = T, ncol = 2)
    ecart_moyen = unlist(lapply(1:length(index),function(j) unlist(map(map(res[[j]],i),2)))) %>% matrix(byrow = T, ncol = 2)
    sous_2 = unlist(lapply(1:length(index),function(j) unlist(map(map(res[[j]],i),3)))) %>% matrix(byrow = T, ncol = 2)
    largeur_inf = unlist(lapply(1:length(index),function(j) unlist(map(map(res[[j]],i),4)))) %>% matrix(byrow = T, ncol = 2)
    largeur_sup = unlist(lapply(1:length(index),function(j) unlist(map(map(res[[j]],i),5)))) %>% matrix(byrow = T, ncol = 2)
    list("95" = matrix(c(erreur_point[,1],ecart_moyen[,1],sous_2[,1],largeur_inf[,1],largeur_sup[,1]), byrow = F, ncol = 5),
         "99_5" = matrix(c(erreur_point[,2],ecart_moyen[,2],sous_2[,2],largeur_inf[,2],largeur_sup[,2]), byrow = F, ncol = 5))
  }
  names(erreur) = label_modeles
  
  
  erreur_95 = map(erreur,1)
  erreur_99_5 = map(erreur,2)
  
  nbre_erreur_95 = unlist(lapply(1:length(modeles), function(i) erreur_95[[i]][,1])) %>% matrix(byrow = F, ncol = length(modeles)) %>%
    as.data.table %>% setnames(label_modeles)
  
  nbre_erreur_99_5 = unlist(lapply(1:length(modeles), function(i) erreur_99_5[[i]][,1])) %>% matrix(byrow = F, ncol = length(modeles)) %>%
    as.data.table %>% setnames(label_modeles)
  
  ecart_moyen_95 = unlist(lapply(1:length(modeles), function(i) erreur_95[[i]][,2])) %>% matrix(byrow = F, ncol = length(modeles)) %>%
    as.data.table %>% setnames(label_modeles)
  
  ecart_moyen_99_5 = unlist(lapply(1:length(modeles), function(i) erreur_99_5[[i]][,2])) %>% matrix(byrow = F, ncol = length(modeles)) %>%
    as.data.table %>% setnames(label_modeles)
  
  # sous_2 = foreach(i=1:(length(modeles)), .combine = cbind) %do% {
  #   erreur_95[[i]][,3]
  # }
  # sous_2 = cbind(sous_2, 1:length(index)) ; colnames(sous_2) = c(label_modeles,"index")
  # 
  # largeur = foreach(i=1:(length(modeles)), .combine = cbind) %do% {
  #   matrix(c(erreur_95[[i]][,4],erreur_95[[i]][,5]), byrow = F, ncol = 2)
  # }
  
  # yield_last = unlist(map(res[[1]],length(res[[1]][[1]])))
  # largeur = cbind(largeur,yield_last)
  # names_largeur = matrix(c(rep(label_modeles,each = 2), rep(c("min","max"),length(label_modeles))), ncol = 2) %>%
  # {unlist(lapply(1:nrow(.), function(i) paste0(.[i,1],"_",.[i,2])))} %>% c("yield_last")
  # colnames(largeur) = names_largeur
  
  # g_erreur = nbre_erreur_95 %>% as.data.table %>% melt(id.vars = "index") %>% ggplot(aes(index,value)) +
  #   geom_line(aes(col=variable)) + xlab("") + ylab("Nombre erreur")
  # g_ecart_moyen = ecart_moyen_95 %>% as.data.table %>% melt(id.vars = "index") %>% ggplot(aes(index,value)) +
  #   geom_line(aes(col=variable)) + xlab("") + ylab("Ecart moyen")
  # g_sous_2 = sous_2 %>% as.data.table %>% melt(id.vars = "index") %>% ggplot(aes(index,value)) +
  #   geom_line(aes(col=variable)) + xlab("") + ylab("Point sous 2")

  
  # names_largeur = names_largeur[-length(names_largeur)]
  # g_largeur_geom = mapply(function(u,v)  paste0("geom_line(aes(y =", paste0("c(",paste0(largeur[,u], 
  #                                                                                  collapse = ","),")"),
  #                              " ,col =",paste0("\"",
  #                                               v,"\""),"))"), 
  #        u = 1:(ncol(largeur)-1), v = map(strsplit(names_largeur,split = "_"),1), SIMPLIFY = F) %>% 
  #   unlist %>% paste(collapse = "+")
  # 
  # g_largeur = quote(ggplot(largeur,aes(1:nrow(largeur),yield_last)) + geom_line() + xlab("") + 
  #                     ylab("Largeur intervalle") + labs(title = "Intervalle(s) de confiance terminaux")) %>% 
  #   deparse %>% paste(collapse = " ") %>% paste("+", g_largeur_geom)%>% {parse(text = .)} %>% eval
  
  
  # plots = list("Erreur" = g_erreur, "Ecart_moyen" = g_ecart_moyen, "Sous_2" = g_sous_2, "Largeur" =
  #                g_largeur, "data_erreur_95" = nbre_erreur_95, "data_ecart_moyen_95" = ecart_moyen_95,
  #              "data_erreur_99_5" = nbre_erreur_99_5, "data_ecart_moyen_99_5" = ecart_moyen_99_5,
  #              "data_sous_2" = sous_2, "data_largeur" = largeur)
  
  res = list("Nbre_erreur_95" = nbre_erreur_95, "Nbre_erreur_99_5" = nbre_erreur_99_5,
             "Ecart_moyen_95" = ecart_moyen_95, "Ecart_moyen_99_5" = ecart_moyen_99_5)
  
  return(res)
}





# Fonctions projections chocs S2 ---------------------------------------------------

# A dé-commenter pour cette section :

# residus_NS = foreach(i=1:30) %dopar% {
#   unlist(lapply(1:3664, function(j) yield_all_years[j,i] - NS_YC(t = i, par = beta_optim[j,])))
# }

reg_BETA_EIOPA = function(theta){
  t = 1:20 ; lambda = 0.598
  beta_1 = theta[1] ; beta_2 = theta[2] ; beta_3 = theta[3]
  
  y_th = beta_1 + beta_2*((1-exp(-lambda*t))/(lambda*t)) + beta_3*((1-exp(-lambda*t))/(lambda*t) - 
                                                                     exp(-lambda*t))
  return( sum((y_th-y)^2))
}

estim_par_EIOPA = function(yc){
  loadings = loadings_DNS(t = c(2,5,10)) %>% inv
  theta = c(yc[2],yc[5],yc[10]) %*% loadings 
  
  y = yc
  beta_initiaux = optim(par = theta, fn = reg_BETA_EIOPA)$par
  
  residus_initiaux = y - NS_YC(t = 1:20, par = beta_initiaux)
  
  return(list("beta" = beta_initiaux, "residus" = residus_initiaux))
}
  
calibration_DNS_original = function(x,taux){
  
  loadings = as.numeric(loadings_DNS(t = taux))
  
  beta = beta_optim[1:length(x),]
  residus_10_ans = residus_NS[[taux]]
  
  AR_beta = lapply(1:3, function(i) lm(beta[2:nrow(beta),i]~beta[1:(nrow(beta)-1),i]) %>%
  {list("intercept" = .$coef[1], "phi" = .$coef[2], "sd" = sd(.$residuals))}) 
  
  fit_residus = list("mean" = mean(residus_10_ans), "sd" = sd(residus_10_ans))
  
  
  return(list("AR_beta" = AR_beta, "fit_residus" = fit_residus))
}

choc_S2_DNS_original = function(param_EIOPA, param_DNS, m = 10^4, n = 260, taux){
  
  beta_initiaux = param_EIOPA$beta ; residu_initial = param_EIOPA$residus[taux]
  AR_beta = param_DNS$AR_beta ; fit_residus = param_DNS$fit_residus
  
  loadings = loadings_DNS(t = taux)
  beta_th = list("B1" = matrix(rep(0,m*(n+1)), nrow = m), "B2" = matrix(rep(0,m*(n+1)), nrow = m),
                 "B3" = matrix(rep(0,m*(n+1)), nrow = m))
  
  for(i in 1:3){
    beta_th[[i]][,1] = rep(beta_initiaux[i],m)
  }
  
  res_beta_th = lapply(AR_beta, function(x) matrix(replicate(m,rnorm(n, mean = 0, sd = x$sd)), byrow = T, 
                                                         nrow = m))
  
  residus_th = matrix(replicate(m,rnorm(n, mean = fit_residus$mean, sd = fit_residus$sd)),byrow = T, nrow = m) 
  
  for(i in 1:3){
    for(j in 2:(n+1)){
      beta_th[[i]][,j] = AR_beta[[i]]$intercept + AR_beta[[i]]$phi*beta_th[[i]][,j-1] + res_beta_th[[i]][,j-1]
    }
  }
  
  residus_th = cbind(rep(residu_initial,m),residus_th)
  
  yield_th = beta_th[[1]]*loadings[1] + beta_th[[2]]*loadings[2] + beta_th[[3]]*loadings[3] + residus_th
  
  quant = quantile(yield_th[,n+1], probs = c(0.005,0.995))
  
  return(quant)
  
}

calibration_DNS_cor = function(x,taux){
  
  loadings = as.numeric(loadings_DNS(t = taux))
  
  beta = beta_optim[1:length(x),]
  residus_10_ans = residus_NS[[taux]]
  
  beta_incr = beta[2:nrow(beta),] - beta[1:(nrow(beta)-1),]
  resi_incr = residus_10_ans[2:length(residus_10_ans)] - residus_10_ans[1:(length(residus_10_ans)-1)]
  
  B1_incr_JSUo = gamlssML(beta_incr[,1], family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                    "tau" = .$tau,"nu" = .$nu)}
  B2_incr_JSUo = gamlssML(beta_incr[,2], family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                    "tau" = .$tau,"nu" = .$nu)}
  B3_incr_JSUo = gamlssML(beta_incr[,3], family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                    "tau" = .$tau,"nu" = .$nu)}
  fit_resi_incr = gamlssML(resi_incr,family = "JSUo") %>% {list("mu" = .$mu,"sigma" = .$sigma,
                                                                "tau" = .$tau,"nu" = .$nu)}
  
  
  a_B1 = alea_JSUo(par_JSUo = B1_incr_JSUo, indice = beta_incr[,1])
  a_B2 = alea_JSUo(par_JSUo = B2_incr_JSUo, indice = beta_incr[,2])
  a_B3 = alea_JSUo(par_JSUo = B3_incr_JSUo, indice = beta_incr[,3])
  a_R = alea_JSUo(par_JSUo = fit_resi_incr, indice = resi_incr)
  
  cor_mat = cbind(a_B1,a_B2,a_B3,a_R) %>% cor
  
  
  return(list("B1_incr_JSUo" = B1_incr_JSUo, "B2_incr_JSUo" = B2_incr_JSUo, "B3_incr_JSUo" = B3_incr_JSUo, "fit_resi_incr" = fit_resi_incr,
              "cor_mat" = cor_mat))
}

choc_S2_DNS_cor = function(param_EIOPA, param_DNS, m = 10^4, n = 260, taux){
  
  beta_initiaux = param_EIOPA$beta ; residu_initial = param_EIOPA$residus[taux]
  
  cor_mat = param_DNS$cor_mat ; B1_incr_JSUo = param_DNS$B1_incr_JSUo ; B2_incr_JSUo = param_DNS$B2_incr_JSUo 
  B3_incr_JSUo = param_DNS$B3_incr_JSUo ; fit_resi_incr = param_DNS$fit_resi_incr
  
  loadings = loadings_DNS(t = taux)
  
  res_th = lapply(1:m, function(i) cor_epsilon(matrice_cor = cor_mat, n=n) %>%
  {list("B1" = c(beta_initiaux[1],simu_JSUo(par = B1_incr_JSUo, epsilon = .[,1])) %>% cumsum,
        "B2" = c(beta_initiaux[2],simu_JSUo(par = B2_incr_JSUo, epsilon = .[,2])) %>% cumsum,
        "B3" = c(beta_initiaux[3],simu_JSUo(par = B3_incr_JSUo, epsilon = .[,3])) %>% cumsum,
        "R" = c(residu_initial,simu_JSUo(par = fit_resi_incr, epsilon = .[,4])) %>% cumsum)})
  
  yield_th = matrix(unlist(lapply(1:m, function(i) res_th[[i]][[1]] + loadings[2]*res_th[[i]][[2]] + 
                                    loadings[3]*res_th[[i]][[3]] + res_th[[i]][[4]])), byrow = T, nrow = m)
  
  quant = quantile(yield_th[,n+1], probs = c(0.005,0.995))
  
  return(quant)
}


# 
# calibration_DNS_cor_TS = function(x){
#   
#   n = length(x)
#   
#   beta = beta_optim[1:n,] %>% {.[2:n,] - .[1:(n-1),]}
#   residus = residus_NS %>% unlist %>% matrix(ncol = 20) %>% .[1:(n-1),]
#   
#   series = cbind(beta,residus)
#                 
#   fit_JSUo = lapply(1:ncol(series), function(i) series[,i] %>% gamlssML(family = "JSUo") %>% 
#                            {list("mu" = .$mu,"sigma" = .$sigma,"tau" = .$tau,"nu" = .$nu)})
#   
#   a_R = lapply(1:23, function(i) alea_JSUo(par_JSUo = fit_JSUo[[i]], indice = series[,i])) %>% 
#     unlist %>% matrix(ncol = 23)
#   
#   cor_mat = a_R %>% cor
#   
#   
#   return(list("fit_JSUo" = fit_JSUo, "cor_mat" = cor_mat))
# }
# 
# choc_S2_DNS_cor_TS = function(param_EIOPA, param_DNS, m = 10^4, n = 260){
#   
#   beta_initiaux = param_EIOPA$beta ; residus_initiaux = param_EIOPA$residus
#   
#   cor_mat = param_DNS$cor_mat ; fit_JSUo = param_DNS$fit_JSUo 
#   
#   loadings = loadings_DNS(t = 1:20)
#   
#   valeurs_initiales = c(beta_initiaux,residus_initiaux)
#   
#   res_th = lapply(1:m, function(j) lapply(1:length(valeurs_initiales), function(i) cor_epsilon(matrice_cor = cor_mat, n = n) %>%
#                                                                     {c(valeurs_initiales[i], simu_JSUo(par = fit_JSUo[[i]], 
#                                                                                                       epsilon = .[,i]))} %>%
#                                                                     cumsum) %>% unlist %>% matrix(ncol = 23))
#   
#   # yield_th = lapply(1:20, function(j) (matrix(unlist(lapply(1:m, function(i) res_th[[i]][[1]] + loadings[2]*res_th[[i]][[2]] + 
#   #                                   loadings[3]*res_th[[i]][[3]] + res_th[[i]][[4]])), byrow = T, nrow = m)))
#   
#   yield_th = lapply(1:m, function(j) lapply(1:20, function(i) res_th[[j]][,1] + loadings[2,i]*res_th[[j]][,2] + loadings[3,i]*res_th[[j]][,3] +
#                                               res_th[[j]][,i+3]) %>% unlist %>% matrix(ncol = 20) %>% t )
#   
#   quantile = lapply(1:20, function(j) lapply(1:m, function(i) yield_th[[i]][j,261]) %>% unlist %>% {.[order(.)]} %>% .[c(0.0025*m,0.9975*m)]) %>%
#     unlist %>% matrix(byrow = T,ncol = 2)
#   
#   
#   yield_low = yield_th[[which_sim[1]]][,261]
#   yield_up = yield_th[[which_sim[2]]][,261]
#   quant = quantile(yield_th[,n+1], probs = c(0.0025,0.9975))
#   
#   
#   return(quant)
# }
# 
# 
# 
# calibration_DNS_original_TS = function(x){
#   
#   
#   beta = beta_optim[1:length(x),]
#   residus = residus_NS %>% unlist %>% matrix(ncol = 20)
#   
#   AR_beta = lapply(1:3, function(i) lm(beta[2:nrow(beta),i]~beta[1:(nrow(beta)-1),i]) %>%
#   {list("intercept" = .$coef[1], "phi" = .$coef[2], "sd" = sd(.$residuals))}) 
#   
#   fit_residus = lapply(1:ncol(residus) , function(j) list("mean" = mean(residus[,j]), "sd" = sd(residus[,j])))
#   
#   
#   return(list("AR_beta" = AR_beta, "fit_residus" = fit_residus))
# }
# 
# choc_S2_DNS_original_TS = function(param_EIOPA, param_DNS, m = 10^4, n = 260){
#   
#   beta_initiaux = param_EIOPA$beta ; residu_initial = param_EIOPA$residus
#   AR_beta = param_DNS$AR_beta ; fit_residus = param_DNS$fit_residus
#   
#   loadings = loadings_DNS(t = 1:20)
#   beta_th = list("B1" = matrix(rep(0,m*(n+1)), nrow = m), "B2" = matrix(rep(0,m*(n+1)), nrow = m),
#                  "B3" = matrix(rep(0,m*(n+1)), nrow = m))
#   
#   for(i in 1:3){
#     beta_th[[i]][,1] = rep(beta_initiaux[i],m)
#   }
#   
#   res_beta_th = lapply(AR_beta, function(x) matrix(replicate(m,rnorm(n, mean = 0, sd = x$sd)), byrow = T, 
#                                                    nrow = m))
#   
#   residus_th = lapply(1:20, function(j) cbind(residu_initial[j],matrix(replicate(m,rnorm(n, mean = fit_residus[[j]]$mean, 
#                                                                                          sd = fit_residus[[j]]$sd)),byrow = T, nrow = m)))
#   
#   for(i in 1:3){
#     for(j in 2:(n+1)){
#       beta_th[[i]][,j] = AR_beta[[i]]$intercept + AR_beta[[i]]$phi*beta_th[[i]][,j-1] + res_beta_th[[i]][,j-1]
#     }
#   }
# 
#   yield_th = lapply(1:20, function(i) beta_th[[1]]*loadings[1,i] + beta_th[[2]]*loadings[2,i] + beta_th[[3]]*loadings[3,i] + residus_th[[i]])
#   
#   quant = quantile(yield_th[,n+1], probs = c(0.005,0.995))
#   
#   return(quant)
#   
# }
# 















