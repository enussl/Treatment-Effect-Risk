rm(list = ls())

library(dplyr)
library(tidyr)
library(cvar)
library(ggplot2)
library(ggsci)
library(gganimate)
library(stats)


setwd("C:/Users/eminu/OneDrive/Desktop/Treatment-Effect-Risk")


simulate.full.data = function(shift, sigma.1, sigma.0,
                              n.obs, mu.1, mu.0, alpha){
  
  # shift: upper limit of uniform r.v. (controls mean shifts)
  # sigma.1: upper limit of std of Y_1
  # sigma.0: upper limit of std of Y_0
  # n.obs: number of observations
  # mu.1: mean Y_1
  # mu.0: mean Y_0
  # alpha: alpha of CVaR
  
  set.seed(42)
  
  max.var = sigma.1^2+sigma.0^2+2*sigma.1*sigma.0
  x = seq(1,10000,1)
  corr = seq(-1, 1, l = 20)
  sd = list(seq(0, sigma.1, l = 20), seq(0, sigma.0, l = 20))

  alpha = 0.05
  
  results.sigma.shift = matrix(data = NA, nrow = length(sd[[1]])*length(corr), ncol = 7)
  
  v = 1
  for (c in corr){
    for(k in 1:length(sd[[1]])){
      
      tau.vec = numeric()
      cdte.mat = matrix(data = NA, nrow = length(x), ncol = 3)
      
      i = 1
      for (j in 1:length(x)){
        
        shift.1 = runif(n = 1, 0, shift)
        shift.0 = runif(n = 1, 0, shift)
        y.1 = rnorm(n = n.obs, mean = mu.1+shift.1, sd = sd[[1]][k])
        y.0 = rnorm(n = n.obs, mean = mu.0+shift.0, sd = sd[[2]][k])
        
        delta = rnorm(n = n.obs, mean = mu.1+shift.1-(mu.0+shift.0),
                      sd = sqrt(sd[[1]][k]^2+sd[[2]][k]^2-2*c*sd[[1]][k]*sd[[2]][k]))
        
        tau.vec[j] = mean(y.1)-mean(y.0)
        
        cdte.mat[i,1] = max(ES(y.1, p_loss = alpha)-ES(y.0, p_loss = alpha),
                            ES(-1*y.0, p_loss = alpha)-ES(-1*y.1, p_loss = alpha))
        cdte.mat[i,2] = ES(delta, p_loss = alpha)
        cdte.mat[i,3] = ES(y.1, p_loss = alpha) + ES(-1*y.0, p_loss = alpha)
        i = i + 1
      }
      
      cvar.kall = ES(tau.vec, p_loss = alpha)
      cvar.us = colMeans(cdte.mat)
      
      results.sigma.shift[v,1] = cvar.us[1]
      results.sigma.shift[v,2] = cvar.us[2]
      results.sigma.shift[v,3] = cvar.us[3]
      results.sigma.shift[v,4] = cvar.kall
      results.sigma.shift[v,5] = c
      results.sigma.shift[v,6] = sd[[1]][k]
      results.sigma.shift[v,7] = sd[[2]][k]
      v = v + 1
    }
  }
  
  sigma.shift.df = data.frame(results.sigma.shift) %>%
    rename("LB" = "X1",
           "True" = "X2",
           "UB" = "X3",
           "Kallus" = "X4",
           "Corr" = "X5",
           "Se.1" = "X6",
           "Se.0" = "X7") %>%
    mutate(var.delta.x = Se.1^2+Se.0^2-2*Corr*Se.1*Se.0,
           var.delta.x.group = cut(var.delta.x, breaks = 30,
                                   include.lowest = T),
           var.delta.x.group = round((as.numeric(var.delta.x.group)-1)*max.var/30, 1),
           var.delta.x.group = factor(var.delta.x.group)) %>%
    group_by(var.delta.x.group) %>%
    summarise(LB = mean(LB),
              True = mean(True),
              UB = mean(UB),
              Kallus = mean(Kallus)) %>%
    pivot_longer(cols = c("LB", "True", "UB", "Kallus")) %>%
    arrange(name, value)
  
  return(sigma.shift.df)
}



simulate.corr.static = function(shift, sigma.1, sigma.0,
                         n.obs, mu.1, mu.0,
                         rho){
  
  # shift: upper limit of uniform r.v. (controls mean shifts)
  # sigma.1: std of Y_1
  # sigma.0: std of Y_0
  # n.obs: number of observations
  # mu.1: mean Y_1
  # mu.0: mean Y_0
  # rho: correlation pot. outcomes
  
  
  set.seed(42)
  n.sample = n.obs
  
  y.1 = rnorm(n = n.sample, mean = mu.1, sd = sigma.1)
  y.0 = rnorm(n = n.sample, mean = mu.0, sd = sigma.0)
  
  delta = rnorm(n = n.sample, mean = mu.1-mu.0,
                sd = sqrt(sigma.1^2+sigma.0^2-2*rho*sigma.1*sigma.0))
  
  alpha = seq(0, 1, 0.05)
  
  results = matrix(data = NA, nrow = length(alpha), ncol = 3)
  
  i = 1
  for(j in alpha){
    results[i,1] = max(ES(y.1, p_loss = j)-ES(y.0, p_loss = j),
                       ES(-1*y.0, p_loss = j)-ES(-1*y.1, p_loss = j))
    results[i,2] = ES(delta, p_loss = j)
    results[i,3] = ES(y.1, p_loss = j) + ES(-1*y.0, p_loss = j)
    i = i + 1
  }
  
  
  results.df = data.frame(results) %>%
    rename("LB" = "X1",
           "True" = "X2",
           "UB" = "X3") %>%
    mutate(alpha = alpha) %>%
    pivot_longer(cols = c("LB", "True", "UB"))
  
  return(results.df)
}


simulate.corr.dynamic = function(shift, sigma.1, sigma.0,
                                 n.obs, mu.1, mu.0){
  
  set.seed(42)
  n.sample = n.obs
  
  corr = seq(-1, 1, l = 20)
  results.corr = list()
  
  v = 1
  for (c in corr){
    
    y.1 = rnorm(n = n.sample, mean = mu.1, sd = sigma.1)
    y.0 = rnorm(n = n.sample, mean = mu.0, sd = sigma.0)
    
    delta = rnorm(n = n.sample, mean = mu.1-mu.0,
                  sd = sqrt(sigma.1^2+sigma.0^2-2*c*sigma.1*sigma.0))
    
    alpha = seq(0, 1, 0.05)
    results = matrix(data = NA, nrow = length(alpha), ncol = 4)
    
    i = 1
    for(j in alpha){
      results[i,1] = max(ES(y.1, p_loss = j)-ES(y.0, p_loss = j),
                         ES(-1*y.0, p_loss = j)-ES(-1*y.1, p_loss = j))
      results[i,2] = ES(delta, p_loss = j)
      results[i,3] = ES(y.1, p_loss = j) + ES(-1*y.0, p_loss = j)
      results[i,4] = j
      i = i + 1
    }
    
    results.corr[[v]] = data.frame(results) %>%
      mutate(corr = c)
    v = v + 1
  }
  
  results = do.call(rbind, results.corr) %>%
    rename("LB" = "X1", "True" = "X2", "UB" = "X3", "alpha" = "X4") %>%
    pivot_longer(cols = c("LB", "True", "UB")) %>%
    mutate(corr = round(corr, 2))
  return(results)
}


# simulate.corr.dynamic = function(shift, sigma.1, sigma.0,
#                                  n.obs, mu.1, mu.0){
#   
#   set.seed(42)
#   n.sample = n.obs
#   
#   corr = seq(-1, 1, l = 20)
#   results.corr = list()
#   
#   v = 1
#   for (c in corr){
#     
#     y.1 = rnorm(n = n.sample, mean = mu.1, sd = sigma.1)
#     y.0 = rnorm(n = n.sample, mean = mu.0, sd = sigma.0)
#     
#     delta.conv = stats::convolve(x = y.1, y = rev(y.0))
#     
#     delta = rnorm(n = n.sample, mean = mu.1-mu.0,
#                   sd = sqrt(sigma.1^2+sigma.0^2-2*c*sigma.1*sigma.0))
#     
#     alpha = seq(0, 1, 0.05)
#     results = matrix(data = NA, nrow = length(alpha), ncol = 5)
#     
#     i = 1
#     for(j in alpha){
#       results[i,1] = max(ES(y.1, p_loss = j)-ES(y.0, p_loss = j),
#                          ES(-1*y.0, p_loss = j)-ES(-1*y.1, p_loss = j))
#       results[i,2] = ES(delta, p_loss = j)
#       results[i,3] = ES(y.1, p_loss = j) + ES(-1*y.0, p_loss = j)
#       results[i,4] = ES(delta.conv, p_loss = j)
#       results[i,5] = j
#       i = i + 1
#     }
#     
#     results.corr[[v]] = data.frame(results) %>%
#       mutate(corr = c)
#     v = v + 1
#   }
#   
#   results = do.call(rbind, results.corr) %>%
#     rename("LB" = "X1", "True" = "X2", "UB" = "X3", "Ind" = "X4" ,"alpha" = "X5") %>%
#     pivot_longer(cols = c("LB", "True", "UB", "Ind")) %>%
#     mutate(corr = round(corr, 2))
#   return(results)
# }


simulate.corr.full = function(shift, sigma.1, sigma.0,
                              n.obs, mu.1, mu.0, alpha){
  
  # shift: upper limit of uniform r.v. (controls mean shifts)
  # sigma.1: upper limit of std of Y_1
  # sigma.0: upper limit of std of Y_0
  # n.obs: number of observations
  # mu.1: mean Y_1
  # mu.0: mean Y_0
  # alpha: alpha of CVaR
  
  set.seed(42)
  x = seq(1,10000,1)
  n.sample = n.obs
  
  corr = seq(-1, 1, l = 20)
  alpha = 0.05
  
  results.corr = matrix(data = NA, nrow = length(corr), ncol = 4)
  
  v = 1
  for (c in corr){
    
    tau.vec = numeric()
    cdte.mat = matrix(data = NA, nrow = length(x), ncol = 3)
    
    i = 1
    for (k in 1:length(x)){
      
      shift.1 = runif(n = 1, 0, shift)
      shift.0 = runif(n = 1, 0, shift)
      y.1 = rnorm(n = n.sample, mean = mu.1 + shift.1, sd = sigma.1)
      y.0 = rnorm(n = n.sample, mean = mu.0 + shift.0, sd = sigma.0)
      
      delta = rnorm(n = n.sample, mean = mu.1+shift.1-(mu.0+shift.0),
                    sd = sqrt(sigma.1^2+sigma.0^2-2*c*sigma.1*sigma.0))
      
      tau.vec[k] = mean(y.1)-mean(y.0)
      
      cdte.mat[i,1] = max(ES(y.1, p_loss = alpha)-ES(y.0, p_loss = alpha),
                          ES(-1*y.0, p_loss = alpha)-ES(-1*y.1, p_loss = alpha))
      cdte.mat[i,2] = ES(delta, p_loss = alpha)
      cdte.mat[i,3] = ES(y.1, p_loss = alpha) + ES(-1*y.0, p_loss = alpha)
      i = i + 1
    }
    
    cvar.kall = ES(tau.vec, p_loss = alpha)
    cvar.us = colMeans(cdte.mat)
    
    results.corr[v,1] = cvar.us[1]
    results.corr[v,2] = cvar.us[2]
    results.corr[v,3] = cvar.us[3]
    results.corr[v,4] = cvar.kall
    v = v + 1
  }
  
  
  corr.df = data.frame(results.corr) %>%
    rename("LB" = "X1",
           "True" = "X2",
           "UB" = "X3",
           "Kallus" = "X4") %>%
    mutate(corr = corr) %>%
    pivot_longer(cols = c("LB", "True", "UB", "Kallus")) %>%
    arrange(name, corr)
  return(corr.df)
}


simulate.shift.dynamic = function(shift, sigma.1, sigma.0,
                                  n.obs, mu.1, mu.0, alpha){
  
  # shift: upper limit of uniform r.v. (controls mean shifts)
  # sigma.1: std of Y_1
  # sigma.0: std of Y_0
  # n.obs: number of observations
  # mu.1: mean Y_1
  # mu.0: mean Y_0
  # alpha: alpha of CVaR
  
  set.seed(42)
  
  x = seq(1,1000,1)
  corr = seq(-1, 1, l = 20)
  shift = seq(0, shift, l = 20)
  
  alpha = 0.05
  
  results.corr.shift = matrix(data = NA, nrow = length(shift)*length(corr), ncol = 6)
  
  v = 1
  for (c in corr){
    for(k in shift){
      
      tau.vec = numeric()
      cdte.mat = matrix(data = NA, nrow = length(x), ncol = 3)
      
      i = 1
      for (j in 1:length(x)){
        
        shift.1 = runif(n = 1, 0, k)
        shift.0 = runif(n = 1, 0, k)
        y.1 = rnorm(n = n.obs, mean = mu.1+shift.1, sd = sigma.1)
        y.0 = rnorm(n = n.obs, mean = mu.0+shift.0, sd = sigma.0)
        
        delta = rnorm(n = n.obs, mean = mu.1+shift.1-(mu.0+shift.0),
                      sd = sqrt(sigma.1^2+sigma.0^2-2*c*sigma.1*sigma.0))
        
        tau.vec[j] = mean(y.1)-mean(y.0)
        
        cdte.mat[i,1] = max(ES(y.1, p_loss = alpha)-ES(y.0, p_loss = alpha),
                            ES(-1*y.0, p_loss = alpha)-ES(-1*y.1, p_loss = alpha))
        cdte.mat[i,2] = ES(delta, p_loss = alpha)
        cdte.mat[i,3] = ES(y.1, p_loss = alpha) + ES(-1*y.0, p_loss = alpha)
        i = i + 1
      }
      
      cvar.kall = ES(tau.vec, p_loss = alpha)
      cvar.us = colMeans(cdte.mat)
      
      results.corr.shift[v,1] = cvar.us[1]
      results.corr.shift[v,2] = cvar.us[2]
      results.corr.shift[v,3] = cvar.us[3]
      results.corr.shift[v,4] = cvar.kall
      results.corr.shift[v,5] = c
      results.corr.shift[v,6] = k
      v = v + 1
    }
  }
  
  corr.shift.df = data.frame(results.corr.shift) %>%
    rename("LB" = "X1",
           "True" = "X2",
           "UB" = "X3",
           "Kallus" = "X4",
           "Corr" = "X5",
           "Shift" = "X6") %>%
    mutate(Shift = round(Shift,2)) %>%
    pivot_longer(cols = c("LB", "True", "UB", "Kallus")) %>%
    arrange(name, value)
  
  return(corr.shift.df)
}


plot.results.full = function(df){
  
  # df: data set of results
  
  plot = ggplot(data = df, aes(x = var.delta.x.group, y = value, group = name, color = name)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    geom_point(size = 1.5) +
    theme_bw() +
    scale_color_jama() +
    theme(legend.position = "top",
          text = element_text(size = 12)) +
    labs(color = "", group = "", x = expression(max(rho) ~ and ~ min(sigma) ~ symbol('\254') ~ Var(delta) ~ given ~ X ~ symbol('\256') ~ min(rho) ~ and ~ max(sigma)),
         y = expression(widehat(CVaR[alpha](delta)))) 
  return(plot)
}


plot.results.corr.static = function(df){
  
  # df: data set of results from static
  
  plot = ggplot(data = df, aes(x = alpha, y = value, group = name,
                                              color = name)) +
    geom_line(alpha = 0.6, position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
    geom_point(position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
    theme_bw() +
    scale_color_jama() +
    theme(legend.position = "top") +
    labs(color = "", group = "", x = expression(alpha),
         y = expression(widehat(CVaR[alpha](delta))))
  return(plot)
}


plot.results.corr.dynamic = function(df){
  
  # df: data set of results from dynamic
  
  plot = ggplot(data = df, aes(x = alpha, y = value, group = name,
                                 color = name)) +
    geom_line(alpha = 0.6, position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
    geom_point(position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
    theme_bw() +
    scale_color_manual(values = pal_jama("default")(4)[c(1,2,3,4)]) +
    theme(legend.position = "top") +
    labs(color = "", group = "", x = expression(alpha),
         y = expression(widehat(CVaR[alpha](delta)))) +
    facet_wrap(~corr)
  return(plot)
}


plot.results.corr.full = function(df){
  
  # df: data set of results from corr full
  
  plot = ggplot(data = df, aes(x = corr, y = value, group = name,
                               color = name)) +
    geom_line(alpha = 0.6, position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
    geom_point(position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
    theme_bw() +
    scale_color_jama() +
    theme(legend.position = "top") +
    labs(color = "", group = "", x = expression(rho),
         y = expression(widehat(CVaR[alpha](delta))))
  return(plot)
}


plot.results.shift.dynamic = function(df){
  
  # df: data set of results from shift dynamic
  
  plot = ggplot(data = df, aes(x = Corr, y = value, group = name,
                               color = name)) +
    geom_line(alpha = 0.6, position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
    geom_point(position = position_jitter(width = 0.01, height = 0.01, seed = 1),
               size = 0.85) +
    theme_bw() +
    scale_color_jama() +
    theme(legend.position = "top") +
    labs(color = "", group = "", x = expression(rho),
         y = expression(widehat(CVaR[alpha](delta)))) +
    facet_wrap(~Shift)
  return(plot)
}
