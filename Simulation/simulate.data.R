rm(list=ls())

library(dplyr)
library(cvar)
library(tidyr)
library(ggpubr)
library(ggsci)

setwd("C:/Users/eminu/OneDrive/Desktop/Treatment-Effect-Risk")


# create plots illustrating our bounds. Assumption: Y_1, Y_0 jointly Gaussian.
set.seed(42)
n.sample = 1000000
mu.1 = 0
mu.0 = 0

sd.1 = 1
sd.0 = 1

rho = -1

y.1 = rnorm(n = n.sample, mean = mu.1, sd = sd.1)
y.0 = rnorm(n = n.sample, mean = mu.0, sd = sd.0)

delta = rnorm(n = n.sample, mean = mu.1-mu.0,
              sd = sqrt(sd.1^2+sd.0^2-2*rho*sd.1*sd.0))

alpha = seq(0, 1, 0.05)

results = matrix(data = NA, nrow = length(alpha), ncol = 3)

i = 1
for(j in alpha){
  results[i,1] = min(ES(y.1, p_loss = j)-ES(y.0, p_loss = j),
                     ES(y.0, p_loss = j)-ES(y.1, p_loss = j))
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


plot.bounds = ggplot(data = results.df, aes(x = alpha, y = value, group = name,
                                            color = name)) +
  geom_line(alpha = 0.6, position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
  theme_bw() +
  scale_color_jama() +
  theme(legend.position = "top") +
  labs(color = "", group = "", x = expression(alpha),
       y = expression(CVaR[alpha](delta)))
plot.bounds
ggsave("./Plots/bounds_neg_corr.png", plot.bounds, width = 15, height = 15, units = "cm")



# compare with Kallus. Do it for varying degrees of sigma and for the uncorrelated case,
# perf. negatively correlated and perf. positively. Here tau(x) is the same for all x, thus the CVaR 
# with Kallus method is zero.
x = seq(1,500,1)

n.sample = 10000
mu.1 = 1
mu.0 = 1

sd = seq(0,2, 0.1)
rho = 0
alpha = 0.05

results.sigma = matrix(data = NA, nrow = length(sd), ncol = 4)

v = 1
for (p in sd){
  
  tau.vec = numeric()
  cdte.mat = matrix(data = NA, nrow = length(x), ncol = 3)
  
  i = 1
  for (k in 1:length(x)){
    
    y.1 = rnorm(n = n.sample, mean = mu.1, sd = p)
    y.0 = rnorm(n = n.sample, mean = mu.0, sd = p)
    
    delta = rnorm(n = n.sample, mean = mu.1-mu.0,
                  sd = sqrt(p^2+p^2-2*rho*p*p))
    
    tau.vec[k] = mean(y.1)-mean(y.0)
    
    cdte.mat[i,1] = min(ES(y.1, p_loss = alpha)-ES(y.0, p_loss = alpha),
                        ES(y.0, p_loss = alpha)-ES(y.1, p_loss = alpha))
    cdte.mat[i,2] = ES(delta, p_loss = alpha)
    cdte.mat[i,3] = ES(y.1, p_loss = alpha) + ES(-1*y.0, p_loss = alpha)
    i = i + 1
  }
  
  cvar.kall = ES(tau.vec, p_loss = alpha)
  cvar.us = colMeans(cdte.mat)
  
  results.sigma[v,1] = cvar.us[1]
  results.sigma[v,2] = cvar.us[2]
  results.sigma[v,3] = cvar.us[3]
  results.sigma[v,4] = cvar.kall
  v = v + 1
}


sigma.df = data.frame(results.sigma) %>%
  rename("LB" = "X1",
         "True" = "X2",
         "UB" = "X3",
         "Kallus" = "X4") %>%
  mutate(sigma = sd) %>%
  pivot_longer(cols = c("LB", "True", "UB", "Kallus")) %>%
  arrange(name, sigma)
  
  
plot.kallus.sigma = ggplot(data = sigma.df %>% filter(sigma <= 1),
                     aes(x = sigma, y = value, group = name, color = name)) +
  geom_line(alpha = 0.6, linewidth = 0.5,
            position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.01, seed = 1),
             size = 1.5) +
  theme_bw() +
  scale_color_jama() +
  theme(legend.position = "top") +
  labs(color = "", group = "", x = expression(Var(delta) ~ given ~ X),
       y = expression(CVaR[alpha](delta)))
plot.kallus.sigma
ggsave("./Plots/kallus_no_corr.png", plot.kallus.sigma, width = 15, height = 15, units = "cm")


# vary the level of heterogeneity in tau(x)
x = seq(1,500,1)

n.sample = 10000
mu.1 = 1
mu.0 = 1

mu.shift = seq(0,2,0.1)

sd.0 = 0.5
sd.1 = 0.5
rho = 0
alpha = 0.05

results.shift = matrix(data = NA, nrow = length(mu.shift), ncol = 4)

v = 1
for (p in mu.shift){
  
  tau.vec = numeric()
  cdte.mat = matrix(data = NA, nrow = length(x), ncol = 3)
  
  i = 1
  for (k in 1:length(x)){
    
    shift.1 = runif(n = 1, 0, p)
    shift.0 = runif(n = 1, 0, p)
    y.1 = rnorm(n = n.sample, mean = mu.1 + shift.1, sd = sd.1)
    y.0 = rnorm(n = n.sample, mean = mu.0 + shift.0, sd = sd.0)
    
    delta = rnorm(n = n.sample, mean = mu.1+shift.1-(mu.0+shift.0),
                  sd = sqrt(sd.1^2+sd.0^2-2*rho*sd.1*sd.0))
    
    tau.vec[k] = mean(y.1)-mean(y.0)
    
    cdte.mat[i,1] = min(ES(y.1, p_loss = alpha)-ES(y.0, p_loss = alpha),
                        ES(y.0, p_loss = alpha)-ES(y.1, p_loss = alpha))
    cdte.mat[i,2] = ES(delta, p_loss = alpha)
    cdte.mat[i,3] = ES(y.1, p_loss = alpha) + ES(-1*y.0, p_loss = alpha)
    i = i + 1
  }
  
  cvar.kall = ES(tau.vec, p_loss = alpha)
  cvar.us = colMeans(cdte.mat)
  
  results.shift[v,1] = cvar.us[1]
  results.shift[v,2] = cvar.us[2]
  results.shift[v,3] = cvar.us[3]
  results.shift[v,4] = cvar.kall
  v = v + 1
}


hetero.df = data.frame(results.shift) %>%
  rename("LB" = "X1",
         "True" = "X2",
         "UB" = "X3",
         "Kallus" = "X4") %>%
  mutate(hetero = mu.shift) %>%
  pivot_longer(cols = c("LB", "True", "UB", "Kallus")) %>%
  arrange(name, hetero)


plot.kallus.shift = ggplot(data = hetero.df,
                     aes(x = hetero, y = value, group = name, color = name)) +
  geom_line(alpha = 0.6, linewidth = 0.5,
            position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.01, seed = 1),
             size = 1.5) +
  theme_bw() +
  scale_color_jama() +
  theme(legend.position = "top") +
  labs(color = "", group = "", x = "Heterogeneity in potential outcome distributions",
       y = expression(CVaR[alpha](delta)))
plot.kallus.shift
ggsave("./Plots/kallus_no_corr_hetero.png", plot.kallus.shift, width = 15, height = 15, units = "cm")



# compare with Kallus. Do it for varying degrees of sigma and for the uncorrelated case,
# perf. negatively correlated and perf. positively. Here we let the distribution of tau(x) no
# non-degenerate by varying the pot. outcome distributions.
n.sample = 10000
mu.1 = 1
mu.0 = 1
p.shift = 1

sd = seq(0,2, 0.1)
rho = 1
alpha = 0.05

results.sigma.shift = matrix(data = NA, nrow = length(sd), ncol = 4)

v = 1
for (p in sd){
  
  tau.vec = numeric()
  cdte.mat = matrix(data = NA, nrow = length(x), ncol = 3)
  
  i = 1
  for (k in 1:length(x)){
    
    shift.1 = runif(n = 1, 0, p.shift)
    shift.0 = runif(n = 1, 0, p.shift)
    y.1 = rnorm(n = n.sample, mean = mu.1+shift.1, sd = p)
    y.0 = rnorm(n = n.sample, mean = mu.0+shift.0, sd = p)
    
    delta = rnorm(n = n.sample, mean = mu.1+shift.1-(mu.0+shift.0),
                  sd = sqrt(p^2+p^2-2*rho*p*p))
    
    tau.vec[k] = mean(y.1)-mean(y.0)
    
    cdte.mat[i,1] = min(ES(y.1, p_loss = alpha)-ES(y.0, p_loss = alpha),
                        ES(y.0, p_loss = alpha)-ES(y.1, p_loss = alpha))
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
  v = v + 1
}


sigma.shift.df = data.frame(results.sigma.shift) %>%
  rename("LB" = "X1",
         "True" = "X2",
         "UB" = "X3",
         "Kallus" = "X4") %>%
  mutate(sigma = sd) %>%
  pivot_longer(cols = c("LB", "True", "UB", "Kallus")) %>%
  arrange(name, sigma)


plot.kallus.sigma.shift = ggplot(data = sigma.shift.df %>% filter(sigma <= 1),
                           aes(x = sigma, y = value, group = name, color = name)) +
  geom_line(alpha = 0.6, linewidth = 0.5,
            position = position_jitter(width = 0.01, height = 0.01, seed = 1)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.01, seed = 1),
             size = 1.5) +
  theme_bw() +
  scale_color_jama() +
  theme(legend.position = "top") +
  labs(color = "", group = "", x = expression(Var(delta) ~ given ~ X),
       y = expression(CVaR[alpha](delta)))
plot.kallus.sigma.shift
ggsave("./Plots/kallus_pos_corr_shift.png", plot.kallus.sigma.shift, width = 15, height = 15, units = "cm")

