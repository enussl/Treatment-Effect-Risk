rm(list = ls())

library(dplyr)
library(drf)
library(cvar)
library(tidyr)
library(ggplot2)


set.seed(1)
setwd("C:/Users/eminu/OneDrive/Desktop/Treatment-Effect-Risk")

data = read.csv("./Application/data.r.csv")


vars.demean = data %>% select(-c("X.Canton_3", "X.oneweek"))
for (k in 1:length(vars.demean)) {
  
  fit = lm(data[,k] ~ factor(data$X.oneweek) + factor(data$X.Canton_3))
  res = fit$residuals 
  data[,k] = fit$residuals
} 

write.csv(data, "./data_demeaned.csv", row.names = F)


X = data %>% select(-Y, -X.oneweek, -X.Canton_3, -W)
X.0 = data %>% filter(W == 0) %>% select(-Y, -X.oneweek, -X.Canton_3, -W) 
X.1 = data %>% filter(W == 1) %>% select(-Y, -X.oneweek, -X.Canton_3, -W) 

Y.0 = data %>% filter(W == 0) %>% select(Y) 
Y.1 = data %>% filter(W == 1) %>% select(Y)

X.t = data %>% 
  select(-X.oneweek, -X.Canton_3, -W, -Y) %>%
  summarise(across(everything(), mean))


drf.0 = drf(Y = Y.0, X = X.0)
drf.1 = drf(Y = Y.1, X = X.1)

mean.pred.0 = predict(drf.0, newdata = X, functional = "mean")
mean.pred.1 = predict(drf.1, newdata = X, functional = "mean")

cate.drf = mean.pred.1[["mean"]] - mean.pred.0[["mean"]]
hist(cate.drf)
mean(cate.drf)
cvar.cate.drf = ES(-1*cate.drf, p_loss = 0.25)


weights.0 = drf::get_sample_weights(drf.0, newdata = X)
weights.1 = drf::get_sample_weights(drf.1, newdata = X)

dist.0 = matrix(weights.0)*as.vector(Y.0)
dist.1 = sweep(as.matrix(weights.1), 1, as.vector(X.1), "*")


alpha.x = seq(0, 1, b = 0.05)
res = matrix(data = NA, nrow = length(alpha.x), ncol = 3)

k = 1
for (alpha in alpha.x) {
  
  n.samples = 5000
  samples.0 = matrix(data = NA, nrow = nrow(X), ncol = n.samples)
  cvar.0.min = matrix(data = NA, nrow = nrow(X), ncol = 1)
  cvar.0.plus = matrix(data = NA, nrow = nrow(X), ncol = 1)
  
  
  samples.1 = matrix(data = NA, nrow = nrow(X), ncol = n.samples)
  cvar.1.min = matrix(data = NA, nrow = nrow(X), ncol = 1)
  cvar.1.plus = matrix(data = NA, nrow = nrow(X), ncol = 1)
  
  
  for (i in 1:nrow(X)){
    samples.0[i,] = sample(as.numeric(unlist(Y.0)), size = n.samples, replace = T, prob = weights.0[i,])
    cvar.0.plus[i,1] = ES(-1*samples.0[i,], p_loss = alpha)
    cvar.0.min[i,1] = ES(samples.0[i,], p_loss = alpha)
    
    
    samples.1[i,] = sample(as.numeric(unlist(Y.1)), size = n.samples, replace = T, prob = weights.1[i,])
    cvar.1.plus[i,1] = ES(-1*samples.1[i,], p_loss = alpha)
    cvar.1.min[i,1] = ES(samples.1[i,], p_loss = alpha)
  }
  
  cvar.delta.lb = mean(cvar.1.plus) - mean(cvar.0.plus)
  cvar.delta.ub = mean(cvar.1.plus) + mean(cvar.0.plus)
  cvar.cate = ES(-1*cate.drf, p_loss = alpha)
  
  res[k, 1] = cvar.delta.lb
  res[k, 2] = cvar.delta.ub
  res[k, 3] = cvar.cate
  
  k = k + 1
}

res.df = data.frame(res) %>% rename("LB" = "X1", "UB" = "X2", "CATE" = "X3") %>%
  mutate(alpha = alpha.x) %>% 
  pivot_longer(cols = c("LB", "UB", "CATE")) 


plot.results = ggplot(data = res.df, aes(x = alpha, y = value, group = name, color = name)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "top") +
  labs(x = expression(alpha), y = expression(widehat(CVaR)[alpha](delta)),
       group = "", color = "")
plot.results
ggsave("./plots/results_mask.png", width = 10, height = 10, units = "cm")


# causal forest
cf = causal_forest(X = as.matrix(data %>% select(-Y, -X.oneweek, -X.Canton_3, -W)),
                   Y = as.matrix(data %>% select(Y)),
                   W = as.matrix(data %>% select(W)))

tau.hat.oob = predict(cf)
hist(tau.hat.oob$predictions, breaks = 20)
mean(tau.hat.oob$predictions)

cvar.cate.grf = cvar::ES(-1*cate, p_loss = 0.25)



