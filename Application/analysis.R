rm(list = ls())

library(dplyr)
library(drf)
library(cvar)
library(tidyr)
library(ggplot2)
library(bayesmeta)


cvar = function(sample, alpha) {
  # sample: a vector of observations
  # alpha: level of cvar
  
  quantile = quantile(sample, probs = 1-alpha)
  tail.sample = sample[sample >= as.numeric(quantile)]
  return(1/length(tail.sample)*sum(tail.sample))
}

var.custom = function(sample, alpha) {
  # sample: a vector of observations
  # alpha: level of cvar
  
  quantile = quantile(sample, probs = 1-alpha)
  tail.sample = sample[sample >= as.numeric(quantile)]
  return(var(tail.sample))
}


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


alpha.x = seq(0, 1, b = 0.05)
res = matrix(data = NA, nrow = length(alpha.x), ncol = 6)

k = 1
for (alpha in alpha.x) {
  
  n.samples = 5000
  samples.0 = matrix(data = NA, nrow = nrow(X), ncol = n.samples)
  cvar.0.min = matrix(data = NA, nrow = nrow(X), ncol = 1)
  cvar.0.plus = matrix(data = NA, nrow = nrow(X), ncol = 1)
  var.0 = matrix(data = NA, nrow = nrow(X), ncol = 1)
  
  
  samples.1 = matrix(data = NA, nrow = nrow(X), ncol = n.samples)
  cvar.1.min = matrix(data = NA, nrow = nrow(X), ncol = 1)
  cvar.1.plus = matrix(data = NA, nrow = nrow(X), ncol = 1)
  var.1 = matrix(data = NA, nrow = nrow(X), ncol = 1)
  
  cvar.ind = matrix(data = NA, nrow = nrow(X), ncol = 1)
  
  
  for (i in 1:nrow(X)){
    #upper and lower bounds
    samples.0[i,] = sample(as.numeric(unlist(Y.0)), size = n.samples, replace = T, prob = weights.0[i,])
    cvar.0.plus[i,1] = cvar(samples.0[i,], alpha)
    cvar.0.min[i,1] = cvar(-1*samples.0[i,], alpha)
    var.0[i,1] = var.custom(samples.0[i,], alpha)
    
    
    samples.1[i,] = sample(as.numeric(unlist(Y.1)), size = n.samples, replace = T, prob = weights.1[i,])
    cvar.1.plus[i,1] = cvar(samples.1[i,], alpha)
    cvar.1.min[i,1] = cvar(-1*samples.1[i,], alpha)
    var.1[i,1] = var.custom(samples.1[i,], alpha)
    
    #independent case
    mean.1 = mean(samples.1[i,])
    mean.0 = mean(samples.0[i,])
    ecdf.1 = ecdf(samples.1[i,])
    ecdf.0 = ecdf(samples.0[i,])
    
    z = seq(min(min(samples.0[i,]), min(samples.1[i,])), max(max(samples.0[i,]), max(samples.1[i,])), length.out = 3000)
    convolution_result = stats::convolve(x = ecdf.1(z),
                                         y = rev(ecdf.0(z)), type = "open")
    normalized_cdf = cumsum(convolution_result)/sum(convolution_result)
    num_samples = 5000
    random_samples = quantile(runif(num_samples), normalized_cdf) -0.23*(mean.1+mean.0)
    cvar.ind[i,1] = cvar(random_samples, alpha)
  }
  
  cvar.delta.lb = max(mean(cvar.1.plus) - mean(cvar.0.plus), mean(cvar.0.min) - mean(cvar.1.min))
  cvar.delta.ub = mean(cvar.1.plus) + mean(cvar.0.plus)
  cvar.ind = mean(cvar.ind)
  cvar.cate = cvar(cate.drf, alpha)
  var.0.fin = mean(var.0)
  var.1.fin = mean(var.1)
  
  
  res[k, 1] = cvar.delta.lb
  res[k, 2] = cvar.delta.ub
  res[k, 3] = cvar.cate
  res[k, 4] = cvar.ind
  res[k, 5] = var.0.fin
  res[k, 6] = var.1.fin
  
  k = k + 1
}

res.df = data.frame(res) %>%
  select(X1, X2, X3, X4) %>%
  rename("LB" = "X1", "UB" = "X2", "Kallus" = "X3", "IND" = "X4") %>%
  mutate(alpha = alpha.x) %>% 
  pivot_longer(cols = c("LB", "UB", "Kallus", "IND")) 


plot.results = ggplot(data = res.df, aes(x = alpha, y = value, group = name, color = name)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor.y = element_blank()) +
  labs(x = expression(alpha), y = expression(widehat(CVaR)[alpha](delta)),
       group = "", color = "") +
  geom_hline(yintercept = as.numeric(res.df[83,3]), linetype = "dashed", size = 0.75) +
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.75, linetype = "dashed") +
  scale_y_continuous(breaks = c(round(as.numeric(res.df[83,3]),2), 0, 1, 2, 4, 6))
plot.results
ggsave("./plots/results_mask.pdf", width = 25, height = 10, units = "cm")


# causal forest
cf = causal_forest(X = as.matrix(data %>% select(-Y, -X.oneweek, -X.Canton_3, -W)),
                   Y = as.matrix(data %>% select(Y)),
                   W = as.matrix(data %>% select(W)))

tau.hat.oob = predict(cf)
hist(tau.hat.oob$predictions, breaks = 20)
mean(tau.hat.oob$predictions)

cvar.cate.grf = cvar::ES(-1*cate, p_loss = 0.25)



