rm(list=ls())

library(dplyr)
library(cvar)


mu.x = 2
sd.x = 1
n.obs = 1000
x = rnorm(n = n.obs, mean = mu.x, sd = sd.x)

mu.delta = 1
sd.delta = 1


samples = list()
for(j in 1:length(x)) {
  delta.given.x = rnorm(n = n.obs, mean = mu.delta*x[j], sd = sd.delta)
  samples[[j]] = delta.given.x
}

res = do.call(rbind, samples)

cvar.delta.x = apply(-1*res, 1, cvar::ES)
estim.our = mean(cvar.delta.x)

tau = apply(res, 1, mean)
estim.kall = mean(tau)









