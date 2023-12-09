rm(list = ls())

library(dplyr)
library(tidyr)
library(cvar)
library(ggplot2)
library(ggsci)
library(gganimate)


setwd("C:/Users/eminu/OneDrive/Desktop/Treatment-Effect-Risk")
source("./Simulation/helperfunctions.R")


# First look at our bounds vs the truth depending on the correlation
data.corr = simulate.corr.dynamic(shift = 2, sigma.1 = 2, sigma.0 = 2, n.obs = 100000, mu.1 = 1, mu.0 = 1)
(plot.corr = plot.results.corr.dynamic(data.corr))
ggsave("./Plots/bounds_corr.png", plot.corr, width = 20, height = 20, units = "cm")


# Now compare our results with Kallus depending on Var(delta\mid X)
data.comp = simulate.full.data(shift = 2, sigma.1 = 2, sigma.0 = 2, n.obs = 10000, mu.1 = 1, mu.0 = 1, alpha = 0.05)
(plot.comp = plot.results.full(data.comp))
ggsave("./Plots/bounds_kallus_vardelta.png", plot.comp, width = 15, height = 15, units = "cm")


# Now compare our results with Kallus depending only on the correlation
data.comp.corr = simulate.corr.full(shift = 2, sigma.1 = 2, sigma.0 = 2, n.obs = 1000, mu.1 = 1, mu.0 = 1, alpha = 0.05)
(plot.comp.corr = plot.results.corr.full(data.comp.corr))
ggsave("./Plots/bounds_kallus_corr.png", plot.comp.corr, width = 15, height = 15, units = "cm")


# Now compare our results with Kallus depending only on the correlation and the level of mean-shift in Y_1,Y_0
data.comp.corr.shift = simulate.shift.dynamic(shift = 3, sigma.1 = 2, sigma.0 = 2, n.obs = 1000, mu.1 = 1, mu.0 = 1, alpha = 0.05)
(plot.comp.corr.shift = plot.results.shift.dynamic(data.comp.corr.shift))
ggsave("./Plots/bounds_kallus_corr_shift.png", plot.comp.corr.shift, width = 20, height = 20, units = "cm")


# Kallus is not always a lower bound. How?
