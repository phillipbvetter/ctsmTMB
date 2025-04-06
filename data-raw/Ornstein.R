## code to prepare `Ornstein` dataset goes here:

library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(ctsmTMB)

############################################################
# Data simulation
############################################################

# Simulate data using Euler Maruyama
set.seed(20)
pars = c(theta=5, mu=3, sigma_x=1, sigma_y=0.1)
# 
dt.sim = 1e-3
t.sim = seq(0,20,by=dt.sim)
dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
u.sim = cumsum(rnorm(length(t.sim),sd=0.05))
x = 3
for(i in 1:(length(t.sim)-1)) {
  x[i+1] = x[i] + pars[1]*(pars[2]-x[i]+u.sim[i])*dt.sim + pars[3]*dw[i]
}

# Extract observations and add noise
dt.obs = 1e-1
ids = seq(1,length(t.sim),by=round(dt.obs / dt.sim))
t.obs = t.sim[ids]
y = x[ids] + pars[4] * rnorm(length(t.obs))
# forcing input
u = u.sim[ids]

# Create data
.data = data.frame(
  t = t.obs,
  y = y,
  u = u
)

# create data
Ornstein = data.frame(
  t = t.obs,
  y = y,
  u = u
)

# plot(Ornstein$t,Ornstein$y,type="l")

usethis::use_data(Ornstein, overwrite = TRUE)
