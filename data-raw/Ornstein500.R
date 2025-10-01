## Code to generate Ornstein-Uhlenbeck process with 5000 rows data:
# dx = theta * (mu - x) * dt + sigma * dw
# y = x + e
## The mean value mu is driven by a random input

############################################################
# Data simulation
############################################################

# Simulate data using Euler Maruyama
set.seed(20)
pars = c(theta=10, mu=1, sigma_x=1, sigma_y=0.1)
# 
dt.sim = 1e-3
t.sim = seq(0,5,by=dt.sim)
dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
u.sim = cumsum(rnorm(length(t.sim),sd=0.05))
x = 3
for(i in 1:(length(t.sim)-1)) {
  x[i+1] = x[i] + pars[1]*(pars[2]-x[i]+u.sim[i])*dt.sim + pars[3]*dw[i]
}

# Extract observations and add noise
dt.obs = 1e-2
ids = seq(1,length(t.sim),by=round(dt.obs / dt.sim))
t.obs = t.sim[ids]
y = x[ids] + pars[4] * rnorm(length(t.obs))
# forcing input
u = u.sim[ids]

# create data
Ornstein500 = data.frame(
  t = t.obs,
  y = y,
  u = u
)

# plot(Ornstein$t,Ornstein$y,type="l")

# uncomment below to update
usethis::use_data(Ornstein500, overwrite = TRUE)
