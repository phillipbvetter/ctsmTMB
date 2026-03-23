## code to prepare `Ornstein` dataset goes here:

# Create model object
model <- newModel()
model$addSystem(dx ~ theta * (mu-x+u) * dt + sigma_x*dw)
model$addObs(y ~ x)
model$setVariance(y ~ sigma_y^2)
model$addInput(u)
model$setParameter(
  theta   = c(initial=2, lower=1e-5, upper=50),
  mu      = c(initial=2, lower=0, upper=5),
  sigma_x = c(initial=1e-2, lower=1e-10, upper=30),
  sigma_y = c(initial=1e-2, lower=1e-10, upper=30)
)
model$setInitialState(list(x0=3, p0=0.01 * diag(1)))

# Create initial data
pars = c(theta=5, mu=3, sigma_x=1, sigma_y=0.1)
dt.sim = 1e-3
t.sim = seq(0,20,by=dt.sim)
dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
u.sim = cumsum(rnorm(length(t.sim),sd=0.05))
df.sim <- data.frame(
  t = t.sim,
  u = u.sim,
  y = NA 
)

# Simulate
cpp.seeds <- c(20,20)
sim <- model$simulate(data=df.sim,
                      pars = pars,
                      method="ekf", 
                      n.sims = 1,
                      cpp.seeds = cpp.seeds)
y.sim <- sim$observations$y$i0[,1]

# Extract observed data
dt.obs = 1e-1
ids = seq(1,length(t.sim),by=round(dt.obs / dt.sim))
# Create data
Ornstein = data.frame(
  t = t.sim[ids],
  y = y.sim[ids],
  u = u.sim[ids]
)


# plot(Ornstein$t,Ornstein$y,type="l")

# uncomment below to update
# usethis::proj_set("~/github/ctsmTMB")
# usethis::use_data(Ornstein, overwrite = TRUE)
