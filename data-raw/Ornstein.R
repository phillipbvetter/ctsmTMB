## code to prepare `Ornstein` dataset goes here:

# One-state one-observation Ornstein-Uhlenbeck system.
#
# State equations:
#   dx = theta * (mu + u - x) * dt + sigma_x * dw
#
# Observation equations:
#   y ~ x

# Create model object
# model <- newModel()
# model$addSystem(dx ~ theta * (mu * u - x) * dt + sigma_x * dw)
# model$addObs(y ~ x)
# model$setVariance(y ~ sigma_y^2)
# model$setAlgebraics(
#   theta   ~ exp(logtheta),
#   sigma_x ~ exp(logsigma_x),
#   sigma_y ~ exp(logsigma_y)
# )
# model$addInput("u")
# model$setParameter(
#   logtheta   = log(c(initial = 1,    lower = 1e-5,  upper = 50)),
#   mu         =     c(initial = 1.5,  lower = 0,     upper = 5),
#   logsigma_x = log(c(initial = 1,    lower = 1e-10, upper = 2)),
#   logsigma_y = log(c(initial = 0.1,  lower = 1e-10, upper = 2))
# )
# model$setInitialState(list(1, 0.1))
model <- create.Ornstein1D.model()

# Create initial data
set.seed(5)
dt.sim <- 1e-3
t.sim <- seq(0, 20, by=dt.sim)
dw <- rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
u.sim <- cumsum(rnorm(length(t.sim),sd=0.05))
df.sim <- data.frame(
  t = t.sim,
  u = u.sim,
  y = NA
)

# Simulate
cpp.seeds <- c(20,20)
sim <- model$simulate(data=df.sim,
                      method="ekf",
                      n.sims = 1,
                      cpp.seeds = cpp.seeds)

# Extract observed data
dt.obs <- 1e-1
ids <- seq(1,length(t.sim),by=round(dt.obs / dt.sim))

# Create data
Ornstein <- data.frame(
  t = t.sim[ids],
  u = u.sim[ids],
  y = sim$observations$y$i0[,1][ids]
)

ggplot(data=Ornstein) +
  geom_line(aes(x=t,y=y)) +
  ctsmTMB:::getggplot2theme()

# uncomment below to update
usethis::proj_set("~/github/ctsmTMB")
usethis::use_data(Ornstein, overwrite = TRUE)
