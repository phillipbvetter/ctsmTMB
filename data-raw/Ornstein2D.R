## code to prepare `Ornstein2D` dataset goes here:

# -----------------------------------------------------------------------
# Description
# -----------------------------------------------------------------------

# Two-state, two-observation Ornstein-Uhlenbeck system.
#
# State equations:
#   dx1 = theta * (mu + u - x1) * dt + sigma_x1 * dw1
#   dx2 = alpha * (x1 - x2)    * dt + sigma_x2 * dw2
#
# Observation equations:
#   y1 ~ x1   (direct observation of x1)
#   y2 ~ x2   (direct observation of x2)

# -----------------------------------------------------------------------
# Create Model
# -----------------------------------------------------------------------
model <- create.Ornstein2D.model()

# -----------------------------------------------------------------------
# Create data for simulation
# -----------------------------------------------------------------------
set.seed(13)
dt.sim <- 1e-3
t.sim = seq(0, 20, by=dt.sim)
u.sim = cumsum(rnorm(length(t.sim), sd = 0.05))
df.sim <- data.frame(
  t = t.sim,
  u = u.sim,
  y1 = NA,
  y2 = NA
)

# -----------------------------------------------------------------------
# Simulate
# -----------------------------------------------------------------------
sim.out <- model$simulate(data=df.sim,
                          cpp.seeds = c(1,10),
                          n.sim=1)

# -----------------------------------------------------------------------
# Create final data
# -----------------------------------------------------------------------
# Extract observed data
dt.obs <- 1e-1
ids <- seq(1, length(t.sim), by=round(dt.obs / dt.sim))

# Create data
Ornstein2D <- data.frame(
  t = t.sim[ids],
  u = u.sim[ids],
  y1 = sim.out$observations$y1$i0[ids,1],
  y2 = sim.out$observations$y2$i0[ids,1]
)

ggplot(data=Ornstein2D) +
  geom_line(aes(x=t,y=y1)) +
  geom_line(aes(x=t,y=y2),color="red") +
  ctsmTMB:::getggplot2theme()

# uncomment below to update
usethis::proj_set("~/github/ctsmTMB")
usethis::use_data(Ornstein2D, overwrite = TRUE)
