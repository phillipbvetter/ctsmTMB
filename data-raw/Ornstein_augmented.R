## code to prepare `Ornstein_augmented` dataset goes here:

############################################################
# Data simulation
############################################################

# Two-state, two-observation Ornstein-Uhlenbeck system.
#
# State equations:
#   dx1 = theta * (mu + u - x1) * dt + sigma_x1 * dw1
#   dx2 = alpha * (x1 - x2)    * dt + sigma_x2 * dw2
#
# Observation equations:
#   y1 ~ x1   (direct observation of x1)
#   y2 ~ x2   (direct observation of x2)
#
# x2 mean-reverts toward x1, so it acts as a lagged/smoothed
# version of x1. This makes the system identifiable and keeps
# the augmentation physically meaningful while being trivially
# simple to verify.

set.seed(42)
pars = c(theta=5, mu=3, alpha=2,
         sigma_x1=1, sigma_x2=0.5,
         sigma_y1=0.1, sigma_y2=0.1)

dt.sim = 1e-3
t.sim  = seq(0, 20, by=dt.sim)
n.sim  = length(t.sim)

dw1 = rnorm(n.sim - 1, sd = sqrt(dt.sim))
dw2 = rnorm(n.sim - 1, sd = sqrt(dt.sim))
u.sim = cumsum(rnorm(n.sim, sd = 0.05))

x1 = numeric(n.sim); x1[1] = 3
x2 = numeric(n.sim); x2[1] = 3

for (i in 1:(n.sim - 1)) {
  x1[i+1] = x1[i] + pars["theta"] * (pars["mu"] + u.sim[i] - x1[i]) * dt.sim +
             pars["sigma_x1"] * dw1[i]
  x2[i+1] = x2[i] + pars["alpha"] * (x1[i] - x2[i]) * dt.sim +
             pars["sigma_x2"] * dw2[i]
}

# Extract observations at coarser time-step and add noise
dt.obs = 1e-1
ids   = seq(1, n.sim, by = round(dt.obs / dt.sim))
t.obs = t.sim[ids]
y1    = x1[ids] + pars["sigma_y1"] * rnorm(length(ids))
y2    = x2[ids] + pars["sigma_y2"] * rnorm(length(ids))
u     = u.sim[ids]

Ornstein_augmented = data.frame(
  t  = t.obs,
  y1 = y1,
  y2 = y2,
  u  = u
)

# uncomment below to update
# usethis::use_data(Ornstein_augmented, overwrite = TRUE)
