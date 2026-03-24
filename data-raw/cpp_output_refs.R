## Generate reference outputs for C++ filter / predict / simulate regression tests.
##
## Uses the 2D Ornstein_augmented model at its initial parameter values
## (no estimation) so the reference is fully deterministic and optimizer-independent.
## Simulate uses set.seed(123) for reproducibility.

library(ctsmTMB)

# -----------------------------------------------------------------------
# Input data: first 10 rows of Ornstein_augmented with a few NA values
# -----------------------------------------------------------------------

data(Ornstein_augmented)
Ornstein_augmented_NA        = Ornstein_augmented[1:10, ]
Ornstein_augmented_NA[3, "y1"] = NA   # missing first observation
Ornstein_augmented_NA[5, "y2"] = NA   # missing second observation
Ornstein_augmented_NA[7, "y1"] = NA
Ornstein_augmented_NA[7, "y2"] = NA   # both observations missing at row 7

# -----------------------------------------------------------------------
# 2D model (identical to the one used in test-estimate.R / test-filter.R)
# -----------------------------------------------------------------------

model2 = ctsmTMB$new()
model2$addSystem(
  dx1 ~ theta * (mu + u - x1) * dt + sigma_x1 * dw1,
  dx2 ~ alpha * (x1 - x2)    * dt + sigma_x2 * dw2
)
model2$addObs(y1 ~ x1, y2 ~ x2)
model2$setVariance(y1 ~ sigma_y1^2, y2 ~ sigma_y2^2)
model2$setAlgebraics(
  theta    ~ exp(logtheta),
  alpha    ~ exp(logalpha),
  sigma_x1 ~ exp(logsigma_x1),
  sigma_x2 ~ exp(logsigma_x2),
  sigma_y1 ~ exp(logsigma_y1),
  sigma_y2 ~ exp(logsigma_y2)
)
model2$addInput(u)
model2$setParameter(
  logtheta    = log(c(initial = 5,   lower = 1e-5,  upper = 50)),
  mu          =     c(initial = 3,   lower = -10,   upper = 10),
  logalpha    = log(c(initial = 2,   lower = 1e-5,  upper = 50)),
  logsigma_x1 = log(c(initial = 1,   lower = 1e-10, upper = 5)),
  logsigma_x2 = log(c(initial = 0.5, lower = 1e-10, upper = 5)),
  logsigma_y1 = log(c(initial = 0.1, lower = 1e-10, upper = 2)),
  logsigma_y2 = log(c(initial = 0.1, lower = 1e-10, upper = 2))
)
model2$setInitialState(list(c(3, 3), diag(0.1, 2)))

# -----------------------------------------------------------------------
# Reference outputs (initial parameters, no estimation)
# -----------------------------------------------------------------------

ref_filter_cpp   = model2$filter(Ornstein_augmented_NA,
                                  method  = "ekf",
                                  use.cpp = TRUE,
                                  silent  = TRUE)

ref_predict_cpp  = model2$predict(Ornstein_augmented_NA,
                                   method  = "ekf",
                                   use.cpp = TRUE,
                                   silent  = TRUE)

ref_simulate_cpp = model2$simulate(Ornstein_augmented_NA,
                                    method    = "ekf",
                                    use.cpp   = TRUE,
                                    cpp.seeds = c(123, 456),
                                    silent    = TRUE)

# -----------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------

usethis::use_data(Ornstein_augmented_NA, overwrite = TRUE)
usethis::use_data(ref_filter_cpp,        overwrite = TRUE)
usethis::use_data(ref_predict_cpp,       overwrite = TRUE)
usethis::use_data(ref_simulate_cpp,      overwrite = TRUE)
