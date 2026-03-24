## Generate reference outputs for estimate() regression tests.
##
## Uses the 2D Ornstein_augmented_NA dataset (10 rows, some NAs).
## One reference object per method; the 'private' environment is excluded.

library(ctsmTMB)

data(Ornstein_augmented_NA)

# Build a fresh model each time to avoid accumulated fit state affecting the
# starting parameters of subsequent estimates.
make_model2 = function() {
  m = ctsmTMB$new()
  m$addSystem(
    dx1 ~ theta * (mu + u - x1) * dt + sigma_x1 * dw1,
    dx2 ~ alpha * (x1 - x2)    * dt + sigma_x2 * dw2
  )
  m$addObs(y1 ~ x1, y2 ~ x2)
  m$setVariance(y1 ~ sigma_y1^2, y2 ~ sigma_y2^2)
  m$setAlgebraics(
    theta    ~ exp(logtheta),
    alpha    ~ exp(logalpha),
    sigma_x1 ~ exp(logsigma_x1),
    sigma_x2 ~ exp(logsigma_x2),
    sigma_y1 ~ exp(logsigma_y1),
    sigma_y2 ~ exp(logsigma_y2)
  )
  m$addInput(u)
  m$setParameter(
    logtheta    = log(c(initial = 5,   lower = 1e-5,  upper = 50)),
    mu          =     c(initial = 3,   lower = -10,   upper = 10),
    logalpha    = log(c(initial = 2,   lower = 1e-5,  upper = 50)),
    logsigma_x1 = log(c(initial = 1,   lower = 1e-10, upper = 5)),
    logsigma_x2 = log(c(initial = 0.5, lower = 1e-10, upper = 5)),
    logsigma_y1 = log(c(initial = 0.1, lower = 1e-10, upper = 2)),
    logsigma_y2 = log(c(initial = 0.1, lower = 1e-10, upper = 2))
  )
  m$setInitialState(list(c(3, 3), diag(0.1, 2)))
  m
}

# Helper: strip the non-serialisable private environment from a fit object
strip_fit = function(fit) {
  out = fit[setdiff(names(fit), "private")]
  class(out) = class(fit)
  out
}

methods = c("ekf", "lkf", "ukf", "laplace", "laplace.thygesen")

for (m in methods) {
  fit      = make_model2()$estimate(Ornstein_augmented_NA, method = m, silent = TRUE, trace = 0)
  ref      = strip_fit(fit)
  obj_name = paste0("ref_estimate_", gsub("\\.", "_", m))
  assign(obj_name, ref)
  save(list = obj_name, file = file.path("data", paste0(obj_name, ".rda")),
       envir = environment(), compress = "bzip2")
  message("Saved ", obj_name)
}
