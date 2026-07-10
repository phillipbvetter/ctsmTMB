## Generate reference outputs for estimate() regression tests.
##
## Uses the 2D Ornstein_augmented_NA dataset (10 rows, some NAs).
## One reference object per method; the 'private' environment is excluded.

library(ctsmTMB)
data(Ornstein2D)
how.many.rows <- 50
df <- Ornstein2D[1:how.many.rows,]

# Helper: strip the non-serialisable private environment from a fit object
strip_fit = function(fit) {
  out = fit[setdiff(names(fit), "private")]
  class(out) = class(fit)
  out
}

methods = c("ekf", "lkf", "ukf", "laplace", "laplace.thygesen")

EstimateReferenceData = list()

for (m in methods) {
  tempmodel <- create.Ornstein2D.model()
  fit   = tempmodel$estimate(df,
                             method = m,
                             ode.solver="rk4",
                             first.order.input.hold = FALSE,
                             silent = TRUE,
                             trace = 0)
  ref      = strip_fit(fit)
  EstimateReferenceData[[m]] <- ref
}

for (m in methods) {
  tempmodel <- create.Ornstein2D.model()
  fit   = tempmodel$estimate(df,
                             method = m,
                             ode.solver="rk4",
                             first.order.input.hold = TRUE,
                             silent = TRUE,
                             trace = 0)
  ref      = strip_fit(fit)
  EstimateReferenceData[[paste0(m,"_foh")]] <- ref
}

usethis::proj_set("~/github/ctsmTMB")
usethis::use_data(EstimateReferenceData, overwrite = TRUE)
