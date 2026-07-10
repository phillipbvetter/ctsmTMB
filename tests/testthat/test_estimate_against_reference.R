#### Setup ####
library(ctsmTMB)
data(Ornstein2D)
how.many.rows <- 50
df <- Ornstein2D[1:how.many.rows,]
data(EstimateReferenceData)

# Helper: strip the non-serialisable private environment from a fit object
strip_fit = function(fit) {
  out = fit[setdiff(names(fit), "private")]
  class(out) = class(fit)
  out
}

methods = c("ekf", "lkf", "ukf", "laplace", "laplace.thygesen")


# ZERO ORDER HOLD
for (m in methods) {
  tempmodel <- create.Ornstein2D.model()
  fit      = tempmodel$estimate(df,
                                method = m,
                                silent = TRUE,
                                ode.solver="rk4",
                                first.order.input.hold = FALSE,
                                trace = 0)
  fit <- strip_fit(fit)
  testthat::test_that("Estimate output matches reference",{
    testthat::expect_equal(fit, EstimateReferenceData[[m]], tolerance=1e-5)
  })
}

for (m in methods) {
  tempmodel <- create.Ornstein2D.model()
  fit      = tempmodel$estimate(df,
                                method = m,
                                silent = TRUE,
                                ode.solver="rk4",
                                first.order.input.hold = TRUE,
                                trace = 0)
  fit <- strip_fit(fit)
  testthat::test_that("Estimate output matches reference",{
    testthat::expect_equal(fit, EstimateReferenceData[[paste0(m,"_foh")]],
                           tolerance=1e-5)
  })
}

