#### Setup ####
data(Ornstein2D)
how.many.rows <- 50
df <- Ornstein2D[1:how.many.rows,]
data(EstimateReferenceData)
ed <- EstimateReferenceData

# Helper: strip the non-serialisable private environment from a fit object
strip_fit = function(fit) {
  out = fit[setdiff(names(fit), "private")]
  class(out) = class(fit)
  out
}

model <- create.Ornstein2D.model()
methods = c("ekf", "lkf", "ukf", "laplace", "laplace.thygesen")

testthat::test_that("Check if estimation matches references for zero order hold",{
  # Author only test
  testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - skipping due to minor numerical inconsistencies across operating systems and versions")

  # ZERO ORDER HOLD
  for (m in methods) {
    fit      = model$estimate(df, method = m, silent = TRUE, trace=0)
    testthat::expect_equal(strip_fit(fit), ed[[m]])
  }

  # FIRST ORDER HOLD
  for (m in methods) {
    fit      = model$estimate(df, method = m, silent = TRUE, trace=0, first.order.input.hold = TRUE)
    testthat::expect_equal(strip_fit(fit), ed[[paste0(m,"_foh")]])
  }
})

