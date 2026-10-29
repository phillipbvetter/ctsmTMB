# Data
# ------------------
data(Ornstein2D)
how.many.rows <- 50
df <- Ornstein2D[1:how.many.rows,]
# add missing data to check those kind of errors
df[1,c("y1","y2")] <- NA
df[2,"y1"] <- NA
df[3,"y2"] <- NA

# Load reference results
# ------------------
data(EstimateReferenceData)
ed <- EstimateReferenceData

tol <- 1e-6

# Helper: strip the non-serialisable private environment from a fit object
strip_fit = function(fit) {
  out = fit[setdiff(names(fit), "private")]
  class(out) = class(fit)
  out
}

model <- create.Ornstein2D.model()
# methods = c("ekf", "lkf", "ukf", "laplace", "laplace.thygesen")
methods = c("ekf", "lkf", "laplace", "laplace.thygesen")

testthat::test_that("Check estimate vs reference (ZOH)",{
  # Author only test
  testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - skipping due to minor numerical inconsistencies across operating systems and versions")

  # ZERO ORDER HOLD
  for (m in methods) {
    fit      = model$estimate(df, method = m, silent = TRUE, trace=0)
    testthat::expect_equal(strip_fit(fit), ed[[m]],tolerance=tol)
  }

})

testthat::test_that("Check estimate vs reference (FOH)",{
  # Author only test
  testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - skipping due to minor numerical inconsistencies across operating systems and versions")

  # FIRST ORDER HOLD
  for (m in methods) {
    fit      = model$estimate(df, method = m, silent = TRUE, trace=0, first.order.input.hold = TRUE)
    testthat::expect_equal(strip_fit(fit), ed[[paste0(m,"_foh")]], tolerance=tol)
  }
})

