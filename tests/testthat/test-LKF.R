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
report.fields <- c("states","observations","residuals")

time.reduction <- 10
dt <- diff(df$t)[1]/time.reduction

tol <- 1e-3

testthat::test_that("Check if LKF is immune to ode.timestep refinement (ZOH)",{
  # Author only test
  testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - skipping due to minor numerical inconsistencies across operating systems and versions")
  fit1 <- model$estimate(df, method = "lkf", silent = TRUE, trace=0)
  fit2 <- model$estimate(df, method = "lkf", ode.timestep = dt, silent = TRUE, trace=0)
  expect_equal(fit1$nll, fit2$nll, tolerance=tol)
  expect_equal(fit1[report.fields], fit2[report.fields], tolerance=tol)
})

testthat::test_that("Check if LKF is immune to ode.timestep refinement (FOH)",{
  # Author only test
  testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - skipping due to minor numerical inconsistencies across operating systems and versions")
  fit1 <- model$estimate(df, method = "lkf", first.order.input.hold = TRUE, silent = TRUE, trace=0)
  fit2 <- model$estimate(df, method = "lkf", ode.timestep = dt, first.order.input.hold = TRUE, silent = TRUE, trace=0)
  expect_equal(fit1$nll, fit2$nll,tolerance=tol)
  expect_equal(fit1[report.fields], fit2[report.fields], tolerance=tol)
})

