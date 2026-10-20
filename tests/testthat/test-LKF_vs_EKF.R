#### Setup ####
data(Ornstein2D)
how.many.rows <- 50
df <- Ornstein2D[1:how.many.rows,]

# Helper: strip the non-serialisable private environment from a fit object
strip_fit = function(fit) {
  out = fit[setdiff(names(fit), "private")]
  class(out) = class(fit)
  out
}

model <- create.Ornstein2D.model()

time.reduction <- 10
dt <- diff(df$t)[1]/time.reduction
tol <- 1e-3

testthat::test_that("Check if LKF is identical to EKF for small time-step",{
  # Author only test
  report.fields <- c("states","observations","residuals")
  testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - skipping due to minor numerical inconsistencies across operating systems and versions")
  fit1 <- model$estimate(df, method = "lkf", silent = TRUE, trace=0)
  fit2 <- model$estimate(df, method = "ekf", ode.timestep=dt, silent = TRUE, trace=0)
  expect_equal(fit1$nll, fit2$nll, tolerance=tol)
  expect_equal(fit1[report.fields], fit2[report.fields], tolerance=tol)
})

testthat::test_that("Check if LKF and EKF smoothers given similar results",{
  # Author only test
  testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - skipping due to minor numerical inconsistencies across operating systems and versions")
  smooth1 <- model$smooth(df, method = "lkf", silent = TRUE, trace=0)
  smooth2 <- model$smooth(df, method = "ekf", silent = TRUE, trace=0, ode.timestep=dt)
  expect_equal(smooth1, smooth2, tolerance=tol)
})

