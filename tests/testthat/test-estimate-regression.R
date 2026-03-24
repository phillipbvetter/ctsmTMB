#### Setup ####

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

# Fields to compare (excludes the non-serialisable 'private' environment)
compare_fields = c("convergence", "nll", "nll.gradient", "nll.hessian",
                   "par.fixed", "sd.fixed", "cov.fixed",
                   "tvalue", "Pr.tvalue")

check_estimate = function(method, ref) {
  fit = make_model2()$estimate(Ornstein_augmented_NA, method = method,
                               silent = TRUE, trace = 0)
  for (f in compare_fields) {
    testthat::expect_equal(fit[[f]], ref[[f]], label = paste0(method, "$", f))
  }
}

#### ekf ####

testthat::test_that("estimate (ekf) output matches reference", {
  data(ref_estimate_ekf)
  check_estimate("ekf", ref_estimate_ekf)
})

#### lkf ####

testthat::test_that("estimate (lkf) output matches reference", {
  data(ref_estimate_lkf)
  check_estimate("lkf", ref_estimate_lkf)
})

#### ukf ####

testthat::test_that("estimate (ukf) output matches reference", {
  data(ref_estimate_ukf)
  check_estimate("ukf", ref_estimate_ukf)
})

#### laplace ####

testthat::test_that("estimate (laplace) output matches reference", {
  data(ref_estimate_laplace)
  check_estimate("laplace", ref_estimate_laplace)
})

#### laplace.thygesen ####

testthat::test_that("estimate (laplace.thygesen) output matches reference", {
  data(ref_estimate_laplace_thygesen)
  check_estimate("laplace.thygesen", ref_estimate_laplace_thygesen)
})
