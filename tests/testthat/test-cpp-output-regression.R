#### Setup ####
# Uses the 2D Ornstein_augmented model at initial parameter values
# (no estimation) so results are fully deterministic.
# Ornstein_augmented_NA has 10 rows with NA values at:
#   row 3 -> y1, row 5 -> y2, row 7 -> y1 and y2

data(Ornstein_augmented_NA)

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

#### filter (use.cpp = TRUE) ####

testthat::test_that("filter (use.cpp=TRUE) output matches reference", {
  data(ref_filter_cpp)
  out = model2$filter(Ornstein_augmented_NA, method = "ekf",
                      use.cpp = TRUE, silent = TRUE)

  testthat::expect_equal(out$states$mean$prior,       ref_filter_cpp$states$mean$prior)
  testthat::expect_equal(out$states$mean$posterior,   ref_filter_cpp$states$mean$posterior)
  testthat::expect_equal(out$states$sd$prior,         ref_filter_cpp$states$sd$prior)
  testthat::expect_equal(out$states$sd$posterior,     ref_filter_cpp$states$sd$posterior)
  testthat::expect_equal(out$states$cov$prior,        ref_filter_cpp$states$cov$prior)
  testthat::expect_equal(out$states$cov$posterior,    ref_filter_cpp$states$cov$posterior)
  testthat::expect_equal(out$observations$mean$prior,     ref_filter_cpp$observations$mean$prior)
  testthat::expect_equal(out$observations$mean$posterior, ref_filter_cpp$observations$mean$posterior)
  testthat::expect_equal(out$residuals,               ref_filter_cpp$residuals)
})

#### predict (use.cpp = TRUE) ####

testthat::test_that("predict (use.cpp=TRUE) output matches reference", {
  data(ref_predict_cpp)
  out = model2$predict(Ornstein_augmented_NA, method = "ekf",
                       use.cpp = TRUE, silent = TRUE)

  testthat::expect_equal(out$states,       ref_predict_cpp$states)
  testthat::expect_equal(out$observations, ref_predict_cpp$observations)
})

#### simulate (use.cpp = TRUE) ####
# cpp.seeds = c(123, 456) seeds the state and observation RNG draws respectively.

testthat::test_that("simulate (use.cpp=TRUE) output matches reference", {
  data(ref_simulate_cpp)
  out = model2$simulate(Ornstein_augmented_NA, method = "ekf",
                        use.cpp = TRUE, cpp.seeds = c(123, 456), silent = TRUE)

  testthat::expect_equal(out$states,       ref_simulate_cpp$states)
  testthat::expect_equal(out$observations, ref_simulate_cpp$observations)
  testthat::expect_equal(out$times,        ref_simulate_cpp$times)
})
