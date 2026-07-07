# Tests if the outputs from filter, predict and simulate against the Ornstein2D data, matches
# the results previously obtained and stored in the OutputReferenceData.

data(Ornstein2D)
how.many.rows <- 25
df <- Ornstein2D[1:how.many.rows,]
data(OutputReferenceData)

model <- create.Ornstein2D.model()

# -----------------------------------------------------------------------
# filter tests
# -----------------------------------------------------------------------

testthat::test_that("filter (use.cpp=TRUE) output matches reference", {
  out = model$filter(df,
                     method = "ekf",
                     use.cpp = TRUE,
                     silent = TRUE)
  testthat::expect_equal(out, OutputReferenceData$filters$ekf)
})

testthat::test_that("filter (use.cpp=TRUE) output matches reference", {
  out = model$filter(df,
                     method = "lkf",
                     use.cpp = TRUE,
                     silent = TRUE)
  testthat::expect_equal(out, OutputReferenceData$filters$lkf)
})

testthat::test_that("filter (use.cpp=TRUE) output matches reference", {
  out = model$filter(df,
                     method = "ukf",
                     use.cpp = TRUE,
                     silent = TRUE)
  testthat::expect_equal(out, OutputReferenceData$filters$ukf)
})

# -----------------------------------------------------------------------
# predict tests
# -----------------------------------------------------------------------

testthat::test_that("predict (use.cpp=TRUE) output matches reference", {
  out = model$predict(df,
                      method = "ekf",
                      use.cpp = TRUE,
                      silent = TRUE)
  testthat::expect_equal(out, OutputReferenceData$predicts$ekf)
})

testthat::test_that("predict (use.cpp=TRUE) output matches reference", {
  out = model$predict(df,
                      method = "lkf",
                      use.cpp = TRUE,
                      silent = TRUE)
  testthat::expect_equal(out, OutputReferenceData$predicts$lkf)
})

testthat::test_that("predict (use.cpp=TRUE) output matches reference", {
  out = model$predict(df,
                      method = "ukf",
                      use.cpp = TRUE,
                      silent = TRUE)
  testthat::expect_equal(out, OutputReferenceData$predicts$ukf)
})

# -----------------------------------------------------------------------
# simulate tests
# -----------------------------------------------------------------------

testthat::test_that("simulate (use.cpp=TRUE) output matches reference", {
  out = model$simulate(df,
                       method = "ekf",
                       use.cpp = TRUE,
                       cpp.seeds = c(123, 456),
                       silent = TRUE)
  testthat::expect_equal(out, OutputReferenceData$simulate$ekf)
})

testthat::test_that("simulate (use.cpp=TRUE) output matches reference", {
  out = model$simulate(df,
                       method = "lkf",
                       use.cpp = TRUE,
                       cpp.seeds = c(123, 456),
                       silent = TRUE)
  testthat::expect_equal(out, OutputReferenceData$simulate$lkf)
})

testthat::test_that("simulate (use.cpp=TRUE) output matches reference", {
  out = model$simulate(df,
                       method = "ukf",
                       use.cpp = TRUE,
                       cpp.seeds = c(123, 456),
                       silent = TRUE)
  testthat::expect_equal(out, OutputReferenceData$simulate$ukf)
})
