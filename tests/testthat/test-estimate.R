#### Setup: 1D Ornstein model ####

data(Ornstein)
df <- Ornstein[1:10, ]

model <- create.Ornstein1D.model()


#### Setup: 2D Ornstein_augmented model ####

data(Ornstein2D)
df2 <- Ornstein2D[1:10, ]

model2 <- create.Ornstein2D.model()


#### Ornstein (1D) ####

testthat::test_that("estimate runs without error (ekf)", {
  testthat::expect_no_error(
    model$estimate(df, method = "ekf", silent = TRUE, trace = 0)
  )
})

testthat::test_that("estimate runs without error (lkf)", {
  testthat::expect_no_error(
    model$estimate(df, method = "lkf", silent = TRUE, trace = 0)
  )
})

testthat::test_that("estimate runs without error (ukf)", {
  testthat::expect_no_error(
    model$estimate(df, method = "ukf", silent = TRUE, trace = 0)
  )
})

testthat::test_that("estimate runs without error (laplace)", {
  testthat::expect_no_error(
    suppressWarnings(model$estimate(df, method = "laplace", silent = TRUE, trace = 0))
  )
})

testthat::test_that("estimate runs without error (laplace.thygesen)", {
  testthat::expect_no_error(
    suppressWarnings(model$estimate(df, method = "laplace.thygesen", silent = TRUE, trace = 0))
  )
})


#### Ornstein_augmented (2D) ####

testthat::test_that("estimate (2D) runs without error (ekf)", {
  testthat::expect_no_error(
    model2$estimate(df2, method = "ekf", silent = TRUE, trace = 0)
  )
})

testthat::test_that("estimate (2D) runs without error (lkf)", {
  testthat::expect_no_error(
    model2$estimate(df2, method = "lkf", silent = TRUE, trace = 0)
  )
})

testthat::test_that("estimate (2D) runs without error (ukf)", {
  testthat::expect_no_error(
    model2$estimate(df2, method = "ukf", silent = TRUE, trace = 0)
  )
})

testthat::test_that("estimate (2D) runs without error (laplace)", {
  testthat::expect_no_error(
    suppressWarnings(model2$estimate(df2, method = "laplace", silent = TRUE, trace = 0))
  )
})

testthat::test_that("estimate (2D) runs without error (laplace.thygesen)", {
  testthat::expect_no_error(
    suppressWarnings(model2$estimate(df2, method = "laplace.thygesen", silent = TRUE, trace = 0))
  )
})
