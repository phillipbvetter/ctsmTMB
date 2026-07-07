#### Setup: 1D Ornstein model ####

data(Ornstein)
df <- Ornstein[1:10, ]

model <- create.Ornstein1D.model()


#### Setup: 2D Ornstein_augmented model ####

data(Ornstein2D)
df2 <- Ornstein2D[1:10, ]

model2 <- create.Ornstein2D.model()


#### Ornstein (1D) ####

testthat::test_that("predict runs without error using C++ (ekf)", {
  testthat::expect_no_error(
    model$predict(df, method = "ekf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("predict runs without error using C++ (lkf)", {
  testthat::expect_no_error(
    model$predict(df, method = "lkf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("predict runs without error using C++ (ukf)", {
  testthat::expect_no_error(
    model$predict(df, method = "ukf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("predict runs without error using R (ekf)", {
  testthat::expect_no_error(
    model$predict(df, method = "ekf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("predict runs without error using R (lkf)", {
  testthat::expect_no_error(
    model$predict(df, method = "lkf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("predict runs without error using R (ukf)", {
  testthat::expect_no_error(
    model$predict(df, method = "ukf", use.cpp = FALSE, silent = TRUE)
  )
})


#### Ornstein_augmented (2D) ####

testthat::test_that("predict (2D) runs without error using C++ (ekf)", {
  testthat::expect_no_error(
    model2$predict(df2, method = "ekf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("predict (2D) runs without error using C++ (lkf)", {
  testthat::expect_no_error(
    model2$predict(df2, method = "lkf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("predict (2D) runs without error using C++ (ukf)", {
  testthat::expect_no_error(
    model2$predict(df2, method = "ukf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("predict (2D) runs without error using R (ekf)", {
  testthat::expect_no_error(
    model2$predict(df2, method = "ekf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("predict (2D) runs without error using R (lkf)", {
  testthat::expect_no_error(
    model2$predict(df2, method = "lkf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("predict (2D) runs without error using R (ukf)", {
  testthat::expect_no_error(
    model2$predict(df2, method = "ukf", use.cpp = FALSE, silent = TRUE)
  )
})
