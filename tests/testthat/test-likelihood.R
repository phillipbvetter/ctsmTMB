#### Setup: 1D Ornstein model ####

data(Ornstein)
df <- Ornstein[1:10, ]

model <- create.Ornstein1D.model()


#### Setup: 2D Ornstein_augmented model ####

data(Ornstein2D)
df2 <- Ornstein2D[1:10, ]

model2 <- create.Ornstein2D.model()


#### Ornstein (1D) ####

testthat::test_that("likelihood runs without error (ekf)", {
  testthat::expect_no_error(
    model$likelihood(df, method = "ekf", silent = TRUE)
  )
})

testthat::test_that("likelihood runs without error (lkf)", {
  testthat::expect_no_error(
    model$likelihood(df, method = "lkf", silent = TRUE)
  )
})

testthat::test_that("likelihood runs without error (ukf)", {
  testthat::expect_no_error(
    model$likelihood(df, method = "ukf", silent = TRUE)
  )
})

testthat::test_that("likelihood runs without error (laplace)", {
  testthat::expect_no_error(
    model$likelihood(df, method = "laplace", silent = TRUE)
  )
})

testthat::test_that("likelihood runs without error (laplace.thygesen)", {
  testthat::expect_no_error(
    model$likelihood(df, method = "laplace.thygesen", silent = TRUE)
  )
})


#### Ornstein_augmented (2D) ####

testthat::test_that("likelihood (2D) runs without error (ekf)", {
  testthat::expect_no_error(
    model2$likelihood(df2, method = "ekf", silent = TRUE)
  )
})

testthat::test_that("likelihood (2D) runs without error (lkf)", {
  testthat::expect_no_error(
    model2$likelihood(df2, method = "lkf", silent = TRUE)
  )
})

testthat::test_that("likelihood (2D) runs without error (ukf)", {
  testthat::expect_no_error(
    model2$likelihood(df2, method = "ukf", silent = TRUE)
  )
})

testthat::test_that("likelihood (2D) runs without error (laplace)", {
  testthat::expect_no_error(
    model2$likelihood(df2, method = "laplace", silent = TRUE)
  )
})

testthat::test_that("likelihood (2D) runs without error (laplace.thygesen)", {
  testthat::expect_no_error(
    model2$likelihood(df2, method = "laplace.thygesen", silent = TRUE)
  )
})
