#### Setup: 1D Ornstein model ####

data(Ornstein)
df <- Ornstein[1:10, ]

model <- create.Ornstein1D.model()


#### Setup: 2D Ornstein_augmented model ####

data(Ornstein2D)
df2 <- Ornstein2D[1:10, ]

model2 <- create.Ornstein2D.model()


#### Ornstein (1D) - basic functionality ####

testthat::test_that("filter runs without error using C++ (ekf)", {
  testthat::expect_no_error(
    model$filter(df, method = "ekf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter runs without error using C++ (lkf)", {
  testthat::expect_no_error(
    model$filter(df, method = "lkf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter runs without error using C++ (ukf)", {
  testthat::expect_no_error(
    model$filter(df, method = "ukf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter runs without error using R (ekf)", {
  testthat::expect_no_error(
    model$filter(df, method = "ekf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter runs without error using R (lkf)", {
  testthat::expect_no_error(
    model$filter(df, method = "lkf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter runs without error using R (ukf)", {
  testthat::expect_no_error(
    model$filter(df, method = "ukf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using C++ (ekf)", {
  df_na <- df
  df_na$y <- NA
  testthat::expect_no_error(
    model$filter(df_na, method = "ekf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using C++ (lkf)", {
  df_na <- df
  df_na$y <- NA
  testthat::expect_no_error(
    model$filter(df_na, method = "lkf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using C++ (ukf)", {
  df_na <- df
  df_na$y <- NA
  testthat::expect_no_error(
    model$filter(df_na, method = "ukf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using R (ekf)", {
  df_na <- df
  df_na$y <- NA
  testthat::expect_no_error(
    model$filter(df_na, method = "ekf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using R (lkf)", {
  df_na <- df
  df_na$y <- NA
  testthat::expect_no_error(
    model$filter(df_na, method = "lkf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using R (ukf)", {
  df_na <- df
  df_na$y <- NA
  testthat::expect_no_error(
    model$filter(df_na, method = "ukf", use.cpp = FALSE, silent = TRUE)
  )
})


#### Ornstein_augmented (2D) - basic functionality ####

testthat::test_that("filter (2D) runs without error using C++ (ekf)", {
  testthat::expect_no_error(
    model2$filter(df2, method = "ekf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) runs without error using C++ (lkf)", {
  testthat::expect_no_error(
    model2$filter(df2, method = "lkf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) runs without error using C++ (ukf)", {
  testthat::expect_no_error(
    model2$filter(df2, method = "ukf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) runs without error using R (ekf)", {
  testthat::expect_no_error(
    model2$filter(df2, method = "ekf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) runs without error using R (lkf)", {
  testthat::expect_no_error(
    model2$filter(df2, method = "lkf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) runs without error using R (ukf)", {
  testthat::expect_no_error(
    model2$filter(df2, method = "ukf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using C++ (ekf)", {
  df2_na <- df2
  df2_na$y1 <- NA
  df2_na$y2 <- NA
  testthat::expect_no_error(
    model2$filter(df2_na, method = "ekf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using C++ (lkf)", {
  df2_na <- df2
  df2_na$y1 <- NA
  df2_na$y2 <- NA
  testthat::expect_no_error(
    model2$filter(df2_na, method = "lkf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using C++ (ukf)", {
  df2_na <- df2
  df2_na$y1 <- NA
  df2_na$y2 <- NA
  testthat::expect_no_error(
    model2$filter(df2_na, method = "ukf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using R (ekf)", {
  df2_na <- df2
  df2_na$y1 <- NA
  df2_na$y2 <- NA
  testthat::expect_no_error(
    model2$filter(df2_na, method = "ekf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using R (lkf)", {
  df2_na <- df2
  df2_na$y1 <- NA
  df2_na$y2 <- NA
  testthat::expect_no_error(
    model2$filter(df2_na, method = "lkf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using R (ukf)", {
  df2_na <- df2
  df2_na$y1 <- NA
  df2_na$y2 <- NA
  testthat::expect_no_error(
    model2$filter(df2_na, method = "ukf", use.cpp = FALSE, silent = TRUE)
  )
})
