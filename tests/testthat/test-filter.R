#### Setup: 1D Ornstein model ####

data(Ornstein)
df <- Ornstein[1:10, ]

model <- ctsmTMB$new()
model$addSystem(dx ~ theta * (mu * u - x) * dt + sigma_x * dw)
model$addObs(y ~ x)
model$setVariance(y ~ sigma_y^2)
model$setAlgebraics(
  theta   ~ exp(logtheta),
  sigma_x ~ exp(logsigma_x),
  sigma_y ~ exp(logsigma_y)
)
model$addInput(u)
model$setParameter(
  logtheta   = log(c(initial = 1,    lower = 1e-5,  upper = 50)),
  mu         =     c(initial = 1.5,  lower = 0,     upper = 5),
  logsigma_x = log(c(initial = 1,    lower = 1e-10, upper = 2)),
  logsigma_y = log(c(initial = 0.1,  lower = 1e-10, upper = 2))
)
model$setInitialState(list(1, 0.1))


#### Setup: 2D Ornstein_augmented model ####

data(Ornstein_augmented)
df2 <- Ornstein_augmented[1:10, ]

model2 <- ctsmTMB$new()
model2$addSystem(
  dx1 ~ theta * (mu + u - x1) * dt + sigma_x1 * dw1,
  dx2 ~ alpha * (x1 - x2)    * dt + sigma_x2 * dw2
)
model2$addObs(
  y1 ~ x1,
  y2 ~ x2
)
model2$setVariance(
  y1 ~ sigma_y1^2,
  y2 ~ sigma_y2^2
)
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
  df_na$y <- NA_real_
  testthat::expect_no_error(
    model$filter(df_na, method = "ekf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using C++ (lkf)", {
  df_na <- df
  df_na$y <- NA_real_
  testthat::expect_no_error(
    model$filter(df_na, method = "lkf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using C++ (ukf)", {
  df_na <- df
  df_na$y <- NA_real_
  testthat::expect_no_error(
    model$filter(df_na, method = "ukf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using R (ekf)", {
  df_na <- df
  df_na$y <- NA_real_
  testthat::expect_no_error(
    model$filter(df_na, method = "ekf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using R (lkf)", {
  df_na <- df
  df_na$y <- NA_real_
  testthat::expect_no_error(
    model$filter(df_na, method = "lkf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter handles all-NA observations without error using R (ukf)", {
  df_na <- df
  df_na$y <- NA_real_
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
  df2_na$y1 <- NA_real_
  df2_na$y2 <- NA_real_
  testthat::expect_no_error(
    model2$filter(df2_na, method = "ekf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using C++ (lkf)", {
  df2_na <- df2
  df2_na$y1 <- NA_real_
  df2_na$y2 <- NA_real_
  testthat::expect_no_error(
    model2$filter(df2_na, method = "lkf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using C++ (ukf)", {
  df2_na <- df2
  df2_na$y1 <- NA_real_
  df2_na$y2 <- NA_real_
  testthat::expect_no_error(
    model2$filter(df2_na, method = "ukf", use.cpp = TRUE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using R (ekf)", {
  df2_na <- df2
  df2_na$y1 <- NA_real_
  df2_na$y2 <- NA_real_
  testthat::expect_no_error(
    model2$filter(df2_na, method = "ekf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using R (lkf)", {
  df2_na <- df2
  df2_na$y1 <- NA_real_
  df2_na$y2 <- NA_real_
  testthat::expect_no_error(
    model2$filter(df2_na, method = "lkf", use.cpp = FALSE, silent = TRUE)
  )
})

testthat::test_that("filter (2D) handles all-NA observations without error using R (ukf)", {
  df2_na <- df2
  df2_na$y1 <- NA_real_
  df2_na$y2 <- NA_real_
  testthat::expect_no_error(
    model2$filter(df2_na, method = "ukf", use.cpp = FALSE, silent = TRUE)
  )
})
