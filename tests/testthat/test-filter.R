#### Setup: build a minimal model used across all filter tests ####

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


#### Basic functionality ####

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
