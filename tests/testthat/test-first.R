data(Ornstein)
df <- Ornstein[1:10,]
df$z = df$y

model = ctsmTMB$new()
testthat::expect_no_error({
model$addSystem(
  dx ~ theta * (mu*u - x) * dt + sigma_x*dw,
  dx2 ~ sigma_x * dw1
)
model$addObs(
  y ~ x,
  z ~ x^2
)
model$setVariance(
  y ~ sigma_y^2,
  z ~ sigma_y^2 * x^2
)
model$setAlgebraics(
  theta   ~ exp(logtheta),
  sigma_x ~ exp(logsigma_x),
  sigma_y ~ exp(logsigma_y)
)
model$addInput(u)
model$setParameter(
  logtheta    = log(c(initial = 1, lower=1e-5, upper=50)),
  mu          = c(initial=1.5, lower=0, upper=5),
  logsigma_x  = log(c(initial=1, lower=1e-10, upper=2)),
  logsigma_y  = log(c(initial=1, lower=1e-10, upper=2)),
  logsigma_y  = log(0.1)
)
model$setInitialState(list(rep(1,2), 0.656*diag(2)))
})

testthat::expect_no_error({
  model$estimate(df, method="ekf", silent=TRUE, trace=0)
  model$estimate(df, method="lkf", silent=TRUE, trace=0)
  model$estimate(df, method="ukf", silent=TRUE, trace=0)
  model$estimate(df, method="laplace", silent=TRUE, trace=0)
  model$estimate(df, method="laplace.thygesen", silent=TRUE, trace=0)
})

testthat::expect_no_error({
  model$filter(df, method="ekf", use.cpp=TRUE, silent=TRUE)
  model$filter(df, method="lkf", use.cpp=TRUE, silent=TRUE)
  model$filter(df, method="ukf", use.cpp=TRUE, silent=TRUE)
})
testthat::expect_no_error({
  model$filter(df, method="ekf", use.cpp=FALSE, silent=TRUE)
  model$filter(df, method="lkf", use.cpp=FALSE, silent=TRUE)
  model$filter(df, method="ukf", use.cpp=FALSE, silent=TRUE)
})

testthat::expect_no_error({
  model$smoother(df, method="laplace", silent=TRUE)
  model$smoother(df, method="laplace.thygesen", silent=TRUE)
})

testthat::expect_no_error({
  model$predict(df, method="ekf", use.cpp=TRUE, silent=TRUE)
  model$predict(df, method="lkf", use.cpp=TRUE, silent=TRUE)
  model$predict(df, method="ukf", use.cpp=TRUE, silent=TRUE)
})
testthat::expect_no_error({
  model$predict(df, method="ekf", use.cpp=FALSE, silent=TRUE)
  model$predict(df, method="lkf", use.cpp=FALSE, silent=TRUE)
  model$predict(df, method="ukf", use.cpp=FALSE, silent=TRUE)
})

testthat::expect_no_error({
  model$simulate(df, method="ekf", use.cpp=TRUE, silent=TRUE)
  model$simulate(df, method="lkf", use.cpp=TRUE, silent=TRUE)
  model$simulate(df, method="ukf", use.cpp=TRUE, silent=TRUE)
})
testthat::expect_no_error({
  model$simulate(df, method="ekf", use.cpp=FALSE, silent=TRUE)
  model$simulate(df, method="lkf", use.cpp=FALSE, silent=TRUE)
  model$simulate(df, method="ukf", use.cpp=FALSE, silent=TRUE)
})

testthat::expect_no_error({
  model$likelihood(df, method="lkf", silent=TRUE)
  model$likelihood(df, method="ukf", silent=TRUE)
  model$likelihood(df, method="ekf", silent=TRUE)
})
