
# This is a testthat script for automatically testing that the functions in
# ctsmTMB are working as intended.

obj = ctsmTMB$new()
testthat::expect_s3_class(obj, class=c("ctsmTMB","R6"))

# Try to create a model
testthat::expect_no_error({
model = ctsmTMB$new()
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

set.seed(10)
fake.data <- data.frame(
  t = c(0,1,2,3,4),
  u = rnorm(5),
  y = rnorm(5),
  z  = rnorm(5)
)

testthat::expect_no_error(
  model$estimate(fake.data, silent=T, control=list(trace=0))
)
testthat::expect_no_error(
  model$likelihood(fake.data, silent=T)
)
testthat::expect_no_error(
  model$predict(fake.data, silent=T)
)
testthat::expect_no_error(
  model$simulate(fake.data, silent=T)
)

# CPP FUNCTIONS
# testthat::expect_no_error(
#   model$predict(fake.data, silent=T, use.cpp=TRUE)
# )

# testthat::expect_no_error(
#   model$simulate(fake.data, silent=T, use.cpp=TRUE)
# )
