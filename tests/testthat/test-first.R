
# This is a testthat script for automatically testing that the functions in
# ctsmTMB are working as intended.

obj = ctsmTMB$new()
testthat::expect_s3_class(obj, class=c("ctsmTMB","R6"))

# Try to create a model
testthat::expect_no_error({
model = ctsmTMB$new()
model$setModelname("ornstein_uhlenbeck")
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
model$setInitialState(list(rep(1,2), 0.656*diag(2)), estimate=F)
})
