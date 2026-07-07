create.Ornstein1D.model <- function(){
  model <- ctsmTMB$new()
  model$addSystem(dx ~ theta * (mu * u - x) * dt + sigma_x * dw)
  model$addObs(y ~ x)
  model$setVariance(y ~ sigma_y^2)
  model$setAlgebraics(
    theta   ~ exp(logtheta),
    sigma_x ~ exp(logsigma_x),
    sigma_y ~ exp(logsigma_y)
  )
  model$addInput("u")
  model$setParameter(
    logtheta   = log(c(initial = 1,    lower = 1e-5,  upper = 50)),
    mu         =     c(initial = 1.5,  lower = 0,     upper = 5),
    logsigma_x = log(c(initial = 1,    lower = 1e-10, upper = 2)),
    logsigma_y = log(c(initial = 0.1,  lower = 1e-10, upper = 2))
  )
  model$setInitialState(list(1, 0.1))

  return(model)
}

# A function for creating a 2D Ornstein Uhlenbeck process with selceted parameter values.
create.Ornstein2D.model <- function(){
  m = newModel()
  m$addSystem(
    dx1 ~ theta * (mu + u - x1) * dt + sigma_x1 * dw1,
    dx2 ~ alpha * (x1 - x2)    * dt + sigma_x2 * dw2
  )
  m$addObs(y1 ~ x1, y2 ~ x2)
  m$setVariance(y1 ~ sigma_y1^2, y2 ~ sigma_y2^2)
  m$setAlgebraics(
    theta    ~ exp(logtheta),
    alpha    ~ exp(logalpha),
    sigma_x1 ~ exp(logsigma_x1),
    sigma_x2 ~ exp(logsigma_x2),
    sigma_y1 ~ exp(logsigma_y1),
    sigma_y2 ~ exp(logsigma_y2)
  )
  m$addInput("u")
  m$setParameter(
    logtheta    = log(c(initial = 5,   lower = 1e-5,  upper = 50)),
    mu          =     c(initial = 3,   lower = -10,   upper = 10),
    logalpha    = log(c(initial = 2,   lower = 1e-5,  upper = 50)),
    logsigma_x1 = log(c(initial = 1,   lower = 1e-10, upper = 5)),
    logsigma_x2 = log(c(initial = 0.5, lower = 1e-10, upper = 5)),
    logsigma_y1 = log(c(initial = 0.1, lower = 1e-10, upper = 2)),
    logsigma_y2 = log(c(initial = 0.1, lower = 1e-10, upper = 2))
  )
  m$setInitialState(list(c(3, 3), diag(0.1, 2)))

  return(m)
}
