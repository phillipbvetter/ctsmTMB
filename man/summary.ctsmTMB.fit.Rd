% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{summary.ctsmTMB.fit}
\alias{summary.ctsmTMB.fit}
\title{Basic summary of ctsmTMB fit object}
\usage{
\method{summary}{ctsmTMB.fit}(object, correlation = FALSE, ...)
}
\arguments{
\item{object}{a ctsmTMB fit object}

\item{correlation}{boolean indicating whether or not to display the
parameter correlation structure}

\item{...}{additional arguments}
}
\value{
a summary of the estimated ctsmTMB model fit
}
\description{
Basic summary of ctsmTMB fit object
}
\examples{
library(ctsmTMB)
model <- ctsmTMB$new()

# create model
model$addSystem(dx ~ theta * (mu+u-x) * dt + sigma_x*dw)
model$addObs(y ~ x)
model$setVariance(y ~ sigma_y^2)
model$addInput(u)
model$setParameter(
  theta   = c(initial = 1, lower=1e-5, upper=50),
  mu      = c(initial=1.5, lower=0, upper=5),
  sigma_x = c(initial=1, lower=1e-10, upper=30),
  sigma_y = 1e-2
)
model$setInitialState(list(1,1e-1))

# fit model to data
fit <- model$estimate(Ornstein)

# print model summary
summary(fit, correlation=TRUE)
}
