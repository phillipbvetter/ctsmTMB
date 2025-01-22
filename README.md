
<!-- README.md is generated from README.Rmd. -->
<!-- Logo -->

# ctsmTMB <img src='logo/logo.png' align="right" height="200" />

<!-- Badges -->

[![CRAN
status](https://www.r-pkg.org/badges/version/geomtextpath)](https://CRAN.R-project.org/package=geomtextpath)
[![R-CMD-check](https://github.com/AllanCameron/geomtextpath/workflows/R-CMD-check/badge.svg)](https://github.com/AllanCameron/geomtextpath/actions)
[![Codecov test
coverage](https://codecov.io/gh/AllanCameron/geomtextpath/branch/main/graph/badge.svg)](https://app.codecov.io/gh/AllanCameron/geomtextpath?branch=main)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/geomtextpath)](https://cran.r-project.org/package=geomtextpath)

<!-- Package Description -->

## Overview

**ctsmTMB** [(Continuous Time Stochastic Modelling using Template Model
Builder)](https://phillipbvetter.github.io/ctsmTMB/index.html) is an R
package for parameter estimation, state filtration and forecasting in
stochastic state space models, an intended successor of, and heavily
inspired by the **CTSM** package [(Continuous Time Stochastic
Modelling)](https://ctsm.info). The package is essentially a wrapper
around the **TMB**/**RTMB** packages [(Template Model
Builder)](https://github.com/kaskr/adcomp) used to automatically
constructs the underlying *(negative log)* likelihood function behind
the scenes (based on a user-specified model). The model is specified
using the implemented OOP-style **R6** `ctsmTMB` class, with its methods
for e.g. specifying system equations (`addSystem`) and observation
equations (`addObs`).

The primary work-horse method of **ctsmTMB** is `estimate`, used for
estimating parameters (and states) with the `stats::nlminb` quasi-Newton
optimizer due to [D.
Gay](https://dl.acm.org/doi/pdf/10.1145/355958.355965). The available
inference methods are non-linear Kalman filters and the Laplace
approxmation.

The secondary work-horse methods are `predict` and `simulate`. These are
used for integrating the stochastic differential equation forward in
time, either using (first and second order) moment differential
equations or by stochastic (euler-maruyama) simulations. The
implementation of these two methods are based on **C++** code using the
`Rcpp` package universe. The computation speed is in particular aided by
the use of the `RcppXPtrUtils` package which facilities creating and
sending **C++** pointers of the model-specific functions (drift,
diffusion, observation and associated jacobians) rather than sending
(slow) **R** functions to the **C++** side.

## Estimation Methods

The following state reconstruction algorithms are currently available:

1.  The (Continous-Discrete) Extended Kalman Filter, `ekf`.

2.  The (Continous-Discrete) Unscented Kalman Filter, `ukf`.

3.  The (Continuous-Discrete) Laplace Approximation `laplace`.

### Kalman Filters

The package is currently mostly tailored towards the Kalman Filter. The
advantages of the methods are:

1.  The hessian of the likelihood function (w.r.t parameters) is
    available.

2.  The model residuals are easier to compute for e.g. model validation.

3.  Multi-step predictions / simulations with state updates are easier
    to compute.

In these cases **TMB** simply provides an easy framework for automatic
differentiation.

The package is currently mostly tailored towards the Kalman Filter, with
its available methods `predict` and `simulate` for k-step-ahead
predictions and simulations. It also has an `S3 method` implementation
of `plot` to be called on the `ctsmTMB.fit` class object returned from
the `estimate` method, which plots a basic residuals analysis using the
`ggplot2` package.

The Unscented Kalman Filter implementation is based on *Algorithm 4.7*
in [S. Särkkä,
2007](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4303242).

### Laplace Filter

The state-reconstructions based on the `laplace` method are *smoothed*
estimates, meaning that all states are optimized jointly, given all
observations in the data. For further mathematicals details, see
[this](https://phillipbvetter.github.io/ctsmTMB/articles/laplace_approx.html)
article on the package webpage. The Laplace Approximation is natively
built-into and completely handled by **TMB**. A few noteworthy
advantages are:

1.  There is no C++ compilation needed (using **RTMB**). In addition the
    AD-compile time i.e. the call to `RTMB::MakeADFun`, is identical to
    that of pre-compiled **C++** code.

2.  The possibility for non-Gaussian (but unimodal) observation
    densities to accommodate the need for e.g. heavier distribution
    tails.

The method *may* be less useful in the context of model-training towards
forecasting because the likelihood contributions are based on these
smoothed estimates, rather than one-step predictions (as is the case of
the Kalman filters).

## Installation

You can install the package by copying the command below into `R`.

``` r
remotes::install_github(repo="phillipbvetter/ctsmTMB", dependencies=TRUE)
```

The user must have a working C++ compiler. In particular windows users
should install Rtools, and Mac users should install Command Line Tools
to get working C++ compilers. You must make sure that these are added to
the `PATH` vislble to `R`. For further information see the `TMB` GitHub
[here](https://github.com/kaskr/adcomp) and associated installation
instructions [here](https://github.com/kaskr/adcomp/wiki/Download)

Linux users need to make sure that GSL is installed for `RcppZiggurat`
which is necessary for the `simulate` method. You can try the following
command, or google yourself.

``` bash
sudo apt-get install libgsl-dev
```

## Package Dependencies

We note that `ctsmTMB` depends on the following packages: 1. `TMB` and
`RTMB` 2. `Rcpp`, `RcppEigen`, `RcppXPtrUtils` and `RcppZiggurat` 3.
`R6` 4. `Deriv` 5. `stringr`

## Getting Started

You can visit the package
[webpage](https://phillipbvetter.github.io/ctsmTMB/index.html) and
browse the vignettes for example uses, in particular see [Getting
Started](https://phillipbvetter.github.io/ctsmTMB/articles/ctsmTMB.html).

## Help

You can access the documentation for all the available methods with

``` r
?ctsmTMB
```

or individually (for a subset of methods) using
i.e. `?ctsmTMB::addSystem`.

The methods documentation is also available on the
[homepage](https://phillipbvetter.github.io/ctsmTMB/reference/ctsmTMB.html).

## Code Example - Inference in 1D Ornstein-Uhlenbeck

``` r
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(ctsmTMB)

############################################################
# Data simulation
############################################################

# Simulate data using Euler Maruyama
set.seed(20)
pars = c(theta=10, mu=1, sigma_x=1, sigma_y=0.1)
# 
dt.sim = 1e-3
t.sim = seq(0,5,by=dt.sim)
dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
u.sim = cumsum(rnorm(length(t.sim),sd=0.05))
x = 3
for(i in 1:(length(t.sim)-1)) {
  x[i+1] = x[i] + pars[1]*(pars[2]-x[i]+u.sim[i])*dt.sim + pars[3]*dw[i]
}

# Extract observations and add noise
dt.obs = 1e-2
ids = seq(1,length(t.sim),by=round(dt.obs / dt.sim))
t.obs = t.sim[ids]
y = x[ids] + pars[4] * rnorm(length(t.obs))
# forcing input
u = u.sim[ids]

# Create data
.data = data.frame(
  t = t.obs,
  y = y,
  u = u
)

############################################################
# Model creation and estimation
############################################################

# Create model object
model = ctsmTMB$new()

# Set name of model (and the created .cpp file)
model$setModelname("ornstein_uhlenbeck")

# Add system equations
model$addSystem(
  dx ~ theta * (mu-x+u) * dt + sigma_x*dw
)

# Add observation equations
model$addObs(
  y ~ x
)

# Set observation equation variances
model$setVariance(
  y ~ sigma_y^2
)

# Specify algebraic relations
model$setAlgebraics(
  theta   ~ exp(logtheta),
  sigma_x ~ exp(logsigma_x),
  sigma_y ~ exp(logsigma_y)
)

# Add vector input
model$addInput(u)

# Specify parameter initial values and lower/upper bounds in estimation
model$setParameter(
  logtheta    = log(c(initial = 1, lower=1e-5, upper=50)),
  mu          = c(initial=1.5, lower=0, upper=5),
  logsigma_x  = log(c(initial=1, lower=1e-10, upper=30)),
  logsigma_y  = log(c(initial=1e-1, lower=1e-10, upper=30))
)

# Set initial state mean and covariance
model$setInitialState(list(x[1], 1e-1*diag(1)))

# Carry out estimation with default settings (extended kalman filter)
fit <- model$estimate(data=.data, method="ekf")

# Check parameter estimates against truth
p0 = fit$par.fixed
cbind(c(exp(p0[1]),p0[2],exp(p0[3]),exp(p0[4])), pars)

# Create plot of one-step predictions, simulated states and observations
t.est = fit$states$mean$prior$t
x.mean = fit$states$mean$prior$x
x.sd = fit$states$sd$prior$x
plot1 = ggplot() +
  geom_ribbon(aes(x=t.est, ymin=x.mean-2*x.sd, ymax=x.mean+2*x.sd),fill="grey", alpha=0.9) +
  geom_line(aes(x=t.est, x.mean),col="steelblue",lwd=1) +
  geom_line(aes(x=t.sim,y=x)) + 
  geom_point(aes(x=t.obs,y=y),col="tomato",size=1) +
  labs(title="1-Step State Estimates vs Observations", x="Time", y="") +
  theme_minimal()

# Predict to obtain k-step-ahead predictions to see model forecasting ability
pred.list = model$predict(data=.data, 
                        k.ahead=10, 
                        method="ekf",
)

# Create plot all 10-step predictions against data
pred = pred.list$states
pred10step = pred %>% dplyr::filter(k.ahead==10)
plot2 = ggplot() +
  geom_ribbon(aes(x=pred10step$t.j, 
                  ymin=pred10step$x-2*sqrt(pred10step$var.x),
                  ymax=pred10step$x+2*sqrt(pred10step$var.x)),fill="grey", alpha=0.9) +
  geom_line(aes(x=pred10step$t.j,pred10step$x),color="steelblue",lwd=1) +
  geom_point(aes(x=t.obs,y=y),color="tomato",size=1) +
  labs(title="10 Step Predictions vs Observations", x="Time", y="") +
  theme_minimal()

# Perform full prediction without data update
pred.list = model$predict(data=.data, 
                        k.ahead=1e6, 
                        method="ekf",
)

# Perform full simulation without data update
sim.list = model$simulate(data=.data, 
                        k.ahead=1e6, 
                        method="ekf"
)

# Collapse simulation data for easy use with ggplot 
sim.df = sim.list$states$x$i0 %>%
  select(!c("i","j","t.i","k.ahead")) %>%
  reshape2::melt(., id.var="t.j")

# Plot all full simulations and the full prediction against observations
# (full means no data-update at all)
plot3 = ggplot() +
  geom_line(data=sim.df, aes(x=t.j, y=value, group=variable),color="grey") +
  geom_line(aes(x=pred.list$states$t.j,y=pred.list$states$x),color="steelblue") +
  geom_point(aes(x=t.obs,y=y),color="tomato",size=1) +
  labs(title="No Update Prediction and Simulations vs Observations", x="Time", y="") +
  theme_minimal() + theme(legend.position = "none")

# Draw both plots
patchwork::wrap_plots(plot1, plot2, plot3, ncol=1)

# Plot one-step-ahead residual analysis using the command below
# plot(fit)
```
