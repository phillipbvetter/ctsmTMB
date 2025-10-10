
<!-- README.md is generated from README.Rmd. -->

<!-- Logo -->

# ctsmTMB <img src='man/figures/logo.png' align="right" height="150" />

<!-- Badges -->

[![CRAN
status](https://www.r-pkg.org/badges/version/ctsmTMB)](https://CRAN.R-project.org/package=ctsmTMB)
[![R-CMD-check](https://github.com/phillipbvetter/ctsmTMB/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/phillipbvetter/ctsmTMB/actions/workflows/R-CMD-check.yaml)
[![](https://cranlogs.r-pkg.org/badges/ctsmTMB?color=brightgreen)](https://cran.rstudio.com/web/packages/ctsmTMB/index.html)
<!-- Package Description -->

# Overview

Welcome to this GitHub repository which hosts **ctsmTMB** [(Continuous
Time Stochastic Modelling using Template Model
Builder)](https://phillipbvetter.github.io/ctsmTMB/index.html), the
intended successor of, and heavily inspired by, the **CTSM** package
[(Continuous Time Stochastic Modelling)](https://ctsm.info). The purpose
of the package is to facilitate a user-friendly tool for (state and
parameter) inference, and forecasting, in (multi-dimensional)
continuous-discrete stochastic state space systems, i.e. systems on the
form

$$
dx_{t} = f\left( t, x_t, u_t, p \right) dt + g\left( t, x_t, u_t, p \right) dB_{t}
$$ $$
y_{k} = h\left( t_k, x_{t_k}, u_{t_k}, p \right)
$$

Here the latent state $x_t$ evolves continuously in time, governed by a
set of stochastic differential equations, and information about the
system is available at discrete times $t_{k}$ by the observations
$y_{k}$.

The **ctsmTMB** package is essentially wrapper around the **TMB/RTMB**
packages [(Template Model Builder)](https://github.com/kaskr/adcomp)
that provide automatic differentiation of the likelihood function, and
access to other computational tools such as the Laplace approximation.
The likelihood function is constructed based on the (symbolic)
user-provided state space equations, which may be specified using the
implemented OOP-style **R6** `ctsmTMB` class, with methods such as
`addSystem` (for defining system equations), and `addObs` (for defining
observation equations).

The primary work-horse of **ctsmTMB** is the `estimate` method, which
carries out inference by minimizing the (negative log) likelihood using
the `stats::nlminb` quasi-Newton optimizer. The resulting object
contains the maximum likelihood parameter and state estimates, and
associated marginal uncertainties. The available inference methods are
the Linear and Extended Kalman filters in addition to filtration
(actually smoothing) using a Laplace approximation approach.
<!-- The user can extract the likelihood function handles (function, gradient and hessian) with the `likelihood` method if e.g. they want to use another optimizer. -->

The package facilities forecasting through the `predict` (moment
forecasts) and `simulate` (stochastic path simulations) methods. The
calculations for these may be carried out in either **R** (default) or
for additional speed in **C++** using `Rcpp`.
<!-- in particular using the `RcppXPtrUtils` package for sending **C++** pointers of the model-specific functions (drift, diffusion, observation and jacobians) rather than sending (slow to evaluate) **R** functions. -->

<!-- Estimation Methods -->

# Estimation Methods

The following state reconstruction algorithms are currently available:

1.  The Linear Kalman Filter, `lkf`.

2.  The Extended Kalman Filter, `ekf`.

3.  The Unscented Kalman Filter, `ukf`.

4.  The Laplace Smoothers `laplace` and `laplace.thygesen`.

## Kalman Filters

The package is currently mostly tailored towards the Kalman Filter. The
advantages of the methods are:

1.  The hessian of the likelihood function (w.r.t the fixed parameters)
    is available.

2.  The model residuals are easier to compute for e.g. model validation.

3.  Multi-step predictions / simulations with state updates are easier
    to compute.

The Unscented Kalman Filter implementation is based on *Algorithm 4.7*
in [S. Särkkä, 2007](https://ieeexplore.ieee.org/document/4303242).

## Laplace Smoother

The state-reconstructions based on the `laplace` (approximation) method
are *smoothed* estimates, meaning that states are optimized jointly
conditioned on all observations. The Laplace approximation is natively
built-into and completely handled by **TMB**. The additional method
`laplace.thygesen` is an implementation of the the stability-improved
laplace method for systems with state-dependent diffusion and is due to
[Thygesen, 2025](https://arxiv.org/abs/2503.21358).

A particular advantage of the Laplace smoother is:

1.  The possibility for unimodal non-Gaussian observation densities to
    accommodate the need for e.g. heavier distribution tails. *Not yet
    implemented*.

<!-- Installation -->

# Installation (Compilers)

**NOTE: YOU NEED TO MAKE SURE YOU HAVE A WORKING C++ COMPILER IN R TO
INSTALL THE PACKAGE!!**

## Windows

C++ compilation in R requires **Rtools**:

1.  Go to <https://cran.r-project.org/bin/windows/Rtools/> and find the
    latest version.

2.  Go to “Control Panel -\> System -\>”Advanced” (tab) -\> “Environment
    Variables” -\> Highlight “Path” -\> “Edit” -\> Add to the character
    string in “Variable Value” the path to your Rtools folder
    **C:\Rtools\bin;C:\Rtools\MinGW\bin**.

## Mac / Unix

Mac users should install *Command Line Tools*. Run the following command
in the Terminal

``` bash
xcode-select --install
```

# Installation (Package)

The package can be installed from source from CRAN:

``` r
install.packages("ctsmTMB", type="source")
```

The development version is available on GitHub and R-universe:

``` r
remotes::install_github(repo="phillipbvetter/ctsmTMB", dependencies=TRUE)
```

``` r
install.packages('ctsmTMB', repos = c('https://phillipbvetter.r-universe.dev'), type="source")
```

If you encounter problems with a locked folder `00LOCK-ctsmTMB` due to a
previously failed installation remove it before reinstalling

``` r
# Enter the path on your system
enter.your.path <- "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/00LOCK-ctsmTMB"
unlink(enter.your.path, recursive = TRUE)
```

## Test Installation

Once you have installed the package is it a good idea to test whether
**TMB** and C++ compilation works. You should be able to run all
examples without compilation errors:

``` r
library(TMB)
runExample(all=TRUE)
```

For further information see the [TMB
GitHub](https://github.com/kaskr/adcomp) and its associated
[installation
instructions](https://github.com/kaskr/adcomp/wiki/Download).

# Getting Started

You can visit the package
[webpage](https://phillipbvetter.github.io/ctsmTMB/index.html) and
browse the vignettes for example uses, in particular see [Getting
Started](https://phillipbvetter.github.io/ctsmTMB/articles/ctsmTMB.html).

## Documentation

You can access the documentation for all methods with

``` r
?ctsmTMB
```

The methods documentation is also available on the [package
homepage](https://phillipbvetter.github.io/ctsmTMB/reference/ctsmTMB.html).

## Code Example - Inference in 1D Ornstein-Uhlenbeck Process

We consider estimating the parameters of the modified Ornstein-Uhlenbeck
process

$$
dx_{t} = \theta \left( \mu + u_t - x_t \right) dt + \sigma_x dB_{t}
$$

where the stationary mean $\mu$ has been augmented with the addition of
a time-varying input $u_{t}$. The observations remain linear and
Gaussian:

$$
y_{k} = x_{t_k} + \varepsilon_{t} \qquad \varepsilon_{t} \sim \mathcal{N}\left(0, \sigma_{y}^2 \right)
$$

The code chunk below simulates data from this process using an
Euler-Maruyama scheme, generates an appropriate `ctsmTMB` model object,
performs parameter estimation using an Extended Kalman Filter (the
Linear Kalman Filter `method='lkf'` could also be used) and inspects the
resulting residuals, moment predictions and stochastic simulations.

``` r
library(ctsmTMB)

############################################################
# Data simulation
############################################################

# Simulate data using Euler Maruyama
set.seed(20)
true.pars = c(theta=10, mu=1, sigma_x=1, sigma_y=0.1)
# 
dt.sim = 1e-3
t.sim = seq(0,5,by=dt.sim)
dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
u.sim = cumsum(rnorm(length(t.sim),sd=0.05))
x.sim = 3
for(i in 1:(length(t.sim)-1)) {
  x.sim[i+1] = x.sim[i] + true.pars["theta"]*(true.pars["mu"]-x.sim[i]+u.sim[i])*dt.sim + true.pars["sigma_x"]*dw[i]
}
df.sim <- data.frame(
  t = t.sim,
  u = u.sim,
  x = x.sim
)

# Extract simulation at observations time points, and create observations by adding noise
dt.obs = 1e-2
ids = seq(1,length(t.sim),by=round(dt.obs / dt.sim))
df.obs <- df.sim[ids,]
df.obs$y = df.obs$x + true.pars["sigma_y"] * rnorm(nrow(df.obs))

############################################################
# Model creation
############################################################

# Create model object
model = ctsmTMB$new()

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

# Add vector input
model$addInput(u)

# Specify parameter initial values and lower/upper bounds in estimation
model$setParameter(
  theta   = c(initial = 1, lower=1e-5, upper=50),
  mu      = c(initial=1.5, lower=0, upper=5),
  sigma_x = c(initial=1, lower=1e-10, upper=30),
  sigma_y = c(initial=1e-1, lower=1e-10, upper=30)
)

# Set initial state mean and covariance
model$setInitialState(list(x.sim[1], 1e-1*diag(1)))


############################################################
# Model estimation
############################################################

# Carry out estimation with default settings (extended kalman filter)
fit <- model$estimate(df.obs, method="ekf")

# Check parameter estimates against truth
fitted.pars <- fit$par.fixed
cbind(true.pars, fitted.pars, difference=true.pars - fitted.pars)

par(mfrow=c(3,1))
# Plot prior predictions (1-step predictions) against simulation (truth) and observations (data)
df.est <- cbind(fit$states$mean$prior, sd=fit$states$sd$prior$x)
plot(x=df.est$t, y=df.est$x, type="n", main="1-Step State Estimates vs Observations", xlab="Time", ylab="",  ylim=c(-3,3))
polygon(c(df.est$t, rev(df.est$t)), c(df.est$x+1.96*df.est$sd, rev(df.est$x-1.96*df.est$sd)), col="grey70", border=NA)
lines(df.est$t, df.est$x, col="steelblue", lwd=2)
points(df.obs$t, df.obs$y, col="tomato", pch=16, cex=0.7)

# Predict to obtain k-step-ahead predictions to see model forecasting ability
pred.10step <- model$predict(df.obs, k.ahead=10, method="ekf", return.k.ahead = 10)

# Plot 10 step predictions vs data
dfp <- pred.10step$states[c("t.j","x","var.x")]
dfp[,4] <- sqrt(dfp["var.x"])
names(dfp) <- c("t","x","var","sd")
plot(x=dfp$t, y=dfp$x, type="n", main="10 Step Predictions vs Observations", xlab="Time", ylab="", ylim=c(-3,3))
polygon(c(dfp$t, rev(dfp$t)), c(dfp$x+1.96*dfp$sd, rev(dfp$x-1.96*dfp$sd)), col="grey70", border=NA)
lines(dfp$t, dfp$x, col="steelblue", lwd=2)
points(df.obs$t, df.obs$y, col="tomato", pch=16, cex=0.7)

# Perform a full prediction i.e. without updating to data along the way
dfp <- model$predict(df.obs, method="ekf")$states

# Perform full simulations - 10 sample trajectories
sdf <- model$simulate(df.obs, method="ekf", n.sims=10)$states$x$i0
sdf.sim <- sdf[,6:ncol(sdf)]
matplot(sdf$t.j, sdf.sim, type="l", lty="solid", col="grey70", main="No Update Prediction and Simulations vs Observations", xlab="Time")
lines(dfp$t.j, dfp$x, col="steelblue", lwd=2)
points(df.obs$t, df.obs$y, col="tomato", pch=16, cex=0.7)

# Perform residual analysis
p2 <- plot(fit)
```

## Bibliography

- U. H. Thygesen and K. Kristensen (2025), *“Inference in stochastic
  differential equations using the Laplace approximation: Demonstration
  and examples”*. In:
  [arXiv:2503.21358v2](https://arxiv.org/abs/2503.21358).

- S. Särkkä, *“On Unscented Kalman Filtering for State Estimation of
  Continuous-Time Nonlinear Systems”*. In: [IEEE Transactions on
  Automatic Control, 52(9),
  1631-1641](https://ieeexplore.ieee.org/document/4303242).
