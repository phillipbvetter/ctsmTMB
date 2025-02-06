
<!-- README.md is generated from README.Rmd. -->
<!-- Logo -->

# ctsmTMB <img src='man/figures/logo.png' align="right" height="150" />

<!-- Badges -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/geomtextpath)](https://CRAN.R-project.org/package=geomtextpath) -->

[![R-CMD-check](https://github.com/phillipbvetter/ctsmTMB/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/phillipbvetter/ctsmTMB/actions/workflows/R-CMD-check.yaml)
<!-- [![metacran downloads](https://cranlogs.r-pkg.org/badges/geomtextpath)](https://cran.r-project.org/package=geomtextpath) -->
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
form $$
\begin{align}
dx_{t} & = f(t, x_t, u_t, \theta) \, dt + g(t, x_t, u_t, \theta) \, dB_{t} \\
y_{t_k} & = h(t, x_t, u_t, \theta)
\end{align}
$$ Here the latent state $x_t$ evolves continuously in time, governed by
a set of stochastic differential equations, and information about the
system is available at discrete times through the observations
$y_{t_k}$.

The **ctsmTMB** package is essentially “just” a convience wrapper for
the **TMB**/**RTMB** packages [(Template Model
Builder)](https://github.com/kaskr/adcomp) that provide automatic
differentiation of the likelihood function, and access to other
computational tools such as the Laplace approximation. The likelihood
function is constructed based on the (symbolic) user-provided state
space equations, which may be specified using the implemented OOP-style
**R6** `ctsmTMB` class, with methods such as `addSystem` (for defining
system equations), and `addObs` (for defining observation equations).

The primary work-horse of **ctsmTMB** is the `estimate` method, which
carries out inference by minimizing the (negative log) likelihood using
the `stats::nlminb` quasi-Newton optimizer due to [D.
Gay](https://dl.acm.org/doi/pdf/10.1145/355958.355965). The resulting
object contains the maximum likelihood parameter and state estimates,
and associated marginal uncertainties. The available inference methods
are non-linear Kalman filters (EKF and UKF) and filtration by the
Laplace approxmation in a random-effects setting. The user can extract
the likelihood function handles (function, gradient and hessian) with
the `likelihood` method if e.g. they want to use another optimizer.

The package facilities forecasting through the `predict` and `simulate`
methods (state updates for k-step ahead forecasts are currently only
available using Kalman filters). The difference between the two is that
the former produces moment predictions (mean and covariance) while the
latter produces stochastic path simulations (distributions). The
calculations are carried out in **C++** through the `Rcpp` package
universe, in particular using the `RcppXPtrUtils` package for sending
**C++** pointers of the model-specific functions (drift, diffusion,
observation and jacobians) rather than sending (slow to evaluate) **R**
functions.

# Estimation Methods

The following state reconstruction algorithms are currently available:

1.  The (Continous-Discrete) Extended Kalman Filter, `ekf`.

2.  The (Continous-Discrete) Unscented Kalman Filter, `ukf`.

3.  The (Continuous-Discrete) Laplace Approximation `laplace`.

## Kalman Filters

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

## Laplace Filter

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

# Installation

You can install the package by copying the command below into `R`.

``` r
remotes::install_github(repo="phillipbvetter/ctsmTMB", dependencies=TRUE)
```

It is important to note that users must have working C++ compilers to
install and use **ctsmTMB**.

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

<!---
Linux also need to make sure that GSL is installed for `RcppZiggurat` which is necessary for the `simulate` method. You can try the following command, or google yourself.
&#10;
``` bash
sudo apt-get install libgsl-dev
```
--->

## Test the Installation

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

or invidually using the standard syntax

``` r
?ctsmTMB::addSystem
?ctsmTMB::estimate
```

The methods documentation is also available on the [package
homepage](https://phillipbvetter.github.io/ctsmTMB/reference/ctsmTMB.html).

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
