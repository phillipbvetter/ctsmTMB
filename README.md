
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
form

<!-- $$ -->
<!-- \begin{align} -->
<!-- dx_{t} & = f(t, x_t, u_t, \theta) \, dt + g(t, x_t, u_t, \theta) \, dB_{t} \\ -->
<!-- y_{t_k} & = h(t, x_t, u_t, \theta) -->
<!-- \end{align} -->
<!-- $$ -->

$$
dx_{t} = f(t, x_t, u_t, \theta) \, dt + g(t, x_t, u_t, \theta) \, dB_{t}
$$ $$
y_{t_k} = h(t_k, x_{t_k}, u_{t_k}, \theta)
$$

Here the latent state $x_t$ evolves continuously in time, governed by a
set of stochastic differential equations, and information about the
system is available at discrete times through the observations
$y_{t_k}$.

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

1.  The (Continuous-Discrete) Linear Kalman Filter, `lkf`.

2.  The (Continuous-Discrete) Extended Kalman Filter, `ekf`.

3.  The (Continuous-Discrete) Laplace “Filter” `laplace`.

## Kalman Filters

The package is currently mostly tailored towards the Kalman Filter. The
advantages of the methods are:

1.  The hessian of the likelihood function (w.r.t parameters) is
    available.

2.  The model residuals are easier to compute for e.g. model validation.

3.  Multi-step predictions / simulations with state updates are easier
    to compute.

In these cases **TMB** simply provides the automatic differentiation.

<!-- The Unscented Kalman Filter implementation is based on *Algorithm 4.7* in [S. Särkkä, 2007](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4303242). -->

## Laplace “Filter”

The state-reconstructions based on the `laplace` (approximation) method
are *smoothed* estimates, meaning that all states are optimized jointly,
given all observations in the data. The Laplace Approximation is
natively built-into and completely handled by **TMB**. The package
implements the stability-improved method due to [Thygesen,
2025](https://arxiv.org/abs/2503.21358).

A particular advantage of Laplace filter is:

1.  The possibility for unimodal non-Gaussian observation densities to
    accommodate the need for e.g. heavier distribution tails. *Not yet
    implemented*.

The method is typically not useful for model-training with the goal of
forecasting because the likelihood contributions are based on smoothed
estimates, rather than the one-step predictions of Kalman filters.

<!-- Installation -->

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

The methods documentation is also available on the [package
homepage](https://phillipbvetter.github.io/ctsmTMB/reference/ctsmTMB.html).

## Code Example - Inference in 1D Ornstein-Uhlenbeck Process

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
model$setInitialState(list(x[1], 1e-1*diag(1)))

# Carry out estimation with default settings (extended kalman filter)
fit <- model$estimate(data=.data, method="ekf")

# Check parameter estimates against truth
p0 = fit$par.fixed
cbind(p0,pars)

# Create plot of one-step predictions, simulated states and observations
t.est = fit$states$mean$prior$t
x.mean = fit$states$mean$prior$x
x.sd = fit$states$sd$prior$x
plot1 = ggplot() +
  geom_ribbon(aes(x=t.est, ymin=x.mean-2*x.sd, ymax=x.mean+2*x.sd),fill="grey", alpha=0.9) +
  geom_line(aes(x=t.est, x.mean),col="steelblue",lwd=1) +
  geom_line(aes(x=t.sim,y=x)) + 
  geom_point(aes(x=t.obs,y=y),col="tomato",size=0.5) +
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
  geom_point(aes(x=t.obs,y=y),color="tomato",size=0.5) +
  labs(title="10 Step Predictions vs Observations", x="Time", y="") +
  theme_minimal()

# Perform prediction ignoring all data
pred.list = model$predict(data=.data,method="ekf")

# Perform simulation ignoring all data
sim.list = model$simulate(data=.data, method="ekf", n.sims=10)

# Collapse simulation data for easy use with ggplot 
sim.df = sim.list$states$x$i0 %>%
  select(!c("i","j","t.i","k.ahead")) %>%
  reshape2::melt(., id.var="t.j")

# Plot all simulations and the prediction against observations
plot3 = ggplot() +
  geom_line(data=sim.df, aes(x=t.j, y=value, group=variable),color="grey") +
  geom_line(aes(x=pred.list$states$t.j,y=pred.list$states$x),color="steelblue") +
  geom_point(aes(x=t.obs,y=y),color="tomato",size=0.5) +
  labs(title="No Update Prediction and Simulations vs Observations", x="Time", y="") +
  theme_minimal() + theme(legend.position = "none")

# Create plot
p1 <- patchwork::wrap_plots(plot1, plot2, plot3, ncol=1)

# Create one-step-ahead residual analysis plot
p2 <- plot(fit)

# Wrap both plots
patchwork::wrap_plots(p1,p2[[1]], ncol=2)
```

## Bibliography

- U. H. Thygesen and K. Kristensen (2025), *“Inference in stochastic
  differential equations using the Laplace approximation: Demonstration
  and examples”*. In:
  [arXiv:2503.21358v2](https://arxiv.org/abs/2503.21358).
