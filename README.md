# Stochastic Differential Equations using Template Model Builder (sdeTMB)

`sdeTMB` is an R package for parameter estimation, state filtration and forecasting in stochastic state space models, heavily inspired by [Continuous Time Stochastic Modelling](https://ctsm.info). 
The package is a user-friendly wrapper for [Template Model Builder](https://github.com/kaskr/adcomp) that frees the user from writing the required C++ file containing the (negative log) likelihood function themselves. Instead, the C++ script is generated automatically based on a model specified by the user using the provided R6 `sdeTMB` class object. The package furthermore employs the `Rcpp` package universe to allow faster calculations of model predictions and stochastic simulation paths.

The package implements the following methods 

1. The (Continous-Discrete) Extended Kalman Filter, `ekf`

2. The (Continous-Discrete) Unscented Kalman Filter, `ukf`
 
3. The Laplace-style approach where latent states are considered random effects (see e.g. [this example]( https://github.com/kaskr/adcomp/blob/master/tmb_examples/sde_linear.cpp)), `laplace`

The main advantage of the Kalman Filter implementations are a large increase in the computation speed, and access to the fixed effects hessian for improved convergence of the optimization. In these cases TMB just provides automatic differentiation.

A district advantage of the `laplace`-style implementation is its use of the Laplace approximation for likelihood calculations which allows state space formulations where the density of the observation residuals are non-Gaussian.

The package is currently mostly tailored towards the Kalman Filter, with its available methods `predict` and `simulate`  for k-step-ahead predictions and simulations. It also has an `S3 method` implementation of `plot` to be called on the `sdeTMB.fit` class object returned from the `estimate` method, which plots a basic residuals analysis using the `ggplot2` package.

## Installation

You can install the package by copying the command below into `R`.
``` r
remotes::install_github(repo="phillipbvetter/sdeTMB", dependencies=TRUE)
```

We note that `sdeTMB` depends on the following packages:
1. `TMB`
2. `Rcpp`
3. `RcppEigen`
4. `RcppXPtrUtils`
5. `RcppZiggurat`
6. `R6`
7. `Deriv`
8. `stringr`
9. `stats`

The user must therefore have a working C++ compiler. In particular windows users should install Rtools, and Mac users should install Command Line Tools to get working C++ compilers. You must make sure that these are added to the `PATH` vislble to `R`. For further information see the `TMB` GitHub [here](https://github.com/kaskr/adcomp) and associated installation instructions [here](https://github.com/kaskr/adcomp/wiki/Download)

Linux users need to make sure that GSL is installed for `RcppZiggurat`. You can try the following command, or google yourself.
``` bash
sudo apt-get install libgsl-dev
```

## How to get started
You can visit the package [webpage](https://phillipbvetter.github.io/sdeTMB/index.html) and browse the vignettes for example uses, in particular see [Getting Started](https://phillipbvetter.github.io/sdeTMB/articles/sdeTMB.html).

## Help
You can access the documentation for all the available methods with
``` r
?sdeTMB
```
or individually (for a subset of methods) using i.e. `?sdeTMB::add_systems`. The methods documentation is also available on the [homepage](https://phillipbvetter.github.io/sdeTMB/reference/sdeTMB.html).

## Example Usage

```r
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(sdeTMB)

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
dt.obs = 1e-1
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
obj = sdeTMB$new()

# Set name of model (and the created .cpp file)
obj$set_modelname("ornstein_uhlenbeck")

# Add system equations
obj$add_systems(
  dx ~ theta * (mu-x+u) * dt + sigma_x*dw
)

# Add observation equations
obj$add_observations(
  y ~ x
)

# Set observation equation variances
obj$add_observation_variances(
  y ~ sigma_y^2
)

# Specify algebraic relations
obj$add_algebraics(
  theta   ~ exp(logtheta),
  sigma_x ~ exp(logsigma_x),
  sigma_y ~ exp(logsigma_y)
)

# Add vector input
obj$add_inputs(u)

# Specify parameter initial values and lower/upper bounds in estimation
obj$add_parameters(
  logtheta    = log(c(initial = 1, lower=1e-5, upper=50)),
  mu          = c(initial=1.5, lower=0, upper=5),
  logsigma_x  = log(c(initial=1, lower=1e-10, upper=30)),
  logsigma_y  = log(c(initial=1e-1, lower=1e-10, upper=30))
)

# Set initial state mean and covariance
obj$set_initial_state(list(x[1], 1e-1*diag(1)))

# Carry out estimation using extended kalman filter method with stats::nlminb as optimizer
fit <- obj$estimate(data=.data, 
                    method="ekf", 
                    use.hessian=T,
                    ode.timestep=1e-2
)

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
pred.list = obj$predict(data=.data, 
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
pred.list = obj$predict(data=.data, 
                        k.ahead=1e6, 
                        method="ekf",
)

# Perform full simulation without data update
sim.list = obj$simulate(data=.data, 
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


