# Laplace estimation (Thygesen formulation) - random effects are the discretized
# states AND the brownian increments dB. States are softly constrained to follow
# the euler-maruyama discretization through a "tiny" covariance.

# integrator for time update
# ------------------------------------
laplace.thygesen.construct.ode.integrator <- function(self, private){

  n.states <- private$dims$states
  n.pars   <- private$dims$pars
  n.inputs <- private$dims$inputs
  n.diffusions <- private$dims$diffusions
  funs <- private$model$state.space.funs

  # tiny covariance softly constraining state increments to the euler-maruyama
  # discretization, and identity for the brownian increment covariance
  tiny <- 1e-5 * diag(n.states)
  I0 <- diag(n.diffusions)

  arglist <- list(
    numeric(n.states), # state at start of step
    numeric(n.states), # state at end of step
    numeric(n.pars),
    numeric(n.inputs),
    numeric(n.diffusions), # brownian increment of the step
    1 # IMPORTANT!: dt CAN'T BE INITIALIZED ZERO OTHERWISE THE TAPE
    # IS WRONG BECAUSE RTMB::dmvnorm decomposes matrix difficult
    # it V is zero!
  )
  step.tape <- RTMB::MakeTape(function(z) {
    x0 <- z[[1]]
    x1 <- z[[2]]
    p <- z[[3]]
    u <- z[[4]]
    dB <- z[[5]]
    dt <- z[[6]]
    f <- funs$f(x0, p, u)
    G <- funs$g(x0, p, u)
    # state increments should follow the euler-maruyama prediction f*dt + G%*%dB
    e <- (x1 - x0) - (f * dt + G %*% dB)
    # soft constraint on state increments + brownian increment likelihood
    -RTMB::dmvnorm(e, Sigma = tiny, log = TRUE) -
      RTMB::dmvnorm(dB, Sigma = dt * I0, log = TRUE)
  }, arglist)

  time.update <- function(X, p, u, du, dB, steps, dt){
    nll <- 0
    for (j in seq_len(steps)) {
      nll <- nll + step.tape(list(X[j,], X[j+1,], p, u, dB[j,], dt))
      u <- u + du
    }
    return(nll)
  }
  return(time.update)
}

laplace.thygesen.construct.data.updater <- function(self, private){

  funs <- private$model$state.space.funs
  dims <- private$dims

  arglist <- list(
    numeric(dims$states),
    numeric(dims$pars),
    numeric(dims$inputs),
    numeric(dims$observations),
    numeric(dims$observations)
  )

  obs.tape <- RTMB::MakeTape(function(z) {
    y <- z[[4]]
    m <- z[[5]]
    h <- do.call(funs$h,z[1:3])
    sd <- sqrt(do.call(funs$hvar, z[1:3]) + (1-m))
    -sum(m * RTMB::dnorm(y, mean = h, sd = sd, log = TRUE))
  }, arglist)


  update.fun <- function(x,p,u,y,m){
    nll <- obs.tape(list(x,p,u,y,m))
    return(nll)
  }

  return(update.fun)
}

# setup function
# ------------------------------------
laplace.thygesen.setup <- function(self, private){

  settings <- list(
    first.order.hold = private$algo.settings$first.order.input.hold,
    estimate.initial.state = private$algo.settings$estimate.initial
  )

  advanced.settings <- list(
    nllreport = private$algo.settings$advanced.settings$nllreport,
    force.ad.on.input = private$algo.settings$advanced.settings$force.ad.on.input
  )

  dims <- list(
    n.pars = private$dims$pars,
    n.obs = private$dims$observations,
    n.diffusions = private$dims$diffusions,
    n.dbs = nrow(private$algo.settings$tmb.initial.state) - 1
  )

  data <- list(
    x = private$algo.settings$initial.state$x0,
    P = private$algo.settings$initial.state$p0,
    U = as.matrix(private$data[private$names$inputs]),
    Y = as.matrix(private$data[private$names$obs]),
    dts = private$algo.settings$ode.stepsizes,
    nsteps = private$algo.settings$ode.number.of.steps,
    nsteps.cumsum = private$algo.settings$ode.number.of.steps.cumsum
  )
  data$dU <- with(c(data,settings),
                  if(first.order.hold) rbind(diff(U)/nsteps,0) else array(0,dim=dim(U)))
  data$M <- with(data,1*!is.na(Y))
  data$Y0 <- with(data, {Y[is.na(Y)]<-0;Y})
  data$K <- with(data, unname(which(rowSums(M)>0)))

  algorithms <- list(
    time.update = laplace.thygesen.construct.ode.integrator(self, private),
    data.update = laplace.thygesen.construct.data.updater(self, private),
    stationary.solver = construct.initial.state.estimators(self, private)
  )

  # return
  return(named_list(dims, data,settings,advanced.settings,algorithms))
}

MakeADFun_Laplace_thygesen <- function(self, private)
{

  # Settings
  configure_ad_tape("RTMB", self, private)
  setup <- laplace.thygesen.setup(self, private)
  list2env(setup$data, envir=environment())
  list2env(setup$dims, envir=environment())
  list2env(setup$settings, envir=environment())
  list2env(setup$advanced.settings, envir=environment())
  list2env(setup$algorithms, envir=environment())
  funs <- private$model$state.space.funs

  # likelihood function --------------------------------------
  laplace.nll = function(p_fixed_and_random){

    ####### Sometimes necessary to avoid rtmb errors #######
    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")
    "c" <- RTMB::ADoverload("c")

    # Convert random effects lists to matrices
    # (fixed effects are the first n.pars components, the brownian
    # increments dB the last component, and the states in between)
    X <- do.call(cbind, p_fixed_and_random[(n.pars+1):(length(p_fixed_and_random)-1)])
    dB <- p_fixed_and_random$dB
    # Convert fixed effects list to vector
    p <- do.call(c, p_fixed_and_random[1:n.pars])

    # Set likelihood to zero
    nll <- 0

    # Cast initial state and inputs to AD
    # (NOTE: Doesnt appear to be necessary here - we leave this option in case of issues)
    if(force.ad.on.input) {
      x <- RTMB::AD(x)
      P <- RTMB::AD(P)
      U <- RTMB::AD(U)
    }

    # Estimate initial conditions as stationary solution with first input
    # This estimate is, as with the user provided initial state, assumed to be
    # the prior mean and covariance x_{0 | -1} and P_{0 | -1}.
    if (any(estimate.initial.state)) {
      get.cov <- estimate.initial.state[2]
      sol <- stationary.solver(p,U[1,],get.cov=get.cov,get.mean=estimate.initial.state[1],x0=x)
      x <- sol[,1]
      if(get.cov) P <- sol[,-1, drop=FALSE]
    }

    # prior contribution
    nll <- nll - RTMB::dmvnorm(X[1,] - x, Sigma=P, log=TRUE)
    for(i in 1:(nrow(Y)-1)) {
      # time update
      X0 <- X[nsteps.cumsum[i]:nsteps.cumsum[i+1] + 1, , drop=FALSE]
      dB0 <- dB[(nsteps.cumsum[i]+1):nsteps.cumsum[i+1], , drop=FALSE]
      nll0 <- time.update(X0, p, U[i,], dU[i,], dB0, nsteps[i], dts[i])
      nll <- nll + nll0
    }
    # data update
    # Y = RTMB::OBS(Y)
    for(i in K){
      nll0 <- data.update(X[nsteps.cumsum[i]+1,,drop=FALSE],p,U[i,],Y0[i,],M[i,])
      nll <- nll + nll0
    }

    # report
    if (nllreport) {
      RTMB::REPORT(X)
    }

    # return
    return(nll)
  }

  ################################################
  # Construct Neg. Log-Likelihood
  ################################################

  # initial brownian increments (one row per discretization step)
  dB <- matrix(0, nrow=n.dbs, ncol=n.diffusions)
  parameters <- c(
    lapply(private$model$parameters, function(x) x$initial),
    private$algo.settings$tmb.initial.state,
    dB = list(dB)
  )
  # parameters <- c(
  #   lapply(private$model$parameters, function(x) x$initial),
  #   X = list(as.matrix(unname(private$algo.settings$tmb.initial.state))),
  #   dB = list(dB)
  # )
  map <- lapply(private$model$fixed.pars, function(x) x$factor)
  nll = RTMB::MakeADFun(func = laplace.nll,
                        parameters=parameters,
                        random=c(private$names$states, "dB"),
                        map = map,
                        silent=TRUE)

  # save objective function
  private$nll = nll

  # return
  return(invisible(self))

}
