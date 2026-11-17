# UKF implementation with taped time-update integrator and taped data-update,
# following the same design as algo_estimate_EKF.R.
# The UKF sigma-point ODE (cov.ode.1step from the RTMB implementation in
# algo_estimate_RTMB.R) is taped once with RTMB::MakeTape, and the data update
# is a taped joint update (ukf.construct.data.updater, copied from
# kalman.data.update.with.nll in get_ukf_update) or a taped scalar update
# processing one observation at a time (ukf.construct.data.updater2).
# This keeps the AD tape small, which speeds up compilation a lot.

# UKF weights from hyperparameters
# ------------------------------------
ukf.compute.weights <- function(self, private){
  n.states <- private$dims$states
  n.sigmapoints <- 2*n.states + 1

  alpha <- private$algo.settings$ukf.hyperpars[1]
  beta  <- private$algo.settings$ukf.hyperpars[2]
  kappa <- private$algo.settings$ukf.hyperpars[3]

  sqrt_c <- sqrt(alpha^2*(n.states + kappa))
  lambda <- sqrt_c^2 - n.states

  # weights for mean and covariance
  W.m <- W.c <- rep(1/(2*(n.states+lambda)), n.sigmapoints)
  W.m[1] <- lambda/(n.states+lambda)
  W.c[1] <- lambda/(n.states+lambda) + (1 - alpha^2 + beta)
  W.m.mat <- replicate(n.sigmapoints, W.m)
  W <- (diag(n.sigmapoints) - W.m.mat) %*% diag(W.c) %*% t(diag(n.sigmapoints) - W.m.mat)

  list(n.states=n.states, n.sigmapoints=n.sigmapoints,
       sqrt_c=sqrt_c, lambda=lambda, W.m=W.m, W.c=W.c, W=W)
}

# time update: taped sigma-point ODE integrator
# ------------------------------------
ukf.construct.time.updater <- function(self, private){
  solver <- private$algo.settings$ode.solver
  n.states <- private$dims$states
  n.pars   <- private$dims$pars
  n.inputs <- private$dims$inputs
  funs <- private$model$state.space.funs

  # weights
  w <- ukf.compute.weights(self, private)
  n.sigmapoints <- w$n.sigmapoints
  sqrt_c <- w$sqrt_c
  W <- w$W
  W.m <- w$W.m

  # get chol(P) from sigma points
  sigma2chol <- function(X){
    (X[,2:(n.states+1)] - X[,1])/sqrt_c
  }

  # rhs of the sigma-point covariance ODE (same as cov.ode.1step in the
  # untaped UKF implementation)
  cov.ode.1step <- function(X, sqrtP, p, u){

    # diffusion at the mean
    G <- funs$g(X[,1], p, u)
    F.sigma <- do.call(cbind, lapply(seq_len(n.sigmapoints), \(i) funs$f(X[,i], p, u)))

    # rhs term 2: [0, chol(P)*Phi(M), -chol(P)*Phi(M)]
    Ainv <- RTMB::solve(sqrtP)
    Phi.M <- Ainv %*% (X %*% W %*% t(F.sigma) + F.sigma %*% W %*% t(X) + G %*% t(G)) %*% t(Ainv)
    Phi.M[!lower.tri(Phi.M, diag=TRUE)] <- 0
    diag(Phi.M) <- diag(Phi.M)/2
    z0 <- sqrtP %*% Phi.M
    
    # dX = rhs1 + sqrt(c) * rhs2
    F.sigma %*% W.m + sqrt_c * cbind(0, z0, -z0)
  }

  # one integration step of the sigma points
  ode.integrator <- switch(solver,
    euler = function(X, sqrtP, p, u, du, dt){
      X + cov.ode.1step(X, sqrtP, p, u) * dt
    },
    rk4 = function(X, sqrtP, p, u, du, dt){
      X0 <- X
      k1 <- cov.ode.1step(X, sqrtP, p, u)
      X <- X0 + k1*dt/2
      sqrtP <- sigma2chol(X)
      k2 <- cov.ode.1step(X, sqrtP, p, u + du/2)
      X <- X0 + k2*dt/2
      sqrtP <- sigma2chol(X)
      k3 <- cov.ode.1step(X, sqrtP, p, u + du/2)
      X <- X0 + k3*dt
      sqrtP <- sigma2chol(X)
      k4 <- cov.ode.1step(X, sqrtP, p, u + du)
      X0 + (k1 + 2*k2 + 2*k3 + k4)/6 * dt
    },
    stop("Unknown solver - choose 'euler' or 'rk4'.")
  )

  # tape the one-step integrator (much faster AD graph compilation)
  arglist <- list(
    X = array(0, dim=c(n.states, n.sigmapoints)),
    sqrtP = array(0, dim=c(n.states, n.states)),
    p = numeric(n.pars),
    u = numeric(n.inputs),
    du = numeric(n.inputs),
    dt = numeric(1)
  )
  taped.integrator <- RTMB::MakeTape(function(arglist) do.call(ode.integrator, arglist), arglist)

  # integrate sigma points forward 'steps' steps
  time.update <- function(X, sqrtP, p, u, du, steps, dt){
    for (j in seq_len(steps)) {
      X <- taped.integrator(list(X, sqrtP, p, u, du, dt))
      sqrtP <- sigma2chol(X)
      u <- u + du
    }
    list(X=X, sqrtP=sqrtP)
  }

  return(time.update)
}

# data update: taped joint UKF update (all available observations at once),
# copied from kalman.data.update.with.nll in get_ukf_update (old implementation)
# ------------------------------------
ukf.construct.data.updater.joint <- function(self, private){

  # Unpack dimensions
  n.states <- private$dims$states
  n.pars   <- private$dims$pars
  n.inputs <- private$dims$inputs
  n.obs    <- private$dims$observations
  funs <- private$model$state.space.funs

  # weights
  w <- ukf.compute.weights(self, private)
  n.sigmapoints <- w$n.sigmapoints
  W <- w$W
  W.m <- w$W.m

  # Loss function (multivariate, as in the old joint update)
  loss.type <- private$algo.settings$loss$loss
  loss.c <- private$algo.settings$loss$c
  loss.function <- construct.loss.function(dim="joint", type=loss.type, c=loss.c)

  # Taped h of sigma points together with the observation variance matrix at the mean
  # output: n.obs x (n.sigmapoints + n.obs) matrix [Y.sigma | V]
  obs.sigma.tape <- RTMB::MakeTape(function(z) {
    X <- z[[1]]
    x <- z[[2]]
    p <- z[[3]]
    u <- z[[4]]
    Y.sigma <- do.call(cbind, lapply(seq_len(n.sigmapoints), \(i) funs$h(X[,i], p, u)))
    cbind(Y.sigma, funs$hvar_matrix(x, p, u))
  }, list(array(0, dim=c(n.states, n.sigmapoints)), numeric(n.states), numeric(n.pars), numeric(n.inputs)))

  # Joint update with all available observations, copied from
  # kalman.data.update.with.nll in get_ukf_update (old implementation).
  # NOTE: the subsetting to available observations (which(y_bool)) must happen
  # outside a MakeTape - tape arguments become AD objects, so which() would fail
  # inside a tape. Here the function runs inside the outer AD tape, where
  # subsetting AD matrices by a plain logical vector works (as in the old code).
  # Only the evaluation of h on the sigma points is taped (obs.sigma.tape).
  joint.update <- function(x, P, X, p, u, y, y_bool){
    obs <- obs.sigma.tape(list(X, x, p, u))
    idx <- which(y_bool)
    H <- obs[idx, seq_len(n.sigmapoints), drop=FALSE]
    V <- obs[idx, n.sigmapoints + idx, drop=FALSE]
    # Innovation
    e <- y[idx] - H %*% W.m
    # Measurement variance
    R <- H %*% W %*% t(H) + V
    # Cross covariance between states and observations
    P.xy <- X %*% W %*% t(H)
    # Kalman gain
    K <- t(RTMB::solve(R, t(P.xy)))
    # State update
    x <- x + K %*% e
    # Cov update
    # P <- P - K %*% R %*% t(K)
    # Generalized Joseph Form
    P <- P - P.xy %*% t(K) - K %*% t(P.xy) + K %*% R %*% t(K)
    # Resymmetrize
    P <- 0.5 * (P + t(P))
    # Likelihood
    nll <- loss.function(e, R)
    return(list(x=x, P=P, nll=nll))
  }

  # Joint likelihood contribution without state update
  # (used when training against full prediction; same as kalman.no.update.with.nll)
  joint.no.update <- function(x, P, X, p, u, y, y_bool){
    obs <- obs.sigma.tape(list(X, x, p, u))
    idx <- which(y_bool)
    H <- obs[idx, seq_len(n.sigmapoints), drop=FALSE]
    V <- obs[idx, n.sigmapoints + idx, drop=FALSE]
    e <- y[idx] - H %*% W.m
    R <- V
    nll <- loss.function(e, R)
    return(list(x=x, P=P, nll=nll))
  }

  create_ukf_update_function <- function(data.update.fun){
    f <- function(x, P, X, p, u, y){
      y_bool = !is.na(y)
      if (any(y_bool)) {
        out <- data.update.fun(x, P, X, p, u, y, y_bool)
      } else {
        out <- list(x=x, P=P, nll=0)
      }
      return(out)
    }
    return(f)
  }

  if(private$algo.settings$train.against.full.prediction){
    update.fun <- create_ukf_update_function(joint.no.update)
  } else {
    update.fun <- create_ukf_update_function(joint.update)
  }

  return(update.fun)
}

# data update: taped scalar (one observation at a time) UKF update
# ------------------------------------
ukf.construct.data.updater.scalar <- function(self, private){

  # Unpack dimensions
  n.states <- private$dims$states
  n.pars   <- private$dims$pars
  n.inputs <- private$dims$inputs
  funs <- private$model$state.space.funs

  # weights
  w <- ukf.compute.weights(self, private)
  n.sigmapoints <- w$n.sigmapoints
  sqrt_c <- w$sqrt_c
  W <- w$W
  W.m <- w$W.m

  # Loss function
  loss.type <- private$algo.settings$loss$loss
  loss.c <- private$algo.settings$loss$c
  loss.function <- construct.loss.function(dim="scalar", type=loss.type, c=loss.c)

  # jitter for the cholesky factorization
  epsI0 <- 1e-8 * diag(n.states)

  # Sample sigma-points from mean and sqrt-covariance
  create.sigmaPoints <- function(x, sqrtP){
    xs <- rep(x, n.sigmapoints)
    dim(xs) <- c(n.states, n.sigmapoints)
    z <- sqrt_c * c(RTMB::AD(numeric(n.states)), sqrtP, -sqrtP)
    dim(z) <- c(n.states, n.sigmapoints)
    xs + z
  }

  # Taped h of sigma points together with observation variances at the mean
  obs.sigma.tape <- RTMB::MakeTape(function(z) {
    X <- z[[1]]
    x <- z[[2]]
    p <- z[[3]]
    u <- z[[4]]
    Y.sigma <- do.call(cbind, lapply(seq_len(n.sigmapoints), \(i) funs$h(X[,i], p, u)))
    cbind(Y.sigma, funs$hvar(x, p, u))
  }, list(array(0, dim=c(n.states, n.sigmapoints)), numeric(n.states), numeric(n.pars), numeric(n.inputs)))

  # Taped scalar UKF update for a single observation j
  scalar.update.tape <- RTMB::MakeTape(function(z) {
    # Unpack list
    x <- z[[1]]
    P <- z[[2]]
    X <- z[[3]]
    Hj <- RTMB::matrix(z[[4]], nrow = 1) # h of sigma points, 1 x n.sigmapoints
    vj <- z[[5]]                          # observation variance
    yj <- z[[6]]                          # observation value
    # Innovation
    e <- yj - drop(Hj %*% W.m)
    # Measurement variance (scalar)
    R <- drop(Hj %*% W %*% t(Hj)) + vj
    # Cross covariance between states and the j'th observation
    P.xy <- X %*% W %*% t(Hj) # n.states x 1
    # Kalman gain
    K <- P.xy / R
    # State update
    x <- x + K * e
    # Cov update - generalized Joseph form, algebraically equal to P - R*K%*%t(K)
    P <- P - P.xy %*% t(K) - K %*% t(P.xy) + R * K %*% t(K)
    # Resymmetrize
    P <- 0.5 * (P + t(P))
    # Likelihood
    nll <- loss.function(e, R)
    # Return
    # It is faster to return a matrix for subsetting than to return a vector with c()
    # so we just repeat nll to match the rows of x/P
    cbind(x, P, rep(nll, n.states))
  }, list(numeric(n.states), array(0, dim=c(n.states,n.states)), array(0, dim=c(n.states,n.sigmapoints)), numeric(n.sigmapoints), 1, 1))

  # Sequentially update with each available observation.
  # The sigma points are regenerated from the current (x, P) before each
  # scalar update, so that the sigma points stay consistent with the state
  # distribution and P remains positive definite.
  sequential.update <- function(x, P, X, p, u, y, y_bool){
    nll <- 0
    idx <- which(y_bool)
    for (j in idx) {
      sqrtP <- t(Matrix::chol(P + epsI0))
      X <- create.sigmaPoints(x, sqrtP)
      obs <- obs.sigma.tape(list(X, x, p, u))
      Y.sigma <- obs[, seq_len(n.sigmapoints), drop=FALSE]
      V <- obs[, n.sigmapoints+1]
      out <- scalar.update.tape(list(x, P, X, Y.sigma[j,], V[j], y[j]))
      x <- out[,1]
      # cbind(x, P, rep(nll, n.states)) flattens P column-major into n.states columns
      P <- out[, 1 + 1:n.states, drop=FALSE]
      nll <- nll + out[1, 2 + n.states]
    }
    return(list(x=x, P=P, nll=nll))
  }

  # Taped scalar likelihood contribution without state update
  # (used when training against full prediction)
  scalar.no.update.tape <- RTMB::MakeTape(function(z) {
    Hj <- RTMB::matrix(z[[1]], nrow = 1)
    vj <- z[[2]]
    yj <- z[[3]]
    e <- yj - drop(Hj %*% W.m)
    R <- vj
    loss.function(e, R)
  }, list(numeric(n.sigmapoints), 1, 1))

  sequential.no.update <- function(x, P, X, p, u, y, y_bool){
    nll <- 0
    idx <- which(y_bool)
    obs <- obs.sigma.tape(list(X, x, p, u))
    Y.sigma <- obs[, seq_len(n.sigmapoints), drop=FALSE]
    V <- obs[, n.sigmapoints+1]
    for (j in idx) {
      nll <- nll + scalar.no.update.tape(list(Y.sigma[j,], V[j], y[j]))
    }
    return(list(x=x, P=P, nll=nll))
  }

  create_ukf_update_function <- function(data.update.fun){
    f <- function(x, P, X, p, u, y){
      y_bool = !is.na(y)
      if (any(y_bool)) {
        out <- data.update.fun(x, P, X, p, u, y, y_bool)
      } else {
        out <- list(x=x, P=P, nll=0)
      }
      return(out)
    }
    return(f)
  }

  if(private$algo.settings$train.against.full.prediction){
    update.fun <- create_ukf_update_function(sequential.no.update)
  } else {
    update.fun <- create_ukf_update_function(sequential.update)
  }

  return(update.fun)
}

ukf.construct.data.updater <- function(self, private, scalar.or.joint=c("joint","scalar")){
  scalar.or.joint <- match.arg(scalar.or.joint)
  if(scalar.or.joint=="joint"){
    return(ukf.construct.data.updater.joint(self, private))
  } else {
    return(ukf.construct.data.updater.scalar(self, private))
  }
}

# setup function
# ------------------------------------
ukf.setup <- function(self, private){

  w <- ukf.compute.weights(self, private)

  settings <- list(
    first.order.hold = private$algo.settings$first.order.input.hold,
    estimate.initial.state = private$algo.settings$estimate.initial,
    n.states = w$n.states,
    n.sigmapoints = w$n.sigmapoints,
    sqrt_c = w$sqrt_c,
    W.m = w$W.m,
    W = w$W
  )

  advanced.settings <- list(
    nllreport = private$algo.settings$advanced.settings$nllreport,
    force.ad.on.input = private$algo.settings$advanced.settings$force.ad.on.input
  )

  data <- list(
    x = private$algo.settings$initial.state$x0,
    P = private$algo.settings$initial.state$p0,
    U = as.matrix(private$data[private$names$inputs]),
    Y = as.matrix(private$data[private$names$obs]),
    ode.stepsizes = private$algo.settings$ode.stepsizes,
    ode.number.of.steps = private$algo.settings$ode.number.of.steps
  )
  data$dU <- with(c(data,settings),
                  if(first.order.hold) rbind(diff(U)/ode.number.of.steps,0) else array(0,dim=dim(U)))

  algorithms <- list(
    time.update = ukf.construct.time.updater(self, private),
    # TO DO: add option to choose between joint and scalar update
    data.update = ukf.construct.data.updater(self, private, scalar.or.joint="joint"),
    stationary.solver = construct.initial.state.estimators(self, private)
  )

  # return
  return(named_list(data,settings,advanced.settings,algorithms))
}

# likelihood constructor
# ------------------------------------
MakeADFun_UKF = function(self, private) {

  # Settings
  configure_ad_tape("RTMB", self, private)
  setup <- ukf.setup(self, private)
  list2env(setup$data, envir=environment())
  list2env(setup$settings, envir=environment())
  list2env(setup$advanced.settings, envir=environment())
  list2env(setup$algorithms, envir=environment())

  # Sample sigma-points from mean and sqrt-covariance
  # [x, x, ..., x] + sqrt(c) [0, chol(P), -chol(P)]
  create.sigmaPoints <- function(x, sqrtP){
    xs <- rep(x, n.sigmapoints)
    dim(xs) <- c(n.states, n.sigmapoints)
    z <- sqrt_c * c(RTMB::AD(numeric(n.states)), sqrtP, -sqrtP)
    dim(z) <- c(n.states, n.sigmapoints)
    xs + z
  }

  # likelihood function --------------------------------------
  ukf.nll = function(p){

    ####### Sometimes necessary to avoid rtmb errors #######
    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")
    "c" <- RTMB::ADoverload("c")

    # Convert parameter list to vector
    p <- do.call(c, p[1:private$dims$pars])

    # Set likelihood to zero
    nll <- 0

    # Cast initial state and inputs to AD
    x <- RTMB::AD(x)
    P <- RTMB::AD(P)
    if (force.ad.on.input) U <- RTMB::AD(U)

    epsI0 <- 1e-8 * diag(n.states)

    # Estimate initial conditions as stationary solution with first input
    # This estimate is, as with the user provided initial state, assumed to be
    # the prior mean and covariance x_{0 | -1} and P_{0 | -1}
    if (any(estimate.initial.state)) {
      get.cov <- estimate.initial.state[2]
      sol <- stationary.solver(p,U[1,],get.cov=get.cov,get.mean=estimate.initial.state[1],x0=x)
      x <- sol[,1]
      if(get.cov) P <- sol[,-1, drop=FALSE]
    }

    # First data update
    sqrtP <- t(Matrix::chol(P + epsI0))
    X <- create.sigmaPoints(x, sqrtP)
    post <- data.update(x,P,X,p,U[1,],Y[1,])
    x <- post$x
    P <- post$P
    nll <- nll + post$nll

    for(i in seq_len(nrow(Y)-1)){
      ###### time update #######
      sqrtP <- t(Matrix::chol(P + epsI0))
      X <- create.sigmaPoints(x, sqrtP)
      sol <- time.update(X, sqrtP, p, U[i,], dU[i,], ode.number.of.steps[i], ode.stepsizes[i])
      X <- sol$X
      x <- X[,1]
      P <- sol$sqrtP %*% t(sol$sqrtP)
      ######## data update ########
      post <- data.update(x,P,X,p,U[i+1,],Y[i+1,])
      x <- post$x
      P <- post$P
      nll <- nll + post$nll
    }

    # return
    return(nll)
  }

  # construct AD-likelihood function ----------------------------------------
  map <- lapply(private$model$fixed.pars, function(x) x$factor)
  parameters <- lapply(private$model$parameters, function(x) x$initial)
  nll <- RTMB::MakeADFun(func = ukf.nll,
                         parameters = parameters,
                         map = map,
                         silent=TRUE)

  # save objective function
  private$nll <- nll

  # return
  return(invisible(self))
}
