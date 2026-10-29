# integrator for time update
# ------------------------------------
ekf.construct.ode.integrator <- function(self, private){
  solver <- private$algo.settings$ode.solver
  n.states <- private$dims$states
  n.pars   <- private$dims$pars
  n.inputs <- private$dims$inputs
  funs <- private$model$state.space.funs

  #moment odes
  dm_dt <- function(x, P, p, u){
    f <- funs$f(x, p, u)
    G <- funs$g(x, p, u)
    AP <- funs$dfdx(x, p, u) %*% P
    list(
      dx = f,
      dP = AP + t(AP) + G %*% t(G)
    )
  }

  ode.integrator <- switch(
    solver,
    euler = function(x, P, p, u, du, dt){
      k <- dm_dt(x, P, p, u)
      list(
        x = x + k$dx * dt,
        P = P + k$dP * dt
      )
    },
    rk4 = function(x, P, p, u, du, dt){
      k1 <- dm_dt(x, P, p, u)
      k2 <- dm_dt(x + k1$dx * dt/2, P + k1$dP * dt/2, p, u + du/2)
      k3 <- dm_dt(x + k2$dx * dt/2, P + k2$dP * dt/2, p, u + du/2)
      k4 <- dm_dt(x + k3$dx * dt, P + k3$dP * dt, p, u + du)
      list(
        x = x + (k1$dx + 2*k2$dx + 2*k3$dx + k4$dx)/6 * dt,
        P = P + (k1$dP + 2*k2$dP + 2*k3$dP + k4$dP)/6 * dt
      )
    },
    stop("Unknown solver - choose 'euler' or 'rk4'.")
  )

  # # tape the ode integrator (3 times faster adfun compilation)
  arglist <- list(
    x = numeric(n.states),
    P = array(0, dim=rep(n.states,2)),
    p = numeric(n.pars),
    u = numeric(n.inputs),
    du = numeric(n.inputs),
    dt = numeric(1)
  )

  taped.integrator <- RTMB::MakeTape(function(arglist) {
    do.call(cbind, do.call(ode.integrator, arglist))
  }, arglist)

  time.update <- function(x,P,p,u,du,steps,dt){
    for (j in seq_len(steps)) {
      sol <- taped.integrator(list(x, P, p, u, du, dt))
      x = sol[,1]
      P = sol[,-1,drop=FALSE]
      u <- u + du
    }
    list(x=x,P=P)
  }

  return(time.update)
}

# data update step
# ------------------------------------
ekf.construct.data.updater <- function(self, private){

  # Unpack dimensions
  n.states <- private$dims$states
  n.pars   <- private$dims$pars
  n.inputs <- private$dims$inputs

  # Loss function
  loss.type <- private$algo.settings$loss$loss
  loss.c <- private$algo.settings$loss$c
  loss.function <- construct.loss.function(dim="scalar",type=loss.type, c=loss.c)

  # Taped observation functions
  funs <- private$model$state.space.funs
  obs.fun.tape <- RTMB::MakeTape(function(x) {
    do.call(cbind, lapply(funs[c("h","hvar","dhdx")], \(f) do.call(f,x)))
  }, list(numeric(n.states), numeric(n.pars), numeric(n.inputs)))


  I0 <- diag(n.states)
  # Taped sequential EKF update
  scalar.update.tape <- RTMB::MakeTape(function(z) {
    # Unpack list
    x <- z[[1]]
    P <- z[[2]]
    x0 <- z[[3]]
    C <- RTMB::matrix(z[[4]], nrow = 1)
    Hj <- z[[5]]
    vj <- z[[6]]
    yj <- z[[7]]
    # Linearized innovation (h(x) \approx h(x*) + dh/dx(x*) * (x-x*))
    e <- yj - Hj - drop(C %*% (x - x0))
    PCt <- P %*% t(C)
    R <- drop(C %*% PCt) + vj
    K <- PCt / R
    # State Update
    x <- x + K * e
    # Cov Update
    # Standard Update Form
    # P <- P - K %*% t(PCt)
    # Joseph Form
    IKC <- I0 - K %*% C
    P <- IKC %*% P %*% t(IKC) + vj * K %*% t(K)
    # Resymmetrize
    P <- 0.5 * (P + t(P))
    # Likelihood
    nll <- loss.function(e,R)
    # Return
    # It is faster to return a matrix for subsetting than to return a vector with c()
    # so we just repeat nll to match the rows of stateVec/covMat
    cbind(x, P, rep(nll, n.states))
  }, list(numeric(n.states), diag(n.states), numeric(n.states), numeric(n.states), 0, 1, 0))

  sequential.update <- function(x, P, p, u, y, y_bool){
    nll <- 0
    idx <- which(y_bool)
    x_prior <- x
    obs.funs <- obs.fun.tape(list(x, p, u))
    H <- obs.funs[,1]
    V <- obs.funs[,2]
    C <- obs.funs[, 2+1:n.states, drop=FALSE]
    for (j in idx) {
      out <- scalar.update.tape(list(x, P, x_prior, C[j,,drop=FALSE], H[j], V[j], y[j]))
      x <- out[,1]
      P <- out[,1 + 1:n.states, drop=FALSE]
      nll <- nll + out[1, 2+n.states]
    }
    return(list(x=x, P=P, nll=nll))
  }

  scalar.no.update.tape <- RTMB::MakeTape(function(z) {
    # Unpack list
    Hj <- z[[1]]
    vj <- z[[2]]
    yj <- z[[3]]
    e <- yj - Hj
    R <- vj
    # Likelihood
    nll <- loss.function(e,R)
    # Return
    return(nll)
  }, list(1,1,1))

  sequential.no.update <- function(x, P, p, u, y, y_bool){
    nll <- 0
    idx <- which(y_bool)
    obs.funs <- obs.fun.tape(list(x, p, u))
    H <- obs.funs[,1]
    V <- obs.funs[,2]
    for (j in idx) {
      out <- scalar.no.update.tape(list(H[j], V[j], y[j]))
      nll <- nll + out
    }
    return(list(x=x, P=P, nll=nll))
  }

  create_ekf_update_function <- function(data.update.fun){
    f <- function(x,P,p,u,y){
      y_bool = !is.na(y)
      if (any(y_bool)) {
        out <- data.update.fun(x, P, p, u, y, y_bool)
      } else {
        out <- list(x=x,P=P,nll=0)
      }
      return(out)
    }
    return(f)
  }


  if(private$algo.settings$train.against.full.prediction){
    update.fun <- create_ekf_update_function(sequential.no.update)
  } else {
    update.fun <- create_ekf_update_function(sequential.update)
  }

  return(update.fun)
}

# setup function
# ------------------------------------
ekf.setup <- function(self, private){

  settings <- list(
    first.order.hold = private$algo.settings$first.order.input.hold,
    estimate.initial.state = private$algo.settings$estimate.initial
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
    time.update = ekf.construct.ode.integrator(self, private),
    data.update = ekf.construct.data.updater(self, private),
    stationary.solver = construct.initial.state.estimators(self, private)
  )

  # return
  return(named_list(data,settings,advanced.settings,algorithms))
}

# likelihood constructor
# ------------------------------------
MakeADFun_EKF = function(self, private) {

  # Settings
  configure_ad_tape("RTMB", self, private)
  setup <- ekf.setup(self, private)
  list2env(setup$data, envir=environment())
  list2env(setup$settings, envir=environment())
  list2env(setup$advanced.settings, envir=environment())
  list2env(setup$algorithms, envir=environment())

  # likelihood function --------------------------------------
  ekf.nll = function(p){

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
    post <- data.update(x,P,p,U[1,],Y[1,])
    x <- post$x
    P <- post$P
    nll <- nll + post$nll

    for(i in seq_len(nrow(Y)-1)){
      ###### time update #######
      sol <- time.update(x,P,p,U[i,],dU[i,],ode.number.of.steps[i],ode.stepsizes[i])
      x <- sol$x; P <- sol$P
      ######## data update ########
      post <- data.update(x,P,p,U[i+1,],Y[i+1,])
      x <- post$x; P <- post$P
      nll <- nll + post$nll
    }

    # return
    return(nll)
  }

  # construct AD-likelihood function ----------------------------------------
  map <- lapply(private$model$fixed.pars, function(x) x$factor)
  parameters <- lapply(private$model$parameters, function(x) x$initial)
  nll <- RTMB::MakeADFun(func = ekf.nll,
                         parameters = parameters,
                         map = map,
                         silent=TRUE)

  # save objective function
  private$nll <- nll


  # # MAP Estimation?
  # MAP_bool = 0L
  # if (!is.null(private$algo.settings$map)) {
  #   bool = self$getParameters()[,"type"] == "free"
  #   MAP_bool = 1L
  #   map_mean__ = private$algo.settings$map$mean[bool]
  #   map_cov__ = private$algo.settings$map$cov[bool,bool]
  #   map_ints__ = as.numeric(bool)
  #   sum_map_ints__ = sum(as.numeric(bool))
  # }

  # return
  return(invisible(self))
}
