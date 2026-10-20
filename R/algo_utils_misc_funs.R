configure_ad_tape <- function(context, self, private){

  if (context=="RTMB"){
    if(is.null(private$algo.settings$advanced.settings$rtmb.tapeconfig)){
      # Default settings
      # atomic: disabling yields 7x optimization speed
      x <- list(atomic="disable", comparison="NA", vectorize="NA")
      do.call(RTMB::TapeConfig, x)
    } else {
      # User-provided settings
      x <- as.list(private$algo.settings$advanced.settings$rtmb.tapeconfig)
      nams <- names(formals(RTMB::TapeConfig))
      arg.names <- nams[nams %in% names(x)]
      x <- x[arg.names]
      do.call(RTMB::TapeConfig, x)
    }
  }

  if(context=="TMB"){
    tmb.config.names <- c(
      "trace.atomic",
      "trace.optimize",
      "nthreads",
      "trace.parallel",
      "tmbad.atomic_sparse_log_determinant",
      "tmbad_deterministic_hash",
      "tape.parallel",
      "optimize.parallel",
      "autopar",
      "optimize.instantly",
      "debug.getListElement",
      "tmbad.sparse_hessian_compress"
    )
    if (is.null(private$algo.settings$advanced.settings$tmb.tapeconfig)) {
      # Default settings
      x <- list(
        trace.atomic = 1,
        trace.optimize = 1,
        nthreads = 1,
        trace.parallel = 1,
        tmbad.atomic_sparse_log_determinant = 1,
        tmbad_deterministic_hash = 1,
        tape.parallel = 1,
        optimize.parallel = 0,
        autopar = 0,
        optimize.instantly = 1,
        debug.getListElement = 0,
        tmbad.sparse_hessian_compress = 0
      )
      do.call(TMB::config, c(x, list(DLL=private$modelname.with.method)))
    } else {
      # User-provided settings
      x <- as.list(private$algo.settings$advanced.settings$tmb.tapeconfig)
      arg.names <- tmb.config.names[tmb.config.names %in% names(x)]
      x <- x[arg.names]
      do.call(TMB::config, c(x, list(DLL=private$modelname.with.method)))
    }
  }

  return(invisible(self))
}

unpack_data <- function(self, private, force.ad=TRUE, .envir=parent.frame()){

  e <- new.env()
  local({

    # initial
    stateVec = private$algo.settings$initial.state$x0
    covMat = private$algo.settings$initial.state$p0

    # input and obs matrix
    inputMat = as.matrix(private$data[private$names$inputs])
    obsMat = as.matrix(private$data[private$names$obs])

    # force ad on the data
    if(force.ad){
      stateVec <- RTMB::AD(stateVec, force=TRUE)
      covMat <- RTMB::AD(covMat, force=TRUE)
      if (private$algo.settings$advanced.settings$force.ad.on.input){
        inputMat <- RTMB::AD(inputMat,force=TRUE)
      }
    }

    ode.stepsizes = private$algo.settings$ode.stepsizes
    ode.number.of.steps = private$algo.settings$ode.number.of.steps

  },
  envir=e)

  # Unpack into the caller environment
  list2env(as.list(e, all.names=TRUE), envir=.envir)
  return(invisible(self))
}

unpack_advanced_settings <- function(self, private, .envir=parent.frame()){

  e <- new.env()
  local({
    force.ad <- private$algo.settings$advanced.settings$forceAD
    nllreport <- private$algo.settings$advanced.settings$nllreport
  },
  envir=e)

  # Unpack into the caller environment
  list2env(as.list(e, all.names=TRUE), envir=.envir)
  return(invisible(self))
}

unpack_system_dimensions <- function(self, private, .envir=parent.frame()){

  e <- new.env()
  local({
    n.states <- private$dims$states
    n.obs <- private$dims$observations
    n.pars <- private$dims$pars
    n.diffusions <- private$dims$diffusions
    n.inputs <- private$dims$inputs
    n.sigmapoints <- 2*n.states + 1
  },
  envir=e)

  # Unpack into the caller environment
  list2env(as.list(e, all.names=TRUE), envir=.envir)

  return(invisible(e))
}

unpack_adjoints <- function(.envir=parent.frame()){

  e <- new.env()
  local({

    # log determinant ----------------------------------------
    # f(X(p)) = log det (X(p))
    # df = Tr[X^{-1} dX]
    logdet <- RTMB::ADjoint(
      function(x) {
        dim(x) <- rep(sqrt(length(x)), 2)
        log(det(x))
      },
      function(x, y, dy) {
        dim(x) <- rep(sqrt(length(x)), 2)
        t(RTMB::solve(x)) * dy
      },
      name = "logdet")

  },
  envir=e)

  # Unpack into the caller environment
  list2env(as.list(e, all.names=TRUE), envir=.envir)

  return(invisible(NULL))
}

get_loss_function <- function(type=c("multivariate","scalar"),.envir=parent.frame()){

  type <- match.arg(type)

  list2env(as.list(.envir), envir = environment())

  # Loss function ----------------------------------------
  if(type=="multivariate"){
    # quadratic loss
    if(private$algo.settings$loss$loss == "quadratic"){
      loss.function = function(e,R) -RTMB::dmvnorm(e, Sigma=R, log=TRUE)
    }
    # huber loss
    if(private$algo.settings$loss$loss == "huber"){
      loss.c <- private$algo.settings$loss$c
      k.smooth <- 5
      sigmoid <- function(r_sqr) 1/(1+exp(-k.smooth*(sqrt(r_sqr)-loss.c)))
      huber.loss <- function(r_sqr) {
        s <- sigmoid(r_sqr)
        r_sqr * (1-s) + loss.c * (2*sqrt(r_sqr)-loss.c)*s
      }
      log2pi = log(2*pi)
      loss.function = function(e,R){
        r_squared <- t(e) %*% RTMB::solve(R,e)
        0.5 * logdet(R) + 0.5 * log2pi * length(e) + 0.5*huber.loss(r_squared)
      }
    }
    # tukey loss
    if(private$algo.settings$loss$loss == "tukey"){
      loss.c <- private$algo.settings$loss$c
      k.smooth <- 5
      sigmoid <- function(r_sqr) 1/(1+exp(-k.smooth*(sqrt(r_sqr)-loss.c)))
      tukey.loss <- function(r_sqr) {
        s <- sigmoid(r_sqr)
        r_sqr * (1-s) + loss.c^2*s
      }
      log2pi = log(2*pi)
      loss.function = function(e,R){
        r_squared <- t(e) %*% RTMB::solve(R,e)
        0.5 * logdet(R) + 0.5 * log2pi * length(e) + 0.5*tukey.loss(r_squared)
      }
    }
  } else if (type=="scalar") {
    # quadratic loss
    if(private$algo.settings$loss$loss == "quadratic"){
      loss.function = function(e,R) -RTMB::dnorm(e, sd=sqrt(R), log=TRUE)
    }
  } else {
    stop("Please select a loss function")
  }

  assign("loss.function", loss.function, envir = .envir)

  return(NULL)
}

get_ode_solvers <- function(.envir=parent.frame()){

  list2env(as.list(.envir), envir = environment())

  # 1-step covariance ODE ----------------------------------------
  cov.ode.1step = function(covMat, stateVec, parVec, inputVec){
    G <- g__(stateVec, parVec, inputVec)
    AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }

  # forward euler ----------------------------------------
  if(private$algo.settings$ode.solver=="euler"){
    ode.integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){

      X1 = stateVec + f__(stateVec, parVec, inputVec) * dt
      P1 = covMat + cov.ode.1step(covMat, stateVec, parVec, inputVec) * dt

      return(list(X1,P1))
    }
  } else if (private$algo.settings$ode.solver=="rk4") {
    # rk4 ----------------------------------------
    ode.integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){

      # Initials
      X0 = stateVec
      P0 = covMat

      # Classical 4th Order Runge-Kutta Method
      # 1. Approx Slope at Initial Point
      k1 = f__(stateVec, parVec, inputVec)
      c1 = cov.ode.1step(covMat, stateVec, parVec, inputVec)

      # 2. First Approx Slope at Midpoint
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + 0.5 * dt * k1
      covMat   = P0 + 0.5 * dt * c1
      k2       = f__(stateVec, parVec, inputVec)
      c2       = cov.ode.1step(covMat, stateVec, parVec, inputVec)

      # 3. Second Approx Slope at Midpoint
      stateVec = X0 + 0.5 * dt * k2
      covMat   = P0 + 0.5 * dt * c2
      k3       = f__(stateVec, parVec, inputVec)
      c3       = cov.ode.1step(covMat, stateVec, parVec, inputVec)

      # 4. Approx Slope at End Point
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + dt * k3
      covMat   = P0 + dt * c3
      k4       = f__(stateVec, parVec, inputVec)
      c4       = cov.ode.1step(covMat, stateVec, parVec, inputVec)

      # ODE UPDATE
      X1 = X0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0 * dt
      P1 = P0 + (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0 * dt

      return(list(X1,P1))
    }
  } else if (private$algo.settings$ode.solver=="implicit_euler"){
    #### Implicit Euler ####

    # id helpers
    sizes <- c(n.states, n.states^2, n.pars, n.inputs, n.inputs, 1, n.states, n.states^2)
    id.end <- cumsum(sizes)
    id.start <- c(1, head(id.end,-1) + 1)
    sub.ids <- mapply(function(s,e) s:e, id.start, id.end)

    # function for MakeTape
    maketape_ode_integrator <- function(x){
      # extract elements of x
      stateVec <- x[sub.ids[[1]]]
      covMat <- RTMB::matrix(x[sub.ids[[2]]], nrow=n.states)
      parVec <- x[sub.ids[[3]]]
      inputVec <- x[sub.ids[[4]]]
      dinputVec <- x[sub.ids[[5]]]
      dt <- x[sub.ids[[6]]]
      stateVec_rootfind <- x[sub.ids[[7]]]
      covMat_rootfind <- RTMB::matrix(x[sub.ids[[8]]],nrow=n.states)

      # implicit euler step
      X1 = stateVec_rootfind - (stateVec + f__(stateVec_rootfind, parVec, inputVec) * dt)
      P1 = covMat_rootfind - (covMat + cov.ode.1step(covMat_rootfind, stateVec_rootfind, parVec, inputVec) * dt)

      # squared sum for minimization (newton)
      cost <- sum(X1*X1) + sum(c(P1*P1))
      return(cost)
    }

    # Create the tape
    x.length <- tail(id.end, 1)
    newton.ids <- c(sub.ids[[7]],sub.ids[[8]])
    ode_integrator0 <- RTMB::MakeTape(maketape_ode_integrator, rep(1,x.length))
    ode_integrator1 <- ode_integrator0$newton(newton.ids)

    # tape wrapper function
    ode.integrator <- function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      x <- c(stateVec, covMat, parVec, inputVec, dinputVec, dt)
      out <- ode_integrator1(x)
      X1 <- head(out, n.states)
      P1 <- RTMB::matrix(tail(out, n.states^2), nrow=n.states)
      return(list(X1,P1))
    }

  }

  # # tape the ode integrator (3 times faster adfun compilation)
  ode.integrator.taped <- RTMB::MakeTape(function(x) {
    sol <- ode.integrator(x[[2]], x[[1]], x[[3]], x[[4]], x[[5]], x[[6]])
    return(cbind(sol[[1]],sol[[2]]))
  }, list(numeric(n.states), array(0,dim=rep(n.states,2)), numeric(n.pars), numeric(n.inputs), numeric(n.inputs), numeric(1)))


  assign("cov.ode.1step", cov.ode.1step, envir=.envir)
  assign("ode.integrator", ode.integrator, envir = .envir)
  assign("ode.integrator.taped", ode.integrator.taped, envir = .envir)

  return(NULL)
}

get_initial_state_estimator_fixed_inputs <- function(inputVec, .envir=parent.frame()){

  # unpack objects from parent
  list2env(as.list(.envir), envir = environment())

  ###################################################
  # STATIONARY SOLVER FOR THE MEAN
  ###################################################

  # id helpers
  sizes <- c(n.states, n.pars)
  n.total <- sum(sizes)
  id.end <- cumsum(sizes)
  id.start <- c(1, head(id.end,-1) + 1)
  ids <- mapply(function(s,e) s:e, id.start, id.end)


  # maketape function
  initial.state.newton <- function(x){
    stateVec <- x[ids[[1]]]
    parVec <- x[ids[[2]]]
    y <- f__(stateVec, parVec, inputVec)
    sum(y*y)
  }
  # create the tape and perform newton
  # f.initial.state.newton0 <- RTMB::MakeTape(initial.state.newton, numeric(n.total))
  f.initial.state.newton0 <- RTMB::MakeTape(initial.state.newton, rep(1,n.total))
  f.initial.state.newton <- f.initial.state.newton0$newton(1:n.states, grad_tol=1e-20)

  assign("f.initial.state.newton_fixed_input", f.initial.state.newton, envir=.envir)
}

get_initial_state_estimator <- function(.envir=parent.frame()){

  # unpack objects from parent
  list2env(as.list(.envir), envir = environment())

  ###################################################
  # STATIONARY SOLVER FOR THE MEAN
  ###################################################

  # id helpers
  sizes <- c(n.states, n.pars, n.inputs)
  id.end <- cumsum(sizes)
  id.start <- c(1, head(id.end,-1) + 1)
  sub.ids <- mapply(function(s,e) s:e, id.start, id.end)

  # maketape function
  initial.state.newton <- function(x){
    stateVec_rootfind <- x[sub.ids[[1]]]
    parVec <- x[sub.ids[[2]]]
    inputVec <- x[sub.ids[[3]]]
    dx <- f__(stateVec_rootfind, parVec, inputVec)
    sum(dx*dx)
  }

  # create the tape and perform newton
  f.initial.state.newton0 <- RTMB::MakeTape(initial.state.newton, rep(1,sum(sizes)))
  f.initial.state.newton <- f.initial.state.newton0$newton(1:n.states, grad_tol=1e-20)

  ###################################################
  # STATIONARY SOLVER FOR THE COVARIANCE
  ###################################################

  ######## Function 2
  ######## Covariance Stationary Solve
  # TODO: Can we implement Bartels–Stewart algorithm with eigen?
  # f.initial.covar.solve <- function(stateVec, parVec, inputVec){
  #   A <- dfdx__(stateVec, parVec, inputVec)
  #   G <- g__(stateVec, parVec, inputVec)
  #   Q <- G %*% t(G)
  #   # this line below causes hessian to NaN because of the RTMB::ADjoints
  #   # kron.left and kron.right may be wrongly specified? not sure
  #   P <- kron.left(A) + kron.right(A)
  #   X <- -RTMB::solve(P, as.numeric(Q))
  #   covMat <- RTMB::matrix(X, nrow=n.states)
  #   return(covMat)
  # }


  ###################################################
  # ASSIGN TO CALLING ENVIRONMENT
  ###################################################

  assign("f.initial.state.newton", f.initial.state.newton, envir=.envir)
  # assign("f.initial.covar.solve", f.initial.covar.solve, envir=.envir)

  return(NULL)
}

get_ukf_utils <- function(.envir=parent.frame()){

  list2env(as.list(.envir), envir = environment())

  # Sample sigma-points from mean and sqrt-covariance
  # [X0, X0,...,X0] + sqrt(c) [0, chol(P), -chol(P)]
  create.sigmaPoints <- function(stateVec, chol.covMat){
    # Copy state vector
    x <- rep(stateVec, n.sigmapoints)
    dim(x) <- c(n.states, n.sigmapoints)

    # Compute sqrt(c) * [0, A, -A], with A = chol(P)
    n.zeros <- RTMB::AD(numeric(n.states))
    y <- sqrt_c * c(n.zeros, chol.covMat, -chol.covMat)
    dim(y) <- c(n.states, n.sigmapoints)

    X.sigma <- x+y
    return(X.sigma)
  }

  # get chol(P) from sigma points
  sigma2chol <- function(X.sigma){
    (X.sigma[,2:(n.states+1)] - X.sigma[,1])/sqrt_c
  }

  # create f of sigma points
  f.sigma <- function(X, parVec, inputVec){
    for(i in 1:n.sigmapoints){
      Fsigma[[1,i]] <- f__(X[,i], parVec, inputVec)
    }
    x <- do.call("c", Fsigma)
    dim(x) <- c(n.states, n.sigmapoints)
    return(x)
  }

  h.sigma <- function(X, parVec, inputVec){
    for(i in 1:n.sigmapoints){
      Hsigma[[1,i]] <- h__(X[,i], parVec, inputVec)
    }
    x <- do.call("c", Hsigma)
    dim(x) <- c(n.obs, n.sigmapoints)
    return(x)
  }

  assign("create.sigmaPoints", create.sigmaPoints, envir=.envir)
  assign("sigma2chol", sigma2chol, envir=.envir)
  assign("f.sigma", f.sigma, envir=.envir)
  assign("h.sigma", h.sigma, envir=.envir)

}

get_ukf_ode_solvers <- function(.envir=parent.frame()){

  list2env(as.list(.envir), envir = environment())

  # # Sample sigma-points from mean and sqrt-covariance
  # # [X0, X0,...,X0] + sqrt(c) [0, chol(P), -chol(P)]
  # create.sigmaPoints <- function(stateVec, chol.covMat){
  #   # Copy state vector
  #   x <- rep(stateVec, n.sigmapoints)
  #   dim(x) <- c(n.states, n.sigmapoints)
  #
  #   # Compute sqrt(c) * [0, A, -A], with A = chol(P)
  #   n.zeros <- RTMB::AD(numeric(n.states))
  #   y <- sqrt_c * c(n.zeros, chol.covMat, -chol.covMat)
  #   dim(y) <- c(n.states, n.sigmapoints)
  #
  #   X.sigma <- x+y
  #   return(X.sigma)
  # }
  #
  # # get chol(P) from sigma points
  # sigma2chol <- function(X.sigma){
  #   (X.sigma[,2:(n.states+1)] - X.sigma[,1])/sqrt_c
  # }
  #
  # # create f of sigma points
  # f.sigma <- function(X, parVec, inputVec){
  #   for(i in 1:n.sigmapoints){
  #     Fsigma[[1,i]] <- f__(X[,i], parVec, inputVec)
  #   }
  #   x <- do.call("c", Fsigma)
  #   dim(x) <- c(n.states, n.sigmapoints)
  #   return(x)
  # }
  #
  # h.sigma <- function(X, parVec, inputVec){
  #   for(i in 1:n.sigmapoints){
  #     Hsigma[[1,i]] <- h__(X[,i], parVec, inputVec)
  #   }
  #   x <- do.call("c", Hsigma)
  #   dim(x) <- c(n.obs, n.sigmapoints)
  #   return(x)
  #   # return(Hsigma)
  # }

  # 1-step UKF ODE ----------------------------------------
  cov.ode.1step = function(X.sigma, chol.covMat, parVec, inputVec){

    n.zerosAD <- RTMB::AD(numeric(n.states))

    # Create G without sigma points (just mean value)
    G <- g__(X.sigma[,1], parVec, inputVec)

    # Send sigma points through f
    F.sigma <- f.sigma(X.sigma, parVec, inputVec)

    # Rhs term 1:
    y <- rep(F.sigma %*% W.m, n.sigmapoints)
    dim(y) <- c(n.states, n.sigmapoints)

    # Rhs term 2: Compute [0 , chol(P) * Phi(M) , -chol(P)*Phi(M)]

    # Compute M and Phi(M)
    Ainv <- RTMB::solve(chol.covMat)
    Phi.M <- Ainv %*% (X.sigma %*% W %*% t(F.sigma) + F.sigma %*% W %*% t(X.sigma) + G %*% t(G) ) %*% t(Ainv)

    # Compute Phi(M) - set upper tri to zero, and divide diagonal by 2
    Phi.M[!lower.tri(Phi.M, diag=TRUE)] <- 0
    diag(Phi.M) <- diag(Phi.M)/2

    # Create [0, chol(P) * Phi(M), -chol(P) * Phi(M)]
    z0 <- chol.covMat %*% Phi.M
    z <- c(n.zerosAD, z0, -z0)
    dim(z) <- c(n.states, n.sigmapoints)

    # Create dX = RHS1 + sqrt(c) * RHS2
    dx <- y + sqrt_c * z

    # return
    return(dx)
  }

  # forward euler ----------------------------------------
  if(private$algo.settings$ode.solver=="euler"){
    ode.integrator = function(X.sigma, chol.covMat, parVec, inputVec, dinputVec, dt){

      X1 <- X.sigma + cov.ode.1step(X.sigma, chol.covMat, parVec, inputVec) * dt

      return(X1)
    }
  }

  # rk4 ----------------------------------------
  if(private$algo.settings$ode.solver=="rk4"){
    ode.integrator = function(X.sigma, chol.covMat, parVec, inputVec, dinputVec, dt){

      # Initials
      X0 <- X.sigma

      # Classical 4th Order Runge-Kutta Method
      # 1. Approx Slope at Initial Point
      k1 <- cov.ode.1step(X.sigma, chol.covMat, parVec, inputVec)

      # 2. First Approx Slope at Midpoint
      inputVec = inputVec + 0.5 * dinputVec
      X.sigma <- X0 + 0.5 * dt * k1
      chol.covMat <- sigma2chol(X.sigma)
      k2 <- cov.ode.1step(X.sigma, chol.covMat, parVec, inputVec)

      # 3. Second Approx Slope at Midpoint
      X.sigma = X0 + 0.5 * dt * k2
      chol.covMat <- sigma2chol(X.sigma)
      k3 <- cov.ode.1step(X.sigma, chol.covMat, parVec, inputVec)

      # 4. Approx Slope at End Point
      inputVec = inputVec + 0.5 * dinputVec
      X.sigma = X0 + dt * k3
      chol.covMat <- sigma2chol(X.sigma)
      k4 <- cov.ode.1step(X.sigma, chol.covMat, parVec, inputVec)

      # ODE UPDATE
      X1 = X0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0 * dt

      return(X1)
    }
  }
  # implicit euler
  if(private$algo.settings$ode.solver=="implicit_euler"){

    # # id helpers
    sizes <- c(n.states*n.sigmapoints, n.pars, n.inputs, n.inputs, 1, n.states*n.sigmapoints)
    id.end <- cumsum(sizes)
    id.start <- c(1, head(id.end,-1) + 1)
    sub.ids <- mapply(function(s,e) s:e, id.start, id.end)

    # function for MakeTape
    maketape_ode_integrator <- function(x){

      # extract elements
      Xsigma <- RTMB::matrix(x[sub.ids[[1]]], nrow=n.states)
      parVec <- x[sub.ids[[2]]]
      inputVec <- x[sub.ids[[3]]]
      dinputVec <- x[sub.ids[[4]]]
      dt <- x[sub.ids[[5]]]
      Xsigma.rootfind <- RTMB::matrix(x[sub.ids[[6]]],nrow=n.states)
      chol.covMat.rootfind <- sigma2chol(Xsigma.rootfind)

      # implicit euler step
      X1 <- Xsigma.rootfind - (Xsigma + cov.ode.1step(Xsigma.rootfind, chol.covMat.rootfind, parVec, inputVec) * dt)

      # squared sum for minimization (newton)
      cost <- sum(X1*X1)

      return(cost)
    }

    # Create the tape
    x.length <- tail(id.end, 1)
    newton.ids <- sub.ids[[6]]

    # The initial guess is important becauase ufe_ode_1step solves the chol.covMat.
    # The chol.covMat via sigma2chol therefore must be invertible
    initials <- rep(0.1, x.length - n.states*n.sigmapoints)
    initial.guess.sigmapoints <- c(cbind(0.1, diag(n.states), diag(n.states)))
    initials <- c(initials, initial.guess.sigmapoints)
    ode_integrator0 <- RTMB::MakeTape(maketape_ode_integrator, initials)

    # perform newton
    ode_integrator1 <- ode_integrator0$newton(newton.ids)

    # tape wrapper function
    ode.integrator = function(X.sigma, chol.covMat, parVec, inputVec, dinputVec, dt){
      x <- c(X.sigma, parVec, inputVec, dinputVec, dt)
      out <- ode_integrator1(x)
      X1 <- RTMB::matrix(out, nrow=n.states)
      return(X1)
    }
  }

  assign("cov.ode.1step", cov.ode.1step, envir=.envir)
  assign("ode.integrator", ode.integrator, envir = .envir)

  return(NULL)
}

get_ukf_weights <- function(.envir=parent.frame()){

  list2env(as.list(.envir), envir = environment())

  # grab hyperparameters
  ukf.alpha <- private$algo.settings$ukf.hyperpars[1]
  ukf.beta <- private$algo.settings$ukf.hyperpars[2]
  ukf.kappa <- private$algo.settings$ukf.hyperpars[3]

  sqrt_c <- sqrt(ukf.alpha^2*(n.states + ukf.kappa))
  ukf.lambda <- sqrt_c^2 - n.states

  # weights
  W.m <- W.c <- rep(1/(2*(n.states+ukf.lambda)), n.sigmapoints)
  W.m[1] <- ukf.lambda/(n.states+ukf.lambda)
  W.c[1] <- ukf.lambda/(n.states+ukf.lambda)+(1-ukf.alpha^2+ukf.beta)
  W.m.mat <- replicate(n.sigmapoints, W.m)
  W <- (diag(n.sigmapoints) - W.m.mat) %*% diag(W.c) %*% t(diag(n.sigmapoints) - W.m.mat)

  assign("sqrt_c", sqrt_c, envir = .envir)
  assign("W", W, envir = .envir)
  assign("W.m", W.m, envir = .envir)
  assign("W.c", W.c, envir = .envir)

}

# getEstimateInitialState

# getTemplate <- function(self, private, .envir=parent.frame()){
#   # unpack objects from parent
#   list2env(as.list(.envir), envir = environment())
#
#   # Define function here
#
#   assign("test",test, envir=.envir)
#
#   return(NULL)
# }

compute_initial_state_estimate <- function(pars, inputvec, self, private) {
  unpack_system_dimensions(self, private)
  apply_forceAD_to_some_variables(force.ad=FALSE)
  create_state_space_functions_for_estimation(force.ad=FALSE)
  get_initial_state_estimator()
  state.vec <- f.initial.state.newton(c(pars, inputvec))

  # stateVec <- try_withWarningRecovery(
  #   f.initial.state.newton(c(pars, inputvec))
  # )
  # if (inherits(stateVec, "try-error")) {
  #   message("The intial state estimation failed - using initial guess.")
  #   stateVec <- private$algo.settings$initial.state$x0
  # }

  state.vec
}

