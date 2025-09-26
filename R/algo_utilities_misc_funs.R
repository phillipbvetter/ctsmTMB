ConfigureTape <- function(context, self, private){
  
  if(context=="RTMB"){
    if(is.null(private$rtmb.tapeconfig)){
      # Default settings
      # atomic: disabling yields 7x optimization speed
      x <- list(atomic="disable", comparison="NA", vectorize="NA")
      do.call(RTMB::TapeConfig, x)
    } else {
      # User-provided settings
      x <- as.list(private$rtmb.tapeconfig)
      nams <- names(formals(RTMB::TapeConfig))
      arg.names <- nams[nams %in% names(x)]
      x <- x[arg.names]
      do.call(RTMB::TapeConfig, x)
    }
  }
  
  if(context=="TMB"){
    tmb.config.names <- c("trace.atomic",
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
                          "tmbad.sparse_hessian_compress")
    if(is.null(private$tmb.tapeconfig)){
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
      x <- as.list(private$tmb.tapeconfig)
      arg.names <- tmb.config.names[tmb.config.names %in% names(x)]
      x <- x[arg.names]
      do.call(TMB::config, c(x, list(DLL=private$modelname.with.method)))
    }
  }
  
  return(invisible(self))
}
getSystemDimensions <- function(private, .envir=parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  n.sigmapoints <- 2*n.states + 1
  
  assign("n.states", n.states, envir = .envir)
  assign("n.obs", n.obs, envir = .envir)
  assign("n.pars", n.pars, envir = .envir)
  assign("n.diffusions", n.diffusions, envir = .envir)
  assign("n.inputs", n.inputs, envir = .envir)
  assign("n.sigmapoints", n.sigmapoints, envir = .envir)
  
  return(invisible(NULL))
}

getAdjoints <- function(.envir=parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  # adjoints ----------------------------------------
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
  
  kron.left <- RTMB::ADjoint(
    function(x) {
      dim(x) <- rep(sqrt(length(x)), 2)
      i <- diag(sqrt(length(x)))
      kronecker(x,i)
    },
    function(x, y, dy) {
      n <- sqrt(length(x))
      out <- matrix(0,nrow=n,ncol=n)
      for(i in 1:n){
        for(j in 1:n){
          id.seq <- 1+(j-1)*n + (n^2+1) * (1:n-1)
          id.seq <- id.seq + (i-1) * n^3
          out[[i,j]] <- sum(dy[id.seq])
        }
      }
      return(out)
    },
    name = "kron.left")
  
  kron.right <- RTMB::ADjoint(
    function(x) {
      dim(x) <- rep(sqrt(length(x)), 2)
      i <- diag(sqrt(length(x)))
      kronecker(i,x)
    },
    function(x, y, dy) {
      n <- sqrt(length(x))
      out <- matrix(0,nrow=n,ncol=n)
      for(i in 1:n){
        for(j in 1:n){
          id.seq <- j + (n^3+n) * (1:n-1)
          id.seq <- id.seq + (i-1) * n^2
          out[[i,j]] <- sum(dy[id.seq])
        }
      }
      return(out)
    },
    name = "kron.right")
  
  assign("logdet", logdet, envir = .envir)
  assign("kron.left", kron.left, envir = .envir)
  assign("kron.right", kron.right, envir = .envir)
  
  return(invisible(NULL))
}

getLossFunction <- function(.envir=parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  # Loss function ----------------------------------------a
  # quadratic loss
  if(private$loss$loss == "quadratic"){
    loss.function = function(e,R) -RTMB::dmvnorm(e, Sigma=R, log=TRUE)
  }
  # huber loss
  if(private$loss$loss == "huber"){
    loss.c <- private$loss$c
    k.smooth <- 5
    sigmoid <- function(r_sqr) 1/(1+exp(-k.smooth*(sqrt(r_sqr)-loss.c)))
    huber.loss <- function(r_sqr) {
      s <- sigmoid(r_sqr)
      r_sqr * (1-s) + loss.c * (2*sqrt(r_sqr)-loss.c)*s
    }
    log2pi = log(2*pi)
    loss.function = function(e,R){
      r_squared <- t(e) %*% RTMB::solve(R) %*% e
      0.5 * logdet(R) + 0.5 * log2pi * length(e) + 0.5*huber.loss(r_squared)
    }
  }
  # tukey loss
  if(private$loss$loss == "tukey"){
    loss.c <- private$loss$c
    k.smooth <- 5
    sigmoid <- function(r_sqr) 1/(1+exp(-k.smooth*(sqrt(r_sqr)-loss.c)))
    tukey.loss <- function(r_sqr) {
      s <- sigmoid(r_sqr)
      r_sqr * (1-s) + loss.c^2*s
    }
    log2pi = log(2*pi)
    loss.function = function(e,R){
      r_squared <- t(e) %*% RTMB::solve(R) %*% e
      0.5 * logdet(R) + 0.5 * log2pi * length(e) + 0.5*tukey.loss(r_squared)
    }
  }
  
  assign("loss.function", loss.function, envir = .envir)
  
  return(NULL)
}

getOdeSolvers <- function(.envir=parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  # 1-step covariance ODE ----------------------------------------
  cov.ode.1step = function(covMat, stateVec, parVec, inputVec){
    G <- g__(stateVec, parVec, inputVec)
    AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  # forward euler ----------------------------------------
  if(private$ode.solver==1){
    ode.integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      X1 = stateVec + f__(stateVec, parVec, inputVec) * dt
      P1 = covMat + cov.ode.1step(covMat, stateVec, parVec, inputVec) * dt
      
      return(list(X1,P1))
    }
  } else if (private$ode.solver==2) {
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
  } else if (private$ode.solver==3){
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
  # else {
  
  #   # desolve ODE solver ----------------------------------------
  #   
  #   # Step 1: construct input interpolators
  #   
  #   # Expand time-domain slightly to avoid NAs in RTMB::interpol1Dfun
  #   # if ODE is evaluated slightly outside domain
  #   time <- inputMat[,1]
  #   time.range <- range(time)
  #   time.range.diff <- diff(time.range)
  #   new.range = time.range + rep(time.range.diff/10, 2) * c(-1,1)
  #   new.time <- seq(new.range[1], new.range[2], by=min(diff(time)))
  #   input.interp.funs <- vector("list",length=n.inputs)
  #   # Create interpolation on uniform grid
  #   for(i in 1:length(input.interp.funs)){
  #     # create uniform time grid
  #     out <- stats::approx(x=time, y=inputMat[,i], xout=new.time, rule=2)
  #     # interpol1Dfun assumes uniform distances in x-coordinates
  #     input.interp.funs[[i]] <- RTMB::interpol1Dfun(z=out$y, xlim=range(new.time), R=1)
  #   }
  #   
  #   # Step 2: construct ode fun for DeSolve
  #   ode.fun <- function(time, stateVec_an d_covMat, parVec){
  #     inputVec <- RTMB::sapply(input.interp.funs, function(f) f(time))
  #     # 
  #     stateVec <- head(stateVec_and_covMat, n.states)
  #     covMat <- RTMB::matrix(tail(stateVec_and_covMat, -n.states),nrow=n.states)
  #     # 
  #     G <- g__(stateVec, parVec, inputVec)
  #     AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
  #     #
  #     dX <- f__(stateVec, parVec, inputVec)
  #     dP <- AcovMat + t(AcovMat) + G %*% t(G)
  #     return(list(c(dX,dP)))
  #   }
  #   
  #   # Step 3: construct function to call in likelihood functon
  #   ode_integrator <- function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
  #     out <- RTMBode::ode(y = c(stateVec, covMat),
  #                         times = c(inputVec[1], inputVec[1]+dt),
  #                         func = ode.fun,
  #                         parms = parVec,
  #                         method=ode.solver)[2,-1]
  #     
  #     return(
  #       list(head(out,n.states),
  #            RTMB::matrix(tail(out,-n.states),nrow=n.states)
  #       )
  #     )
  #   }
  # }
  
  assign("cov.ode.1step", cov.ode.1step, envir=.envir)
  assign("ode.integrator", ode.integrator, envir = .envir)
  
  return(NULL)
}

getInitialStateEstimator <- function(.envir=parent.frame()){
  
  # unpack objects from parent
  list2env(as.list(.envir), envir = environment())
  
  ######## Function 1
  ######## Mean Stationary Solver
  
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
    
    # Root-find drift function (stationary 1st moment ODE solution)
    X1 <- f__(stateVec_rootfind, parVec, inputVec)
    
    cost <- sum(X1*X1)
    return(cost)
  }
  
  # create the tape and perform newton
  f.initial.state.newton0 <- RTMB::MakeTape(initial.state.newton, rep(1,sum(sizes)))
  f.initial.state.newton <- f.initial.state.newton0$newton(1:n.states)
  
  ######## Function 2
  ######## Covariance Stationary Solve
  # TODO: Can we implement Bartels–Stewart algorithm with eigen?
  f.initial.covar.solve <- function(stateVec, parVec, inputVec){
    A <- dfdx__(stateVec, parVec, inputVec)
    G <- g__(stateVec, parVec, inputVec)
    Q <- G %*% t(G)
    # this line below causes hessian to crash because of the RTMB::ADjoints
    # kron.left and kron.right
    P <- kron.left(A) + kron.right(A)
    X <- -RTMB::solve(P, as.numeric(Q))
    covMat <- RTMB::matrix(X, nrow=n.states)
    return(covMat)
  }
  
  assign("f.initial.state.newton", f.initial.state.newton, envir=.envir)
  assign("f.initial.covar.solve", f.initial.covar.solve, envir=.envir)
  
  return(NULL)
}

getUkfOdeSolvers <- function(.envir=parent.frame()){
  
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
  f.sigma <- function(sigmaPoints, parVec, inputVec){
    # Fsigma <- RTMB::matrix(0,nrow=n.states,ncol=nn)
    for(i in 1:n.sigmapoints){
      # Fsigma[,i] <- f__(sigmaPoints[,i], parVec, inputVec)
      Fsigma[[1,i]] <- f__(sigmaPoints[,i], parVec, inputVec)
    }
    x <- do.call("c", Fsigma)
    dim(x) <- c(n.states, n.sigmapoints)
    return(x)
    # return(Fsigma)
  }
  
  h.sigma <- function(sigmaPoints, parVec, inputVec){
    # Hsigma <- RTMB::matrix(0,nrow=n.obs,ncol=nn)
    for(i in 1:n.sigmapoints){
      # Hsigma[,i] <- h__(sigmaPoints[,i], parVec, inputVec)
      Hsigma[[1,i]] <- h__(sigmaPoints[,i], parVec, inputVec)
    }
    x <- do.call("c", Hsigma)
    dim(x) <- c(n.obs, n.sigmapoints)
    return(x)
    # return(Hsigma)
  }
  
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
  if(private$ode.solver==1){
    ode.integrator = function(X.sigma, chol.covMat, parVec, inputVec, dinputVec, dt){
      
      X1 <- X.sigma + cov.ode.1step(X.sigma, chol.covMat, parVec, inputVec) * dt
      
      return(X1)
    }
  }
  
  # rk4 ----------------------------------------
  if(private$ode.solver==2){
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
  if(private$ode.solver==3){
    
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
  
  assign("create.sigmaPoints", create.sigmaPoints, envir=.envir)
  assign("sigma2chol", sigma2chol, envir=.envir)
  assign("f.sigma", f.sigma, envir=.envir)
  assign("h.sigma", h.sigma, envir=.envir)
  assign("cov.ode.1step", cov.ode.1step, envir=.envir)
  assign("ode.integrator", ode.integrator, envir = .envir)
  
  return(NULL)
}

getUkfSigmaWeights <- function(.envir=parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  # grab hyperparameters
  ukf.alpha <- private$ukf_hyperpars[1]
  ukf.beta <- private$ukf_hyperpars[2]
  ukf.kappa <- private$ukf_hyperpars[3]
  
  sqrt_c <- sqrt(ukf.alpha^2*(n.states + ukf.kappa))
  ukf.lambda <- sqrt_c^2 - n.states
  
  # weights
  W.m <- W.c <- rep(1/(2*(n.states+ukf.lambda)), n.sigmapoints)
  W.m[1] <- ukf.lambda/(n.states+ukf.lambda)
  W.c[1] <- ukf.lambda/((n.states+ukf.lambda)+(1-ukf.alpha^2+ukf.beta))
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



