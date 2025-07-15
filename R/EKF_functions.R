#######################################################
#######################################################
# EKF RTMB-IMPLEMENTATION (FOR OPTIMIZATION)
#######################################################
#######################################################

makeADFun_ekf_rtmb = function(self, private)
{
  
  # Tape Configration ----------------------
  # The best options for tape configuration was seen to be disabling atomic 
  # (x7 speed improvement of optimization) enabled by default, and keeping 
  # vectorized disabled
  RTMB::TapeConfig(atomic="disable")
  RTMB::TapeConfig(vectorize="disable")
  
  # Data ----------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # methods and purpose
  ode.solver = private$ode.solver
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # parameters ----------------------------------------
  parameters = lapply(private$parameters, function(x) x[["initial"]])
  
  # Loss function ----------------------------------------a
  # quadratic loss
  if(private$loss$loss == "quadratic"){
    loss_fun = function(e,R) -RTMB::dmvnorm(e, Sigma=R, log=TRUE)
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
    loss_fun = function(e,R){
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
    loss_fun = function(e,R){
      r_squared <- t(e) %*% RTMB::solve(R) %*% e
      0.5 * logdet(R) + 0.5 * log2pi * length(e) + 0.5*tukey.loss(r_squared)
    }
  }
  
  # # MAP Estimation?
  # MAP_bool = 0L
  # if (!is.null(private$map)) {
  #   bool = self$getParameters()[,"type"] == "free"
  #   MAP_bool = 1L
  #   map_mean__ = private$map$mean[bool]
  #   map_cov__ = private$map$cov[bool,bool]
  #   map_ints__ = as.numeric(bool)
  #   sum_map_ints__ = sum(as.numeric(bool))
  # }
  
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
  
  kron_left <- RTMB::ADjoint(
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
    name = "kron_left")
  
  kron_right <- RTMB::ADjoint(
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
    name = "kron_right")
  
  # 1-step covariance ODE ----------------------------------------
  cov_ode_1step = function(covMat, stateVec, parVec, inputVec){
    G <- g__(stateVec, parVec, inputVec)
    AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  # Implicit Euler Solver
  # ids.states <- 1:n.states
  # ids.cov <- max(ids.states) + 1:n.states^2
  # ids.pars <- max(ids.cov) + 1:n.pars
  # ids.inputs <- max(ids.pars) + 1:n.inputs
  # 
  # makeTape.implicit.euler <- function(x){
  #   stateVec <- x[ids.states]
  #   covMat <- RTMB::matrix(x[ids.cov],ncol=n.states)
  #   parVec <- x[ids.pars]
  #   inputVec <- x[ids.inputs]
  #   
  #   a <- stateVec - stateVec0 - f__(stateVec, parVec, inputVec) * dt
  #   b <- covMat - covMat0 - cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
  #   
  #   return(c(a,b))
  # }
  # RTMB::MakeTape(makeTape.implicit.euler, numeric())
  
  # forward euler ----------------------------------------
  if(ode.solver==1){
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      X1 = stateVec + f__(stateVec, parVec, inputVec) * dt
      P1 = covMat + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
      
      return(list(X1,P1))
    }
  } else if (ode.solver==2) {
    # rk4 ----------------------------------------
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      # Initials
      X0 = stateVec
      P0 = covMat
      
      # Classical 4th Order Runge-Kutta Method
      # 1. Approx Slope at Initial Point
      k1 = f__(stateVec, parVec, inputVec)
      c1 = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 2. First Approx Slope at Midpoint
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + 0.5 * dt * k1
      covMat   = P0 + 0.5 * dt * c1
      k2       = f__(stateVec, parVec, inputVec)
      c2       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 3. Second Approx Slope at Midpoint
      stateVec = X0 + 0.5 * dt * k2
      covMat   = P0 + 0.5 * dt * c2
      k3       = f__(stateVec, parVec, inputVec)
      c3       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 4. Approx Slope at End Point
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + dt * k3
      covMat   = P0 + dt * c3
      k4       = f__(stateVec, parVec, inputVec)
      c4       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # ODE UPDATE
      X1 = X0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0 * dt
      P1 = P0 + (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0 * dt
      
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
  
  # user-defined functions ---------------------------
  for(i in seq_along(private$rtmb.function.strings.indexed)){
    eval(parse(text=private$rtmb.function.strings.indexed[[i]]))
  }
  
  # new functions ----------------------------------------
  #error function
  erf = function(x){
    y <- sqrt(2) * x
    2*RTMB::pnorm(y)-1
  }
  
  # AD overwrites ----------------------------------------
  f_vec <- RTMB::AD(numeric(n.states),force=TRUE)
  dfdx_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.states, ncol=n.states),force=TRUE)
  g_mat <- RTMB::AD(RTMB::matrix(0,nrow=n.states, ncol=n.diffusions),force=TRUE)
  h_vec <- RTMB::AD(numeric(n.obs),force=TRUE)
  dhdx_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.obs, ncol=n.states),force=TRUE)
  hvar_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.obs, ncol=n.obs),force=TRUE)
  
  stateVec <- RTMB::AD(stateVec,force=TRUE)
  covMat <- RTMB::AD(covMat,force=TRUE)
  
  # obsMat <- RTMB::AD(obsMat,force=T)
  inputMat <- RTMB::AD(inputMat,force=T)
  
  # likelihood function --------------------------------------
  ekf.nll = function(p){
    
    ####### Sometimes necessary to avoid rtmb errors #######
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    
    ####### Parameters into vector #######
    parVec <- do.call(c, p[1:n.pars])
    # stateVec <- do.call(c, p[(n.pars+1):(n.pars+n.states)])
    
    ####### Neg. LogLikelihood #######
    nll <- 0
    
    ####### Pre-Allocated Object #######
    I0 <- RTMB::diag(n.states)
    E0 <- RTMB::diag(n.obs)
    
    ####### Stationary Solution #######
    inputVec = inputMat[1,]
    if(estimate.initial){
      # 1. Mean (root-find)
      .F <- RTMB::MakeTape(function(x) sum(f__(x, parVec, inputVec)^2), numeric(n.states))
      stateVec <- .F$newton(1:n.states)(numeric(0))
      # 2. Covar (Lyapunov via vectorization)
      A <- dfdx__(stateVec, parVec, inputVec)
      G <- g__(stateVec, parVec, inputVec)
      Q <- G %*% t(G)
      P <- kron_left(A) + kron_right(A)
      X <- -RTMB::solve(P, as.numeric(Q))
      covMat <- RTMB::matrix(X, nrow=n.states)
    }
    
    ######## Data update ########
    obsVec = obsMat[1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      y = obsVec[obsVec_bool]
      E = E0[obsVec_bool,, drop=FALSE]
      C = E %*% dhdx__(stateVec, parVec, inputVec)
      e = y - E %*% h__(stateVec, parVec, inputVec)
      V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
      R = C %*% covMat %*% t(C) + V
      K = covMat %*% t(C) %*% RTMB::solve(R)
      # Likelihood Contribution
      nll <- nll + loss_fun(e,R)
      # Update State/Cov
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    }
    
    ###### Main Loop #######
    for(i in 1:(nrow(obsMat)-1)){
      
      # load input vector
      inputVec = inputMat[i,]
      dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
      
      ###### time update - ode solve moments #######
      for(j in 1:ode_timesteps[i]){
        sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i])
        stateVec = sol[[1]]
        covMat = sol[[2]]
        # update input vec (first order interpolation ish)
        inputVec = inputVec + dinputVec
      }
      
      ######## data update ########
      inputVec = inputMat[i+1,]
      obsVec = obsMat[i+1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        y = obsVec[obsVec_bool]
        E = E0[obsVec_bool,, drop=FALSE]
        C = E %*% dhdx__(stateVec, parVec, inputVec)
        e = y - E %*% h__(stateVec, parVec, inputVec)
        V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
        R = C %*% covMat %*% t(C) + V
        K = covMat %*% t(C) %*% RTMB::solve(R)
        # Likelihood Contribution
        nll <- nll + loss_fun(e,R)
        # Update State/Cov
        stateVec = stateVec + K %*% e
        covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
      }
      # end of main loop
    }
    
    # return
    return(nll)
  }
  
  # construct AD-likelihood function ----------------------------------------
  
  # add states to parameters
  # state.list <- as.list(stateVec)
  # names(state.list) <- private$state.names
  # parameters <- c(parameters, state.list)
  
  map <- lapply(private$fixed.pars, function(x) x$factor)
  
  nll <- RTMB::MakeADFun(func = ekf.nll,
                         parameters = parameters,
                         map = map,
                         silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}

#######################################################
#######################################################
# EKF TMB-IMPLEMENTATION (FOR OPTIMIZATION)
#######################################################
#######################################################
makeADFun_ekf_tmb = function(self, private){
  
  # Data ----------------------------------------
  
  # add mandatory entries to data
  tmb.data = list(
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names]),
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # ode
    ode_solver = private$ode.solver,
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    
    # loss function
    loss_type = private$loss$loss,
    loss_c = private$loss$c,
    
    # system size
    n_states = private$number.of.states,
    n_obs = private$number.of.observations,
    n_inputs = private$number.of.inputs,
    
    # estimate stationary levels
    estimate_stationary_initials = as.numeric(private$estimate.initial)
  )
  
  # MAP Estimation?
  tmb.map.data = list(
    MAP_bool = 0L
  )
  if (!is.null(private$map)) {
    bool = self$getParameters()[,"type"] == "free"
    tmb.map.data = list(
      MAP_bool = 1L,
      map_mean__ = private$map$mean[bool],
      map_cov__ = private$map$cov[bool,bool],
      map_ints__ = as.numeric(bool),
      sum_map_ints__ = sum(as.numeric(bool))
    )
  }
  
  # construct final data list
  data = c(tmb.data, tmb.map.data)
  
  
  # Parameters ----------------------------------------
  parVec = sapply(private$parameters, function(x) x[["initial"]])
  parameters = list(parVec = parVec)
  
  # Create map for fixed parameters ----------------------------------------
  pseq <- 1:private$number.of.pars
  id.fixed.pars <- private$parameter.names %in% names(private$fixed.pars)
  pseq[id.fixed.pars] <- NA
  map <- list(parVec = factor(pseq))
  
  # Create AD-likelihood function ---------------------------------------
  nll <- TMB::MakeADFun(data = data,
                        parameters = parameters,
                        map = map,
                        DLL = private$modelname.with.method,
                        silent = TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}


#######################################################
#######################################################
# EKF FILTER (FOR REPORTING)
#######################################################
#######################################################

ekf_filter_r = function(parVec, self, private)
{
  
  # parameters ----------------------------------------
  if(missing(parVec)){
    parVec = sapply(private$parameters, function(x) x[["initial"]])
  }
  
  # Data ----------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # methods and purpose
  ode.solver = private$ode.solver
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  
  ode_timesteps = private$ode.timesteps
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  
  # 1-step covariance ODE ----------------------------------------
  cov_ode_1step = function(covMat, stateVec, parVec, inputVec){
    G <- g__(stateVec, parVec, inputVec)
    AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  
  # forward euler ----------------------------------------
  if(ode.solver==1){
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      X1 = stateVec + f__(stateVec, parVec, inputVec) * dt
      P1 = covMat + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
      
      return(list(X1,P1))
    }
  } else if (ode.solver==2) {
    # rk4 ----------------------------------------
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      # Initials
      X0 = stateVec
      P0 = covMat
      
      # Classical 4th Order Runge-Kutta Method
      # 1. Approx Slope at Initial Point
      k1 = f__(stateVec, parVec, inputVec)
      c1 = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 2. First Approx Slope at Midpoint
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + 0.5 * dt * k1
      covMat   = P0 + 0.5 * dt * c1
      k2       = f__(stateVec, parVec, inputVec)
      c2       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 3. Second Approx Slope at Midpoint
      stateVec = X0 + 0.5 * dt * k2
      covMat   = P0 + 0.5 * dt * c2
      k3       = f__(stateVec, parVec, inputVec)
      c3       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 4. Approx Slope at End Point
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + dt * k3
      covMat   = P0 + dt * c3
      k4       = f__(stateVec, parVec, inputVec)
      c4       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # ODE UPDATE
      X1 = X0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0 * dt
      P1 = P0 + (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0 * dt
      
      return(list(X1,P1))
    }
  } else {
    # Construct input interpolators
    #---------------------------------
    input.interp.funs <- vector("list",length=n.inputs)
    for(i in 1:length(input.interp.funs)){
      input.interp.funs[[i]] <- stats::approxfun(x=inputMat[,1], y=inputMat[,i], rule=2)
    }
    
    # Construct ode fun for DeSolve
    #---------------------------------
    ode.fun <- function(time, stateVec_and_covMat, parVec){
      inputVec <- sapply(input.interp.funs, function(f) f(time))
      # 
      stateVec <- head(stateVec_and_covMat, n.states)
      covMat <- matrix(tail(stateVec_and_covMat, -n.states),nrow=n.states)
      # 
      G <- g__(stateVec, parVec, inputVec)
      AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
      #
      dX <- f__(stateVec, parVec, inputVec)
      dP <- AcovMat + t(AcovMat) + G %*% t(G)
      return(list(c(dX,dP)))
    }
    
    # Construct function to call in likelihood functon
    #---------------------------------
    ode_integrator <- function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      out <- deSolve::ode(y = c(stateVec, covMat),
                          times = c(inputVec[1], inputVec[1]+dt),
                          func = ode.fun,
                          parms = parVec,
                          method=ode.solver)[2,-1]
      
      return(
        list(head(out,n.states),
             matrix(tail(out,-n.states),nrow=n.states)
        )
      )
    }
  }
  
  # user-defined functions ---------------------------
  for(i in seq_along(private$rekf.function.strings)){
    eval(parse(text=private$rekf.function.strings[[i]]))
  }
  
  # new functions ----------------------------------------
  # error function ----------------------------------------
  erf = function(x){
    y <- sqrt(2) * x
    2*RTMB::pnorm(y)-1
  }
  
  ####### STORAGE #######
  xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
  
  ####### Neg. LogLikelihood #######
  nll <- 0
  
  ####### Pre-Allocated Object #######
  I0 <- diag(n.states)
  E0 <- diag(n.obs)
  
  ####### INITIAL STATE / COVARIANCE #######
  # The state/covariance is either given by user or obtained from solving the
  # stationary mean, and then solving for the covariance.
  # In principle these are coupled equations, but believe that root-finding
  # both simultaneously can lead to likelihood blow-up.
  inputVec = inputMat[1,]
  if(estimate.initial){
    # 1. Root-find stationary mean
    opt <- stats::nlminb(stateVec, function(x) sum(f__(x, parVec, inputVec)^2),control=list(trace=0))
    stateVec <- opt$par
    # 2. Use stationary mean to solve lyapunov eq. for associated covariance
    A <- dfdx__(stateVec, parVec, inputVec)
    G <- g__(stateVec, parVec, inputVec)
    Q <- G %*% t(G)
    P <- kronecker(A,I0) + kronecker(I0,A)
    X <- -solve(P, as.numeric(Q))
    covMat <- matrix(X, nrow=n.states)
  }
  xPrior[[1]] <- stateVec
  pPrior[[1]] <- covMat
  
  ######## (PRE) DATA UPDATE ########
  # This is done to include the first measurements in the provided data
  # We update the state and covariance based on the "new" measurement
  obsVec = obsMat[1,]
  obsVec_bool = !is.na(obsVec)
  if(any(obsVec_bool)){
    y = obsVec[obsVec_bool]
    E = E0[obsVec_bool,, drop=FALSE]
    C = E %*% dhdx__(stateVec, parVec, inputVec)
    e = y - E %*% h__(stateVec, parVec, inputVec)
    V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
    R = C %*% covMat %*% t(C) + V
    K = covMat %*% t(C) %*% solve(R)
    # Update State/Cov
    stateVec = stateVec + K %*% e
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    # Store innovation and covariance
    Innovation[[1]] = e
    InnovationCovariance[[1]] = R
  }
  xPost[[1]] <- stateVec
  pPost[[1]] <- covMat
  
  ###### MAIN LOOP START #######
  for(i in 1:(nrow(obsMat)-1)){
    
    inputVec = inputMat[i,]
    dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
    
    ###### TIME UPDATE #######
    # We solve the first two moments forward in time
    for(j in 1:ode_timesteps[i]){
      sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i])
      stateVec = sol[[1]]
      covMat = sol[[2]]
      inputVec = inputVec + dinputVec
    }
    xPrior[[i+1]] = stateVec
    pPrior[[i+1]] = covMat
    
    ######## DATA UPDATE ########
    # We update the state and covariance based on the "new" measurement
    inputVec = inputMat[i+1,]
    obsVec = obsMat[i+1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      y = obsVec[obsVec_bool]
      E = E0[obsVec_bool,, drop=FALSE] #permutation matrix with rows removed
      C = E %*% dhdx__(stateVec, parVec, inputVec)
      e = y - E %*% h__(stateVec, parVec, inputVec)
      V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
      R = C %*% covMat %*% t(C) + V
      K = covMat %*% t(C) %*% solve(R)
      # Update State/Cov
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
      # Store innovation and covariance
      Innovation[[i+1]] <- e
      InnovationCovariance[[i+1]] <- R
    }
    xPost[[i+1]] = stateVec
    pPost[[i+1]] = covMat
  }
  ###### MAIN LOOP END #######
  
  ####### RETURN #######
  returnlist <- list(xPost = xPost,
                     pPost = pPost,
                     xPrior = xPrior,
                     pPrior = pPrior,
                     Innovation = Innovation,
                     InnovationCovariance = InnovationCovariance)
  
  return(invisible(returnlist))
}

#######################################################
#######################################################
# RETURN FIT FOR EKF RTMB
#######################################################
#######################################################

calculate_fit_statistics_ekf <- function(self, private){
  
  # Initialization and Clearing -----------------------------------
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  # clear fit
  private$fit = NULL
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  
  # Parameters and Uncertainties -----------------------------------
  
  # objective value
  private$fit$nll = private$opt$objective
  
  # gradient
  private$fit$nll.gradient = try_withWarningRecovery(
    {
      nll.grad = as.vector(private$nll$gr(private$opt$par))
      names(nll.grad) = names(private$free.pars)
      nll.grad
    }
  )
  if (inherits(private$fit$nll.gradient,"try-error")) {
    private$fit$nll.gradient = NULL
  }
  
  # # hessian
  private$fit$nll.hessian = try_withWarningRecovery(
    {
      nll.hess = private$nll$he(private$opt$par)
      rownames(nll.hess) = names(private$free.pars)
      colnames(nll.hess) = names(private$free.pars)
      nll.hess
    }
  )
  if (inherits(private$fit$nll.hessian, "try-error")) {
    private$fit$nll.hessian = NULL
  }
  
  # parameter estimates and standard deviation
  npars <- length(private$opt$par)
  private$fit$par.fixed = private$opt$par
  private$fit$sd.fixed = rep(NA,npars)
  private$fit$cov.fixed = array(NA,dim=c(npars,npars))
  
  # parameter std. error and full covariance by hessian inversion
  if(!is.null(private$fit$nll.hessian)){
    
    # OPTION 0 -----------------------------------
    # Invert full hessian
    temp.hessian = private$fit$nll.hessian
    covariance = try_withWarningRecovery(solve(temp.hessian))
    
    # OPTION 1 -----------------------------------
    # Remove all row/cols where the diagonal elements are smalller than threshold
    min.diag = 1e-8
    remove.ids <- which(diag(temp.hessian) < min.diag)
    if(inherits(covariance,"try-error") && any(remove.ids)){
      
      # try to invert reduced hessian
      covariance = try_withWarningRecovery(solve(temp.hessian[-remove.ids, -remove.ids]))
      std.dev = try_withWarningRecovery(sqrt(diag(covariance)))
      
      if(!inherits(covariance,"try-error")){
        covtemp = array(NA, dim=dim(temp.hessian))
        covtemp[-remove.ids, -remove.ids] = covariance
        covariance <- covtemp
      }
      if(!inherits(std.dev,"try-error")){
        stdtemp = rep(NA, length(private$fit$par.fixed))
        stdtemp[-remove.ids] <- std.dev
        std.dev <- stdtemp
      }
    }
    
    # OPTION 2 -----------------------------------
    # Remove small diagonal element one by one until solve is succesful
    failed.to.invert.hessian = TRUE
    id.diag.hess <- order(diag(temp.hessian))
    i = 1
    if(inherits(covariance,"try-error")){
      while(failed.to.invert.hessian){
        remove.ids <- id.diag.hess[1:i]
        covariance <- try_withWarningRecovery(solve(temp.hessian[-remove.ids,-remove.ids]))
        std.dev = try_withWarningRecovery(sqrt(diag(covariance)))
        i = i + 1
        
        # if unsuccesful break while loop
        if(i == nrow(temp.hessian)){
          break
        }
        
        # if succesful update results
        if(!inherits(covariance,"try-error")){
          failed.to.invert.hessian <- FALSE
          # 
          covtemp = array(NA, dim=dim(temp.hessian))
          covtemp[-remove.ids, -remove.ids] = covariance
          covariance <- covtemp
          # 
          if(!inherits(std.dev,"try-error")){
            stdtemp = rep(NA, length(private$fit$par.fixed))
            stdtemp[-remove.ids] <- std.dev
            std.dev <- stdtemp
          }
        }
      }
    }
    
    if(!inherits(covariance,"try-error")){
      private$fit$cov.fixed <- covariance
      std.dev <- try_withWarningRecovery(sqrt(diag(covariance)))
    }
    if(!inherits(std.dev,"try-error")){
      private$fit$sd.fixed <- std.dev
    }
    
  }
  
  # States -----------------------------------
  
  # Extract reported items from nll
  estimated_pars <- self$getParameters()[,"estimate"]
  rep <- ekf_filter_r(estimated_pars, self, private)
  
  # Prior States
  temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPrior))))
  temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPrior,diag)))))
  
  colnames(temp.states) = c("t", private$state.names)
  colnames(temp.sd) = c("t",private$state.names)
  private$fit$states$mean$prior = as.data.frame(temp.states)
  private$fit$states$sd$prior = as.data.frame(temp.sd)
  private$fit$states$cov$prior = rep$pPrior
  names(private$fit$states$cov$prior) = paste("t = ",private$data$t,sep="")
  
  # Posterior States
  temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPost))))
  temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPost,diag)))))
  colnames(temp.states) = c("t",private$state.names)
  colnames(temp.sd) = c("t",private$state.names)
  private$fit$states$mean$posterior = as.data.frame(temp.states)
  private$fit$states$sd$posterior = as.data.frame(temp.sd)
  private$fit$states$cov$posterior = rep$pPost
  names(private$fit$states$cov$posterior) = paste("t = ",private$data$t,sep="")
  
  # Residuals -----------------------------------
  rowNAs = as.matrix(!is.na(private$data[private$obs.names]))
  sumrowNAs = rowSums(rowNAs)
  
  nr <- nrow(private$data)
  temp.res = matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
  temp.var =  matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
  
  nr <- nrow(private$data)
  innovation = rep$Innovation
  innovation.cov = rep$InnovationCovariance
  names(innovation.cov) = paste("t = ", private$data$t, sep="")
  for (i in 1:nr) {
    if (sumrowNAs[i] > 0) {
      temp.res[i,rowNAs[i,]] = innovation[[i]]
      temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
    }
  }
  temp.res <- data.frame(private$data$t, temp.res)
  temp.sd = data.frame(private$data$t, try_withWarningRecovery(sqrt(temp.var)))
  names(temp.res) = c("t", private$obs.names)
  names(temp.sd) = c("t", private$obs.names)
  # should we remove the empty matrices?
  # innovation.cov = innovation.cov[sumrowNAs!=0]
  
  private$fit$residuals$residuals = temp.res
  private$fit$residuals$sd = temp.sd
  private$fit$residuals$normalized = temp.res
  private$fit$residuals$normalized[,-1] = temp.res[,-1]/temp.sd[,-1]
  private$fit$residuals$cov = innovation.cov
  
  
  # Observations -----------------------------------
  
  # We need all states, inputs and parameter values to evaluate the observation
  # put them in a list
  listofvariables.prior = c(
    # states
    as.list(private$fit$states$mean$prior[-1]),
    # estimated free parameters
    as.list(private$fit$par.fixed),
    # fixed parameters
    lapply(private$fixed.pars, function(x) x$initial),
    # inputs
    as.list(private$data)
  )
  
  listofvariables.posterior = c(
    # states
    as.list(private$fit$states$mean$posterior[-1]),
    # estimated free parameters
    as.list(private$fit$par.fixed),
    # fixed parameters
    lapply(private$fixed.pars, function(x) x$initial),
    # inputs
    as.list(private$data)
  )
  obs.df.prior = as.data.frame(
    lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.prior)})
  )
  obs.df.posterior = as.data.frame(
    lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.posterior)})
  )
  
  private$fit$observations$mean$prior = data.frame(t=private$data$t, obs.df.prior)
  private$fit$observations$mean$posterior = data.frame(t=private$data$t, obs.df.posterior)
  
  # Observation variances -----------------------------------
  
  # The observation variance (to first order) is:
  #
  # y = h(x) + e -> var(y) = dhdx var(x) dhdx^T + var(e)
  
  n <- private$number.of.states
  m <- private$number.of.observations
  
  # create expression for observation covariance to be 'eval'-uated
  jac.h = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      jac.h = c(jac.h, private$diff.terms.obs[[i]][[j]])
    }
  }
  dhdx <- parse(text=sprintf("matrix(c(%s),nrow=%s,ncol=%s)", paste(jac.h,collapse=","), m, n))[[1]]
  obs.var <- c()
  for(i in seq_along(private$obs.var.trans)){
    obs.var = c(obs.var, private$obs.var.trans[[i]]$rhs)
  }
  obsvar <- parse(text=sprintf("diag(%s)*c(%s)", m, paste(obs.var,collapse=",")))[[1]]
  yCov <- substitute(dhdx %*% xCov %*% t(dhdx) + eCov, list(dhdx=dhdx, eCov=obsvar))
  
  # evaluate prior and posterior variance
  list.of.parameters <- c(
    as.list(private$fit$par.fixed),
    lapply(private$fixed.pars, function(x) x$initial)
  )
  
  # prior
  obsvar.prior <- list()
  for(i in seq_along(private$data$t)){
    obsvar.prior[[i]] <- eval(expr = yCov,
                              envir = c(list.of.parameters,
                                        as.list(private$fit$states$mean$prior[i,-1]),
                                        list(xCov = private$fit$states$cov$prior[[i]]),
                                        as.list(private$data[i,-1])
                              ))
  }
  names(obsvar.prior) <- names(private$fit$states$cov$prior)
  private$fit$observations$cov$prior <- obsvar.prior
  obs.sd.prior <- cbind(private$fit$states$mean$prior["t"], do.call(rbind, lapply(obsvar.prior, diag)))
  rownames(obs.sd.prior) <- NULL
  names(obs.sd.prior) <- c("t",private$obs.names)
  private$fit$observations$sd$prior <- obs.sd.prior
  
  # posterior
  obsvar.post <- list()
  for(i in seq_along(private$data$t)){
    obsvar.post[[i]] <- eval(expr = yCov,
                             envir = c(list.of.parameters,
                                       as.list(private$fit$states$mean$posterior[i,-1]),
                                       list(xCov = private$fit$states$cov$posterior[[i]]),
                                       as.list(private$data[i,-1])
                             ))
  }
  names(obsvar.post) <- names(private$fit$states$cov$posterior)
  private$fit$observations$cov$posterior <- obsvar.post
  obs.sd.post <- cbind(private$fit$states$mean$posterior["t"], do.call(rbind, lapply(obsvar.post, diag)))
  rownames(obs.sd.post) <- NULL
  names(obs.sd.post) <- c("t",private$obs.names)
  private$fit$observations$sd$posterior <- obs.sd.post
  
  # t-values and Pr( t > t_test ) -----------------------------------
  
  private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
  private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
  
  # clone private into fit -----------------------------------
  private$fit$private <- self$clone()$.__enclos_env__$private
  
  # set s3 class -----------------------------------
  class(private$fit) = "ctsmTMB.fit"
  
  # return -----------------------------------
  return(invisible(self))
  
}

#######################################################
#######################################################
# PREDICTION EKF R IMPLEMENTATION
#######################################################
#######################################################

ekf_r_prediction = function(self, private)
{
  
  # parameters ----------------------------------------
  parVec <- private$pars
  
  # Data ----------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # methods and purpose
  ode.solver = private$ode.solver
  k.ahead <- private$n.ahead
  last.pred.index <- private$last.pred.index
  
  # initial
  stateVec = private$pred.initial.state$x0
  covMat = private$pred.initial.state$p0
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  
  # 1-step covariance ODE ----------------------------------------
  cov_ode_1step = function(covMat, stateVec, parVec, inputVec){
    G <- g__(stateVec, parVec, inputVec)
    AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  
  # forward euler ----------------------------------------
  if(ode.solver==1){
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      X1 = stateVec + f__(stateVec, parVec, inputVec) * dt
      P1 = covMat + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
      
      return(list(X1,P1))
    }
  } else if (ode.solver==2) {
    # rk4 ----------------------------------------
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      # Initials
      X0 = stateVec
      P0 = covMat
      
      # Classical 4th Order Runge-Kutta Method
      # 1. Approx Slope at Initial Point
      k1 = f__(stateVec, parVec, inputVec)
      c1 = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 2. First Approx Slope at Midpoint
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + 0.5 * dt * k1
      covMat   = P0 + 0.5 * dt * c1
      k2       = f__(stateVec, parVec, inputVec)
      c2       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 3. Second Approx Slope at Midpoint
      stateVec = X0 + 0.5 * dt * k2
      covMat   = P0 + 0.5 * dt * c2
      k3       = f__(stateVec, parVec, inputVec)
      c3       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 4. Approx Slope at End Point
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + dt * k3
      covMat   = P0 + dt * c3
      k4       = f__(stateVec, parVec, inputVec)
      c4       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # ODE UPDATE
      X1 = X0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0 * dt
      P1 = P0 + (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0 * dt
      
      return(list(X1,P1))
    }
  } else {
    # Construct input interpolators
    #---------------------------------
    input.interp.funs <- vector("list",length=n.inputs)
    for(i in 1:length(input.interp.funs)){
      input.interp.funs[[i]] <- stats::approxfun(x=inputMat[,1], y=inputMat[,i], rule=2)
    }
    
    # Construct ode fun for DeSolve
    #---------------------------------
    ode.fun <- function(time, stateVec_and_covMat, parVec){
      inputVec <- sapply(input.interp.funs, function(f) f(time))
      # 
      stateVec <- head(stateVec_and_covMat, n.states)
      covMat <- matrix(tail(stateVec_and_covMat, -n.states),nrow=n.states)
      # 
      G <- g__(stateVec, parVec, inputVec)
      AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
      #
      dX <- f__(stateVec, parVec, inputVec)
      dP <- AcovMat + t(AcovMat) + G %*% t(G)
      return(list(c(dX,dP)))
    }
    
    # Construct function to call in likelihood functon
    #---------------------------------
    ode_integrator <- function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      out <- deSolve::ode(y = c(stateVec, covMat),
                          times = c(inputVec[1], inputVec[1]+dt),
                          func = ode.fun,
                          parms = parVec,
                          method=ode.solver)[2,-1]
      
      return(
        list(head(out,n.states),
             matrix(tail(out,-n.states),nrow=n.states)
        )
      )
    }
  }
  
  # user-defined functions ---------------------------
  for(i in seq_along(private$rekf.function.strings)){
    eval(parse(text=private$rekf.function.strings[[i]]))
  }
  
  # new functions ----------------------------------------
  # error function ----------------------------------------
  erf = function(x){
    y <- sqrt(2) * x
    2*RTMB::pnorm(y)-1
  }
  
  ####### STORAGE #######
  predMats <- lapply(1:last.pred.index, function(x) matrix(NA,nrow=k.ahead+1,ncol=n.states+n.states^2))
  
  ####### Pre-Allocated Object #######
  I0 <- diag(n.states)
  E0 <- diag(n.obs)
  
  ####### INITIAL STATE / COVARIANCE #######
  # The state/covariance is either given by user or obtained from solving the
  # stationary mean, and then solving for the covariance.
  # In principle these are coupled equations, but believe that root-finding
  # both simultaneously can lead to likelihood blow-up.
  inputVec = inputMat[1,]
  if(estimate.initial){
    # 1. Root-find stationary mean
    opt <- stats::nlminb(numeric(n.states), function(x) sum(f__(x, parVec, inputVec)^2))
    stateVec <- opt$par
    # 2. Use stationary mean to solve lyapunov eq. for associated covariance
    A <- dfdx__(stateVec, parVec, inputVec)
    G <- g__(stateVec, parVec, inputVec)
    Q <- G %*% t(G)
    P <- kronecker(A,I0) + kronecker(I0,A)
    X <- -solve(P, as.numeric(Q))
    covMat <- matrix(X, nrow=n.states)
  }
  
  ######## (PRE) DATA UPDATE ########
  # This is done to include the first measurements in the provided data
  # We update the state and covariance based on the "new" measurement
  obsVec = obsMat[1,]
  obsVec_bool = !is.na(obsVec)
  if(any(obsVec_bool)){
    y = obsVec[obsVec_bool]
    E = E0[obsVec_bool,, drop=FALSE]
    C = E %*% dhdx__(stateVec, parVec, inputVec)
    e = y - E %*% h__(stateVec, parVec, inputVec)
    V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
    R = C %*% covMat %*% t(C) + V
    K = covMat %*% t(C) %*% solve(R)
    # Update State/Cov
    stateVec = stateVec + K %*% e
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
  }
  
  ###### TIME LOOP ####### 
  for(i in 1:last.pred.index){
    # Start Loop
    
    predMats[[i]][1,] <- c(stateVec, covMat)
    
    ###### K-STEP AHEAD LOOP ####### 
    for(k in 1:k.ahead){
      
      inputVec = inputMat[i+k-1,]
      dinputVec = (inputMat[i+k,] - inputVec)/ode_timesteps[i+k-1]
      
      # Solve moment ODEs
      for(j in 1:ode_timesteps[i+k-1]){
        sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i+k-1])
        stateVec = sol[[1]]
        covMat = sol[[2]]
        inputVec = inputVec + dinputVec
      }
      
      predMats[[i]][k+1,] <- c(stateVec, covMat)
      
    }
    
    stateVec <- head(predMats[[i]][2,], n.states)
    covMat <- matrix(tail(predMats[[i]][2,], n.states^2), nrow=n.states)
    
    ######## DATA UPDATE ########
    # We update the state and covariance based on the "new" measurement
    inputVec = inputMat[i+1,]
    obsVec = obsMat[i+1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      y = obsVec[obsVec_bool]
      E = E0[obsVec_bool,, drop=FALSE] #permutation matrix with rows removed
      C = E %*% dhdx__(stateVec, parVec, inputVec)
      e = y - E %*% h__(stateVec, parVec, inputVec)
      V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
      R = C %*% covMat %*% t(C) + V
      K = covMat %*% t(C) %*% solve(R)
      # Update State/Cov
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    }
    
    # End Loop
  }
  ###### MAIN LOOP END #######
  
  ####### STORE PREDICTION #######
  private$prediction = list(predMats=predMats)
  
  ####### RETURN #######
  return(invisible(self))
}

#######################################################
#######################################################
# PREDICTION EKF C++ IMPLEMENTATION
#######################################################
#######################################################

ekf_rcpp_prediction = function(self, private){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # predict using c++ function
  mylist <- execute_ekf_prediction2(private$rcpp_function_ptr$f,
                                    private$rcpp_function_ptr$g,
                                    private$rcpp_function_ptr$dfdx,
                                    private$rcpp_function_ptr$h,
                                    private$rcpp_function_ptr$dhdx,
                                    private$rcpp_function_ptr$hvar,
                                    obsMat,
                                    inputMat,
                                    private$pars,
                                    private$pred.initial.state$p0,
                                    private$pred.initial.state$x0,
                                    private$ode.timestep.size,
                                    private$ode.timesteps,
                                    numeric_is_not_na_obsMat,
                                    number_of_available_obs,
                                    private$number.of.states,
                                    private$number.of.observations,
                                    private$last.pred.index,
                                    private$n.ahead,
                                    private$ode.solver)
  
  ####### STORE PREDICTION #######
  private$prediction = mylist
  
  ####### RETURN #######
  return(invisible(self))
}

create_ekf_predict_return = function(return.covariance, return.k.ahead, self, private){
  
  # Simlify variable names
  n = private$number.of.states
  n.ahead = private$n.ahead
  state.names = private$state.names
  last.pred.index = private$last.pred.index
  
  # Create return data.frame
  df.out = data.frame(matrix(nrow=last.pred.index*(n.ahead+1), ncol=5+n+n^2))
  disp_names = sprintf(rep("cor.%s.%s",n^2),rep(state.names, each=n),rep(state.names,n))
  disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  if(return.covariance){
    disp_names = sprintf(rep("cov.%s.%s",n^2),rep(state.names,each=n),rep(state.names,n))
    disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  }
  names(df.out) = c("i.","j.","t.i","t.j","k.ahead",state.names,disp_names)
  
  # Fill out data.frame
  ran = 0:(last.pred.index-1)
  df.out["i."] = rep(ran,each=n.ahead+1)
  df.out["j."] = df.out["i."] + rep(0:n.ahead,last.pred.index)
  df.out["t.i"] = rep(private$data$t[ran+1],each=n.ahead+1)
  df.out["t.j"] = private$data$t[df.out[,"i."]+1+rep(0:n.ahead,last.pred.index)]
  df.out["k.ahead"] = rep(0:n.ahead,last.pred.index)
  df.obs = df.out[c("i.","j.","t.i","t.j","k.ahead")]
  
  df.out[, c(state.names, disp_names)] <- do.call(rbind, private$prediction$predMats)
  if(!return.covariance){
    diag.ids <- seq(from=1,to=n^2,by=n+1)
    .seq <- seq(from=1,to=n^2,by=1)
    non.diag.ids <- .seq[!(.seq %in% diag.ids)]
    df.out[,disp_names[non.diag.ids]] <- t(apply(df.out, 1, function(x) as.vector(stats::cov2cor(matrix(tail(x, n^2), nrow=n)))))[,non.diag.ids]
  }
  
  ##### OBSERVATION PREDICTIONS #####
  # calculate observations at every time-step in predict
  inputs.df = private$data[df.out[,"j."]+1,private$input.names]
  
  named.pars.list = as.list(private$pars)
  names(named.pars.list) = private$parameter.names
  # create environment
  env.list = c(
    # states
    as.list(df.out[state.names]),
    # inputs
    as.list(inputs.df),
    # free and fixed parameters
    named.pars.list
  )
  
  # calculate observations
  obs.df.predict = as.data.frame(
    lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = env.list)})
  )
  names(obs.df.predict) = paste(private$obs.names)
  
  # add data observation to output data.frame 
  obs.df.data = private$data[df.out[,"j."]+1, private$obs.names, drop=F]
  names(obs.df.data) = paste(private$obs.names,".data",sep="")
  
  df.obs = cbind(df.obs, obs.df.predict, obs.df.data)
  
  # return only specific n.ahead
  df.out = df.out[df.out[,"k.ahead"] %in% return.k.ahead,]
  df.obs = df.obs[df.obs[,"k.ahead"] %in% return.k.ahead,]
  
  list.out = list(states = df.out, observations = df.obs)
  class(list.out) = c(class(list.out), "ctsmTMB.pred")
  
  private$prediction = list.out
  
  return(list.out)
}

#######################################################
#######################################################
# EKF SIMULATION R
#######################################################
#######################################################

ekf_r_simulation = function(self, private, nsims)
{
  
  if(!private$silent) message("Simulating with R...")
  
  # parameters ----------------------------------------
  parVec <- private$pars
  
  # Data ----------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # methods and purpose
  ode.solver = private$ode.solver
  k.ahead <- private$n.ahead
  last.pred.index <- private$last.pred.index
  
  # initial
  stateVec = private$pred.initial.state$x0
  covMat = private$pred.initial.state$p0
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  sde_timestep_size <- private$simulation.timestep.size
  sde_timesteps <- private$simulation.timesteps
  
  # 1-step covariance ODE ----------------------------------------
  cov_ode_1step = function(covMat, stateVec, parVec, inputVec){
    G <- g__(stateVec, parVec, inputVec)
    AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  # forward euler ----------------------------------------
  if(ode.solver==1){
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      X1 = stateVec + f__(stateVec, parVec, inputVec) * dt
      P1 = covMat + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
      
      return(list(X1,P1))
    }
  } else if (ode.solver==2) {
    # rk4 ----------------------------------------
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      # Initials
      X0 = stateVec
      P0 = covMat
      
      # Classical 4th Order Runge-Kutta Method
      # 1. Approx Slope at Initial Point
      k1 = f__(stateVec, parVec, inputVec)
      c1 = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 2. First Approx Slope at Midpoint
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + 0.5 * dt * k1
      covMat   = P0 + 0.5 * dt * c1
      k2       = f__(stateVec, parVec, inputVec)
      c2       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 3. Second Approx Slope at Midpoint
      stateVec = X0 + 0.5 * dt * k2
      covMat   = P0 + 0.5 * dt * c2
      k3       = f__(stateVec, parVec, inputVec)
      c3       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 4. Approx Slope at End Point
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + dt * k3
      covMat   = P0 + dt * c3
      k4       = f__(stateVec, parVec, inputVec)
      c4       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # ODE UPDATE
      X1 = X0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0 * dt
      P1 = P0 + (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0 * dt
      
      return(list(X1,P1))
    }
  } else {
    # Construct input interpolators
    #---------------------------------
    input.interp.funs <- vector("list",length=n.inputs)
    for(i in 1:length(input.interp.funs)){
      input.interp.funs[[i]] <- stats::approxfun(x=inputMat[,1], y=inputMat[,i], rule=2)
    }
    
    # Construct ode fun for DeSolve
    #---------------------------------
    ode.fun <- function(time, stateVec_and_covMat, parVec){
      inputVec <- sapply(input.interp.funs, function(f) f(time))
      # 
      stateVec <- head(stateVec_and_covMat, n.states)
      covMat <- tail(stateVec_and_covMat, -n.states)
      # 
      G <- g__(stateVec, parVec, inputVec)
      AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
      #
      dX <- f__(stateVec, parVec, inputVec)
      dP <- AcovMat + t(AcovMat) + G %*% t(G)
      return(list(c(dX,dP)))
    }
    
    # Construct function to call in likelihood functon
    #---------------------------------
    ode_integrator <- function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      out <- deSolve::ode(y = c(stateVec, covMat),
                          times = c(inputVec[1], inputVec[1]+dt),
                          func = ode.fun,
                          parms = parVec,
                          method=ode.solver)[2,-1]
      
      return(
        list(head(out,n.states),
             matrix(tail(out,-n.states),nrow=n.states)
        )
      )
    }
  }
  
  euler_maruyama_simulation <- function(stateMat, parVec, inputVec, sim.dt){
    # We need to perform an euler maruyama step for each column in the stateMat,
    
    # Each column is f(stateVec)
    .F <- apply(stateMat, 2, function(x) f__(x, parVec, inputVec))
    # Each column is g(stateVec) as a vector
    .G <- lapply(1:nsims, function(i) g__(stateMat[,i],parVec, inputVec))
    # Draw nsims random variables for each db
    dB <- matrix(stats::rnorm(n.diffusions*nsims), ncol=nsims)
    # Perform multiplication
    GdB <- sapply(1:nsims, function(i) .G[[i]] %*% dB[,i])
    
    # Produce euler-maruyama step for all states simultaneously
    return(stateMat + .F * sim.dt + GdB * sqrt(sim.dt))
  }
  
  # user-defined functions ---------------------------
  for(i in seq_along(private$rekf.function.strings)){
    eval(parse(text=private$rekf.function.strings[[i]]))
  }
  
  # new functions ----------------------------------------
  # error function ----------------------------------------
  erf = function(x){
    y <- sqrt(2) * x
    2*RTMB::pnorm(y)-1
  }
  
  ####### STORAGE #######
  xSim <- lapply(1:last.pred.index, function(x) vector("list", length=k.ahead+1))
  
  ####### Pre-Allocated Object #######
  I0 <- diag(n.states)
  E0 <- diag(n.obs)
  
  ####### INITIAL STATE / COVARIANCE #######
  # The state/covariance is either given by user or obtained from solving the
  # stationary mean, and then solving for the covariance.
  # In principle these are coupled equations, but believe that root-finding
  # both simultaneously can lead to likelihood blow-up.
  inputVec = inputMat[1,]
  if(estimate.initial){
    # 1. Root-find stationary mean
    opt <- stats::nlminb(numeric(n.states), function(x) sum(f__(x, parVec, inputVec)^2))
    stateVec <- opt$par
    # 2. Use stationary mean to solve lyapunov eq. for associated covariance
    A <- dfdx__(stateVec, parVec, inputVec)
    G <- g__(stateVec, parVec, inputVec)
    Q <- G %*% t(G)
    P <- kronecker(A,I0) + kronecker(I0,A)
    X <- -solve(P, as.numeric(Q))
    covMat <- matrix(X, nrow=n.states)
  }
  
  ######## (PRE) DATA UPDATE ########
  # This is done to include the first measurements in the provided data
  # We update the state and covariance based on the "new" measurement
  obsVec = obsMat[1,]
  obsVec_bool = !is.na(obsVec)
  if(any(obsVec_bool)){
    y = obsVec[obsVec_bool]
    E = E0[obsVec_bool,, drop=FALSE]
    C = E %*% dhdx__(stateVec, parVec, inputVec)
    e = y - E %*% h__(stateVec, parVec, inputVec)
    V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
    R = C %*% covMat %*% t(C) + V
    K = covMat %*% t(C) %*% solve(R)
    # Update State/Cov
    stateVec = stateVec + K %*% e
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
  }
  
  ###### TIME LOOP ####### 
  for(i in 1:last.pred.index){
    
    # We draw nsims samples of stateVecs from a multivariate normal; z = u + A * dB
    # where u is the mean (stateVec) and A is cholesky factor of covariance matrix (sqrt(covMat))
    # and dB is i.d.d normal vector
    # Each column in stateMat is a sample of a stateVec
    stateMat <- matrix(rep(stateVec, times=nsims), ncol=nsims) + Matrix::chol(covMat) %*% 
      matrix(stats::rnorm(nsims*n.states),ncol=nsims)
    xSim[[i]][[1]] <- t(stateMat)
    
    ###### K-STEP AHEAD LOOP ####### 
    for(k in 1:k.ahead){
      inputVec = inputMat[i+k-1,]
      dinputVec = (inputMat[i+k,] - inputVec)/sde_timesteps[i+k-1]
      
      # Solve moment ODEs
      for(j in 1:ode_timesteps[i+k-1]){
        stateMat <- euler_maruyama_simulation(stateMat, parVec, inputVec, sde_timestep_size[i+k-1])
        inputVec = inputVec + dinputVec
      }
      
      xSim[[i]][[k+1]] <- t(stateMat)
    }
    
    inputVec = inputMat[i,]
    dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
    
    # Solve ODE 1-Step Forward to get next posterior for simulations
    for(j in 1:ode_timesteps[i]){
      sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i])
      stateVec = sol[[1]]
      covMat = sol[[2]]
      inputVec = inputVec + dinputVec
    }
    
    ######## DATA UPDATE ########
    # We update the state and covariance based on the "new" measurement
    inputVec = inputMat[i+1,]
    obsVec = obsMat[i+1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      y = obsVec[obsVec_bool]
      E = E0[obsVec_bool,, drop=FALSE] #permutation matrix with rows removed
      C = E %*% dhdx__(stateVec, parVec, inputVec)
      e = y - E %*% h__(stateVec, parVec, inputVec)
      V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
      R = C %*% covMat %*% t(C) + V
      K = covMat %*% t(C) %*% solve(R)
      # Update State/Cov
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    }
    
    # End Loop
  }
  ###### MAIN LOOP END #######
  
  ####### RETURN #######
  private$simulation = xSim
  return(invisible(self))
}

#######################################################
#######################################################
# EKF SIMULATION C++
#######################################################
#######################################################

# Stochastic Euler-Maruyama simulation function that calls the underlying Rcpp simulation function
ekf_rcpp_simulation = function(self, private, n.sims){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # Call C++ function to perform simulation
  mylist = execute_ekf_simulation2(private$rcpp_function_ptr$f,
                                   private$rcpp_function_ptr$g,
                                   private$rcpp_function_ptr$dfdx,
                                   private$rcpp_function_ptr$h,
                                   private$rcpp_function_ptr$dhdx,
                                   private$rcpp_function_ptr$hvar,
                                   obsMat,
                                   inputMat,
                                   private$pars,
                                   private$pred.initial.state$p0,
                                   private$pred.initial.state$x0,
                                   private$ode.timestep.size,
                                   private$ode.timesteps,
                                   private$simulation.timestep.size,
                                   private$simulation.timesteps,
                                   numeric_is_not_na_obsMat,
                                   number_of_available_obs,
                                   private$number.of.states,
                                   private$number.of.observations,
                                   private$number.of.diffusions,
                                   private$last.pred.index,
                                   private$n.ahead,
                                   private$ode.solver,
                                   n.sims)
  
  private$simulation = mylist
  
  return(invisible(NULL))
}


create_ekf_simulation_return = function(return.k.ahead, n.sims, self, private){
  
  list.out = vector("list",length=private$number.of.states)
  names(list.out) = private$state.names
  
  # Compute the prediction times for each horizon
  ran = 0:(private$last.pred.index-1)
  t.j = private$data$t[rep(ran,each=private$n.ahead+1)+1+rep(0:private$n.ahead,private$last.pred.index)]
  t.j.splitlist = split(t.j, ceiling(seq_along(t.j)/(private$n.ahead+1)))
  list.of.time.vectors = lapply(t.j.splitlist, function(x) data.frame(t.j=x))
  
  for(i in seq_along(list.out)){
    list.out[[i]] = stats::setNames(
      lapply(private$simulation, function(ls.outer){
        t(do.call(cbind, lapply(ls.outer, function(ls.inner) ls.inner[,i])))
      }),
      paste0("i",ran)
    )
  }
  
  for(i in seq_along(list.out)){
    for(j in seq_along(list.out[[i]])){
      list.out[[i]][[j]] = data.frame(i = j-1, 
                                      j = (j-1):(j+private$n.ahead-1), 
                                      t.i = rep(private$data$t[i],private$n.ahead+1),
                                      t.j = list.of.time.vectors[[j]][,"t.j"], 
                                      k.ahead = 0:private$n.ahead,
                                      list.out[[i]][[j]]
      )
      nams = paste0(private$state.names,1:n.sims)
      names(list.out[[i]][[j]]) = c("i","j","t.i","t.j","k.ahead",nams)
    }
  }
  
  # Observations
  # eval(parse(text=private$rekf.function.strings$h))
  
  # We should take the state vector, and pass it through the observation vector, right?
  # The inputs should come from the inputMat for the j'th row...
  # There are so many though...

  private$simulation = list( states = list.out, observations = list() )
  
  return(invisible(NULL))
}
