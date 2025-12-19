#######################################################
# EXTENDED KALMAN FILTER # EXTENDED KALMAN FILTER
#######################################################

makeADFun_ekf_rtmb = function(self, private)
{
  
  # Tape Configration ----------------------
  ConfigureTape("RTMB", self, private)
  
  # Data ----------------------------------------
  getSystemDimensions()
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  # input and obs matrix
  inputMat = as.matrix(private$data[private$input.names])
  obsMat = as.matrix(private$data[private$obs.names])
  
  # create and load state space functions
  force.ad <- private$advanced.settings$forceAD
  create.state.space.functions.for.estimation(force.ad)
  
  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  getAdjoints()
  getLossFunction()
  getOdeSolvers()
  if(private$estimate.initial) {
    getInitialStateEstimator()
  }
  getKalmanFunctions()
  date.update.fun <- kalman.data.update.with.nll
  if(private$train.against.full.prediction){
    date.update.fun <- kalman.no.update.with.nll
  }
  
  # Timesteps, Observations, Inputs and Parameters ----------------------------
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  
  ####### Pre-Allocated Object #######
  I0 <- RTMB::diag(n.states)
  E0 <- RTMB::diag(n.obs)
  
  # likelihood function --------------------------------------
  ekf.nll = function(p){
    
    ####### Sometimes necessary to avoid rtmb errors #######
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    
    ####### Parameters into vector #######
    parVec <- do.call(c, p[1:n.pars])
    
    ####### Neg. LogLikelihood #######
    nll <- 0
    
    ####### Stationary Solution #######
    inputVec = inputMat[1,]
    if(private$estimate.initial){
      stateVec <- f.initial.state.newton(c(parVec, inputVec))
      # covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
    }
    
    ######## Data Update ########
    obsVec = obsMat[1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      data.update <- date.update.fun(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
      stateVec <- data.update[[1]]
      covMat <- data.update[[2]]
      nll <- nll + data.update[[3]]
    }
    
    # STANDARD KALMAN FILTER TRAINING
    ###### Main Loop #######
    for(i in 1:(nrow(obsMat)-1)){
      
      # load input vector
      inputVec = inputMat[i,]
      dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
      
      ###### time update - ode solve moments #######
      for(j in 1:ode_timesteps[i]){
        sol = ode.integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i])
        stateVec = sol[[1]]
        covMat = sol[[2]]
        inputVec = inputVec + dinputVec
      }
      
      ######## data update ########
      inputVec = inputMat[i+1,]
      obsVec = obsMat[i+1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        data.update <- date.update.fun(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
        stateVec <- data.update[[1]]
        covMat <- data.update[[2]]
        nll <- nll + data.update[[3]]
        # print(nll)
      }
      # end of main loop
    }
    
    # return
    return(nll)
  }
  
  # construct AD-likelihood function ----------------------------------------
  map <- lapply(private$fixed.pars, function(x) x$factor)
  parameters <- lapply(private$parameters, function(x) x$initial)
  
  # parameters <- c(parameters, list(private$initial.state$x0))
  # print(parameters)
  
  nll <- RTMB::MakeADFun(func = ekf.nll,
                         parameters = parameters,
                         map = map,
                         silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  
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
  
  # return
  return(invisible(self))
}


#######################################################
# LINEAR KALMAN FILTER # LINEAR KALMAN FILTER
#######################################################

makeADFun_lkf_rtmb = function(self, private)
{
  
  # Tape Configration ----------------------
  ConfigureTape("RTMB", self, private)
  
  # Data ----------------------------------------
  getSystemDimensions()
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  # input and obs matrix
  inputMat = as.matrix(private$data[private$input.names])
  obsMat = as.matrix(private$data[private$obs.names])
  
  # create and load state space functions
  force.ad <- private$advanced.settings$forceAD
  create.state.space.functions.for.estimation(force.ad)
  
  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  getAdjoints()
  getLossFunction()
  # getOdeSolvers()
  if(private$estimate.initial) {
    getInitialStateEstimator()
  }
  getKalmanFunctions()
  date.update.fun <- kalman.data.update.with.nll
  if(private$train.against.full.prediction){
    date.update.fun <- kalman.no.update.with.nll
  }
  
  # Timesteps, Observations, Inputs and Parameters ----------------------------
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  
  ####### Pre-Allocated Object #######
  I0 <- RTMB::diag(n.states)
  E0 <- RTMB::diag(n.obs)
  
  # detect time-variations etc ----------------------
  constant.time.diff <- FALSE
  if(var(diff(ode_timestep_size)) < 1e-15){
    constant.time.diff <- TRUE
    fixed.timestep.size <- ode_timestep_size[1]
  }
  
  # likelihood function --------------------------------
  lkf.nll = function(p){
    
    "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    
    ####### Parameters into vector #######
    parVec <- do.call(c, p[1:n.pars])
    
    ####### Neg. LogLikelihood #######
    nll <- 0
    
    inputVec <- inputMat[1,]
    ####### Compute Matrix Exponentials for 1-Step Mean and Variance #######
    # dX = A*X + B*U + G*dB
    A <- dfdx__(stateVec, parVec, inputVec) #states
    B <- dfdu__(stateVec, parVec, inputVec) #inputs
    G <- g__(stateVec, parVec, inputVec) #diffusions
    # [A B \\ 0 0]
    Phi1 <- 0.0 * RTMB::diag(n.states+n.inputs+1)
    Phi1[1:n.states,1:n.states] <- A
    Phi1[1:n.states,(n.states+1):ncol(Phi1)] <- B
    # ePhi1 <- Matrix::expm(Phi1 * fixed.timestep.size)
    ePhi1 <- as.matrix(Matrix::expm(Phi1 * fixed.timestep.size))
    # [-A GG^T \\ 0 A^t]
    Phi2 <- rbind(cbind(-A, G %*% t(G)),cbind(0*A,t(A)))
    # ePhi2 <- Matrix::expm(Phi2 * fixed.timestep.size)
    ePhi2 <- as.matrix(Matrix::expm(Phi2 * fixed.timestep.size))
    
    # A and B (for mean calculations)
    Ahat <- ePhi1[1:n.states, 1:n.states]
    Ahat_T <- t(Ahat)
    Bhat <- ePhi1[1:n.states, (n.states+1):ncol(ePhi1)]
    # V (for variance)
    Q22 <- ePhi2[(n.states+1):ncol(ePhi2), (n.states+1):ncol(ePhi2)]
    Q12 <- ePhi2[1:n.states, (n.states+1):ncol(ePhi2)]
    Vhat <- t(Q22) %*% Q12
    
    # if the time difference is not constant, we can still do good by
    # using the eigen-decomposition:
    # if(constant.time.diff){
    #   .dt <- ode.timestep.size[1]
    # then repeat the above code
    # } else {
    # we compute the matrix exponential of the eigen decompositions
    # .dt <- diff(ode.timesteps)
    # out1 <- RTMB::eigen(Phi1)
    # out2 <- RTMB::eigen(Phi2)
    # Lambda1 <- diag(out1$values)
    # Lambda2 <- diag(out2$values)
    # Q1 <- out1$vectors
    # Q1inv <- solve(Q1)
    # Q2 <- out2$vectors
    # Q2inv <- solve(Q2)
    # Then in each iteration we can do
    # exp_Phi1_dt = Q1 %*% Lambda1^dt %*% Q1inv
    # exp_Phi2_dt = Q2 %*% Lambda2^dt %*% Q2inv
    # }
    
    ####### INITIAL STATE / COVARIANCE #######
    inputVec = inputMat[1,]
    if(private$estimate.initial){
      stateVec <- f.initial.state.newton(c(parVec, inputVec))
      # covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
    }
    
    ######## (PRE) DATA UPDATE ########
    # This is done to include the first measurements in the provided data
    # We update the state and covariance based on the "new" measurement
    inputVec = inputMat[1,]
    obsVec = obsMat[1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      data.update <- date.update.fun(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
      stateVec <- data.update[[1]]
      covMat <- data.update[[2]]
      nll <- nll + data.update[[3]]
    }
    
    ###### MAIN LOOP START #######
    for(i in 1:(nrow(obsMat)-1)){
      
      ###### TIME UPDATE #######
      # augment input vector to account for constant terms in B
      inputVec = c(1, inputMat[i,])
      # Perform one-step prediction of mean and covariance
      stateVec <- Ahat %*% stateVec + Bhat %*% inputVec
      covMat <- Ahat %*% covMat %*% Ahat_T + Vhat
      
      ######## DATA UPDATE ########
      # We update the state and covariance based on the "new" measurement
      inputVec = inputMat[i+1,]
      obsVec = obsMat[i+1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        data.update <- date.update.fun(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
        stateVec <- data.update[[1]]
        covMat <- data.update[[2]]
        nll <- nll + data.update[[3]]
      }
    }
    ###### MAIN LOOP END #######
    
    # ###### RETURN #######
    return(nll)
  }
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  parameters <- lapply(private$parameters, function(x) x[["initial"]])
  nll = RTMB::MakeADFun(func = lkf.nll,
                        parameters=parameters,
                        map = lapply(private$fixed.pars, function(x) x$factor),
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}


#######################################################
#######################################################
# UKF RTMB-IMPLEMENTATION (FOR OPTIMIZATION)
#######################################################
#######################################################

makeADFun_ukf_rtmb = function(self, private)
{
  
  # Tape Configration ----------------------
  ConfigureTape("RTMB", self, private)
  
  # Data ----------------------------------------
  getSystemDimensions()
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # State Space Functions
  force.ad <- private$advanced.settings$forceAD
  create.state.space.functions.for.estimation(force.ad)
  
  # Weights
  getUkfSigmaWeights()
  
  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  getAdjoints()
  getLossFunction()
  getUkfOdeSolvers()
  if(private$estimate.initial) {
    getInitialStateEstimator()
  }
  getUkfKalmanFunctions()
  date.update.fun <- kalman.data.update.with.nll
  if(private$train.against.full.prediction){
    date.update.fun <- kalman.no.update.with.nll
  }
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  
  ####### Pre-Allocated Object #######
  I0 <- RTMB::diag(n.states)
  E0 <- RTMB::diag(n.obs)
  
  no.report.fun <- TRUE
  if(no.report.fun){
    # likelihood function --------------------------------------
    ukf.nll = function(p){
      
      # "[<-" <- RTMB::ADoverload("[<-")
      # "diag<-" <- RTMB::ADoverload("diag<-")
      # "c" <- RTMB::ADoverload("c")
      
      ####### Parameters into vector #######
      parVec <- do.call(c, p[1:n.pars])
      
      ####### Neg. LogLikelihood #######
      nll <- 0
      
      ####### INITIAL STATE / COVARIANCE #######
      inputVec = inputMat[1,]
      if(private$estimate.initial){
        stateVec <- f.initial.state.newton(c(parVec, inputVec))
        # covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
      }
      # Compute sigma points for data update
      chol.covMat <- t(Matrix::chol(covMat))
      
      # chol.covMat <- covMat
      X.sigma <- create.sigmaPoints(stateVec, chol.covMat)
      
      ######## (PRE) DATA UPDATE ########
      obsVec = obsMat[1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        data.update <- date.update.fun(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
        stateVec <- data.update[[1]]
        covMat <- data.update[[2]]
        nll <- nll + data.update[[3]]
      }
      
      ###### MAIN LOOP START #######
      for(i in 1:(nrow(obsMat)-1)){
        # Compute sigma points
        chol.covMat <- t(Matrix::chol(covMat))
        X.sigma <- create.sigmaPoints(stateVec, chol.covMat)
        
        # Inputs
        inputVec = inputMat[i,]
        dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
        
        ###### TIME UPDATE #######
        # We solve sigma points forward in time
        for(j in 1:ode_timesteps[i]){
          X.sigma <- ode.integrator(X.sigma, chol.covMat, parVec, inputVec, dinputVec, ode_timestep_size[i])
          chol.covMat <- sigma2chol(X.sigma)
          inputVec = inputVec + dinputVec
        }
        # Extract mean and covariance for data update
        stateVec <- X.sigma[,1]
        covMat <- chol.covMat %*% t(chol.covMat)
        
        ######## DATA UPDATE ########
        # We update the state and covariance based on the "new" measurement
        inputVec = inputMat[i+1,]
        obsVec = obsMat[i+1,]
        obsVec_bool = !is.na(obsVec)
        if(any(obsVec_bool)){
          data.update <- date.update.fun(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
          stateVec <- data.update[[1]]
          covMat <- data.update[[2]]
          nll <- nll + data.update[[3]]
        }
      }
      ###### MAIN LOOP END #######
      
      # ###### RETURN #######
      return(nll)
    }
  } else {
    # likelihood function --------------------------------------
    ukf.nll = function(p){
      
      # "[<-" <- RTMB::ADoverload("[<-")
      # "diag<-" <- RTMB::ADoverload("diag<-")
      # "c" <- RTMB::ADoverload("c")
      
      ####### Parameters into vector #######
      parVec <- do.call(c, p[1:n.pars])
      
      ####### Neg. LogLikelihood #######
      nll <- 0
      
      # cholesky stabilizer
      
      xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- XSigmaPrior <- XSigmaPost <- 
        list()
      # vector("list", length=ncol(obsMat))
      ####### INITIAL STATE / COVARIANCE #######
      inputVec = inputMat[1,]
      if(private$estimate.initial){
        stateVec <- f.initial.state.newton(c(parVec, inputVec))
        covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
      }
      xPrior[[1]] <- stateVec
      pPrior[[1]] <- covMat
      
      # Compute sigma points for data update
      chol.covMat <- t(Matrix::chol(covMat))
      X.sigma <- create.sigmaPoints(stateVec, chol.covMat)
      XSigmaPrior[[1]] <- X.sigma
      
      ######## (PRE) DATA UPDATE ########
      obsVec = obsMat[1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        data.update <- date.update.fun(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
        stateVec <- data.update[[1]]
        covMat <- data.update[[2]]
        nll <- nll + data.update[[3]]
      }
      xPost[[1]] <- stateVec
      pPost[[1]] <- covMat
      chol.covMat <- t(Matrix::chol(covMat))
      X.sigma <- create.sigmaPoints(stateVec, chol.covMat)
      XSigmaPost[[1]] <- X.sigma
      
      ###### MAIN LOOP START #######
      for(i in 1:(nrow(obsMat)-1)){
        
        # Inputs
        inputVec = inputMat[i,]
        dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
        
        ###### TIME UPDATE #######
        # We solve sigma points forward in time
        for(j in 1:ode_timesteps[i]){
          X.sigma <- ode.integrator(X.sigma, chol.covMat, parVec, inputVec, dinputVec, ode_timestep_size[i])
          chol.covMat <- sigma2chol(X.sigma)
          inputVec = inputVec + dinputVec
        }
        # Extract mean and covariance for data update
        stateVec <- X.sigma[,1]
        covMat <- chol.covMat %*% t(chol.covMat)
        xPrior[[i+1]] <- stateVec
        pPrior[[i+1]] <- covMat
        XSigmaPrior[[i+1]] <- X.sigma
        
        ######## DATA UPDATE ########
        # We update the state and covariance based on the "new" measurement
        inputVec = inputMat[i+1,]
        obsVec = obsMat[i+1,]
        obsVec_bool = !is.na(obsVec)
        if(any(obsVec_bool)){
          data.update <- date.update.fun(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
          stateVec <- data.update[[1]]
          covMat <- data.update[[2]]
          nll <- nll + data.update[[3]]
        }
        chol.covMat <- t(Matrix::chol(covMat))
        X.sigma <- create.sigmaPoints(stateVec, chol.covMat)
        
        xPost[[i+1]] <- stateVec
        pPost[[i+1]] <- covMat
        XSigmaPost[[i+1]] <- X.sigma
      }
      ###### MAIN LOOP END #######
      
      RTMB::REPORT(xPrior)
      RTMB::REPORT(xPost)
      RTMB::REPORT(pPrior)
      RTMB::REPORT(pPost)
      RTMB::REPORT(XSigmaPrior)
      RTMB::REPORT(XSigmaPost)
      
      # ###### RETURN #######
      return(nll)
    }
  }
  
  # construct AD-likelihood function ----------------------------------------
  
  # parameters ----------------------------------------
  map <- lapply(private$fixed.pars, function(x) x$factor)
  parameters <- lapply(private$parameters, function(x) x$initial)
  nll = RTMB::MakeADFun(func = ukf.nll,
                        parameters=parameters,
                        map = map,
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}

#######################################################
# CONSTRUCT LAPLACE MAKEADFUN WITH RTMB
#######################################################

makeADFun_laplace_rtmb = function(self, private)
{
  
  # Tape Configration ----------------------
  ConfigureTape("RTMB", self, private)
  
  # Data ----------------------------------------
  getSystemDimensions()
  
  # initial states and covariance
  stateVec <- private$initial.state$x0
  covMat <- private$initial.state$p0
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # create and load state space functions
  force.ad <- private$advanced.settings$forceAD
  create.state.space.functions.for.estimation(force.ad)
  
  
  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  getAdjoints()
  if(private$estimate.initial) {
    getInitialStateEstimator()
  }
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # indices with non-na observations
  iobs <- private$iobs
  
  # likelihood function --------------------------------------
  laplace.nll = function(p){
    
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    
    # set negative log-likelihood
    nll = 0
    
    # small identity matrix
    small_identity = diag(1e-8, nrow=n.states, ncol=n.states)
    
    # fixed effects parameter vector
    parVec <- do.call(c, p[1:n.pars])
    
    inputVec = inputMat[1,]
    if(private$estimate.initial){
      stateVec <- f.initial.state.newton(c(parVec, inputVec))
      # covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
    }
    
    # extract state random effects and fix initial condition
    stateMat <- do.call(cbind, p[(n.pars+1):length(p)])
    
    # prior contribution
    z0 <- stateMat[1,] - stateVec
    nll <- nll -  RTMB::dmvnorm(z0, Sigma=covMat, log=TRUE)
    
    ###### TIME LOOP START #######
    for(i in 1:(nrow(obsMat)-1)) {
      
      # Define inputs and use first order input interpolation
      inputVec = inputMat[i,]
      dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
      
      ###### BETWEEN TIME POINTS LOOP START #######
      for(j in 1:ode_timesteps[i]){
        
        # grab current and next state
        x_now = stateMat[ode_cumsum_timesteps[i]+j,]
        x_next = stateMat[ode_cumsum_timesteps[i]+j+1,]
        
        # compute drift (vector) and diffusion (matrix)
        f = f__(x_now, parVec, inputVec)
        g = g__(x_now, parVec, inputVec)
        inputVec = inputVec + dinputVec
        
        # assume multivariate gauss distribution according to euler-step
        # and calculate the likelihood
        z = x_next - (x_now + f * ode_timestep_size[i])
        v = (g %*% t(g) + small_identity) * ode_timestep_size[i]
        nll = nll - RTMB::dmvnorm(z, Sigma=v, log=TRUE)
      }
      ###### BETWEEN TIME POINTS LOOP END #######
    }
    ###### TIME LOOP END #######
    
    obsMat = RTMB::OBS(obsMat)
    ###### DATA UPDATE START #######
    for(i in 1:n.obs){
      iobs.vec <- iobs[[i]]
      for(j in 1:length(iobs[[i]])){
        
        # Get index where observation is available
        k <- iobs.vec[[j]]
        
        # Get corresponding input, state and observation
        inputVec = inputMat[k,]
        stateVec = stateMat[ode_cumsum_timesteps[k]+1,]
        # obsScalar = obsMat[[k,i]]
        obsScalar = obsMat[k,i]
        
        # Observation equation and variance
        h_x = h__(stateVec, parVec, inputVec)[[i]]
        hvar_x = hvar__(stateVec, parVec, inputVec)[[i]]
        
        # likelihood contribution
        nll = nll - RTMB::dnorm(obsScalar, mean=h_x, sd=sqrt(hvar_x), log=TRUE)
      }
    }
    ###### DATA UPDATE END #######
    
    # return
    return(nll)
  }
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  parameters <- c(
    lapply(private$parameters, function(x) x$initial),
    private$tmb.initial.state
  )
  map <- lapply(private$fixed.pars, function(x) x$factor)
  nll = RTMB::MakeADFun(func = laplace.nll, 
                        parameters=parameters, 
                        random=private$state.names,
                        map = map,
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}


#######################################################
# CONSTRUCT LAPLACE (UFFE TINY FORMULATION) MAKEADFUN WITH RTMB
#######################################################

makeADFun_laplace2_rtmb = function(self, private)
{
  
  # Tape Configration ----------------------
  ConfigureTape("RTMB", self, private)
  
  # Data ----------------------------------------
  getSystemDimensions()
  n.dbs <- nrow(private$tmb.initial.state) - 1
  
  # initial states and covariance
  stateVec <- private$initial.state$x0
  covMat <- private$initial.state$p0
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # create and load state space functions
  force.ad <- private$advanced.settings$forceAD
  create.state.space.functions.for.estimation(force.ad)
  
  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  getAdjoints()
  if(private$estimate.initial) {
    getInitialStateEstimator()
  }
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # indices with non-na observations
  iobs <- private$iobs
  
  # likelihood function --------------------------------------
  laplace2.nll = function(p){
    
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    
    # set negative log-likelihood
    nll = 0
    
    # small identity matrix
    # tiny = diag(1e-8, nrow=n.states, ncol=n.states)
    tiny <- 1e-5 * diag(n.states)
    
    # fixed effects parameter vector
    parVec <- do.call(c, p[1:n.pars])
    
    inputVec = inputMat[1,]
    if(private$estimate.initial){
      stateVec <- f.initial.state.newton(c(parVec, inputVec))
      # covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
    }
    
    # extract state random effects and fix initial condition
    stateMat <- do.call(cbind, p[(n.pars+1):(length(p)-1)])
    dstateMat <- RTMB::apply(stateMat, 2, diff)
    I0 <- diag(n.diffusions)
    
    # prior contribution
    z0 <- stateMat[1,] - stateVec
    nll <- nll -  RTMB::dmvnorm(z0, Sigma=covMat, log=TRUE)
    
    ###### TIME LOOP START #######
    for(i in 1:(nrow(obsMat)-1)) {
      
      # Define inputs and use first order input interpolation
      inputVec = inputMat[i,]
      dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
      
      ###### BETWEEN TIME POINTS LOOP START #######
      for(j in 1:ode_timesteps[i]){
        
        # current index in state matrix
        cur.id <- ode_cumsum_timesteps[i]+j
        
        # compute drift (vector) and diffusion (matrix)
        stateVec <- stateMat[cur.id,]
        f = f__(stateVec, parVec, inputVec)
        g = g__(stateVec, parVec, inputVec)
        inputVec = inputVec + dinputVec
        
        # Compute expected dX from Euler Maruyama
        dstateVecPred <- f * ode_timestep_size[i] + g %*% p$dB[cur.id,]
        
        # Likelihood contribution from state difference (diagonal covariance)
        z <- dstateMat[cur.id,] - dstateVecPred
        nll = nll - RTMB::dmvnorm(z, Sigma=tiny, log=TRUE)
        
        # Likelihood contribution from dBs
        nll = nll - RTMB::dmvnorm(p$dB[cur.id,], Sigma=ode_timestep_size[i]*I0, log=TRUE)
      }
      ###### BETWEEN TIME POINTS LOOP END #######
    }
    ###### TIME LOOP END #######
    
    obsMat = RTMB::OBS(obsMat)
    ##### DATA UPDATE START #######
    for(i in 1:n.obs){
      iobs.vec <- iobs[[i]]
      for(j in 1:length(iobs[[i]])){
        
        # Get index where observation is available
        k <- iobs.vec[[j]]
        
        # Get corresponding input, state and observation
        inputVec = inputMat[k,]
        stateVec = stateMat[ode_cumsum_timesteps[k]+1,]
        # obsScalar = obsMat[[k,i]]
        obsScalar = obsMat[k,i]
        
        # Observation equation and variance
        h_x = h__(stateVec, parVec, inputVec)[[i]]
        hvar_x = hvar__(stateVec, parVec, inputVec)[[i]]
        
        # likelihood contribution
        nll = nll - RTMB::dnorm(obsScalar, mean=h_x, sd=sqrt(hvar_x), log=TRUE)
      }
    }
    ###### DATA UPDATE END #######
    
    # return
    return(nll)
  }
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  # parameters ----------------------------------------
  db.len <- n.diffusions * n.dbs
  dB = matrix(numeric(n.diffusions * n.dbs), nrow=n.dbs, ncol=n.diffusions)
  parameters = c(
    lapply(private$parameters, function(x) x[["initial"]]),
    private$tmb.initial.state,
    dB = list(dB)
  )
  
  nll = RTMB::MakeADFun(func = laplace2.nll, 
                        parameters=parameters, 
                        random=c(private$state.names,"dB"),
                        map = lapply(private$fixed.pars, function(x) x$factor),
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}
