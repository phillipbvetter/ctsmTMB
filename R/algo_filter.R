#######################################################
#######################################################
# EKF FILTER (FOR REPORTING)
#######################################################
#######################################################

ekf_filter_r = function(parVec, self, private)
{
  
  # Data ----------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  create.state.space.functions.for.filtering()

  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  getOdeSolvers()
  if(estimate.initial) {
    getInitialStateEstimator()
  }
  getKalmanFunctions()
  
  # Timesteps, Observations, Inputs and Parameters ----------------------------
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  inputMat = as.matrix(private$data[private$input.names])
  obsMat = as.matrix(private$data[private$obs.names])
  
  ####### STORAGE #######
  xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
  
  ####### Neg. LogLikelihood #######
  nll <- 0
  
  ####### Pre-Allocated Object #######
  I0 <- diag(n.states)
  E0 <- diag(n.obs)
  
  ####### INITIAL STATE / COVARIANCE #######
  inputVec = inputMat[1,]
  if(estimate.initial){
    stateVec <- f.initial.state.newton(c(parVec, inputVec))
    # covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
  }
  xPrior[[1]] <- stateVec
  pPrior[[1]] <- covMat
  
  ######## (PRE) DATA UPDATE ########
  # This is done to include the first measurements in the provided data
  # We update the state and covariance based on the "new" measurement
  obsVec = obsMat[1,]
  obsVec_bool = !is.na(obsVec)
  if(any(obsVec_bool)){
    data.update <- kalman.data.update(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
    stateVec <- data.update[[1]]
    covMat <- data.update[[2]]
    Innovation[[1]] = data.update[[3]]
    InnovationCovariance[[1]] = data.update[[4]]
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
      sol = ode.integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i])
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
      data.update <- kalman.data.update(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
      stateVec <- data.update[[1]]
      covMat <- data.update[[2]]
      Innovation[[i+1]] = data.update[[3]]
      InnovationCovariance[[i+1]] = data.update[[4]]
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
# LKF FILTER R-IMPLEMENTATION
#######################################################
#######################################################
lkf_filter_r = function(parVec, self, private)
{
  
  # Data ----------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  create.state.space.functions.for.filtering()

  if(estimate.initial) {
    getInitialStateEstimator()
  }
  getKalmanFunctions()
  
  # Timesteps, Observations, Inputs and Parameters ----------------------------
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  inputMat = as.matrix(private$data[private$input.names])
  obsMat = as.matrix(private$data[private$obs.names])
  
  # detect time-variations etc ----------------------
  constant.time.diff <- FALSE
  if(var(diff(ode_timestep_size)) < 1e-15){
    constant.time.diff <- TRUE
    fixed.timestep.size <- ode_timestep_size[1]
  }
  
  ####### STORAGE #######
  xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
  
  ####### Neg. LogLikelihood #######
  # nll <- 0
  
  ####### Pre-Allocated Object #######
  I0 <- diag(n.states)
  E0 <- diag(n.obs)
  
  ####### Compute Matrix Exponentials for 1-Step Mean and Variance #######
  # dX = A*X + B*U + G*dB
  inputVec <-  inputMat[1,]
  A <- dfdx__(stateVec, parVec, inputVec) #states
  B <- dfdu__(stateVec, parVec, inputVec) #inputs
  G <- g__(stateVec, parVec, inputVec) #diffusions
  H <- dhdx__(stateVec,parVec, inputVec) #observation
  V0 <- hvar__matrix(stateVec, parVec, inputVec) #observation variance
  
  # [A B \\ 0 0]
  Phi1 <- 0.0 * diag(n.states+n.inputs+1)
  Phi1[1:n.states,1:n.states] <- A
  Phi1[1:n.states,(n.states+1):ncol(Phi1)] <- B
  # [-A GG^T \\ 0 A^t]
  Phi2 <- rbind(cbind(-A,G %*% t(G)),cbind(0*A,t(A)))
  ePhi1 <- as.matrix(Matrix::expm(Phi1 * fixed.timestep.size))
  ePhi2 <- as.matrix(Matrix::expm(Phi2 * fixed.timestep.size))
  
  # A and B (for mean calculations)
  Ahat <- ePhi1[1:n.states,1:n.states]
  Ahat_T <- t(Ahat)
  Bhat <- ePhi1[1:n.states, (n.states+1):ncol(ePhi1)]
  # V (for variance)
  Q22 <- ePhi2[(n.states+1):ncol(ePhi2), (n.states+1):ncol(ePhi2)]
  Q12 <- ePhi2[1:n.states, (n.states+1):ncol(ePhi2)]
  Vhat <- t(Q22) %*% Q12
  
  ####### INITIAL STATE / COVARIANCE #######
  # The state/covariance is either given by user or obtained from solving the
  # stationary mean, and then solving for the covariance.
  # In principle these are coupled equations, but believe that root-finding
  # both simultaneously can lead to likelihood blow-up.
  inputVec = inputMat[1,]
  if(estimate.initial){
    stateVec <- f.initial.state.newton(c(parVec, inputVec))
    covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
  }
  xPrior[[1]] <- stateVec
  pPrior[[1]] <- covMat
  
  ######## (PRE) DATA UPDATE ########
  # This is done to include the first measurements in the provided data
  # We update the state and covariance based on the "new" measurement
  obsVec = obsMat[1,]
  obsVec_bool = !is.na(obsVec)
  if(any(obsVec_bool)){
    data.update <- kalman.data.update(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
    stateVec <- data.update[[1]]
    covMat <- data.update[[2]]
    Innovation[[1]] = data.update[[3]]
    InnovationCovariance[[1]] = data.update[[4]]
  }
  xPost[[1]] <- stateVec
  pPost[[1]] <- covMat
  
  
  ###### MAIN LOOP START #######
  for(i in 1:(nrow(obsMat)-1)){
    
    ###### TIME UPDATE #######
    # augment input vector to account for constant terms in B
    inputVec = c(1,inputMat[i,])
    # Perform one-step prediction of mean and covariance
    stateVec <- Ahat %*% stateVec + Bhat %*% inputVec
    covMat <- Ahat %*% covMat %*% Ahat_T + Vhat
    # Store prior predictions
    xPrior[[i+1]] = stateVec
    pPrior[[i+1]] = covMat
    
    ######## DATA UPDATE ########
    # We update the state and covariance based on the "new" measurement
    inputVec = inputMat[i+1,]
    obsVec = obsMat[i+1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      data.update <- kalman.data.update(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
      stateVec <- data.update[[1]]
      covMat <- data.update[[2]]
      Innovation[[i+1]] = data.update[[3]]
      InnovationCovariance[[i+1]] = data.update[[4]]
    }
    xPost[[i+1]] = stateVec
    pPost[[i+1]] = covMat
  }
  ###### MAIN LOOP END #######
  
  ####### RETURN #######
  returnlist <- list(
    # nll=nll,
    xPost = xPost,
    pPost = pPost,
    xPrior = xPrior,
    pPrior = pPrior,
    Innovation = Innovation,
    InnovationCovariance = InnovationCovariance
  )
  
  return(invisible(returnlist))
}

#######################################################
#######################################################
# EKF R-IMPLEMENTATION (FOR REPORTING)
#######################################################
#######################################################

ukf_filter_r = function(parVec, self, private)
{
  
  # Data ----------------------------------------
  estimate.initial <- private$estimate.initial
  getSystemDimensions()
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  # inputs and observations
  inputMat = as.matrix(private$data[private$input.names])
  obsMat = as.matrix(private$data[private$obs.names])
  
  create.state.space.functions.for.filtering()
  # Fsigma <- array(list(),c(1,n.sigmapoints))
  # Hsigma <- array(list(),c(1,n.sigmapoints))
  getUkfSigmaWeights()
  getUkfOdeSolvers()
  if(private$estimate.initial) {
    getInitialStateEstimator()
  }
  getUkfKalmanFunctions()
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  
  ####### STORAGE #######
  xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
  
  ####### Neg. LogLikelihood #######
  nll <- 0
  
  ####### Pre-Allocated Object #######
  I0 <- diag(n.states)
  E0 <- diag(n.obs)
  
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
  
  ######## (PRE) DATA UPDATE ########
  obsVec = obsMat[1,]
  obsVec_bool = !is.na(obsVec)
  if(any(obsVec_bool)){
    data.update <- kalman.data.update(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
    stateVec <- data.update[[1]]
    covMat <- data.update[[2]]
    Innovation[[1]] = data.update[[3]]
    InnovationCovariance[[1]] = data.update[[4]]
  }
  xPost[[1]] <- stateVec
  pPost[[1]] <- covMat
  
  ###### MAIN LOOP START #######
  for(i in 1:(nrow(obsMat)-1)){
    
    # Compute cholesky factorization
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
    xPrior[[i+1]] = stateVec
    pPrior[[i+1]] = covMat
    
    ######## DATA UPDATE ########
    # We update the state and covariance based on the "new" measurement
    inputVec = inputMat[i+1,]
    obsVec = obsMat[i+1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      data.update <- kalman.data.update(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
      stateVec <- data.update[[1]]
      covMat <- data.update[[2]]
      Innovation[[i+1]] = data.update[[3]]
      InnovationCovariance[[i+1]] = data.update[[4]]
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
                     InnovationCovariance = InnovationCovariance
  )
  
  return(invisible(returnlist))
}
