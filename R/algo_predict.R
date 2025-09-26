#######################################################
#######################################################
# PREDICTION EKF R IMPLEMENTATION
#######################################################
#######################################################

ekf_predict_r = function(parVec, self, private)
{
  
  # Data ----------------------------------------
  getSystemDimensions()
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # Utilities
  create.state.space.functions.for.filtering()
  getOdeSolvers()
  
  # Get filtered states
  filt <- ekf_filter_r(parVec, self, private)
  xPost <- filt$xPost
  pPost <- filt$pPost
  stateVec <- xPost[[1]]
  covMat <- matrix(pPost[[1]],nrow=n.states)
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  
  # prediction settings
  k.ahead <- private$n.ahead
  last.pred.index <- private$last.pred.index
  
  ####### STORAGE #######
  predMats <- lapply(1:last.pred.index, function(x) matrix(NA, nrow=k.ahead+1, ncol=n.states+n.states^2))
  
  for(i in 1:last.pred.index){ 
    stateVec <- xPost[[i]]
    covMat <- pPost[[i]]
    predMats[[i]][1,] <- c(stateVec, covMat)
    for(k in 1:k.ahead){
      inputVec = inputMat[i+k-1,]
      dinputVec = (inputMat[i+k,] - inputVec)/ode_timesteps[i+k-1]
      # Solve moment ODEs
      for(j in 1:ode_timesteps[i+k-1]){
        sol = ode.integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i+k-1])
        stateVec = sol[[1]]
        covMat = sol[[2]]
        inputVec = inputVec + dinputVec
      }
      predMats[[i]][k+1,] <- c(stateVec, covMat)
    }
  } 
  
  ####### RETURN #######
  return(invisible(predMats))
}

#######################################################
#######################################################
# PREDICTION LKF R IMPLEMENTATION
#######################################################
#######################################################

lkf_predict_r = function(parVec, self, private)
{
  
  # Data ----------------------------------------
  getSystemDimensions()
  
  # Timesteps, Observations, Inputs and Parameters ----------------------------
  ode_timestep_size = private$ode.timestep.size
  inputMat = as.matrix(private$data[private$input.names])
  obsMat = as.matrix(private$data[private$obs.names])
  
  # prediction settings
  k.ahead <- private$n.ahead
  last.pred.index <- private$last.pred.index
  
  create.state.space.functions.for.filtering()
  
  # Get filtered states
  filt <- lkf_filter_r(parVec, self, private)
  xPost <- filt$xPost
  pPost <- filt$pPost
  stateVec <- xPost[[1]]
  covMat <- matrix(pPost[[1]],nrow=n.states)
  
  ####### STORAGE #######
  predMats <- lapply(1:last.pred.index, function(x) matrix(NA,nrow=k.ahead+1,ncol=n.states+n.states^2))
  
  ####### Compute Matrix Exponentials for 1-Step Mean and Variance #######
  constant.time.diff <- FALSE
  if(var(diff(ode_timestep_size)) < 1e-15){
    constant.time.diff <- TRUE
    fixed.timestep.size <- ode_timestep_size[1]
  }
  # dX = A*X + B*U + G*dB
  inputVec <- inputMat[1,]
  A <- dfdx__(stateVec, parVec, inputVec) #states
  B <- dfdu__(stateVec, parVec, inputVec) #inputs
  G <- g__(stateVec, parVec, inputVec) #diffusions
  # [A B \\ 0 0]
  Phi1 <- 0.0 * diag(n.states+n.inputs+1)
  Phi1[1:n.states,1:n.states] <- A
  Phi1[1:n.states,(n.states+1):ncol(Phi1)] <- B
  # [-A GG^T \\ 0 A^t]
  Phi2 <- rbind(cbind(-A, G %*% t(G)),cbind(0*A,t(A)))
  ePhi1 <- Matrix::expm(Phi1 * fixed.timestep.size)
  ePhi2 <- Matrix::expm(Phi2 * fixed.timestep.size)
  
  # A and B (for mean calculations)
  Ahat <- ePhi1[1:n.states,1:n.states]
  Ahat_T <- t(Ahat)
  Bhat <- ePhi1[1:n.states, (n.states+1):ncol(ePhi1)]
  # V (for variance)
  Q22 <- ePhi2[(n.states+1):ncol(ePhi2), (n.states+1):ncol(ePhi2)]
  Q12 <- ePhi2[1:n.states, (n.states+1):ncol(ePhi2)]
  Vhat <- t(Q22) %*% Q12
  
  for(i in 1:last.pred.index){
    stateVec <- xPost[[i]]
    covMat <- matrix(pPost[[i]],nrow=n.states)
    predMats[[i]][1,] <- c(stateVec, covMat)
    ###### K-STEP AHEAD LOOP ####### 
    for(k in 1:k.ahead){
      # augment input vector to account for constant terms in B
      inputVec = c(1,inputMat[i+k-1,])
      # Perform one-step prediction of mean and covariance
      stateVec <- Ahat %*% stateVec + Bhat %*% inputVec
      covMat <- Ahat %*% covMat %*% Ahat_T + Vhat
      predMats[[i]][k+1,] <- c(stateVec, covMat)
    }
  }
  
  ####### RETURN #######
  return(invisible(predMats))
}

ukf_predict_r = function(parVec, self, private)
{
  
  # Data ----------------------------------------
  getSystemDimensions()
  
  # Timesteps, Observations, Inputs and Parameters ----------------------------
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  inputMat = as.matrix(private$data[private$input.names])
  obsMat = as.matrix(private$data[private$obs.names])
  
  # prediction settings
  k.ahead <- private$n.ahead
  last.pred.index <- private$last.pred.index
  
  create.state.space.functions.for.filtering()
  getUkfSigmaWeights()
  getUkfOdeSolvers()
  getUkfKalmanFunctions()
  
  # Get filtered states
  filt <- ukf_filter_r(parVec, self, private)
  xPost <- filt$xPost
  pPost <- filt$pPost
  stateVec <- xPost[[1]]
  covMat <- matrix(pPost[[1]],nrow=n.states)
  
  ####### STORAGE #######
  predMats <- lapply(1:last.pred.index, function(x) matrix(NA,nrow=k.ahead+1,ncol=n.states+n.states^2))
  
  # we need this to stabilize cholesky factor
  eps.chol <- 1e-10
  
  ###### TIME LOOP #######
  for(i in 1:last.pred.index){
    
    # Initial storage
    predMats[[i]][1,] <- c(stateVec, covMat)
    
    # Compute cholesky factorization
    covMat <- covMat + eps.chol * diag(n.states)
    chol.covMat <- t(Matrix::chol(covMat))
    X.sigma <- create.sigmaPoints(stateVec, chol.covMat)
    
    for(k in 1:k.ahead){ ###### K-STEP AHEAD LOOP ####### 
      inputVec = inputMat[i+k-1,]
      dinputVec = (inputMat[i+k,] - inputVec)/ode_timesteps[i+k-1]
      
      # Solve moment ODEs
      for(j in 1:ode_timesteps[i+k-1]){
        X.sigma <- ode.integrator(X.sigma, chol.covMat, parVec, inputVec, dinputVec, ode_timestep_size[i+k-1])
        chol.covMat <- sigma2chol(X.sigma)
        inputVec = inputVec + dinputVec
      }
      stateVec <- X.sigma[,1]
      covMat <- chol.covMat %*% t(chol.covMat)
      predMats[[i]][k+1,] <- c(stateVec, covMat)
    }
  }
  
  ####### RETURN #######
  return(invisible(predMats))
}

#######################################################
#######################################################
# FULL PREDICTION LAPLACE
#######################################################
#######################################################

laplace_prediction  <- function(self, private) {
  
  # We can perform a laplace prediction by not providing data to the likelihood function.
  # So we can reconstruct the AD graph with NA-data.
  # Is there a better way to re-use the existing data?
  
  # nll <- private$nll
  
  
  return(invisible(self))
  
}
