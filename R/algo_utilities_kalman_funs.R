
getKalmanFunctions <- function(.envir = parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  kalman.data.update <- function(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    # Remove NA's from obsVec
    y = obsVec[obsVec_bool]
    # Create permutation matrix
    E = E0[obsVec_bool,, drop=FALSE]
    # Obs Jacobian
    C = E %*% dhdx__(stateVec, parVec, inputVec)
    # Innovation
    e = y - E %*% h__(stateVec, parVec, inputVec)
    # Observation variance
    V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
    # Innovation variance
    R = C %*% covMat %*% t(C) + V
    # Kalman Gain
    K = covMat %*% t(C) %*% RTMB::solve(R)
    # Posterior State and Covariance
    stateVec = stateVec + K %*% e
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    # Likelihood
    # Return
    return(list(stateVec, covMat, e, R))
  }
  
  kalman.data.update.with.nll <- function(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    # Remove NA's from obsVec
    y = obsVec[obsVec_bool]
    # Create permutation matrix
    E = E0[obsVec_bool,, drop=FALSE]
    # Obs Jacobian
    C = E %*% dhdx__(stateVec, parVec, inputVec)
    # Innovation
    e = y - E %*% h__(stateVec, parVec, inputVec)
    # Observation variance
    V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
    # Innovation variance
    R = C %*% covMat %*% t(C) + V
    # Kalman Gain
    K = covMat %*% t(C) %*% RTMB::solve(R)
    # Posterior State and Covariance
    stateVec = stateVec + K %*% e
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    # Likelihood
    nll <- loss.function(e,R)
    # Return
    return(list(stateVec, covMat, nll))
  }
  
  assign("kalman.data.update", kalman.data.update, envir = .envir)
  assign("kalman.data.update.with.nll", kalman.data.update.with.nll, envir = .envir)
  
  return(NULL)
}


getUkfKalmanFunctions <- function(.envir = parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  kalman.data.update <- function(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    # observations
    y = obsVec[obsVec_bool]
    # permutation matrix
    E = E0[obsVec_bool,, drop=FALSE]
    # H of sigma points
    H <- h.sigma(X.sigma, parVec, inputVec)
    # Innovation
    e <- y - E %*% (H %*% W.m)
    # Observation variance
    V = hvar__matrix(stateVec, parVec, inputVec)
    # Measurement variance
    R <- E %*% (H %*% W %*% t(H) + V) %*% t(E)
    # Cross covariance X and Y
    Cxy <- X.sigma %*% W %*% t(H) %*% t(E)
    # Kalman gain
    K <- Cxy %*% RTMB::solve(R)
    # Update State/Cov
    stateVec = stateVec + K %*% e
    # Jacobian of H needed for Joseph covariance update
    C <- dhdx__(stateVec, parVec, inputVec)
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    # covMat <- covMat - K %*% R %*% t(K)
    return(list(stateVec, covMat, e, R))
  }
  
  kalman.data.update.with.nll <- function(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    # observations
    y = obsVec[obsVec_bool]
    # permutation matrix
    E = E0[obsVec_bool,, drop=FALSE]
    # H of sigma points
    H <- h.sigma(X.sigma, parVec, inputVec)
    # Innovation
    e <- y - E %*% (H %*% W.m)
    # Observation variance
    V = hvar__matrix(stateVec, parVec, inputVec)
    # Measurement variance
    R <- E %*% (H %*% W %*% t(H) + V) %*% t(E)
    # Cross covariance X and Y
    Cxy <- X.sigma %*% W %*% t(H) %*% t(E)
    # Kalman gain
    K <- Cxy %*% RTMB::solve(R)
    # Update State/Cov
    stateVec = stateVec + K %*% e
    # Jacobian of H needed for Joseph covariance update
    C <- dhdx__(stateVec, parVec, inputVec)
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    # covMat <- covMat - K %*% R %*% t(K)
    # Likelihood
    nll <- loss.function(e,R)
    # Return
    return(list(stateVec, covMat, nll))
  }
  
  assign("kalman.data.update", kalman.data.update, envir = .envir)
  assign("kalman.data.update.with.nll", kalman.data.update.with.nll, envir = .envir)
  
  return(NULL)
}

getSimulationFunctions <- function(.envir = parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  euler.maruyama.simulation <- function(stateMat, parVec, inputVec, sim.dt){
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
  
  assign("euler.maruyama.simulation", euler.maruyama.simulation, envir = .envir)
  
}

