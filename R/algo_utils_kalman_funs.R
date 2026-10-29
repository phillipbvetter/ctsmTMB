get_ekf_update_functions <- function(.envir = parent.frame()){

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
    # K = covMat %*% t(C) %*% RTMB::solve(R)
    K = t(RTMB::solve(R, C %*% covMat))
    # Posterior State and Covariance
    stateVec = stateVec + K %*% e
    IKC = I0 - K %*% C
    covMat = IKC %*% covMat %*% t(IKC) + K %*% V %*% t(K)
    covMat <- 0.5 * (covMat + t(covMat))
    # Return
    return(list(stateVec, covMat, e, R))
  }

  kalman.no.update.with.nll <- function(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    # C = dhdx__(stateVec, parVec, inputVec)[obsVec_bool,,drop=FALSE]
    # Innovation
    e = (obsVec - h__(stateVec, parVec, inputVec))[obsVec_bool]
    # Observation variance
    V = hvar__matrix(stateVec, parVec, inputVec)[obsVec_bool, obsVec_bool, drop=FALSE]
    # Innovation variance
    R <- V
    # Posterior State and Covariance
    stateVec <- stateVec
    covMat <- covMat
    # Likelihood
    nll <- loss.function(e, R)
    # Return
    return(list(stateVec, covMat, nll))
  }

  kalman.data.update.with.nll <- function(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    # Remove NA's from obsVec
    # Obs Jacobian
    C = dhdx__(stateVec, parVec, inputVec)[obsVec_bool,,drop=FALSE]
    # Innovation
    e = (obsVec - h__(stateVec, parVec, inputVec))[obsVec_bool]
    # Observation variance
    V = hvar__matrix(stateVec, parVec, inputVec)[obsVec_bool, obsVec_bool, drop=FALSE]
    # Innovation variance
    C_covMat <- C %*% covMat
    R = C_covMat %*% t(C) + V
    # Kalman Gain
    K = t(RTMB::solve(R, C_covMat))
    # K = P t(C) R^{-1}
    # Posterior State and Covariance
    stateVec = stateVec + K %*% e
    # IKC = I0 - K %*% C
    # covMat = IKC %*% covMat %*% t(IKC) + K %*% V %*% t(K)
    covMat <- covMat - K %*% R %*% t(K)
    covMat <- 0.5 * (covMat + t(covMat))
    # Likelihood
    nll <- loss.function(e,R)
    # Return
    return(list(stateVec, covMat, nll, e, R))
  }

  # Three functions for taped sequential update (fastest update version so far)
  # ------------------------------------------------------------------------
  #
  # Tape 1: observation function, jacobian, obs variance (vectorized)
  obs.fun.tape <- RTMB::MakeTape(function(x) {
    # stateVec, parVec, inputVec
    # cbind(
    #   h__(x[[1]], x[[2]], x[[3]]),
    #   hvar__(x[[1]], x[[2]], x[[3]]),
    #   dhdx__(x[[1]], x[[2]], x[[3]])
    # )
    cbind(
      do.call(h__, x),
      do.call(hvar__, x),
      do.call(dhdx__, x)
    )
  }, list(numeric(n.states), numeric(n.pars), numeric(n.inputs)))

  # Tape 2: one scalar measurement update, identical for every output j
  scalar.update.tape <- RTMB::MakeTape(function(x) {
    # Unpack list
    stateVec <- x[[1]]
    covMat <- x[[2]]
    x0 <- x[[3]]
    C <- RTMB::matrix(x[[4]], nrow = 1)
    Hj <- x[[5]]
    vj <- x[[6]]
    yj <- x[[7]]
    # Linearized innovation (h(x) \approx h(x*) + dh/dx(x*) * (x-x*))
    e <- yj - Hj - drop(C %*% (stateVec - x0))
    PCt <- covMat %*% t(C)
    R <- drop(C %*% PCt) + vj
    K <- PCt / R
    # State and Cov Update
    stateVec <- stateVec + K * e
    covMat <- covMat - K %*% t(PCt)
    covMat <- 0.5 * (covMat + t(covMat))
    # Likelihood
    nll <- loss.function(e,R)
    # Return
    # It is faster to return a matrix for subsetting than to return a vector with c()
    # so we just repeat nll to match the rows of stateVec/covMat
    cbind(stateVec, covMat, rep(nll, n.states))
  }, list(numeric(n.states), diag(n.states), numeric(n.states), numeric(n.states), 0, 1, 0))

  kalman.sequential.update.taped <- function(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    nll <- 0
    idx <- which(obsVec_bool)
    stateVec_prior <- stateVec
    x <- list(stateVec, parVec, inputVec)
    obs.funs <- obs.fun.tape(x)
    H <- obs.funs[,1]
    V <- obs.funs[,2]
    C <- obs.funs[, 2+1:n.states, drop=FALSE]
    for (j in idx) {
      out <- scalar.update.tape(list(stateVec, covMat, stateVec_prior, C[j,], H[j], V[j], obsVec[j]))
      stateVec <- out[,1]
      covMat <- out[,1 + 1:n.states, drop=FALSE]
      nll <- nll + out[1, 2+n.states]
    }
    return(list(stateVec, covMat, nll))
  }
  # ------------------------------------------------------------------------

  assign("kalman.data.update", kalman.data.update, envir = .envir)
  assign("kalman.data.update.with.nll", kalman.data.update.with.nll, envir = .envir)
  assign("kalman.no.update.with.nll", kalman.no.update.with.nll, envir = .envir)
  assign("kalman.sequential.update.taped", kalman.sequential.update.taped, envir = .envir)

  return(NULL)
}


get_ukf_update <- function(.envir = parent.frame()){

  list2env(as.list(.envir), envir = environment())

  kalman.data.update <- function(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    # H of sigma points
    H <- h.sigma(X.sigma, parVec, inputVec)[obsVec_bool,,drop=FALSE]
    # Innovation
    e <- obsVec[obsVec_bool] - H %*% W.m
    # Observation variance
    V = hvar__matrix(stateVec, parVec, inputVec)[obsVec_bool, obsVec_bool, drop=FALSE]
    # Measurement variance
    R <- H %*% W %*% t(H) + V
    Cxy <- X.sigma %*% W %*% t(H)
    # Kalman gain
    K = t(RTMB::solve(R, t(Cxy)))
    # Update State/Cov
    stateVec <- stateVec + K %*% e
    covMat <- covMat - K %*% R %*% t(K)
    covMat <- 0.5 * (covMat + t(covMat))
    # C <- E %*% dhdx__(stateVec, parVec, inputVec)
    # IKC <- (I0 - K %*% C)
    # covMat = IKC %*% covMat %*% t(IKC) + K %*% V %*% t(K)
    return(list(stateVec, covMat, e, R))
  }

  kalman.data.update.with.nll <- function(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    # H of sigma points
    H <- h.sigma(X.sigma, parVec, inputVec)[obsVec_bool,,drop=FALSE]
    # Innovation
    e <- obsVec[obsVec_bool] - H %*% W.m
    # Observation variance
    V = hvar__matrix(stateVec, parVec, inputVec)[obsVec_bool, obsVec_bool, drop=FALSE]
    # Measurement variance
    R <- H %*% W %*% t(H) + V
    # Cross covariance X and Y
    Cxy <- X.sigma %*% W %*% t(H)
    # Kalman gain
    K = t(RTMB::solve(R, t(Cxy)))
    # Update State/Cov
    stateVec = stateVec + K %*% e
    covMat = covMat - K %*% R %*% t(K)
    covMat <- 0.5 * (covMat + t(covMat))
    # C <- dhdx__(stateVec, parVec, inputVec)[obsVec_bool,,drop=FALSE]
    # IKC <- (I0 - K %*% C)
    # covMat = IKC %*% covMat %*% t(IKC) + K %*% V %*% t(K)
    # Likelihood
    nll <- loss.function(e,R)
    # Return
    # return(list(stateVec, covMat, nll))
    return(list(stateVec, covMat, nll, e, R))
  }

  kalman.no.update.with.nll <- function(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0){
    # H of sigma points
    H <- h.sigma(X.sigma, parVec, inputVec)[obsVec_bool,,drop=FALSE]
    # Innovation
    e <- obsVec[obsVec_bool] - H %*% W.m
    # Observation variance
    V = hvar__matrix(stateVec, parVec, inputVec)[obsVec_bool, obsVec_bool, drop=FALSE]
    # Measurement variance
    R <- V
    stateVec <- stateVec
    covMat <- covMat
    # Likelihood
    nll <- loss.function(e, R)
    # Return
    return(list(stateVec, covMat, nll))
  }

  assign("kalman.data.update", kalman.data.update, envir = .envir)
  assign("kalman.data.update.with.nll", kalman.data.update.with.nll, envir = .envir)
  assign("kalman.no.update.with.nll", kalman.no.update.with.nll, envir = .envir)

  return(NULL)
}

get_simulate_functions <- function(.envir = parent.frame()){

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

