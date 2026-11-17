MakeADFun_UKF2 <- function(self, private) {

  # Tape Configration ----------------------
  configure_ad_tape("RTMB", self, private)

  # Data ----------------------------------------
  unpack_system_dimensions(self, private)

  # initial
  stateVec = private$algo.settings$initial.state$x0
  covMat = private$algo.settings$initial.state$p0

  # inputs
  inputMat = as.matrix(private$data[private$names$inputs])
  # observations
  obsMat = as.matrix(private$data[private$names$obs])

  # State Space Functions
  force.ad <- private$algo.settings$advanced.settings$forceAD
  apply_forceAD_to_some_variables(force.ad)
  create_state_space_functions_for_estimation(force.ad)

  # enable report functionality?
  nllreport <- private$algo.settings$advanced.settings$nllreport
  if(nllreport) {
    statePrior <- statePosterior <- covPrior <- covPosterior <-
      nllList <- Innov <- InnovCov <- vector("list",length=nrow(obsMat))
  }

  # Weights
  get_ukf_weights()
  unpack_adjoints()
  get_loss_function()
  get_ukf_utils_knudsen()
  get_initial_state_estimator_fixed_inputs(inputMat[1,])

  lambda <- 2
  core.ids <- 1:n.sigmapoints
  last.core.id <- tail(core.ids,1)

  # time-steps
  ode.stepsizes = private$algo.settings$ode.stepsizes
  ode.number.of.steps = private$algo.settings$ode.number.of.steps

  ####### Pre-Allocated Object #######
  I0 <- RTMB::diag(n.states)

  n.augmented.states <- n.states + ode.number.of.steps[1] * n.diffusions
  n.aug.sigmapoints <- 2 * n.augmented.states + 1
  k <- sqrt(lambda + n.augmented.states)
  w <- c(lambda / (k^2), rep(1 / (2 * k^2), n.aug.sigmapoints - 1))
  X = RTMB::matrix(stateVec, nrow=n.states, ncol=n.aug.sigmapoints)

  ukf_time_update <- function(stateVec, covMat, parVec, inputVec, dinputVec,
                              ode.stepsizes, ode.number.of.steps) {
    # dimensions and indices
    # n.augmented.states <- n.states + ode.number.of.steps * n.diffusions
    # n.aug.sigmapoints <- 2 * n.augmented.states + 1
    # k <- sqrt(lambda + n.augmented.states)
    # w <- c(lambda / (k^2), rep(1 / (2 * k^2), n.aug.sigmapoints - 1))
    # X = RTMB::matrix(stateVec, nrow=n.states, ncol=n.aug.sigmapoints)

    Lambda <- sqrt_covmat_from_eigen(covMat)
    X[,core.ids] <- stateVec
    X[,core.ids[-1]] <- X[,core.ids[-1]] + k * cbind(Lambda, -Lambda)
    idx <- last.core.id + seq_len(2*n.diffusions)

    dt <- ode.stepsizes
    for(i in seq_len(ode.number.of.steps)){
      G <- g__(X[,1], parVec, inputVec)
      # Propagate the state-related sigma points (full ODE)
      X[, core.ids] <- X[, core.ids] + f.sigma(X[,core.ids], parVec, inputVec) * dt
      # We also add diffusion for the active step
      X[,idx] <- X[,1] + k * sqrt(dt) * cbind(G, -G)
      idx <- idx + 2 * n.diffusions
      # Update input
      inputVec = inputVec + dinputVec
    }

    return(ut_mean_and_cov(X, w))
  }

  k.obs <- sqrt(lambda + n.states)
  w.obs <- c(lambda / (k.obs^2), rep(1 / (2 * k.obs^2), n.sigmapoints - 1))
  # W.m.mat <- replicate(n.sigmapoints, w.obs)
  # W <- (diag(n.sigmapoints) - W.m.mat) %*% diag(w.obs) %*% t(diag(n.sigmapoints) - W.m.mat)
  ukf_measurement_update <- function(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool){
    # Calculate sigma points
    Lambda <- sqrt_covmat_from_eigen(covMat)
    X = RTMB::matrix(stateVec, nrow=n.states, ncol=n.sigmapoints)
    X[,-1] <- X[,-1] + k.obs * cbind(Lambda, -Lambda)
    # H of sigma points
    H <- h.sigma(X, parVec, inputVec)[obsVec_bool,,drop=FALSE]
    # Calculate mean, covariance and cross_covariance via Unscented Transform
    out <- ut_measurement(X,H,w.obs)
    # Innovation
    e <- obsVec[obsVec_bool] - out$my
    # Observation variance
    V = hvar__matrix(stateVec, parVec, inputVec)[obsVec_bool, obsVec_bool, drop=FALSE]
    # Measurement variance
    R <- out$Cyy + V
    # Kalman Gain
    K = t(RTMB::solve(R, t(out$Cxy)))
    # Update State/Cov
    stateVec = stateVec + K %*% e
    covMat = covMat - K %*% R %*% t(K)
    # Ensure Positive Definite
    covMat <- 0.5 * (covMat + t(covMat))
    # Negative Log-Likelihood
    nll <- loss.function(e,R)
    # Return
    return(list(stateVec, covMat, nll, e, R))
  }

  # likelihood function --------------------------------------
  ukf.nll = function(p){

    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")
    "c" <- RTMB::ADoverload("c")

    ####### Parameters into vector #######
    parVec <- do.call(c, p[1:n.pars])
    epsI0 <- 1e-8 * I0

    ####### Neg. LogLikelihood #######
    nll <- 0

    ####### INITIAL STATE / COVARIANCE #######
    inputVec = inputMat[1,]
    if(private$algo.settings$estimate.initial[1]){
      stateVec <- f.initial.state.newton_fixed_input(parVec)
      # covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
    }
    if (nllreport) {
      statePrior[[1]] <- stateVec
      covPrior[[1]] <- covMat
    }

    ######## (PRE) DATA UPDATE ########
    obsVec = obsMat[1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      data.update <- ukf_measurement_update(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool)
      stateVec <- data.update[[1]]
      covMat <- data.update[[2]]
      nll <- nll + data.update[[3]]
      if (nllreport) {
        nllList[[1]] <- data.update[[3]]
        Innov[[1]] <- data.update[[4]]
        InnovCov[[1]] <- data.update[[5]]
      }
    }
    if (nllreport) {
      statePosterior[[1]] <- stateVec
      covPosterior[[1]] <- covMat
    }

    ###### MAIN LOOP START #######
    for(i in 1:(nrow(obsMat)-1)){
      # Inputs
      inputVec = inputMat[i,]
      # zero-order or first-order hold
      if(private$algo.settings$first.order.input.hold){
        dinputVec = (inputMat[i+1,] - inputMat[i,])/ode.number.of.steps[i]
      } else {
        dinputVec = 0
      }

      ###### TIME UPDATE #######
      # We solve sigma points forward in time
      out <- ukf_time_update(stateVec, covMat, parVec, inputVec, dinputVec, ode.stepsizes[i], ode.number.of.steps[i])
      stateVec <- out[[1]]
      covMat <- out[[2]]
      if (nllreport) {
        statePrior[[i+1]] <- stateVec
        covPrior[[i+1]] <- covMat
      }

      ######## DATA UPDATE ########
      # We update the state and covariance based on the "new" measurement
      inputVec = inputMat[i+1,]
      obsVec = obsMat[i+1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        data.update <- ukf_measurement_update(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool)
        stateVec <- data.update[[1]]
        covMat <- data.update[[2]]
        nll <- nll + data.update[[3]]
        if (nllreport) {
          nllList[[i+1]] <- data.update[[3]]
          Innov[[i+1]] <- data.update[[4]]
          InnovCov[[i+1]] <- data.update[[5]]
        }
      }
      if (nllreport) {
        statePosterior[[i+1]] <- stateVec
        covPosterior[[i+1]] <- covMat
      }
    }
    ###### MAIN LOOP END #######

    if (nllreport) {
      RTMB::REPORT(statePrior)
      RTMB::REPORT(statePosterior)
      RTMB::REPORT(covPrior)
      RTMB::REPORT(covPosterior)
      RTMB::REPORT(nllList)
      RTMB::REPORT(Innov)
      RTMB::REPORT(InnovCov)
    }

    # ###### RETURN #######
    return(nll)
  }

  # construct AD-likelihood function ----------------------------------------

  # parameters ----------------------------------------
  map <- lapply(private$model$fixed.pars, function(x) x$factor)
  parameters <- lapply(private$model$parameters, function(x) x$initial)
  nll = RTMB::MakeADFun(func = ukf.nll,
                        parameters=parameters,
                        map = map,
                        silent=TRUE)

  # save objective function
  private$nll = nll

  # return
  return(invisible(self))
}

get_ukf_utils_knudsen <- function(.envir=parent.frame()){

  list2env(as.list(.envir), envir = environment())

  f.sigma <- function(X, parVec, inputVec){
    for(i in 1:n.sigmapoints){
      # Fsigma[[1,i]] <- f__(X[,i], parVec, inputVec)
      Fsigma[[i]] <- f__(X[,i], parVec, inputVec)
    }
    x <- do.call("c", Fsigma)
    dim(x) <- c(n.states, n.sigmapoints)
    return(x)
  }

  h.sigma <- function(X, parVec, inputVec){
    for(i in 1:n.sigmapoints){
      # Hsigma[[1,i]] <- h__(X[,i], parVec, inputVec)
      Hsigma[[i]] <- h__(X[,i], parVec, inputVec)
    }
    x <- do.call("c", Hsigma)
    dim(x) <- c(n.obs, n.sigmapoints)
    return(x)
  }

  # Construct sqrt matrix from covariance using eigendecomposition
  sqrt_covmat_from_eigen <- function(covMat){
    E <- RTMB::eigen(covMat, symmetric=TRUE)
    E$vectors %*% RTMB::diag(sqrt(E$values))
  }

  # Weighted mean and covariance of the columns of X
  ut_mean_and_cov <- function(X, w) {
    mx <- drop(X %*% w)
    Dx <- X - mx
    # Cxx <- (Dx * rep(w, each = nrow(X))) %*% t(Dx)
    Cxx <- Dx %*% diag(w) %*% t(Dx)
    list(m=mx, P=Cxx)
  }

  ut_measurement <- function(X, Y, w) {
    mx <- drop(X %*% w)
    my <- drop(Y %*% w)
    Dx <- (X - mx)
    Dy <- (Y - my)
    #
    Cyy <- Dy %*% diag(w) %*% t(Dy)
    Cxy <- Dx %*% diag(w) %*% t(Dy)
    # Cyy <- (Dy * rep(w, each = nrow(X))) %*% t(Dy)
    # Cxy <- (Dx * rep(w, each = nrow(X))) %*% t(Dy)
    list(my=my, Cyy=Cyy, Cxy=Cxy)
  }

  assign("f.sigma", f.sigma, envir=.envir)
  assign("h.sigma", h.sigma, envir=.envir)
  assign("sqrt_covmat_from_eigen", sqrt_covmat_from_eigen, envir=.envir)
  assign("ut_mean_and_cov", ut_mean_and_cov, envir=.envir)
  assign("ut_measurement", ut_measurement, envir=.envir)


}
