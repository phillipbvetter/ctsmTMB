makeadfun_ukf_knudsen_rtmb <- function(self, private)
{

  # Tape Configration ----------------------
  configure_ad_tape("RTMB", self, private)

  # Data ----------------------------------------
  get_sys_dims()

  # initial
  stateVec = private$algo.settings$initial.state$x0
  covMat = private$algo.settings$initial.state$p0

  # inputs
  inputMat = as.matrix(private$data[private$names$inputs])
  # observations
  obsMat = as.matrix(private$data[private$names$obs])

  # State Space Functions
  force.ad <- private$algo.settings$advanced.settings$forceAD
  create_state_space_functions_for_estimation(force.ad)

  # Weights
  get_ukf_weights()

  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  get_adjoints()
  get_loss_function()
  get_ukf_ode_solvers()
  if(private$algo.settings$estimate.initial) {
    get_initial_state_estimator()
  }
  get_ukf_update()
  data.update.fun <- kalman.data.update.with.nll
  if(private$algo.settings$train.against.full.prediction){
    data.update.fun <- kalman.no.update.with.nll
  }

  # time-steps
  ode.stepsizes = private$algo.settings$ode.stepsizes
  ode.number.of.steps = private$algo.settings$ode.number.of.steps

  ####### Pre-Allocated Object #######
  I0 <- RTMB::diag(n.states)
  E0 <- RTMB::diag(n.obs)

  # Sample sigma-points from mean and sqrt-covariance
  # [X0, X0,...,X0] + sqrt(c) [0, chol(P), -chol(P)]
  lambda <- 2
  create.sigmapoints.from.eigen <- function(stateVec, covMat, n.timesteps, dt){

    # dimensions for state augmented with brownian increments
    n.brownians <- n.timesteps * n.diffusions
    n.augmented.states <- n.states + n.brownians
    n.total <- 2*n.augmented.states+1

    # Unscented transform weights
    k <- sqrt(n.total + lambda)
    w <- c(lambda/k^2, rep(0.5/k^2, 2*n.total))

    # Construct augmented mean and covariances
    eigs <- RTMB::eigen(covMat, symmetric = TRUE)
    Cx <- diag(rep(sqrt(dt), n.augmented.states))
    Cx[1:n.states,1:n.states] <- eigs$vectors %*% diag(sqrt(eigs$values))
    Mu <- rep(0,n.augmented.states)
    Mu[1:n.states] <- stateVec




  }

  # likelihood function --------------------------------------
  ukf.nll = function(p){

    ####### Parameters into vector #######
    parVec <- do.call(c, p[1:n.pars])

    ####### Neg. LogLikelihood #######
    nll <- 0

    ####### INITIAL STATE / COVARIANCE #######
    inputVec = inputMat[1,]
    if(private$algo.settings$estimate.initial){
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
      data.update <- data.update.fun(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
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
        X.sigma <- ode.integrator(X.sigma, chol.covMat, parVec, inputVec, dinputVec, ode.stepsizes[i])
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
        data.update <- data.update.fun(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
        stateVec <- data.update[[1]]
        covMat <- data.update[[2]]
        nll <- nll + data.update[[3]]
      }
    }
    ###### MAIN LOOP END #######

    if(!force.ad){
      RTMB::REPORT(postList)
      RTMB::REPORT(priorList)
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
