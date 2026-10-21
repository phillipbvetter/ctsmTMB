#######################################################
# EXTENDED KALMAN FILTER # EXTENDED KALMAN FILTER
#######################################################


MakeADFun_EKF = function(self, private)
{

  # Settings
  configure_ad_tape("RTMB", self, private)
  unpack_advanced_settings(self, private)

  # Unpacking data
  unpack_system_dimensions(self, private)
  unpack_data(self, private)
  unpack_state_space_functions(self, private)
  unpack_adjoints()
  get_ode_solvers()
  get_initial_state_estimator_fixed_inputs(inputMat[1,])

  # Data Update
  # data.update.fun <- kalman.data.update.with.nll
  # We can train on full predictions - this means no posterior update step, only likelihood calculation
  if (private$algo.settings$train.against.full.prediction) {
    get_loss_function(type="multivariate")
    get_ekf_update_functions()
    data.update.fun <- kalman.no.update.with.nll
  } else {
    get_loss_function(type="scalar")
    get_ekf_update_functions()
    data.update.fun <- kalman.sequential.update.taped
  }

  ####### Pre-Allocated Object #######
  I0 <- RTMB::diag(n.states)
  E0 <- RTMB::diag(n.obs)
  dinputVec <- numeric(n.inputs)
  if(nllreport) {
    statePrior <- statePosterior <- covPrior <- covPosterior <-
      nllList <- Innov <- InnovCov <- vector("list",length=nrow(obsMat))
  }

  # likelihood function --------------------------------------
  ekf.nll = function(p){

    stateVec <- RTMB::AD(stateVec)
    covMat <- RTMB::AD(covMat)
    if (private$algo.settings$advanced.settings$force.ad.on.input){
      inputMat <- RTMB::AD(inputMat)
    }

    ####### Sometimes necessary to avoid rtmb errors #######
    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")
    "c" <- RTMB::ADoverload("c")

    ####### Parameters into vector #######
    parVec <- do.call(c, p[1:n.pars])

    ####### Neg. LogLikelihood #######
    nll <- 0

    ####### Stationary Solution #######
    inputVec = inputMat[1,]
    if (private$algo.settings$estimate.initial) {
      stateVec <- f.initial.state.newton_fixed_input(parVec)
    }
    if (nllreport) {
      statePrior[[1]] <- stateVec
      covPrior[[1]] <- covMat
    }

    ######## Data Update ########
    obsVec = obsMat[1,]
    obsVec_bool = !is.na(obsVec)
    if (any(obsVec_bool)) {
      data.update <- data.update.fun(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
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

    ###### Main Loop #######
    for(i in 1:(nrow(obsMat)-1)){

      # load input vector
      inputVec = inputMat[i,]
      # zero-order or first-order hold
      if (private$algo.settings$first.order.input.hold) {
        dinputVec = (inputMat[i+1,] - inputMat[i,])/ode.number.of.steps[i]
      }

      ###### time update - ode solve moments #######
      for (j in 1:ode.number.of.steps[i]) {
        sol <- ode.integrator.taped(list(stateVec, covMat, parVec, inputVec, dinputVec, ode.stepsizes[i]))
        stateVec = sol[,1]
        covMat = sol[,-1,drop=FALSE]
        inputVec = inputVec + dinputVec
      }
      if (nllreport) {
        statePrior[[i+1]] <- stateVec
        covPrior[[i+1]] <- covMat
      }

      ######## data update ########
      inputVec = inputMat[i+1,]
      obsVec = obsMat[i+1,]
      obsVec_bool = !is.na(obsVec)
      if (any(obsVec_bool)) {
        data.update <- data.update.fun(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
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
      # end of main loop
    }

    if (nllreport) {
      RTMB::REPORT(statePrior)
      RTMB::REPORT(statePosterior)
      RTMB::REPORT(covPrior)
      RTMB::REPORT(covPosterior)
      RTMB::REPORT(nllList)
      RTMB::REPORT(Innov)
      RTMB::REPORT(InnovCov)
    }

    # return
    return(nll)
  }

  # construct AD-likelihood function ----------------------------------------
  map <- lapply(private$model$fixed.pars, function(x) x$factor)
  parameters <- lapply(private$model$parameters, function(x) x$initial)
  nll <- RTMB::MakeADFun(func = ekf.nll,
                         parameters = parameters,
                         map = map,
                         silent=TRUE)

  # save objective function
  private$nll <- nll


  # # MAP Estimation?
  # MAP_bool = 0L
  # if (!is.null(private$algo.settings$map)) {
  #   bool = self$getParameters()[,"type"] == "free"
  #   MAP_bool = 1L
  #   map_mean__ = private$algo.settings$map$mean[bool]
  #   map_cov__ = private$algo.settings$map$cov[bool,bool]
  #   map_ints__ = as.numeric(bool)
  #   sum_map_ints__ = sum(as.numeric(bool))
  # }

  # return
  return(invisible(self))
}


#######################################################
# LINEAR KALMAN FILTER # LINEAR KALMAN FILTER
#######################################################

MakeADFun_LKF = function(self, private)
{

  # Tape Configration ----------------------
  configure_ad_tape("RTMB", self, private)

  # Data ----------------------------------------
  unpack_system_dimensions(self, private)

  # initial
  stateVec = private$algo.settings$initial.state$x0
  covMat = private$algo.settings$initial.state$p0

  # input and obs matrix
  inputMat = as.matrix(private$data[private$names$inputs])
  inputMat_aug <- cbind(1, inputMat)
  obsMat = as.matrix(private$data[private$names$obs])
  # detect time-variations etc ----------------------
  data.dt <- diff(inputMat[,1])
  if (diff(range(data.dt)) > 1e-10 * data.dt[1])
    stop("The LKF requires equidistant observation times.")
  fixed.timestep.size <- data.dt[1]
  # augmented diff input matrix
  if(private$algo.settings$first.order.input.hold){
    # We add a row of zeros to match dimensions of inputMat
    dinputMat_aug <- rbind(diff(inputMat_aug),0) / fixed.timestep.size
  }

  # create and load state space functions
  force.ad <- private$algo.settings$advanced.settings$forceAD
  apply_forceAD_to_some_variables(force.ad)
  create_state_space_functions_for_estimation(force.ad)

  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  unpack_adjoints()
  get_loss_function()
  get_initial_state_estimator_fixed_inputs(inputMat[1,])
  get_ekf_update_functions()

  # data.update.fun <- kalman.data.update.with.nll
  data.update.fun <- kalman.sequential.update.taped
  if(private$algo.settings$train.against.full.prediction){
    data.update.fun <- kalman.no.update.with.nll
  }

  ####### Pre-Allocated Object #######
  I0 <- RTMB::diag(n.states)
  E0 <- RTMB::diag(n.obs)

  # AD Force
  stateVec <- RTMB::AD(stateVec, force=TRUE)
  covMat <- RTMB::AD(covMat, force=TRUE)
  # inputMat <- RTMB::AD(inputMat,force=TRUE)

  lkf.step.tape <- RTMB::MakeTape(function(x) {
    # x = list(Ahat, stateVec, inputTerm, Vhat, covMat)
    cbind(
      # Ahat %*% stateVec + B %*% inputVec (+ B1 %*% dinputVec)
      x[[1]] %*% x[[2]] + x[[3]],
      # Ahat %*% covMat %*% t(Ahat) + Vhat
      x[[1]] %*% x[[5]] %*% t(x[[1]]) + x[[4]]
    )
  }, list(diag(n.states), numeric(n.states), numeric(n.states), diag(n.states), diag(n.states))
  )

  # likelihood function --------------------------------
  lkf.nll = function(p){

    ####### Sometimes necessary to avoid rtmb errors #######
    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")
    "c" <- RTMB::ADoverload("c")

    ####### Parameters into vector #######
    parVec <- do.call(c, p[1:n.pars])

    ####### Neg. LogLikelihood #######
    nll <- 0

    inputVec <- inputMat[1,]
    ####### Compute Matrix Exponentials for 1-Step Mean and Variance #######
    # dX = A*X*dt + B*U*dt + G*dB
    A <- dfdx__(stateVec, parVec, inputVec) #states
    B <- dfdu__(stateVec, parVec, inputVec) #inputs
    G <- g__(stateVec, parVec, inputVec) #diffusions

    # first or zero order hold
    if(private$algo.settings$first.order.input.hold){
      # For FOH we construct the augmented matrix with entries
      # [ A B 0
      #   0 0 I
      #   0 0 0 ]
      Phi1 <- 0.0 * RTMB::diag(n.states+2*(n.inputs+1))
      Phi1[1:n.states,1:n.states] <- A
      Phi1[1:n.states,(n.states+1):(n.states+1+n.inputs)] <- B
      Phi1[(n.states+1):(n.states+1+n.inputs),(n.states+n.inputs+2):ncol(Phi1)] <- RTMB::diag(n.inputs+1)
      ePhi1 <- as.matrix(Matrix::expm(Phi1 * fixed.timestep.size))
      # A, B0 and B1
      Ahat <- ePhi1[1:n.states, 1:n.states, drop=FALSE]
      Ahat_T <- t(Ahat)
      Bhat <- ePhi1[1:n.states, (n.states+1):(n.states+n.inputs+1), drop=FALSE]
      B1hat <- ePhi1[1:n.states, (n.states+n.inputs+2):ncol(Phi1), drop=FALSE]
    } else {
      # For ZOH we construct the augmented matrix with entries
      # [ A B
      #   0 0 ]
      # Mean calculations
      Phi1 <- 0.0 * RTMB::diag(n.states+n.inputs+1)
      Phi1[1:n.states,1:n.states] <- A
      Phi1[1:n.states,(n.states+1):ncol(Phi1)] <- B
      ePhi1 <- as.matrix(Matrix::expm(Phi1 * fixed.timestep.size))
      #
      Ahat <- ePhi1[1:n.states, 1:n.states, drop=FALSE]
      Ahat_T <- t(Ahat)
      Bhat <- ePhi1[1:n.states, (n.states+1):ncol(ePhi1), drop=FALSE]
    }
    # V (for variance)
    # [-A GG^T \\ 0 A^t]
    Phi2 <- rbind(cbind(-A, G %*% t(G)),cbind(0*A,t(A)))
    ePhi2 <- as.matrix(Matrix::expm(Phi2 * fixed.timestep.size))
    #
    Q22 <- ePhi2[(n.states+1):ncol(ePhi2), (n.states+1):ncol(ePhi2), drop=FALSE]
    Q12 <- ePhi2[1:n.states, (n.states+1):ncol(ePhi2), drop=FALSE]
    Vhat <- t(Q22) %*% Q12

    # Compute input term for all time indices
    BU <- Bhat %*% t(inputMat_aug) # Zero Order Hold
    if (private$algo.settings$first.order.input.hold) {
      BU <- BU + B1hat %*% t(dinputMat_aug) # First Order Hold
    }

    ####### INITIAL STATE / COVARIANCE #######
    inputVec = inputMat[1,]
    if(private$algo.settings$estimate.initial){
      stateVec <- f.initial.state.newton_fixed_input(parVec)
      # covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
    }

    ######## (PRE) DATA UPDATE ########
    # This is done to include the first measurements in the provided data
    # We update the state and covariance based on the "new" measurement
    inputVec = inputMat[1,]
    obsVec = obsMat[1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      data.update <- data.update.fun(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
      stateVec <- data.update[[1]]
      covMat <- data.update[[2]]
      nll <- nll + data.update[[3]]
    }

    ###### MAIN LOOP START #######
    for(i in seq_len(nrow(obsMat)-1)){

      # TIME UPDATE
      sol <- lkf.step.tape(list(Ahat,stateVec, BU[,i], Vhat, covMat))
      stateVec <- sol[,1]
      covMat <- sol[,-1,drop=FALSE]

      ######## DATA UPDATE ########
      # We update the state and covariance based on the "new" measurement
      inputVec = inputMat[i+1,]
      obsVec = obsMat[i+1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        data.update <- data.update.fun(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
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

  parameters <- lapply(private$model$parameters, function(x) x[["initial"]])
  nll = RTMB::MakeADFun(func = lkf.nll,
                        parameters = parameters,
                        map = lapply(private$model$fixed.pars, function(x) x$factor),
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

MakeADFun_UKF = function(self, private)
{

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
  get_ukf_utils()
  get_ukf_ode_solvers()
  get_initial_state_estimator()
  get_initial_state_estimator_fixed_inputs(inputMat[1,])
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
    if(private$algo.settings$estimate.initial){
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
      chol.covMat <- t(Matrix::chol(covMat + epsI0))
      X.sigma <- create.sigmaPoints(stateVec, chol.covMat)
      data.update <- data.update.fun(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
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
      # Compute sigma points
      chol.covMat <- t(Matrix::chol(covMat + epsI0))
      X.sigma <- create.sigmaPoints(stateVec, chol.covMat)

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
      for(j in 1:ode.number.of.steps[i]){
        X.sigma <- ode.integrator(X.sigma, chol.covMat, parVec, inputVec, dinputVec, ode.stepsizes[i])
        chol.covMat <- sigma2chol(X.sigma)
        inputVec = inputVec + dinputVec
      }
      # Extract mean and covariance for data update
      stateVec <- X.sigma[,1]
      covMat <- chol.covMat %*% t(chol.covMat)
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
        data.update <- data.update.fun(X.sigma, stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
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

#######################################################
# CONSTRUCT LAPLACE MAKEADFUN WITH RTMB
#######################################################

MakeADFun_Laplace = function(self, private)
{

  ####### Sometimes necessary to avoid rtmb errors #######
  "[<-" <- RTMB::ADoverload("[<-")
  "diag<-" <- RTMB::ADoverload("diag<-")
  "c" <- RTMB::ADoverload("c")

  # Tape Configration ----------------------
  configure_ad_tape("RTMB", self, private)

  # Data ----------------------------------------
  unpack_system_dimensions(self, private)

  # initial states and covariance
  stateVec <- private$algo.settings$initial.state$x0
  covMat <- private$algo.settings$initial.state$p0

  # inputs
  inputMat = as.matrix(private$data[private$names$inputs])
  # observations
  obsMat = as.matrix(private$data[private$names$obs])

  # create and load state space functions
  force.ad <- private$algo.settings$advanced.settings$forceAD
  apply_forceAD_to_some_variables(force.ad)
  create_state_space_functions_for_estimation(force.ad)

  # enable report functionality?
  nllreport <- private$algo.settings$advanced.settings$nllreport


  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  unpack_adjoints()
  get_initial_state_estimator_fixed_inputs(inputMat[1,])
  get_initial_state_estimator()

  # time-steps
  ode.stepsizes = private$algo.settings$ode.stepsizes
  ode.number.of.steps = private$algo.settings$ode.number.of.steps
  ode_cumsum_timesteps = private$algo.settings$ode.number.of.steps.cumsum
  iobs <- private$algo.settings$iobs # indices with non-na observations

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
    if(private$algo.settings$estimate.initial){
      # stateVec <- f.initial.state.newton(c(parVec, inputVec))
      stateVec <- f.initial.state.newton_fixed_input(parVec)
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
      # zero-order or first-order hold
      if(private$algo.settings$first.order.input.hold){
        dinputVec = (inputMat[i+1,] - inputMat[i,])/ode.number.of.steps[i]
      } else {
        dinputVec = 0
      }

      ###### BETWEEN TIME POINTS LOOP START #######
      for(j in 1:ode.number.of.steps[i]){

        # grab current and next state
        x_now = stateMat[ode_cumsum_timesteps[i]+j,]
        x_next = stateMat[ode_cumsum_timesteps[i]+j+1,]

        # compute drift (vector) and diffusion (matrix)
        f = f__(x_now, parVec, inputVec)
        g = g__(x_now, parVec, inputVec)
        inputVec = inputVec + dinputVec

        # assume multivariate gauss distribution according to euler-step
        # and calculate the likelihood
        z = x_next - (x_now + f * ode.stepsizes[i])
        v = (g %*% t(g) + small_identity) * ode.stepsizes[i]
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
    lapply(private$model$parameters, function(x) x$initial),
    private$algo.settings$tmb.initial.state
  )
  map <- lapply(private$model$fixed.pars, function(x) x$factor)
  nll = RTMB::MakeADFun(func = laplace.nll,
                        parameters=parameters,
                        random=private$names$states,
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

MakeADFun_Laplace_thygesen = function(self, private)
{

  ####### Sometimes necessary to avoid rtmb errors #######
  "[<-" <- RTMB::ADoverload("[<-")
  "diag<-" <- RTMB::ADoverload("diag<-")
  "c" <- RTMB::ADoverload("c")

  # Tape Configration ----------------------
  configure_ad_tape("RTMB", self, private)

  # Data ----------------------------------------
  unpack_system_dimensions(self, private)
  n.dbs <- nrow(private$algo.settings$tmb.initial.state) - 1

  # initial states and covariance
  stateVec <- private$algo.settings$initial.state$x0
  covMat <- private$algo.settings$initial.state$p0
  # inputs
  inputMat = as.matrix(private$data[private$names$inputs])
  # observations
  obsMat = as.matrix(private$data[private$names$obs])

  # create and load state space functions
  force.ad <- private$algo.settings$advanced.settings$forceAD
  apply_forceAD_to_some_variables(force.ad=force.ad)
  create_state_space_functions_for_estimation(force.ad=force.ad)

  # enable report functionality?
  nllreport <- private$algo.settings$advanced.settings$nllreport

  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  unpack_adjoints()
  get_initial_state_estimator()
  get_initial_state_estimator_fixed_inputs(inputMat[1,])

  # time-steps
  ode.stepsizes = private$algo.settings$ode.stepsizes
  ode.number.of.steps = private$algo.settings$ode.number.of.steps
  ode_cumsum_timesteps = private$algo.settings$ode.number.of.steps.cumsum
  iobs <- private$algo.settings$iobs # indices with non-na observations

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
    if(private$algo.settings$estimate.initial){
      # stateVec <- f.initial.state.newton(c(parVec, inputVec))
      stateVec <- f.initial.state.newton_fixed_input(parVec)
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
      # zero-order or first-order hold
      if(private$algo.settings$first.order.input.hold){
        dinputVec = (inputMat[i+1,] - inputMat[i,])/ode.number.of.steps[i]
      } else {
        dinputVec = 0
      }

      ###### BETWEEN TIME POINTS LOOP START #######
      for(j in 1:ode.number.of.steps[i]){

        # current index in state matrix
        cur.id <- ode_cumsum_timesteps[i]+j

        # compute drift (vector) and diffusion (matrix)
        stateVec <- stateMat[cur.id,]
        f = f__(stateVec, parVec, inputVec)
        g = g__(stateVec, parVec, inputVec)
        inputVec = inputVec + dinputVec

        # Compute expected dX from Euler Maruyama
        dstateVecPred <- f * ode.stepsizes[i] + g %*% p$dB[cur.id,]

        # Likelihood contribution from state difference (diagonal covariance)
        z <- dstateMat[cur.id,] - dstateVecPred
        nll = nll - RTMB::dmvnorm(z, Sigma=tiny, log=TRUE)

        # Likelihood contribution from dBs
        nll = nll - RTMB::dmvnorm(p$dB[cur.id,], Sigma=ode.stepsizes[i]*I0, log=TRUE)
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
        obsScalar = obsMat[k,i]

        # Observation equation and variance
        h_x = h__(stateVec, parVec, inputVec)[[i]]
        hvar_x = hvar__(stateVec, parVec, inputVec)[[i]]

        # likelihood contribution
        nll = nll - RTMB::dnorm(obsScalar, mean=h_x, sd=sqrt(hvar_x), log=TRUE)
      }
    }
    ###### DATA UPDATE END #######

    if (nllreport) {
      RTMB::REPORT(stateMat)
    }

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
    lapply(private$model$parameters, function(x) x[["initial"]]),
    private$algo.settings$tmb.initial.state,
    dB = list(dB)
  )

  nll = RTMB::MakeADFun(func = laplace2.nll,
                        parameters=parameters,
                        random=c(private$names$states,"dB"),
                        map = lapply(private$model$fixed.pars, function(x) x$factor),
                        silent=TRUE)

  # save objective function
  private$nll = nll

  # return
  return(invisible(self))

}
