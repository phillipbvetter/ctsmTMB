#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

construct_makeADFun = function(self, private){
  
  # TMB::openmp(max=TRUE, autopar=TRUE, DLL=private$modelname.with.method)
  
  if(any(private$method == c("ekf_cpp","ukf"))){
    comptime <- system.time(construct_kalman_cpp_makeADFun(self, private))
  }
  
  if(private$method == "ekf"){
    comptime <- system.time(construct_rtmb_ekf_makeADFun(self, private))
  }
  
  if(private$method == "ekf_rcpp"){
    comptime <- system.time(rcpp_ekf_estimate(self, private))
  }
  
  if(private$method=="laplace"){
    comptime <- system.time(construct_rtmb_laplace_makeADFun(self, private))
  }
  
  comptime = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
  if(!private$silent) message("...took: ", comptime, " seconds.")
  
  return(invisible(self))
}

#######################################################
# CONSTRUCT KALMAN MAKEADFUN
#######################################################
construct_kalman_cpp_makeADFun = function(self, private){
  
  ################################################
  # Data
  ################################################
  
  # add mandatory entries to data
  tmb.data = list(
    
    # methods and purpose
    # estimation_method = switch(private$method, ekf = 1, ukf = 2),
    ode_solver = private$ode.solver,
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # time-steps
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    
    # loss function
    loss_function = private$loss$loss,
    loss_threshold_value = private$loss$c,
    tukey_loss_parameters = private$tukey.pars,
    
    # system size
    number_of_state_eqs = private$number.of.states,
    number_of_obs_eqs = private$number.of.observations,
    number_of_diffusions = private$number.of.diffusions,
    
    # estimate stationary levels
    estimate_stationary_initials = as.numeric(private$estimate.initial),
    initial_variance_scaling = private$initial.variance.scaling,
    
    # parameter bounds
    par_lb = sapply(private$parameters, function(x) x$lower),
    par_ub = sapply(private$parameters, function(x) x$lower),
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names])
  )
  
  # unscented parameters
  ukf_hyperpars_list = list()
  if(private$method=="ukf")
  {
    ukf_hyperpars_list = list(
      ukf_alpha = private$ukf_alpha,
      ukf_beta = private$ukf_beta,
      ukf_kappa = private$ukf_kappa
    )
  }
  
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
  data = c(tmb.data, tmb.map.data, ukf_hyperpars_list)
  
  ################################################
  # Parameters
  ################################################
  
  parameters = lapply(private$parameters, function(x) x[["initial"]]) # Initial parameter values
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = TMB::MakeADFun(data = data,
                       parameters = parameters,
                       map = lapply(private$fixed.pars, function(x) x$factor),
                       DLL = private$modelname.with.method,
                       silent = TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}

#######################################################
# CONSTRUCT KALMAN MAKEADFUN WITH RTMB
#######################################################

construct_rtmb_ekf_makeADFun = function(self, private)
{
  
  ################################################
  # Data
  ################################################
  
  # methods and purpose
  ode_solver = private$ode.solver
  
  # initial
  stateVec = cbind(private$initial.state$x0)
  covMat = private$initial.state$p0
  
  # loss function
  loss_function = private$loss$loss
  loss_threshold_value = private$loss$c
  tukey_loss_parameters = private$tukey.pars
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # system size
  number_of_state_eqs = private$number.of.states
  number_of_obs_eqs = private$number.of.observations
  number_of_diffusions = private$number.of.diffusions
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # MAP Estimation?
  MAP_bool = 0L
  if (!is.null(private$map)) {
    bool = self$getParameters()[,"type"] == "free"
    MAP_bool = 1L
    map_mean__ = private$map$mean[bool]
    map_cov__ = private$map$cov[bool,bool]
    map_ints__ = as.numeric(bool)
    sum_map_ints__ = sum(as.numeric(bool))
  }
  
  ################################################
  # Initial Parameters
  ################################################
  
  parameters = lapply(private$parameters, function(x) x[["initial"]])
  
  ################################################
  # Define functions
  ################################################
  
  # user-defined functions ---------------------------
  for(i in seq_along(private$rtmb.function.strings)){
    eval(parse(text=private$rtmb.function.strings[[i]]))
  }
  
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
      is.seq = list()
      dy <- c()
      for(i in 1:n){
        for(j in 1:n){
          id.seq <- 1+(j-1)*n + (n^2+1) * (1:n-1)
          id.seq <- id.seq + (i-1) * n^3
          dy <- c(dy, sum(y[id.seq]))
        }
      }
      return(dy)
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
      is.seq = list()
      dy <- c()
      for(i in 1:n){
        for(j in 1:n){
          id.seq <- j + (n^3+n) * (1:n-1)
          id.seq <- id.seq + (i-1) * n^2
          dy <- c(dy, sum(y[id.seq]))
        }
      }
      return(dy)
    },
    name = "kron_right")
  
  # # ODE Solver
  ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt, ode_solver){
    
    # Initials
    X0 = stateVec
    P0 = covMat
    
    # Explicit Forward Euler
    if(ode_solver==1){
      X1 = X0 + f__(stateVec, parVec, inputVec) * dt
      P1 = P0 + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
    }
    
    # Classical 4th Order Runge-Kutta Method
    if(ode_solver==2){
      
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
    }
    return(invisible(list(X1,P1)))
  }
  
  # Covariance ODE 1-Step
  cov_ode_1step = function(covMat, stateVec, parVec, inputVec){
    A = dfdx__(stateVec, parVec, inputVec)
    G = g__(stateVec, parVec, inputVec)
    AcovMat = A %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  # error function in terms of pnorm
  erf = function(x){
    y <- sqrt(2) * x
    2*RTMB::pnorm(y)-1
  }
  
  NUMBER_OF_STATES <- private$number.of.states
  NUMBER_OF_OBSERVATIONS <- private$number.of.observations
  NUMBER_OF_PARS <- private$number.of.pars
  
  ekf.nll = function(p){
    
    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")
    "c" <- RTMB::ADoverload("c")
    
    ####### INITIALIZATION #######
    xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
    xPrior[[1]] <- stateVec
    pPrior[[1]] <- covMat
    I0 <- diag(NUMBER_OF_STATES)
    E0 <- diag(NUMBER_OF_OBSERVATIONS)
    nll <- 0
    parVec <- do.call(c, p[1:NUMBER_OF_PARS])
    
    ####### INITIAL VALUES FROM STATIONARY SOLUTION #######
    # inputVec = inputMat[1,]
    # F <- RTMB::MakeTape(function(x) f__(x, inputVec, parVec)^2, numeric(NUMBER_OF_STATES))
    # X0 <- F$newton(1:NUMBER_OF_STATES)
    # X1 <- X0(numeric(0))
    # stateVec <- X1
    # A <- dfdx__(stateVec, parVec, inputVec)
    # G <- g__(stateVec, parVec, inputVec)
    # Q <- G %*% t(G)
    # I <- diag(rep(1, nrow(A)))
    # P <- kron_left(A) + kron_right(A)
    # X <- -RTMB::solve(P, as.numeric(Q))
    # covMat <- X
    
    ######## FIRST DATA UPDATE ######## 
    obsVec = obsMat[1,]
    inputVec = inputMat[1,]
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
      nll = nll - RTMB::dmvnorm(e, Sigma=R, log=TRUE)
      # Update State/Cov
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
      # Store innovation and covariance
      Innovation[[1]] = e
      InnovationCovariance[[1]] = R
    }
    xPost[[1]] <- stateVec
    pPost[[1]] <- covMat
    
    # ###### TIME LOOP START #######
    for(i in 1:(nrow(obsMat)-1)){
      
      # # Define inputs and use first order input interpolation
      inputVec = inputMat[i,]
      dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
      
      ###### TIME UPDATE - ODE SOLVER #######
      for(j in 1:ode_timesteps[i]){
        sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i], ode_solver)
        stateVec = sol[[1]]
        covMat = sol[[2]]
        inputVec = inputVec + dinputVec
      }
      xPrior[[i+1]] = stateVec
      pPrior[[i+1]] = covMat
      
      ######## DATA UPDATE ######## 
      obsVec = obsMat[i+1,]
      inputVec = inputMat[i+1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        y = obsVec[obsVec_bool]
        E = E0[obsVec_bool,, drop=FALSE] #permutation matrix with rows removed
        C = E %*% dhdx__(stateVec, parVec, inputVec)
        e = y - E %*% h__(stateVec, parVec, inputVec)
        V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
        R = C %*% covMat %*% t(C) + V
        K = covMat %*% t(C) %*% RTMB::solve(R)
        # Likelihood Contribution
        nll = nll - RTMB::dmvnorm(e, Sigma=R, log=TRUE)
        # Update State/Cov
        stateVec = stateVec + K %*% e
        covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
        # Store innovation and covariance
        Innovation[[i+1]] = e
        InnovationCovariance[[i+1]] = R
      }
      xPost[[i+1]] = stateVec
      pPost[[i+1]] = covMat
    }
    
    # ###### MAXIMUM A POSTERIOR #######
    
    ##### REPORT #######
    RTMB::REPORT(Innovation)
    RTMB::REPORT(InnovationCovariance)
    RTMB::REPORT(xPrior)
    RTMB::REPORT(xPost)
    RTMB::REPORT(pPrior)
    RTMB::REPORT(pPost)
    
    # ###### RETURN #######
    return(nll)
  }
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = RTMB::MakeADFun(func = ekf.nll,
                        parameters=parameters,
                        map = lapply(private$fixed.pars, function(x) x$factor),
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}

#######################################################
# CONSTRUCT LAPLACE C++ MAKEADFUN
#######################################################
construct_laplace_cpp_makeADFun = function(self, private){
  
  ################################################
  # Data
  ################################################
  
  # add mandatory entries to data
  tmb.data = list(
    
    # time-steps
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    ode_timesteps_cumsum = private$ode.timesteps.cumsum,
    
    # system size
    number_of_state_eqs = private$number.of.states,
    number_of_obs_eqs = private$number.of.observations,
    number_of_diffusions = private$number.of.diffusions,
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names])
  )
  
  # construct final data list
  data = c(tmb.data, private$iobs)
  
  ################################################
  # Parameters
  ################################################
  
  parameters = c(
    lapply(private$parameters, function(x) x[["initial"]]),
    list(stateMat = do.call(cbind,private$tmb.initial.state.for.parameters))
  )
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = TMB::MakeADFun(data = data,
                       parameters = parameters,
                       random="stateMat",
                       map = lapply(private$fixed.pars, function(x) x$factor),
                       DLL = private$modelname.with.method,
                       silent = TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}

#######################################################
# CONSTRUCT LAPLACE MAKEADFUN WITH RTMB
#######################################################

construct_rtmb_laplace_makeADFun = function(self, private)
{
  
  ################################################
  # Data
  ################################################
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # system size
  number_of_state_eqs = private$number.of.states
  number_of_obs_eqs = private$number.of.observations
  number_of_diffusions = private$number.of.diffusions
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # indices in state parameter vectors corresponding to indices in observations / inputs
  # add 1 because too lazy to change private$iobs from 0-index to 1-indexed.
  iobs = lapply(private$iobs,function(x) x+1)
  
  ################################################
  # Parameters
  ################################################
  
  parameters = c(
    lapply(private$parameters, function(x) x[["initial"]]),
    private$tmb.initial.state.for.parameters
  )
  
  ################################################
  # Functions
  ################################################
  
  for(i in seq_along(private$rtmb.function.strings)){
    eval(parse(text=private$rtmb.function.strings[[i]]))
  }
  
  # error function in terms of pnorm
  erf = function(x){
    y <- sqrt(2) * x
    2*RTMB::pnorm(y)-1
  }

  NUMBER_OF_STATES <- private$number.of.states
  NUMBER_OF_OBSERVATIONS <- private$number.of.observations
  NUMBER_OF_PARS <- private$number.of.pars

  laplace.nll = function(p){

    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")
    "c" <- RTMB::ADoverload("c")

    # set negative log-likelihood
    nll = 0

    # small identity matrix
    small_identity = diag(1e-8, nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_STATES)

    # fixed effects parameter vector
    parVec <- do.call(c, p[1:NUMBER_OF_PARS])
    
    # extract state random effects
    stateMat <- do.call(cbind, p[(NUMBER_OF_PARS+1):length(p)])


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
    for(i in 1:NUMBER_OF_OBSERVATIONS){
      for(j in 1:length(iobs[[i]])){

        # Get index where observation is available
        k = iobs[[i]][j]

        # Get corresponding input, state and observation
        inputVec = inputMat[k,]
        stateVec = stateMat[ode_cumsum_timesteps[k]+1,]
        obsScalar = obsMat[k,i]

        # Observation equation and variance
        h_x = h__(stateVec, parVec, inputVec)[i]
        hvar_x = hvar__(stateVec, parVec, inputVec)[i]

        # likelihood contribution
        nll = nll - RTMB::dnorm(obsScalar, mean=h_x, sd=sqrt(hvar_x), log=TRUE)
      }}
    ###### DATA UPDATE END #######

    # return
    return(invisible(nll))
  }

  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  
  nll = RTMB::MakeADFun(func = laplace.nll, 
                        parameters=parameters, 
                        random=private$state.names,
                        map = lapply(private$fixed.pars, function(x) x$factor),
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}


#######################################################
# OPTIMISE AD FUN
#######################################################

perform_estimation = function(self, private) {
  
  # Parameter Bounds
  lower.parameter.bound = unlist(lapply(private$free.pars, function(par) par$lower))
  upper.parameter.bound = unlist(lapply(private$free.pars, function(par) par$upper))
  if(private$unconstrained.optim){
    lower.parameter.bound = -Inf
    upper.parameter.bound = Inf
  }
  
  # IF METHOD IS KALMAN FILTER
  if (any(private$method==c("ekf","ukf","ekf_cpp"))) {
    
    # use function, gradient and hessian
    if (private$use.hessian) {
      comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                           objective = private$nll$fn,
                                                                           gradient = private$nll$gr,
                                                                           hessian = private$nll$he,
                                                                           lower = lower.parameter.bound,
                                                                           upper = upper.parameter.bound,
                                                                           control=private$control.nlminb))
      )
      # or just function and gradient
    } else {
      comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                           objective = private$nll$fn,
                                                                           gradient = private$nll$gr,
                                                                           lower = lower.parameter.bound,
                                                                           upper = upper.parameter.bound,
                                                                           control=private$control.nlminb))
      )
    }
    
  }
  
  # IF METHOD IS TMB
  if (any(private$method == c("laplace"))) {
    comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                         objective = private$nll$fn,
                                                                         gradient = private$nll$gr,
                                                                         lower = lower.parameter.bound,
                                                                         upper = upper.parameter.bound,
                                                                         control=private$control.nlminb))
    )
  }
  
  # DID THE OPTIMIZATION FAIL?
  if (inherits(opt,"try-error")) {
    
    message("The optimisation failed due to the following error: \n\n\t",opt)
    
    if(stringr::str_detect(opt,"NA/NaN")){
      
      message("You should consider the following to circumvent the error:
              1. Reduce the ODE step-size via argument 'ode.timestep'.
              2. Run the optimization with / without the hessian via argument 'use.hessian'.
              3. Change / explore parameter initial values.
              4. Extract the function handlers with the 'construct_nll' method and investigate outputs, or try other optimizers.
              5. Change the optimization tolerances for 'nlminb' with the 'control' argument.")
      
    }
    
    private$opt = NULL
    
    return(invisible(self))
  }
  
  # store optimization object
  private$opt = opt
  
  # extract maxmimum gradient component, and format computation time to 5 digits
  outer_mgc = max(abs(private$nll$gr(opt$par)))
  comp.time = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
  
  # print convergence and timing result
  if(!private$silent){
    # if(outer_mgc > 1){
    # message("BEWARE: THE MAXIMUM GRADIENT COMPONENT APPEARS TO BE LARGE ( > 1 ) - THE FOUND OPTIMUM MIGHT BE INVALID.")
    # }
    message("\t Optimization finished!:
            Elapsed time: ", comp.time, " seconds.
            The objective value is: ",format(opt$objective,scientific=T),"
            The maximum gradient component is: ",format(outer_mgc,digits=2,scientific=T),"
            The convergence message is: ", opt$message,"
            Iterations: ",opt$iterations,"
            Evaluations: Fun: ",opt$evaluations["function"]," Grad: ",opt$evaluations[["gradient"]],"
            See stats::nlminb for available tolerance/control arguments."
    )
  }
  
  # For TMB method: run sdreport
  if (any(private$method== c("laplace"))) {
    if(!private$silent) message("Calculating random effects standard deviation...")
    comptime = system.time(
      # private$sdr <- RTMB::sdreport(private$nll, getJointPrecision=T)
      private$sdr <- TMB::sdreport(private$nll, getJointPrecision=T)
      # NOTE: The state covariances can be retrived by inverting sdr$jointPrecision
      # but this takes very long time. Should it be an option?
    )
    comptime = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
    if(!private$silent) message("...took: ", comptime, " seconds.")
  }
  
  # return
  return(invisible(self))
}




#######################################################
# Create return fit after calling estimate 
#######################################################

create_fit = function(self, private, calculate.laplace.onestep.residuals) {
  
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  # clear fit
  private$fit = NULL
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  # store the object in fit - these gives access to e.g.
  # private$fit$.__object__ = self$clone()
  
  ################################################
  # FOR KALMAN FILTERS
  ################################################
  
  if (any(private$method == c("ekf_cpp","ukf"))) {
    
    
    ################################################
    # BASICS
    ################################################
    
    # objective value
    private$fit$nll = private$opt$objective
    
    # objective gradient
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
    
    # compute hessian
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
    
    # parameter estimates
    private$fit$par.fixed = private$opt$par
    
    # parameter std. error and full covariance by hessian inversion
    if (is.null(private$fit$nll.hessian)){
      private$fit$cov.fixed <- NaN * diag(length(private$fit$par.fixed))
      private$fit$sd.fixed <- rep(NaN,length(private$fit$par.fixed))
    } else {
      
      ####### OPTION 0 #######
      # This option tries to invert the full hessain
      
      full.hess = private$fit$nll.hessian
      covar <- try_withWarningRecovery(solve(full.hess))
      if (inherits(covar,"try-error")) {
        private$fit$cov.fixed = NaN * full.hess
        private$fit$sd.fixed = rep(NaN,length(private$fit$par.fixed))
      } else {
        private$fit$cov.fixed <- covar
        private$fit$sd.fixed <- try_withWarningRecovery(sqrt(diag(private$fit$cov.fixed)))
      }
      
      ####### OPTION 1 #######
      # This option filters out all row/cols of the hessian where the diagonal
      # element is smaller than some set threshold
      
      min.diag = 1e-8
      keep.ids = !(diag(full.hess) < min.diag)
      if(inherits(covar,"try-error") && any(keep.ids)){
        reduced.hess <- full.hess[keep.ids, keep.ids]
        cov <- try_withWarningRecovery(solve(reduced.hess))
        private$fit$cov.fixed <- NaN * full.hess
        private$fit$sd.fixed <- rep(NaN,length(private$fit$par.fixed))
        if (!inherits(cov,"try-error")) {
          private$fit$cov.fixed[keep.ids, keep.ids] <- cov
          private$fit$sd.fixed[keep.ids] = try_withWarningRecovery(sqrt(diag(cov)))
        }
      }
      
      ####### OPTION 2 #######
      # This option tries to invert the hessian by recursively removing row/cols
      # of the hessian where the diagonal are smallest
      if (inherits(covar,"try-error")) {
        private$fit$cov.fixed <- NaN * full.hess
        private$fit$sd.fixed <- rep(NaN,length(private$fit$par.fixed))
        ids = sort(diag(full.hess), index.return=T)$ix
        for(i in seq_along(ids)){
          id <- ids[1:i]
          hess <- full.hess[-id,-id]
          cov <- try_withWarningRecovery(solve(hess))
          if (!inherits(cov,"try-error")){
            private$fit$cov.fixed[-id,-id] <- cov
            private$fit$sd.fixed[-id] <- try_withWarningRecovery(sqrt(diag(cov)))
            break
          }
        }
      }
    }
    
    ################################################
    # STATES, RESIDUALS, OBSERVATIONS ETC.
    ################################################
    
    # Extract reported items from nll
    rep = private$nll$report()
    
    # likelihood terms
    private$fit$nll.terms = rep$nll_report
    
    # Prior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, do.call(rbind,rep$xPrior)))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPrior,diag)))))
    
    colnames(temp.states) = c("t", private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$prior = as.data.frame(temp.states)
    private$fit$states$sd$prior = as.data.frame(temp.sd)
    private$fit$states$cov$prior = rep$pPrior
    names(private$fit$states$cov$prior) = paste("t = ",private$data$t,sep="")
    
    # Posterior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, do.call(rbind,rep$xPost)))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPost,diag)))))
    colnames(temp.states) = c("t",private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$posterior = as.data.frame(temp.states)
    private$fit$states$sd$posterior = as.data.frame(temp.sd)
    private$fit$states$cov$posterior = rep$pPost
    names(private$fit$states$cov$posterior) = paste("t = ",private$data$t,sep="")
    
    # Residual
    rowNAs = as.matrix(!is.na(private$data[private$obs.names]))
    sumrowNAs = rowSums(rowNAs)
    
    innovation = rep$Innovation
    innovation.cov = rep$InnovationCovariance
    
    temp.res = matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
    temp.var =  matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
    
    # do.call(rbind, lapply(rep$Innovation, "length<-", private$m))
    for (i in seq_along(private$data$t)) {
      if (sumrowNAs[i] > 0) {
        temp.res[i,rowNAs[i,]] = innovation[[i]]
        temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
      }
    }
    temp.res = cbind(private$data$t, temp.res)
    temp.sd = cbind(private$data$t, try_withWarningRecovery(sqrt(temp.var)))
    
    names(innovation.cov) = paste("t = ",private$data$t[-1],sep="")
    
    # should we remove the empty matrices?
    # innovation.cov = innovation.cov[sumrowNAs!=0]
    
    colnames(temp.res) = c("t",private$obs.names)
    colnames(temp.sd) = c("t",private$obs.names)
    private$fit$residuals$mean = as.data.frame(temp.res)
    private$fit$residuals$sd = as.data.frame(temp.sd)
    private$fit$residuals$normalized = as.data.frame(temp.res)
    private$fit$residuals$normalized[,-1] = private$fit$residuals$normalized[,-1]/temp.sd[,-1]
    private$fit$residuals$cov = innovation.cov
    
    
    ### Observations ###
    
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
    
    obs.df.prior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.prior)})
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
    
    obs.df.posterior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.posterior)})
    )
    
    private$fit$observations$mean$prior = data.frame(t=private$data$t, obs.df.prior)
    private$fit$observations$mean$posterior = data.frame(t=private$data$t, obs.df.posterior)
    
    # observation variances
    # The observation variance (to first order) is: y = h(x) + e -> var(y) = dhdx var(x) dhdx^T + var(e)
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
    
    # Evaluate prior and posterior variance
    list.of.parameters <- c(
      as.list(private$fit$par.fixed),
      lapply(private$fixed.pars, function(x) x$initial)
    )
    # prior
    obsvar.prior <- list()
    for(i in seq_along(private$data$t)){
      list.of.states <- as.list(private$fit$states$mean$prior[i,-1,drop=F])
      state.cov <- list(xCov = private$fit$states$cov$prior[[i]])
      list.of.inputs <- as.list(private$data[i,-1,drop=F])
      obsvar.prior[[i]] <- eval(expr = yCov,
                                envir = c(list.of.parameters, 
                                          list.of.states,
                                          state.cov,
                                          list.of.inputs
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
      list.of.states <- as.list(private$fit$states$mean$posterior[i,-1,drop=F])
      state.cov <- list(xCov = private$fit$states$cov$posterior[[i]])
      list.of.inputs <- as.list(private$data[i,-1,drop=F])
      obsvar.post[[i]] <- eval(expr = yCov,
                               envir = c(list.of.parameters, 
                                         list.of.states,
                                         state.cov,
                                         list.of.inputs
                               ))
    }
    names(obsvar.post) <- names(private$fit$states$cov$posterior)
    private$fit$observations$cov$posterior <- obsvar.post
    obs.sd.post <- cbind(private$fit$states$mean$posterior["t"], do.call(rbind, lapply(obsvar.post, diag)))
    rownames(obs.sd.post) <- NULL
    names(obs.sd.post) <- c("t",private$obs.names)
    private$fit$observations$sd$posterior <- obs.sd.post
    
    # t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
  }
  
  if (any(private$method == c("ekf"))) {
    
    
    ################################################
    # BASICS
    ################################################
    
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
    
    # hessian
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
    private$fit$par.fixed = private$opt$par
    
    # parameter std. error and full covariance by hessian inversion
    if(!is.null(private$fit$nll.hessian)){
      
      # Step 1 - invert full hessian
      temp.hessian = private$fit$nll.hessian
      covariance = try(solve(temp.hessian), silent=T)
      
      private$fit$cov.fixed = covariance
      private$fit$sd.fixed = try_withWarningRecovery(sqrt(diag(covariance)))
      
      if(inherits(private$fit$sd.fixed,"try-error")){
        private$fit$sd.fixed = rep(NA,length(private$fit$par.fixed))
      }
      if(inherits(private$fit$cov.fixed,"try-error")){
        private$fit$cov.fixed = matrix(NA,nrow=length(private$fit$par.fixed),ncol=length(private$fit$par.fixed))
      }
      
      # Options 1 - If the above fails, remove all row/cols where the diagonal
      # elements in small than min.diag
      min.diag = 1e-8
      keep.ids = !(diag(temp.hessian) < min.diag)
      if(inherits(covariance,"try-error") && any(keep.ids)){
        
        covariance = temp.hessian[keep.ids, keep.ids]
        covariance = try(solve(covariance), silent=T)
        
        sd.fixed = rep(NA,length(private$fit$par.fixed))
        sd.fixed[keep.ids] = try_withWarningRecovery(sqrt(diag(covariance)))
        private$fit$sd.fixed = sd.fixed
        
        cov.fixed = temp.hessian * NA
        cov.fixed[keep.ids, keep.ids] = covariance
        private$fit$cov.fixed = cov.fixed
      }
      
      # Option 2 - Recursive remove the smallest parameter
      # ids = sort(diag(temp.hessian), index.return=T)$ix
      # i = 0
      # while(inherits(hess,"try-error") && i < length(private$fit$par.fixed)){
      #   i = i + 1
      #   covariance = try(solve(temp.hessian[-ids[1:i],-ids[1:i]]), silent=T)
      # }
      
    }
    
    ################################################
    # STATES, RESIDUALS, OBSERVATIONS ETC.
    ################################################
    
    # Extract reported items from nll
    rep = private$nll$report()
    
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
    
    # Residual
    # rowNAs = as.matrix(!is.na(do.call(cbind, private$data[private$obs.names]))[-1,])
    rowNAs = as.matrix(!is.na(private$data[private$obs.names])[-1,])
    sumrowNAs = rowSums(rowNAs)
    
    innovation = rep$Innovation
    innovation.cov = rep$InnovationCovariance
    innovation[[1]] = NULL
    innovation.cov[[1]] = NULL
    
    temp.res = matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
    temp.var =  matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
    
    # do.call(rbind, lapply(rep$Innovation, "length<-", private$m))
    for (i in seq_along(private$data$t[-1])) {
      if (sumrowNAs[i] > 0) {
        temp.res[i,rowNAs[i,]] = innovation[[i]]
        temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
      }
    }
    temp.res = cbind(private$data$t[-1], temp.res)
    temp.sd = cbind(private$data$t[-1], try_withWarningRecovery(sqrt(temp.var)))
    
    names(innovation.cov) = paste("t = ",private$data$t[-1],sep="")
    
    # should we remove the empty matrices?
    # innovation.cov = innovation.cov[sumrowNAs!=0]
    
    colnames(temp.res) = c("t",private$obs.names)
    colnames(temp.sd) = c("t",private$obs.names)
    private$fit$residuals$mean = as.data.frame(temp.res)
    private$fit$residuals$sd = as.data.frame(temp.sd)
    private$fit$residuals$normalized = as.data.frame(temp.res)
    private$fit$residuals$normalized[,-1] = private$fit$residuals$normalized[,-1]/temp.sd[,-1]
    private$fit$residuals$cov = innovation.cov
    
    
    # Observations
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
    
    # observation variances
    # The observation variance (to first order) is: y = h(x) + e -> var(y) = dhdx var(x) dhdx^T + var(e)
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
    
    # Evaluate prior and posterior variance
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
    
    # t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
    
  }
  
  ################################################
  # FOR TMB
  ################################################
  
  if (any(private$method == c("laplace"))) {
    
    n <- private$number.of.states
    m <- private$number.of.observations
    
    # Objective and Gradient
    private$fit$nll = private$opt$objective
    private$fit$nll.gradient = private$sdr$gradient.fixed
    names(private$fit$nll.gradient) = names(private$free.pars)
    
    # Parameter Estimate
    private$fit$par.fixed = private$opt$par
    private$fit$sd.fixed = diag(private$sdr$cov.fixed)
    
    # Parameter Covariance
    private$fit$cov.fixed = private$sdr$cov.fixed
    
    # Posterior States (Smoothed)
    temp = cbind(matrix(private$sdr$par.random, ncol=n),
                 matrix(sqrt(private$sdr$diag.cov.random), ncol=n))[private$ode.timesteps.cumsum+1, ]
    temp.states = cbind(private$data$t, matrix(temp[,1:n],nrow=length(private$data$t)))
    temp.sd = cbind(private$data$t, matrix(temp[,(n+1):(2*n)],nrow=length(private$data$t)))
    private$fit$states$mean$smoothed = as.data.frame(temp.states)
    private$fit$states$sd$smoothed = as.data.frame(temp.sd)
    colnames(private$fit$states$sd$smoothed) = c("t",private$state.names)
    colnames(private$fit$states$mean$smoothed) = c("t",private$state.names)
    
    # Residuals
    rowNAs = as.matrix(!is.na(do.call(cbind,private$data[private$obs.names]))[-1,])
    sumrowNAs = rowSums(rowNAs)
    
    # compute one-step residuals
    if(calculate.laplace.onestep.residuals){
      message("Calculating one-step ahead residuls...")
      private$fit$residuals = RTMB::oneStepPredict(private$nll,
                                                   observation.name="obsMat",
                                                   method="oneStepGaussian",
                                                   trace=TRUE)
    }
    
    # Parameter t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
    
    ################################################
    # OBSERVATIONS
    ################################################
    
    listofvariables.smoothed = c(
      # states
      as.list(private$fit$states$mean$smoothed[-1]),
      # estimated free parameters 
      as.list(private$fit$par.fixed),
      # fixed parameters
      lapply(private$fixed.pars, function(x) x$initial),
      # inputs
      as.list(private$data)
    )
    
    # Observations
    obs.df.smoothed = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.smoothed)})
    )
    private$fit$observations$mean$smoothed = data.frame(t=private$data$t, obs.df.smoothed)
    
    # Observation Variance
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
    
    # Evaluate prior and posterior variance
    list.of.parameters <- c(
      as.list(private$fit$par.fixed),
      lapply(private$fixed.pars, function(x) x$initial)
    )
    # prior
    # obsvar.smooth <- list()
    # for(i in seq_along(private$data$t)){
    #   obsvar.smooth[[i]] <- eval(expr = yCov,
    #                             envir = c(list.of.parameters, 
    #                                       as.list(private$fit$states$mean$smoothed[i,-1]),
    #                                       list(xCov = private$fit$states$cov$smoothed[[i]]),
    #                                       as.list(private$data[i,-1])
    #                             ))
    # }
    # names(obsvar.smooth) <- names(private$fit$states$cov$prior)
    # private$fit$observations$cov$prior <- obsvar.prior
    # obs.sd.prior <- cbind(private$fit$states$mean$prior["t"], do.call(rbind, lapply(obsvar.prior, diag)))
    # rownames(obs.sd.prior) <- NULL
    # names(obs.sd.prior) <- c("t",private$obs.names)
    # private$fit$observations$sd$prior <- obs.sd.prior
    
  }
  
  # Set S3 class
  # clone private and return fit
  self.clone <- self$clone()
  private$fit$private = self.clone$.__enclos_env__$private
  
  class(private$fit) = "ctsmTMB.fit"
  
  return(invisible(self))
  
}

