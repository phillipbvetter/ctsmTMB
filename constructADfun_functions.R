# Functions for construction likelihood functions using TMB / RTMB

#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

construct_makeADFun = function(self, private){
  
  # TMB::openmp(max=TRUE, autopar=TRUE, DLL=private$modelname.with.method)
  
  if(any(private$method == c("ekf","ukf"))){
    comptime <- system.time(construct_kalman_cpp_makeADFun(self, private))
  }
  
  if(private$method == "ekf_rtmb"){
    comptime <- system.time(construct_rtmb_ekf_makeADFun(self, private))
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
    estimation_method = switch(private$method, ekf = 1, ukf = 2),
    ode_solver = private$ode.solver,
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # time-steps
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    ode_cumsum_timesteps = private$ode.timesteps.cumsum,
    
    # loss function
    loss_function = private$loss$loss,
    loss_threshold_value = private$loss$c,
    tukey_loss_parameters = private$tukey.pars,
    
    # system size
    number_of_state_eqs = private$number.of.states,
    number_of_obs_eqs = private$number.of.observations,
    number_of_diffusions = private$number.of.diffusions,
    
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
  data = c(tmb.data, private$iobs, tmb.map.data, ukf_hyperpars_list)
  
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
  
  # ODE Solver
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
  
  for(i in seq_along(private$rtmb.function.strings)){
    eval(parse(text=private$rtmb.function.strings[[i]]))
  }
  eval(parse(text=private$rtmb.nll.strings$ekf))

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
  eval(parse(text=private$rtmb.nll.strings$laplace))
  
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
