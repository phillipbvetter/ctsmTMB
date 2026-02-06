#######################################################
#######################################################
# EKF TMB-IMPLEMENTATION (FOR OPTIMIZATION)
#######################################################
#######################################################

makeADFun_ekf_tmb = function(self, private){
  
  # Tape Configration ----------------------
  configure_ad_tape("TMB", self, private)
  
  # Data ----------------------------------------
  
  # add mandatory entries to data
  tmb.data = list(
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names]),
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # ode
    ode_solver = private$ode.solver,
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    
    # loss function
    loss_type = private$loss$loss,
    loss_c = private$loss$c,
    
    # system size
    n_states = private$number.of.states,
    n_obs = private$number.of.observations,
    n_inputs = private$number.of.inputs,
    
    # estimate stationary levels
    estimate_stationary_initials = as.numeric(private$estimate.initial)
  )
  
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
  data = c(tmb.data, tmb.map.data)
  
  
  # Parameters ----------------------------------------
  parVec = sapply(private$parameters, function(x) x[["initial"]])
  parameters = list(parVec = parVec)
  
  # Create map for fixed parameters ----------------------------------------
  pseq <- 1:private$number.of.pars
  id.fixed.pars <- private$parameter.names %in% names(private$fixed.pars)
  pseq[id.fixed.pars] <- NA
  map <- list(parVec = factor(pseq))
  
  # Create AD-likelihood function ---------------------------------------
  nll <- TMB::MakeADFun(data = data,
                        parameters = parameters,
                        map = map,
                        DLL = private$modelname.with.method,
                        silent = TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}

#######################################################
#######################################################
# LKF TMB-IMPLEMENTATION (FOR OPTIMIZATION)
#######################################################
#######################################################

makeADFun_lkf_tmb = function(self, private){
  
  # Tape Configration ----------------------
  configure_ad_tape("TMB", self, private)
  
  # Data ----------------------------------------
  
  # add mandatory entries to data
  tmb.data = list(
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names]),
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # loss function
    loss_type = private$loss$loss,
    loss_c = private$loss$c,
    
    # system size
    n_states = private$number.of.states,
    n_obs = private$number.of.observations,
    n_inputs = private$number.of.inputs,
    
    # estimate stationary levels
    estimate_stationary_initials = as.numeric(private$estimate.initial)
  )
  
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
  data = c(tmb.data, tmb.map.data)
  
  
  # Parameters ----------------------------------------
  parVec = sapply(private$parameters, function(x) x[["initial"]])
  parameters = list(parVec = parVec)
  
  # Create map for fixed parameters ----------------------------------------
  pseq <- 1:private$number.of.pars
  id.fixed.pars <- private$parameter.names %in% names(private$fixed.pars)
  pseq[id.fixed.pars] <- NA
  map <- list(parVec = factor(pseq))
  
  # Create AD-likelihood function ---------------------------------------
  nll <- TMB::MakeADFun(data = data,
                        parameters = parameters,
                        map = map,
                        DLL = private$modelname.with.method,
                        silent = TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}

#######################################################
#######################################################
# UKF TMB-IMPLEMENTATION (FOR OPTIMIZATION)
#######################################################
#######################################################

makeADFun_ukf_tmb = function(self, private){
  
  # Tape Configration ----------------------
  configure_ad_tape("TMB", self, private)
  
  # Data ----------------------------------------
  
  # add mandatory entries to data
  tmb.data = list(
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names]),
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # ode
    ode_solver = private$ode.solver,
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    
    # loss function
    loss_type = private$loss$loss,
    loss_c = private$loss$c,
    
    # ukf parameters
    ukf_pars = private$ukf.hyperpars,
    
    # system size
    n_states = private$number.of.states,
    n_obs = private$number.of.observations,
    n_inputs = private$number.of.inputs,
    
    # estimate stationary levels
    estimate_stationary_initials = as.numeric(private$estimate.initial)
  )
  
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
  data = c(tmb.data, tmb.map.data)
  
  
  # Parameters ----------------------------------------
  parVec = sapply(private$parameters, function(x) x[["initial"]])
  parameters = list(parVec = parVec)
  
  # Create map for fixed parameters ----------------------------------------
  pseq <- 1:private$number.of.pars
  id.fixed.pars <- private$parameter.names %in% names(private$fixed.pars)
  pseq[id.fixed.pars] <- NA
  map <- list(parVec = factor(pseq))
  
  # Create AD-likelihood function ---------------------------------------
  nll <- TMB::MakeADFun(data = data,
                        parameters = parameters,
                        map = map,
                        DLL = private$modelname.with.method,
                        silent = TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}
