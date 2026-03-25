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
    obsMat = as.matrix(private$data[private$names$obs]),
    
    # inputs
    inputMat = as.matrix(private$data[private$names$inputs]),
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # ode
    ode_solver = switch(private$algo.settings$ode.solver, euler=1, rk4=2),
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    
    # loss function
    loss_type = private$algo.settings$loss$loss,
    loss_c = private$algo.settings$loss$c,
    
    # system size
    n_states = private$dims$states,
    n_obs = private$dims$observations,
    n_inputs = private$dims$inputs,
    
    # estimate stationary levels
    estimate_stationary_initials = as.numeric(private$algo.settings$estimate.initial)
  )
  
  # MAP Estimation?
  tmb.map.data = list(
    MAP_bool = 0L
  )
  if (!is.null(private$algo.settings$map)) {
    bool = self$getParameters()[,"type"] == "free"
    tmb.map.data = list(
      MAP_bool = 1L,
      map_mean__ = private$algo.settings$map$mean[bool],
      map_cov__ = private$algo.settings$map$cov[bool,bool],
      map_ints__ = as.numeric(bool),
      sum_map_ints__ = sum(as.numeric(bool))
    )
  }
  
  # construct final data list
  data = c(tmb.data, tmb.map.data)
  
  
  # Parameters ----------------------------------------
  parVec = sapply(private$model$parameters, function(x) x[["initial"]])
  parameters = list(parVec = parVec)
  
  # Create map for fixed parameters ----------------------------------------
  pseq <- 1:private$dims$pars
  id.fixed.pars <- private$names$parameters %in% names(private$model$fixed.pars)
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
    obsMat = as.matrix(private$data[private$names$obs]),
    
    # inputs
    inputMat = as.matrix(private$data[private$names$inputs]),
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # loss function
    loss_type = private$algo.settings$loss$loss,
    loss_c = private$algo.settings$loss$c,
    
    # system size
    n_states = private$dims$states,
    n_obs = private$dims$observations,
    n_inputs = private$dims$inputs,
    
    # estimate stationary levels
    estimate_stationary_initials = as.numeric(private$algo.settings$estimate.initial)
  )
  
  # MAP Estimation?
  tmb.map.data = list(
    MAP_bool = 0L
  )
  if (!is.null(private$algo.settings$map)) {
    bool = self$getParameters()[,"type"] == "free"
    tmb.map.data = list(
      MAP_bool = 1L,
      map_mean__ = private$algo.settings$map$mean[bool],
      map_cov__ = private$algo.settings$map$cov[bool,bool],
      map_ints__ = as.numeric(bool),
      sum_map_ints__ = sum(as.numeric(bool))
    )
  }
  
  # construct final data list
  data = c(tmb.data, tmb.map.data)
  
  
  # Parameters ----------------------------------------
  parVec = sapply(private$model$parameters, function(x) x[["initial"]])
  parameters = list(parVec = parVec)
  
  # Create map for fixed parameters ----------------------------------------
  pseq <- 1:private$dims$pars
  id.fixed.pars <- private$names$parameters %in% names(private$model$fixed.pars)
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
    obsMat = as.matrix(private$data[private$names$obs]),
    
    # inputs
    inputMat = as.matrix(private$data[private$names$inputs]),
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # ode
    ode_solver = switch(private$algo.settings$ode.solver, euler=1, rk4=2),
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    
    # loss function
    loss_type = private$algo.settings$loss$loss,
    loss_c = private$algo.settings$loss$c,
    
    # ukf parameters
    ukf_pars = private$algo.settings$ukf.hyperpars,
    
    # system size
    n_states = private$dims$states,
    n_obs = private$dims$observations,
    n_inputs = private$dims$inputs,
    
    # estimate stationary levels
    estimate_stationary_initials = as.numeric(private$algo.settings$estimate.initial)
  )
  
  # MAP Estimation?
  tmb.map.data = list(
    MAP_bool = 0L
  )
  if (!is.null(private$algo.settings$map)) {
    bool = self$getParameters()[,"type"] == "free"
    tmb.map.data = list(
      MAP_bool = 1L,
      map_mean__ = private$algo.settings$map$mean[bool],
      map_cov__ = private$algo.settings$map$cov[bool,bool],
      map_ints__ = as.numeric(bool),
      sum_map_ints__ = sum(as.numeric(bool))
    )
  }
  
  # construct final data list
  data = c(tmb.data, tmb.map.data)
  
  
  # Parameters ----------------------------------------
  parVec = sapply(private$model$parameters, function(x) x[["initial"]])
  parameters = list(parVec = parVec)
  
  # Create map for fixed parameters ----------------------------------------
  pseq <- 1:private$dims$pars
  id.fixed.pars <- private$names$parameters %in% names(private$model$fixed.pars)
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
