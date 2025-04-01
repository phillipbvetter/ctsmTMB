# These data checking and setting functions are called automatically when the user
# requests an estimation, prediction or construct_nll.


###############################################################
# TOP-LAYER FUNCTION CALLING ALL OTHERS DEFINED IN THIS SCRIPT
###############################################################

check_and_set_data = function(data, self, private) {
  
  # convert to data.frame
  data = as.data.frame(data)
  
  # calculate "complex" right-hand side observation equations
  data = calculate_complex_observation_lefthandsides(data, self, private)
  
  # Check that inputs, and observations are there
  basic_data_check(data, self, private)
  
  # save data
  # only store the obs.names, not the parsed data 
  # example: if we have obs eq log(y) ~ x with name log_y, then we store log_y, but not y itself.
  private$data = data[c(private$obs.names, private$input.names)]
  
  # set timestep
  set_ode_timestep(data, self, private)
  set_simulation_timestep(data, self, private)
  
  # various calculations for tmb's laplace method
  if(private$method=="laplace"){
    set_data_for_laplace_method(data, self, private)
  }
  if(private$method=="laplace2"){
    set_data_for_laplace_method(data, self, private)
  }
  
  # Return
  return(invisible(self))
}

#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################

#######################################################
# CHECK AND SET DATA BEFORE OPTIMIZATION
#######################################################

basic_data_check = function(data, self, private) {
  
  # check if data is list or data.frame
  if (!(is.list(data) | is.data.frame(data))) {
    stop("The data should be a data.frame or a list")
  }
  
  # check if all inputs are in the data
  bool = private$input.names %in% names(data)
  if (any(!bool)){
    stop("The following required inputs were not provided in the data: 
         ", paste(private$input.names[!bool],collapse=", "))
  }
  
  # check if all (basic) observations are in the data.
  # The basic observations are not obs.names, but variables on lhs of obs equations, 
  # due to the ability to say e.g. log(y) ~ ..., obsname=log(y)
  required.obs = unique(unlist(sapply(private$obs.eqs.trans, function(ls) all.vars(ls$lhs))))
  bool = required.obs %in% names(data)
  if (any(!bool)){
    stop("The following required observations were not provided in the data: 
         ", required.obs[!bool])
  }
  
  # time vector must be increasing
  if (any(diff(data$t)<=0)) {
    ids = which(diff(data$t)<=0)
    stop(sprintf("The time-vector is non-increasing at the following indice(s) %s",paste(ids,collapse=", ")))
  }
  
  # time vector must only contain numerics
  if (any(is.na(data$t))) {
    ids = which(is.na(data$t))
    stop(sprintf("The time-vector is NA at the following indice(s) %s",paste(ids,collapse=", ")))
  }
  
  return(invisible(self))
}

#######################################################
# CALCULATE THE COMPLEX OBSERVATION NAMES
#######################################################

calculate_complex_observation_lefthandsides = function(data, self, private){
  
  # The class of the quote(log(x)) is 'call' whereas quote(x) is 'name', so complex
  # observation equations can be identified as being class 'call'
  
  # detect the complex obs lhs
  bool = as.vector(unlist(lapply(private$obs.eqs.trans, function(ls) inherits(ls$lhs, "call"))))
  
  # if there are none, return
  if(!any(bool)){
    return(data)
  }
  
  # otherwise, calculate these variables using data variables
  temp.data = list()
  for(i in seq_along(private$obs.eqs.trans)[bool]){
    
    # get name and lhs
    lhs = private$obs.eqs.trans[[i]]$lhs
    name = private$obs.eqs.trans[[i]]$name
    
    # Check if the variables are available in data
    bool = all.vars(lhs) %in% names(data)
    if(!all(bool)){
      stop("Unable to compute the observation ",name," because the following variable(s) are not in the provided data: 
           ",all.vars(lhs)[!bool])
    }
    
    # Compute the complex observation
    new.data.entry = with(data, eval(lhs))
    
    # append to the data
    temp.data[[name]] = new.data.entry
  }
  
  # concatenate data frames
  newdata = data.frame(data, temp.data)
  
  # return
  return(newdata)
}

#######################################################
# SETTINGS FOR ODE TIMESTEP 
#######################################################

set_ode_timestep = function(data, self, private){
  
  # If the required number of steps is N + epsilon or larger (e.g. 3+0.01) then increase step by 1, and reduce timestep there.
  # :::::EXAMPLE:::::
  # data$t = [0 , 1 , 2, 4.5], so data.dt = [1,1,2.5]
  # ode.timestep = 1. There are therefore [1, 1, 2.5] steps required. The last (2.5) has residual larger than 0.05 so (2.5 %% 1 = 0.5 > 0.05)
  # so we round up the number of steps there i.e. ode.N = [1 , 1 , 3]. The last entry is the important one. 
  # We take 3 steps, so for last entry, we must reduce the step-size to data.dt[3] / ode.N[3] = 2.5 / 3  = 0.88883333 
  # down from the set ode.timestep = 1
  
  
  # check for no provided time.step
  # if(is.null(private$ode.timestep)){
  # private$ode.timestep = diff(data$t)
  # }
  
  # check that ode.timestep has length 1 or at least nrow(data)-1.
  if (length(private$ode.timestep) == 1) {
    
    # Recycle to correct length
    private$ode.timestep = rep(private$ode.timestep, nrow(data)-1)
    
  }  else if (length(private$ode.timestep) == nrow(data) -1) {
    
    # do nothing
    
  } else if (length(private$ode.timestep) > nrow(data) - 1 ) {
    
    # trim if it is larger than nrow(data)-1
    private$ode.timestep = head(private$ode.timestep, nrow(data)-1)
    warning("The provided ode.timestep was longer than nrow(data) - 1, only using first nrow(data)-1 entries.")
    
  } else {
    
    # throw error
    stop("Error: The provided ode.timestep must have length 1 or nrow(data)-1")
  }
  
  
  # Data time-difference, and ode.timesteps (number of steps to take)
  data.dt = diff(data$t)
  ode.timesteps = rep(1,length(data.dt))
  
  # For the data time-gaps larger than ode.timestep, we set the time-step to the requested ode.timestep
  # For all others just take a timestep equal to the time-gap in the data (which will be smaller than ode.timestep)
  ode_timestep_size = data.dt
  bool = data.dt > private$ode.timestep
  ode_timestep_size[bool] = private$ode.timestep[bool]
  
  # where data.dt > ode.timestep, calculate the required number of steps using ode.timestep
  ode.timesteps[bool] = data.dt[bool] / private$ode.timestep[bool]
  
  # Find those indices where the number of steps must increase
  epsilon.step  = 1e-3
  residual.step.bool = (ode.timesteps %% 1) > epsilon.step 
  
  # increment these steps from N to N+1
  ode.timesteps[residual.step.bool] = ceiling(ode.timesteps[residual.step.bool]) 
  
  # round other steps down (e.g. any number less than N+epsilon becomes N)
  ode.timesteps[!residual.step.bool] = floor(ode.timesteps[!residual.step.bool])
  
  # now reduce the timestep at these locations to match the integer number of steps such that ode.timesteps * ode.timestep = data.dt
  # e.g. ode.timestep = data.dt / ode.timesteps
  ode_timestep_size[residual.step.bool] = data.dt[residual.step.bool] / ode.timesteps[residual.step.bool]
  
  # store the calculated step-size, number of steps, and cumulative number of steps
  private$ode.timestep.size = ode_timestep_size
  private$ode.timesteps = ode.timesteps
  private$ode.timesteps.cumsum = c(0,cumsum(private$ode.timesteps)) #this is used in the laplace method
  
  # return
  return(invisible(self))
}


#######################################################
# IOBS VECTOR FOR LAPLACE TO AVOID USING IS.NA
#######################################################

# The reason we want to avoid using is.na is probably that it enables us
# to use the one-step-residual function from TMB...???

set_data_for_laplace_method = function(data, self, private){
  
  # create iobs vector  ------------------------------------
  iobs = list()
  for (i in seq_along(private$obs.names)) {
    iobs[[i]] = seq_along(data$t)[!is.na(data[[private$obs.names[i]]])]
  }
  names(iobs) = paste("iobs_",private$obs.names,sep="")
  private$iobs = iobs
  
  # initial guess on random effects  ------------------------------------
  
  # set state values using only initial guess
  tempdata <- as.data.frame(matrix(0, nrow=nrow(data), ncol=private$number.of.states))
  names(tempdata) <- private$state.names
  for(i in seq_along(private$state.names)){
    tempdata[i] <- rep(private$initial.state$x0[i], nrow(data)) 
  }
  # now overwrite if initial guesses were provided in the data
  bool <- private$state.names %in% names(data)
  tempdata[private$state.names[bool]] <- data[private$state.names[bool]]
  
  # next we need to repeat these each of these state values to create 
  #intermediate points determined by the user-selected ode.timestep variable
  private$tmb.initial.state <- vector("list",length=private$number.of.states)
  for(i in seq_along(private$state.names)){
    private$tmb.initial.state[[i]] <- rep(tempdata[[i]], times=c(private$ode.timesteps,1))
  }
  names(private$tmb.initial.state) = private$state.names
  private$tmb.initial.state <- as.data.frame(private$tmb.initial.state)
  
  
  # return
  return(invisible(self))
}

#######################################################
# SETTINGS FOR ODE TIMESTEP 
#######################################################

set_simulation_timestep = function(data, self, private){
  
  # If the required number of steps is N + epsilon or larger (e.g. 3+0.01) then increase step by 1, and reduce timestep there.
  # :::::EXAMPLE:::::
  # data$t = [0 , 1 , 2, 4.5], so data.dt = [1,1,2.5]
  # ode.timestep = 1. There are therefore [1, 1, 2.5] steps required. The last (2.5) has residual larger than 0.05 so (2.5 %% 1 = 0.5 > 0.05)
  # so we round up the number of steps there i.e. ode.N = [1 , 1 , 3]. The last entry is the important one. 
  # We take 3 steps, so for last entry, we must reduce the step-size to data.dt[3] / ode.N[3] = 2.5 / 3  = 0.88883333 
  # down from the set ode.timestep = 1
  
  # check that simulation.timestep has length 1 or at least nrow(data)-1.
  if (length(private$simulation.timestep) == 1) {
    
    # Recycle to correct length
    private$simulation.timestep = rep(private$simulation.timestep, nrow(data)-1)
    
  }  else if (length(private$simulation.timestep) == nrow(data) - 1) {
    
    # do nothing
    
  } else if (length(private$simulation.timestep) > nrow(data) - 1 ) {
    
    # trim if it is larger than nrow(data)-1
    private$simulation.timestep = head(private$simulation.timestep, nrow(data)-1)
    # warning("The provided simulation.timestep was longer than nrow(data) - 1, only using first nrow(data)-1 entries.")
    
  } else {
    
    # throw error
    stop("Error: The provided simulation.timestep must have length 1 or nrow(data)-1")
  }
  
  
  # Data time-difference, and simulation.timesteps (number of steps to take)
  data.dt = diff(data$t)
  simulation.timesteps = rep(1,length(data.dt))
  
  # For the data time-gaps larger than ode.timestep, we set the time-step to the requested ode.timestep
  # For all others just take a timestep equal to the time-gap in the data (which will be smaller than ode.timestep)
  simulation_timestep_size = data.dt
  bool = data.dt > private$simulation.timestep
  simulation_timestep_size[bool] = private$simulation.timestep[bool]
  
  # where data.dt > ode.timestep, calculate the required number of steps using ode.timestep
  simulation.timesteps[bool] = data.dt[bool] / private$simulation.timestep[bool]
  
  # Find those indices where the number of steps must increase
  epsilon.step  = 1e-2
  residual.step.bool = (simulation.timesteps %% 1) > epsilon.step 
  
  # increment these steps from N to N+1
  simulation.timesteps[residual.step.bool] = ceiling(simulation.timesteps[residual.step.bool]) 
  
  # round other steps down (e.g. any number less than N+epsilon becomes N)
  simulation.timesteps[!residual.step.bool] = floor(simulation.timesteps[!residual.step.bool])
  
  # now reduce the timestep at these locations to match the integer number of steps such that simulation.timesteps * ode.timestep = data.dt
  # e.g. ode.timestep = data.dt / simulation.timesteps
  simulation_timestep_size[residual.step.bool] = data.dt[residual.step.bool] / simulation.timesteps[residual.step.bool]
  
  # store the calculated step-size, number of steps, and cumulative number of steps
  private$simulation.timestep.size = simulation_timestep_size
  private$simulation.timesteps = simulation.timesteps
  
  # return
  return(invisible(self))
}


#######################################################
# SET PARAMETERS
#######################################################
set_parameters = function(pars, silent, self, private){
  
  ###### PARAMETERS #######
  if(is.null(pars)){
    # if the estimation has been run, then use these parameters
    if(!is.null(private$fit)){
      # if(!silent) message("Using estimated parameter values")
      pars = self$getParameters(value="estimate")
    } else {
      # if(!silent) message("Using initial parameter values")
      pars = self$getParameters(value="initial")
    }
  } else {
    # check if parameters is with or without fixed parameters
    lp = length(private$parameter.names)
    fp = length(private$fixed.pars)
    
    # if not correct length give error
    if(!any(length(pars) == c(lp, lp-fp))){
      stop("Incorrect number of parameters supplied (",length(pars),"). ", "Please supply either ",lp," or ", lp-fp, ", i.e. with or without fixed parameters.")
    }
    # if parameters does not contain fixed parameters - add these
    if(length(pars)==lp-fp){
      
      # Get free and fixes ids
      par.type.free <- self$getParameters()[,"type"] == "free"
      par.type.fixed <- !par.type.free
      
      # Create new par-vector
      full.parVec <- rep(NA,lp)
      
      # Assign free and fixed pars 
      full.parVec[par.type.free] <- pars
      full.parVec[par.type.fixed] <- self$getParameters(type="fixed",value="initial")
      
      # return 
      pars <- full.parVec
    }
  }
  
  private$pars = pars
  
  return(invisible(NULL))
}

########################################################################
# SET K STEP AHEAD (DEPENDS ON PRIVATE$DATA)
########################################################################
set_k_ahead = function(k.ahead, self, private) {
  
  # check if n.ahead is positive with length 1
  if (!(is.numeric(k.ahead)) | !(length(k.ahead==1)) | !(k.ahead >= 1)) {
    stop("k.ahead must be a non-negative numeric integer")
  }
  
  # Find last prediction index to avoid exciting boundary
  last.pred.index = nrow(private$data) - k.ahead
  if(last.pred.index < 1){
    k.ahead = nrow(private$data) - 1
    last.pred.index = 1
  }
  
  # set values
  private$n.ahead = k.ahead
  private$last.pred.index = last.pred.index
  
  # return values
  return(invisible(self))
}
