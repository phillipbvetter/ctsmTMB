###############################################################
# TOP-LAYER FUNCTION CALLING ALL OTHERS DEFINED IN THIS SCRIPT
###############################################################

check_and_set_data = function(data, self, private) {

  # AD-Check: Did data change or can we re-use old nll
  if(private$procedure %in% c("estimation", "likelihood")){

    # check if we need to rebuild the ad graph due to data changes
    check_for_data_rebuild(data, self, private)

    # Exit if data is unchanged
    if(!private$rebuild$data){
      return(invisible(self))
    }

    # If new data: revert flags and store for check next time
    private$rebuild$data <- FALSE
    private$rebuild$ad <- TRUE
    private$old.data$entry.data <- data
  }

  if(!private$algo.settings$silent) message("Checking data...")

  # Check that inputs, and observations are there
  basic_data_check(data, self, private)

  # calculate "complex" right-hand side observation equations
  data <- calculate_complex_observation_lefthandsides(data, self, private)

  # save data
  # only store the obs.names, not the parsed data
  # example: if we have obs eq log(y) ~ x with name log_y, then we store log_y, but not y itself.
  private$data = data[c(private$names$obs, private$names$inputs)]

  # set timesteps
  compute_timestep("ode", data, self, private)
  compute_timestep("simulation", data, self, private)

  # various calculations for laplace method
  set_data_for_laplace_method(data, self, private)

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

  # check if data is a data.frame
  if (!(is.data.frame(data))) {
    stop("The data should be a data.frame.")
  }

  # check if all inputs are in the data
  bool = private$names$inputs %in% names(data)
  if (any(!bool)){
    stop("The following required inputs were not provided in the data:
         ", paste(private$names$inputs[!bool],collapse=", "))
  }

  # check if all (basic) observations are in the data.
  # The basic observations are not obs.names, but variables on lhs of obs equations,
  # due to the ability to say e.g. log(y) ~ ..., obsname=log(y)
  required.obs = unique(unlist(sapply(private$model$obs.eqs.trans, function(ls) all.vars(ls$lhs))))
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
  bool = as.vector(unlist(lapply(private$model$obs.eqs.trans, function(ls) inherits(ls$lhs, "call"))))

  # if there are none, return
  if(!any(bool)){
    return(data)
  }

  # otherwise, calculate these variables using data variables
  temp.data = list()
  for(i in seq_along(private$model$obs.eqs.trans)[bool]){

    # get name and lhs
    lhs = private$model$obs.eqs.trans[[i]]$lhs
    name = private$model$obs.eqs.trans[[i]]$name

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

compute_timestep = function(type, data, self, private, epsilon.step = 1e-3){

  # If the required number of steps is N + epsilon or larger (e.g. 3+0.01) then increase step by 1, and reduce timestep there.
  # :::::EXAMPLE:::::
  # data$t = [0 , 1 , 2, 4.5], so data.dt = [1,1,2.5]
  # timestep = 1. There are therefore [1, 1, 2.5] steps required. The last (2.5) has residual larger than epsilon so (2.5 %% 1 = 0.5 > epsilon)
  # so we round up the number of steps there i.e. N = [1, 1, 3]. The last entry is the important one.
  # We take 3 steps, so for last entry, we must reduce the step-size to data.dt[3] / N[3] = 2.5 / 3 = 0.88883333

  n <- nrow(data) - 1
  dt <- private$algo.settings[[paste0(type, ".timestep")]]

  # check that dt has length 1 or nrow(data)-1
  if (length(dt) == 1) {
    dt <- rep(dt, n)
  } else if (length(dt) > n) {
    warning(sprintf("The provided %s.timestep was longer than nrow(data) - 1, only using first nrow(data)-1 entries.", type))
    dt <- head(dt, n)
  } else if (length(dt) < n) {
    warning(sprintf("The provided %s.timestep was shorter than nrow(data) - 1, only using first entry.", type))
    dt <- rep(dt[1],n)
  }

  # Data time-differences and number of steps to take
  data.dt <- diff(data$t)
  timesteps <- rep(1, n)
  timestep.size <- data.dt

  # For gaps larger than dt, use dt as step-size; otherwise use given dt and modify
  bool <- data.dt > dt
  timestep.size[bool] <- dt[bool]
  timesteps[bool] <- data.dt[bool] / dt[bool]
  # now round up where residual exceeds epsilon, round down otherwise
  residual.bool <- (timesteps %% 1) > epsilon.step
  timesteps[residual.bool]  <- ceiling(timesteps[residual.bool])
  timesteps[!residual.bool] <- floor(timesteps[!residual.bool])

  # Adjust step-size so that timesteps * timestep.size == data.dt exactly
  timestep.size[residual.bool] <- data.dt[residual.bool] / timesteps[residual.bool]

  # Store results
  private$algo.settings[[paste0(type, ".timestep.size")]] <- timestep.size
  private$algo.settings[[paste0(type, ".timesteps")]]     <- timesteps
  if (type == "ode") {
    private$algo.settings$ode.timesteps.cumsum <- c(0, cumsum(timesteps))
  }

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
  for (i in seq_along(private$names$obs)) {
    iobs[[i]] = seq_along(data$t)[!is.na(data[[private$names$obs[i]]])]
  }
  names(iobs) = paste("iobs_",private$names$obs,sep="")
  private$algo.settings$iobs = iobs

  # initial guess on random effects  ------------------------------------

  # set state values using only initial guess
  tempdata <- as.data.frame(matrix(0, nrow=nrow(data), ncol=private$dims$states))
  names(tempdata) <- private$names$states
  for(i in seq_along(private$names$states)){
    tempdata[i] <- rep(private$algo.settings$initial.state$x0[i], nrow(data))
  }
  # now overwrite if initial guesses were provided in the data
  bool <- private$names$states %in% names(data)
  tempdata[private$names$states[bool]] <- data[private$names$states[bool]]

  # next we need to repeat these each of these state values to create
  #intermediate points determined by the user-selected ode.timestep variable
  private$algo.settings$tmb.initial.state <- vector("list",length=private$dims$states)
  for(i in seq_along(private$names$states)){
    private$algo.settings$tmb.initial.state[[i]] <- rep(tempdata[[i]], times=c(private$algo.settings$ode.timesteps,1))
  }
  names(private$algo.settings$tmb.initial.state) = private$names$states
  private$algo.settings$tmb.initial.state <- as.data.frame(private$algo.settings$tmb.initial.state)


  # return
  return(invisible(self))
}

########################################################################
# SET PARAMETERS (NEW VERSION - FOR TESTING BEFORE REPLACING ABOVE)
########################################################################
set_parameters = function(pars, self, private){

  # This function sets the parameters used by estimations, filters etc.

  lp = length(private$names$parameters)
  fp = length(private$model$fixed.pars)

  # Build base vector: per-parameter best available value ----------------
  base.pars = self$getParameters(value = "initial")
  if (!is.null(private$results$fit$par.fixed)) {
    estimated = self$getParameters(value = "estimate")
    has.estimate = !is.na(estimated)
    base.pars[has.estimate] = estimated[has.estimate]
  }

  # Apply user-supplied overrides ----------------------------------------
  if (is.null(pars)) {
    # Nothing supplied: use base vector as-is
    pars = base.pars
  } else if (!is.null(names(pars))) {
    # Named vector: replace only the named entries
    base.pars[names(pars)] <- pars
    pars = base.pars
  } else {
    # Unnamed vector: must be length lp or lp-fp
    if (!any(length(pars) == c(lp, lp - fp))) {
      stop("Incorrect number of parameters supplied (", length(pars), "). ",
           "Please supply either ", lp, " (all) or ", lp - fp,
           " (free only), i.e. with or without fixed parameters.")
    }
    if (length(pars) == lp - fp) {
      # Free parameters only: slot them in at the free-parameter positions
      par.type.free = self$getParameters()[, "type"] == "free"
      base.pars[par.type.free] = pars
      pars = base.pars
    }
    # else length == lp: use pars directly (full vector in natural order)
  }

  # set names and leave
  names(pars) = private$names$parameters
  private$algo.settings$argument.parameters = pars

  return(invisible(NULL))
}

########################################################################
# SET K STEP AHEAD (DEPENDS ON PRIVATE$DATA)
########################################################################
set_k_ahead = function(k.ahead, self, private) {

  # check if k.ahead is positive with length 1
  if (!(is.numeric(k.ahead)) | !(length(k.ahead==1)) | !(k.ahead >= 0)) {
    stop("k.ahead must be a non-negative numeric scalar")
  }
  # make sure its an integer
  k.ahead <- round(k.ahead)


  if(!is.finite(k.ahead)){
    k.ahead <- nrow(private$data) - 1
  }

  # Find last prediction index to avoid exciting boundary
  last.pred.index = nrow(private$data) - k.ahead
  if(last.pred.index < 1){
    k.ahead = nrow(private$data) - 1
    last.pred.index = 1
  }

  # set values
  private$algo.settings$k.ahead = k.ahead
  private$algo.settings$last.pred.index = last.pred.index

  # return values
  return(invisible(self))
}

