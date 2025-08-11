
# These functions are used to check whether or not it is necessary to e.g. 
# rebuild the symbolic model, set the data again, recompile the cpp function,
# or construct the AD functions.

save_settings_for_check <- function(self, private){
  
  cloned.self <- self$clone(deep=TRUE)
  cloned.private <- cloned.self$.private()
  
  # settings
  private$old.data$method = cloned.private$method
  private$old.data$ode.solver = cloned.private$ode.solver
  private$old.data$loss = cloned.private$loss
  private$old.data$estimate.initial = cloned.private$estimate.initial
  
  # ukf hyperpars
  private$old.data$ukf.hyperpars = cloned.private$ukf.hyperpars
  
  # Used for predictions
  private$old.data$ode.timestep = cloned.private$ode.timestep
  private$old.data$simulation.timestep = cloned.private$simulation.timestep
  
  return(invisible(self))
}

check_for_data_rebuild <- function(data, self, private){
  
  # We check if the data needs to reset
  bool <- c(
    private$rebuild.data,
    !identical(private$old.data$entry.data, data),
    !identical(private$old.data$ode.timestep, private$ode.timestep),
    !identical(private$old.data$simulation.timestep, private$simulation.timestep)
  )
  
  private$rebuild.data <- any(bool)
  
  return(invisible(self))
}

check_for_ADfun_rebuild <- function(self, private){
  
  # We perform checks against the old data on the entries that would
  # require a new call to RTMB::MakeADFun (i.e. entries that affect the
  # calculations in the likelihood function)
  bool <- c(
    private$rebuild.ad,
    !identical(private$old.data$method, private$method),
    !identical(private$old.data$ode.solver, private$ode.solver),
    !identical(private$old.data$loss, private$loss),
    !identical(private$old.data$estimate.initial, private$estimate.initial),
    !identical(private$old.data$ukf.hyperpars, private$ukf.hyperpars)
  )
  
  private$rebuild.ad <- any(bool)
  
  return(invisible(self))
}
