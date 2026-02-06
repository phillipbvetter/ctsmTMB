
# These functions are used to check whether or not it is necessary to e.g. 
# rebuild the symbolic model, set the data again, recompile the cpp function,
# or construct the AD functions.

save_settings_for_check <- function(self, private){
  
  # Extract only necessary entries from the private fields
  entries.to.extract <- c("method", 
                          "ode.solver",
                          "loss",
                          "estimate.initial",
                          "ukf.hyperpars",
                          "ode.timestep",
                          "simulation.timestep")
  cloned.private <- self$clone(deep=TRUE)$getPrivateFields()
  private$old.data[entries.to.extract] <- mget(entries.to.extract, envir = cloned.private)
  
  return(invisible(self))
}

check_for_data_rebuild <- function(data, self, private){
  
  # Check if the data, or the requested ode/sde time-steps has changed
  # since the last call
  bool <- c(
    private$rebuild.data,
    !identical(private$old.data$entry.data, data),
    !identical(private$old.data$ode.timestep, private$ode.timestep),
    !identical(private$old.data$simulation.timestep, private$simulation.timestep)
  )
  
  private$rebuild.data <- any(bool)
  
  return(invisible(self))
}

check_for_ad_rebuild <- function(self, private){
  
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
