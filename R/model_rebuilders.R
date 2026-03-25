
# These functions are used to check whether or not it is necessary to e.g.
# rebuild the symbolic model, set the data again, recompile the cpp function,
# or construct the AD functions.

save_settings_for_check <- function(self, private){

  private$old.data$method              <- private$algo.settings$method
  private$old.data$ode.solver          <- private$algo.settings$ode.solver
  private$old.data$loss                <- private$algo.settings$loss
  private$old.data$estimate.initial    <- private$algo.settings$estimate.initial
  private$old.data$ukf.hyperpars       <- private$algo.settings$ukf.hyperpars
  private$old.data$ode.timestep        <- private$ode.timestep
  private$old.data$simulation.timestep <- private$simulation.timestep

  return(invisible(self))
}

check_for_data_rebuild <- function(data, self, private){

  # Check if the data, or the requested ode/sde time-steps has changed
  # since the last call
  bool <- c(
    private$rebuild$data,
    !identical(private$old.data$entry.data, data),
    !identical(private$old.data$ode.timestep, private$ode.timestep),
    !identical(private$old.data$simulation.timestep, private$simulation.timestep)
  )
  private$rebuild$data <- any(bool)

  return(invisible(self))
}

check_for_ad_rebuild <- function(self, private){

  # We perform checks against the old data on the entries that would
  # require a new call to RTMB::MakeADFun (i.e. entries that affect the
  # calculations in the likelihood function)
  bool <- c(
    private$rebuild$ad,
    !identical(private$old.data$method, private$algo.settings$method),
    !identical(private$old.data$ode.solver, private$algo.settings$ode.solver),
    !identical(private$old.data$loss, private$algo.settings$loss),
    !identical(private$old.data$estimate.initial, private$algo.settings$estimate.initial),
    !identical(private$old.data$ukf.hyperpars, private$algo.settings$ukf.hyperpars)
  )

  private$rebuild$ad <- any(bool)

  return(invisible(self))
}
