
# These functions are used to check whether or not it is necessary to e.g.
# rebuild the symbolic model, set the data again, recompile the cpp function,
# or construct the AD functions.

save_settings_for_ad_construct_check <- function(self, private){

  private$old.data$method              <- private$algo.settings$method
  private$old.data$ode.solver          <- private$algo.settings$ode.solver
  private$old.data$first.order.input.hold <- private$algo.settings$first.order.input.hold
  private$old.data$loss                <- private$algo.settings$loss
  private$old.data$estimate.initial    <- private$algo.settings$estimate.initial
  private$old.data$ukf.hyperpars       <- private$algo.settings$ukf.hyperpars
  private$old.data$ode.timestep        <- private$algo.settings$ode.timestep
  private$old.data$simulation.timestep <- private$algo.settings$simulation.timestep

  return(invisible(self))
}

check_for_data_rebuild <- function(data, self, private){

  # Check if the data, or the requested ode/sde time-steps has changed since last call
  bool <- c(
    private$rebuild$data,
    !identical(private$old.data$entry.data, data),
    !identical(private$old.data$ode.timestep, private$algo.settings$ode.timestep),
    !identical(private$old.data$simulation.timestep, private$algo.settings$simulation.timestep)
  )
  private$rebuild$data <- any(bool)

  return(invisible(self))
}

check_for_ad_rebuild <- function(self, private){

  fields <- c("method", "ode.solver", "loss", "estimate.initial", "ukf.hyperpars", "first.order.input.hold")
  bool <- unlist(lapply(fields, function(s) !identical(private$old.data[[s]], private$algo.settings[[s]])))
  private$rebuild$ad <- any(private$rebuild$ad, bool)

  return(invisible(self))
}
