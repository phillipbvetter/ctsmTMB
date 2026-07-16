
# These functions are used to check whether or not it is necessary to e.g.
# rebuild the symbolic model, set the data again, recompile the cpp function,
# or construct the AD functions.

save_settings_for_ad_construct_check <- function(self, private){

  private$old.data$method              <- private$algo.settings$method
  private$old.data$ode.solver          <- private$algo.settings$ode.solver
  private$old.data$loss                <- private$algo.settings$loss
  private$old.data$estimate.initial    <- private$algo.settings$estimate.initial
  private$old.data$ukf.hyperpars       <- private$algo.settings$ukf.hyperpars

  return(invisible(self))
}

# Helper function for triggering re-computation of intermediates and the RTMB AD graph
flick_data_rebuild_switches <- function(str, self, private){

  # change flags in each case and update old data field for next time
  if(str=="data"){
    private$rebuild$data <- FALSE
    private$rebuild$ode.dt <- TRUE
    private$rebuild$sim.dt <- TRUE
  }

  if(str=="ode.dt"){
    private$rebuild$ode.dt <- FALSE
  }

  if(str=="sim.dt"){
    private$rebuild$sim.dt <- FALSE
  }

  # We must always rebuild the ad graph when these change
  private$rebuild$ad <- TRUE

  return(invisible(self))
}

# This function checks new entries against old to see if changes require rebuilding the AD graph
check_for_ad_rebuild <- function(self, private){

  # Which fields to check
  check.these.entries <- c(
    "method",
    "ode.solver",
    "loss",
    "estimate.initial",
    "ukf.hyperpars"
  )

  # perform check
  bools <- do.call(c,lapply(check.these.entries, function(s) !identical(private$old.data[[s]],private$algo.settings[[s]])))
  bool <- any(private$rebuild$ad, bools)

  # bool <- any(c(
  #   private$rebuild$ad,
  #   !identical(private$old.data$method, private$algo.settings$method),
  #   !identical(private$old.data$ode.solver, private$algo.settings$ode.solver),
  #   !identical(private$old.data$loss, private$algo.settings$loss),
  #   !identical(private$old.data$estimate.initial, private$algo.settings$estimate.initial),
  #   !identical(private$old.data$ukf.hyperpars, private$algo.settings$ukf.hyperpars)
  # ))

  # update ad logical
  private$rebuild$ad <- bool

  return(invisible(self))
}
