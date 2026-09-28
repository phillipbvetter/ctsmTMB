##########################################################################################
# AD REBUILD SETTINGS
##########################################################################################

# These are the constant fields that we must watch for changes
# If they change we must recompile the AD graph via RTMB::MakeADFun.
ad.recompile.fields <- c(
  "method",
  "ode.solver",
  "loss",
  "estimate.initial",
  "ukf.hyperpars",
  "first.order.input.hold",
  NULL
)

# This function saves the settings that when changed requires rebuildind the AD graph via RTMB::MakeADFun.
save_ad_rebuild_fields <- function(self, private) {
  private$old.data[ad.recompile.fields] <- private$algo.settings[ad.recompile.fields]
  private$rebuild$ad <- FALSE
  return(invisible(self))
}

# This function checks the newest parsed data / settings against the ones used in the previous call.
# If they have changed we must rebuild the AD likelihood function via RTMB::MakeADFun.
check_for_ad_rebuild <- function(self, private) {
  bool <- unlist(lapply(ad.recompile.fields, function(s) !identical(private$old.data[[s]], private$algo.settings[[s]])))
  private$rebuild$ad <- any(private$rebuild$ad, bool)
  return(invisible(self))
}



##########################################################################################
# DATA REBUILD SETTINGS (ALSO CAUSES NEED FOR AD REBUILD)
##########################################################################################

# Helper function for triggering re-computation of data entries
flick_data_rebuild_switches <- function(str, self, private){

  # If the data is changed then we must also rebuild timesteps
  if(str=="data"){
    private$rebuild$data <- FALSE
    private$rebuild$ode.timestep <- TRUE
    private$rebuild$sim.timestep <- TRUE
  }

  if(str=="ode"){
    private$rebuild$ode.timestep <- FALSE
  }

  if(str=="sim"){
    private$rebuild$sim.timestep <- FALSE
  }

  # We must always rebuild the ad graph when these change
  private$rebuild$ad <- TRUE

  return(invisible(self))
}
