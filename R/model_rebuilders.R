
# This function either checks the newest parsed data / settings against the ones used last call.
# If they have changed we must rebuild the AD likelihood function via RTMB::MakeADFun.
check_or_save_for_ad_rebuild <- function(type, self, private){

  # These are the fields that require an AD-fun rebuild
  fields <- c(
    "method",
    "ode.solver",
    "loss",
    "estimate.initial",
    "ukf.hyperpars",
    "first.order.input.hold"
    )

  # CHECK MODE:
  if (type=="check") {
    bool <- unlist(lapply(fields, function(s) !identical(private$old.data[[s]], private$algo.settings[[s]])))
    private$rebuild$ad <- any(private$rebuild$ad, bool)
  }

  # SAVE MODE:
  if (type == "save") {
    private$old.data[fields] <- private$algo.settings[fields]
    private$rebuild$ad <- FALSE
  }

  return(invisible(self))

}

# Helper function for triggering re-computation of intermediates and the RTMB AD graph
flick_data_rebuild_switches <- function(str, self, private){

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
