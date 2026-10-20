#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

laplace_smooth <- function(self, private){

  # Get parameters
  par.type.free <- self$getParameters()[,"type"] == "free"
  free.pars <- private$algo.settings$argument.parameters[par.type.free]

  # Integrate out random effect states (perform laplace approx)
  private$nll$fn(free.pars)

  # Calculate random effects standard deviations
  sdr <- RTMB::sdreport(private$nll,
                        ignore.parm.uncertainty=T,
                        skip.delta.method=T)

  # Store result
  private$sdr <- sdr

  # return
  return(invisible(self))
}
