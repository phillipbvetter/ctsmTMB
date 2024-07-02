#######################################################
# OPTIMISE AD FUN
#######################################################

optimize_negative_loglikelihood = function(self, private) {
  
  lower.parameter.bound = unlist(lapply(private$free.pars, function(par) par$lower))
  upper.parameter.bound = unlist(lapply(private$free.pars, function(par) par$upper))
  # If unconstrained optimization was requested remove these bounds
  if(private$unconstrained.optim){
    lower.parameter.bound = -Inf
    upper.parameter.bound = Inf
  }
  
  # IF METHOD IS KALMAN FILTER
  if (any(private$method==c("ekf","ukf","ekf_rtmb"))) {
    
    # use function, gradient and hessian
    if (private$use.hessian) {
      comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                           objective = private$nll$fn,
                                                                           gradient = private$nll$gr,
                                                                           hessian = private$nll$he,
                                                                           lower = lower.parameter.bound,
                                                                           upper = upper.parameter.bound,
                                                                           control=private$control.nlminb))
      )
      # or just function and gradient
    } else {
      comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                           objective = private$nll$fn,
                                                                           gradient = private$nll$gr,
                                                                           lower = lower.parameter.bound,
                                                                           upper = upper.parameter.bound,
                                                                           control=private$control.nlminb))
      )
    }
    
  }
  
  # IF METHOD IS TMB
  if (private$method =="laplace") {
    comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                         objective = private$nll$fn,
                                                                         gradient = private$nll$gr,
                                                                         lower = lower.parameter.bound,
                                                                         upper = upper.parameter.bound,
                                                                         control=private$control.nlminb))
    )
  }
  
  # DID THE OPTIMIZATION FAIL?
  if (inherits(opt,"try-error")) {
    
    message("The optimisation failed due to the following error: \n\n\t",opt)
    
    if(stringr::str_detect(opt,"NA/NaN")){
      
      message("You should consider the following to circumvent the error:
              1. Reduce the ODE step-size via argument 'ode.timestep'.
              2. Run the optimization with / without the hessian via argument 'use.hessian'.
              3. Change / explore parameter initial values.
              4. Extract the function handlers with the 'construct_nll' method and investigate outputs, or try other optimizers.
              5. Change the optimization tolerances for 'nlminb' with the 'control' argument.")
      
    }
    
    private$opt = NULL
    
    return(invisible(self))
  }
  
  # store optimization object
  private$opt = opt
  
  # extract maxmimum gradient component, and format computation time to 5 digits
  outer_mgc = max(abs(private$nll$gr(opt$par)))
  comp.time = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
  
  # print convergence and timing result
  if(!private$silent){
    if(outer_mgc > 1){
      message("BEWARE: THE MAXIMUM GRADIENT COMPONENT APPEARS TO BE LARGE ( > 1 ) - THE FOUND OPTIMUM MIGHT BE INVALID.")
    }
    message("\t Optimization finished!:
            Elapsed time: ", comp.time, " seconds.
            The objective value is: ",format(opt$objective,scientific=T),"
            The maximum gradient component is: ",format(outer_mgc,digits=2,scientific=T),"
            The convergence message is: ", opt$message,"
            Iterations: ",opt$iterations,"
            Evaluations: Fun: ",opt$evaluations["function"]," Grad: ",opt$evaluations[["gradient"]],"
            See stats::nlminb for available tolerance/control arguments."
    )
  }
  
  # For TMB method: run sdreport
  if (private$method=="laplace") {
    if(!private$silent) message("Calculating random effects standard deviation...")
    comptime = system.time(private$sdr <- RTMB::sdreport(private$nll))
    comptime = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
    if(!private$silent) message("...took: ", comptime, " seconds.")
  }
  
  # return
  return(invisible(self))
}
