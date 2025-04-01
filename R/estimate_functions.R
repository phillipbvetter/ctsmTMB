#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

construct_makeADFun = function(self, private){
  
  # TMB::openmp(max=TRUE, autopar=TRUE, DLL=private$modelname.with.method)
  
  if(private$method == "lkf"){
    comptime <- system.time(
      makeADFun_lkf_rtmb(self, private)
    )
  }
  
  if(private$method == "ekf"){
    makeADFun_ekf_rtmb(self, private)
  }
  
  if(private$method=="laplace"){
    comptime <- system.time(
      makeADFun_laplace_rtmb(self, private)
    )
  }
  
  if(private$method=="laplace2"){
    comptime <- system.time(
      makeADFun_laplace2_rtmb(self, private)
    )
  }
  
  # experimental
  # if(private$method == "ekf_rcpp"){
  #   comptime <- system.time(
  #     rcpp_ekf_estimate(self, private)
  #     )
  # }
  
  # comptime = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
  # if(!private$silent) message("...took: ", comptime, " seconds.")
  
  return(invisible(self))
}

#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

create_fit = function(self, private, laplace.residuals){
  
  # TMB::openmp(max=TRUE, autopar=TRUE, DLL=private$modelname.with.method)
  if(private$method == "lkf"){
    calculate_fit_statistics_lkf(self, private)
  }
  
  if(private$method == "ekf"){
    calculate_fit_statistics_ekf(self, private)
  }
  
  # if(private$method == "ekf_rcpp"){
  # calculate_fit_statistics_ekf(self, private)
  # }
  
  if(private$method=="laplace"){
    calculcate_fit_statistics_laplace(self, private, laplace.residuals)
  }
  
  if(private$method=="laplace2"){
    calculcate_fit_statistics_laplace2(self, private, laplace.residuals)
  }
  
  return(invisible(self))
}

#######################################################
# OPTIMISE AD FUN
#######################################################

perform_estimation = function(self, private) {
  
  # Parameter Bounds
  lower.parameter.bound = unlist(lapply(private$free.pars, function(par) par$lower))
  upper.parameter.bound = unlist(lapply(private$free.pars, function(par) par$upper))
  if(private$unconstrained.optim){
    lower.parameter.bound = -Inf
    upper.parameter.bound = Inf
  }
  
  # IF METHOD IS KALMAN FILTER
  if (any(private$method==c("lkf","ekf","ukf_cpp","ekf_cpp"))) {
    
    if(private$method=="ekf_cpp"){
      # TMB::openmp(max=TRUE, DLL=private$modelname.with.method, autopar=TRUE)
      # TMB::config(trace.atomic=FALSE, DLL=private$modelname.with.method)
    }
    
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
  
  # IF METHOD IS LAPLACE
  if (any(private$method == c("laplace","laplace2"))) {
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
              1. Explore other parameter initial values - watch out of boundaries.
              2. Consider parameter transformations that ensure appropriate domains.
              3. Consider relative parameter values - they should ideally be similar.
              4. Consider reducing the 'ode.timestep' (also reduces the SDE timestep for the laplace method).
              5. The Kalman filters may benefit from optimization with the hessian i.e. 'use.hessian'
              6. Try other optimizations using the function handlers from the 'likelihood' method.
              7. Change the optimization tolerances for 'nlminb' with the 'control' argument.")
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
    # if(outer_mgc > 1){
    # message("BEWARE: THE MAXIMUM GRADIENT COMPONENT APPEARS TO BE LARGE ( > 1 ) - THE FOUND OPTIMUM MIGHT BE INVALID.")
    # }
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
  if (any(private$method== c("laplace","laplace2"))) {
    if(!private$silent) message("Calculating standard deviations...")
    comptime = system.time(
      private$sdr <- TMB::sdreport(private$nll, getJointPrecision=T)
      
      # NOTE: The state covariances can be retrived by inverting sdr$jointPrecision
      # but this takes very long time. Should it be an option?
    )
    # comptime = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
    # if(!private$silent) message("...took: ", comptime, " seconds.")
  }
  
  # return
  return(invisible(self))
}

