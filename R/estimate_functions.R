#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

construct_makeADFun = function(self, private){
  
  # TMB::openmp(n=1, autopar=TRUE, DLL=private$modelname.with.method)
  
  # Check for AD rebuild
  check_for_ADfun_rebuild(self, private)
  if(!private$rebuild.ad) return(invisible(self))
  save_settings_for_check(self, private)
  private$rebuild.ad <- FALSE
  
  if(!private$silent) message("Constructing objective function and derivative tables...")
  
  comptime <- system.time({
    
    if(private$method == "lkf"){
      makeADFun_lkf_rtmb(self, private)
    }
    
    if(private$method == "lkf_cpp"){
      makeADFun_lkf_tmb(self, private)
    }
    
    if(private$method == "ekf"){
      makeADFun_ekf_rtmb(self, private)
    }
    
    if(private$method == "ekf_cpp"){
      makeADFun_ekf_tmb(self, private)
    }
    
    if(private$method == "ukf"){
      makeADFun_ukf_rtmb(self, private)
    }
    
    if(private$method == "ukf_cpp"){
      makeADFun_ukf_tmb(self, private)
    }
    
    if(private$method=="laplace"){
      makeADFun_laplace_rtmb(self, private)
    }
    
    if(private$method=="laplace2"){
      makeADFun_laplace2_rtmb(self, private)
    }
    
    # experimental
    # if(private$method == "ekf_rcpp"){
    #     rcpp_ekf_estimate(self, private)
    # }
    
  }, gcFirst = FALSE)
  
  private$timer_construct_adfun <- comptime
  
  return(invisible(self))
}

#######################################################
# MAIN RETURN FIT FUNCTION THAT CALL OTHERS
#######################################################

create_return_fit = function(self, private, laplace.residuals){
  
  if(!private$silent) message("Returning results...")
  
  if(private$method == "lkf"){
    calculate_fit_statistics_lkf(self, private)
  }
  
  if(private$method == "lkf_cpp"){
    calculate_fit_statistics_lkf(self, private)
  }
  
  if(private$method == "ekf"){
    calculate_fit_statistics_ekf(self, private)
  }
  
  if(private$method == "ekf_cpp"){
    calculate_fit_statistics_ekf(self, private)
  }
  
  if(private$method == "ukf"){
    calculate_fit_statistics_ukf(self, private)
  }
  
  if(private$method == "ukf_cpp"){
    calculate_fit_statistics_ukf(self, private)
  }
  
  if(private$method=="laplace"){
    calculate_fit_statistics_laplace(self, private, laplace.residuals)
  }
  
  if(private$method=="laplace2"){
    calculate_fit_statistics_laplace2(self, private, laplace.residuals)
  }
  
  return(invisible(self))
}

#######################################################
# OPTIMISE AD FUN
#######################################################

perform_estimation = function(self, private) {
  
  if(!private$silent) message("Minimizing the negative log-likelihood...")
  
  # Parameter Bounds
  # lower.parameter.bound = unlist(lapply(private$free.pars, function(par) par$lower))
  # upper.parameter.bound = unlist(lapply(private$free.pars, function(par) par$upper))
  lower.parameter.bound = sapply(private$free.pars, function(par) par$lower)
  upper.parameter.bound = sapply(private$free.pars, function(par) par$upper)
  if(private$unconstrained.optim){
    lower.parameter.bound = -Inf
    upper.parameter.bound = Inf
  }
  
  comptime <- system.time(
    {
      
      # IF METHOD IS KALMAN FILTER
      if ( any(private$method == c("lkf","lkf_cpp","ekf","ekf_cpp","ukf","ukf_cpp")) ) {
        
        # use function, gradient and hessian
        if (private$use.hessian) {
          opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                       objective = private$nll$fn,
                                                       gradient = private$nll$gr,
                                                       hessian = private$nll$he,
                                                       lower = lower.parameter.bound,
                                                       upper = upper.parameter.bound,
                                                       control=private$control.nlminb))
          # or just function and gradient
        } else {
          opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                       objective = private$nll$fn,
                                                       gradient = private$nll$gr,
                                                       lower = lower.parameter.bound,
                                                       upper = upper.parameter.bound,
                                                       control=private$control.nlminb))
        }
      }
      
      # IF METHOD IS LAPLACE
      if ( any(private$method == c("laplace","laplace2")) ) {
        opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                     objective = private$nll$fn,
                                                     gradient = private$nll$gr,
                                                     lower = lower.parameter.bound,
                                                     upper = upper.parameter.bound,
                                                     control=private$control.nlminb))
      }
      
    }, gcFirst = FALSE)
  
  # add timer to estimation
  private$timer_estimation <- comptime
  
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
    
    # exit if optimization failed
    private$opt = NULL
    return(invisible(self))
  }
  
  # store optimization object
  names(opt$par) <- names(private$free.pars)
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
    
    # NOTE: The state covariances can be retrived by inverting sdr$jointPrecision
    # but this takes very long time. Should it be an option?
    comptime <- system.time(
      {
        
      private$sdr <- TMB::sdreport(private$nll, getJointPrecision=T)
      
      }, gcFirst = FALSE)
  }
  
  # return
  return(invisible(self))
}

