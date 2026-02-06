#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

construct_makeADFun = function(self, private){
  
  # TMB::openmp(n=1, autopar=TRUE, DLL=private$modelname.with.method)
  
  # Check for AD rebuild
  check_for_ad_rebuild(self, private)
  if(!private$rebuild.ad) return(invisible(self))
  save_settings_for_check(self, private)
  private$rebuild.ad <- FALSE
  
  if(!private$silent) message("Constructing objective function and derivative tables...")
  
  comptime <- system.time({
    
    switch(private$method,
           ekf = makeADFun_ekf_rtmb(self, private),
           ekf.cpp = makeADFun_ekf_tmb(self, private),
           #
           lkf = makeADFun_lkf_rtmb(self, private),
           lkf.cpp = makeADFun_lkf_tmb(self, private),
           #
           ukf = {
             if(!private$silent) message("The RTMB UKF implementation may be unstable. You can try the TMB version instead 'method=ukf.cpp'.")
             makeADFun_ukf_rtmb(self, private)
           },
           ukf.cpp = makeADFun_ukf_tmb(self, private),
           #
           laplace = makeADFun_laplace_rtmb(self, private),
           laplace.thygesen = makeADFun_laplace2_rtmb(self, private)
    )
    
  }, gcFirst = FALSE)
  
  private$timer_construct_adfun <- comptime
  
  return(invisible(self))
}

#######################################################
# OPTIMISE AD FUN
#######################################################

perform_estimation = function(self, private) {
  
  if(!private$silent) message("Minimizing the negative log-likelihood...")
  
  # Parameter Bounds
  initial.parameters <- sapply(private$parameters[names(private$free.pars)], 
                               function(par) par$initial)
  lower.parameter.bound <- sapply(private$free.pars, function(par) par$lower)
  upper.parameter.bound <- sapply(private$free.pars, function(par) par$upper)
  
  # unconstrained optimization?
  if(private$unconstrained.optim){
    lower.parameter.bound = -Inf
    upper.parameter.bound = Inf
  }
  
  comptime <- system.time(
    {
      
      # IF METHOD IS KALMAN FILTER
      if ( any(private$method == c("lkf","lkf.cpp","ekf","ekf.cpp","ukf","ukf.cpp")) ) {
        
        # use function, gradient and hessian
        if (private$use.hessian) {
          opt <- try_with_warning_recovery(stats::nlminb(start = initial.parameters,
                                                       objective = private$nll$fn,
                                                       gradient = private$nll$gr,
                                                       hessian = private$nll$he,
                                                       lower = lower.parameter.bound,
                                                       upper = upper.parameter.bound,
                                                       control=private$control.nlminb))
          # or just function and gradient
        } else {
          opt <- try_with_warning_recovery(stats::nlminb(start = initial.parameters,
                                                       objective = private$nll$fn,
                                                       gradient = private$nll$gr,
                                                       lower = lower.parameter.bound,
                                                       upper = upper.parameter.bound,
                                                       control=private$control.nlminb))
        }
      }
      
      # IF METHOD IS LAPLACE
      if ( any(private$method == c("laplace","laplace.thygesen")) ) {
        opt <- try_with_warning_recovery(stats::nlminb(start = initial.parameters,
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
  
  # 
  # For TMB method: run sdreport
  if (any(private$method== c("laplace","laplace.thygesen"))) {
    if(!private$silent) message("Calculating standard deviations...")
    # NOTE: The state covariances can be retrived by inverting sdr$jointPrecision
    # but this takes very long time. Should it be an option?
    private$sdr <- TMB::sdreport(private$nll, getJointPrecision=FALSE)
  }
  
  # return
  return(invisible(self))
}

#######################################################
# MAIN RETURN FIT FUNCTION THAT CALL OTHERS
#######################################################

create_estimation_return_fit = function(self, private, report, laplace.residuals){
  
  if(!private$silent) message("Returning results...")
  
  # Initialization and Clearing -----------------------------------
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  # clear fit
  private$fit = NULL
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  # Fit Info -----------------------------------
  compute_mle_gradient_and_hessian(self, private)
  
  # Parameters and Uncertainties -----------------------------------
  compute_mle_parameters_and_std_errors(self, private)
  
  if(report){
    
    if(private$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")) {
      
      # Call filter with silenced settings
      silent.setting <- private$silent
      on.exit(private$silent<-silent.setting, add=TRUE)
      self$filter(private$data, silent=TRUE)
      private$fit = c(private$fit, private$filtration)
      
    } 
    
    if(private$method %in% c("laplace","laplace.thygesen")) {
      
      laplace_report(self, private, laplace.residuals)
      
    }
    
  }
  
  
  # clone private into fit -----------------------------------
  # no need for deep copy??
  private$fit$private <- self$clone()$.__enclos_env__$private
  
  # set s3 class -----------------------------------
  class(private$fit) = "ctsmTMB.fit"
  
  # return -----------------------------------
  return(invisible(self))
}
