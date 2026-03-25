#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

construct_makeADFun = function(self, private){

  # TMB::openmp(n=1, autopar=TRUE, DLL=private$modelname.with.method)

  # Check for AD rebuild
  check_for_ad_rebuild(self, private)
  if(!private$rebuild$ad) return(invisible(self))
  save_settings_for_ad_construct_check(self, private)
  private$rebuild$ad <- FALSE

  if(!private$silent) message("Constructing objective function and derivative tables...")

  comptime <- system.time({

    switch(private$algo.settings$method,
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

  private$timers$construct_adfun <- comptime

  return(invisible(self))
}

#######################################################
# OPTIMISE AD FUN
#######################################################

perform_estimation = function(self, private) {

  if(!private$silent) message("Minimizing the negative log-likelihood...")

  # Parameter Bounds
  initial.parameters <- sapply(private$model$parameters[names(private$model$free.pars)],
                               function(par) par$initial)
  lower.parameter.bound <- sapply(private$model$free.pars, function(par) par$lower)
  upper.parameter.bound <- sapply(private$model$free.pars, function(par) par$upper)

  # unconstrained optimization?
  if(private$optim.settings$unconstrained.optim){
    lower.parameter.bound = -Inf
    upper.parameter.bound = Inf
  }

  comptime <- system.time(
    {

      # IF METHOD IS KALMAN FILTER
      if ( any(private$algo.settings$method == c("lkf","lkf.cpp","ekf","ekf.cpp","ukf","ukf.cpp")) ) {

        # use function, gradient and hessian
        if (private$optim.settings$use.hessian) {
          opt <- try_with_warning_recovery(stats::nlminb(start = initial.parameters,
                                                       objective = private$nll$fn,
                                                       gradient = private$nll$gr,
                                                       hessian = private$nll$he,
                                                       lower = lower.parameter.bound,
                                                       upper = upper.parameter.bound,
                                                       control=private$optim.settings$control.nlminb))
          # or just function and gradient
        } else {
          opt <- try_with_warning_recovery(stats::nlminb(start = initial.parameters,
                                                       objective = private$nll$fn,
                                                       gradient = private$nll$gr,
                                                       lower = lower.parameter.bound,
                                                       upper = upper.parameter.bound,
                                                       control=private$optim.settings$control.nlminb))
        }
      }

      # IF METHOD IS LAPLACE
      if ( any(private$algo.settings$method == c("laplace","laplace.thygesen")) ) {
        opt <- try_with_warning_recovery(stats::nlminb(start = initial.parameters,
                                                     objective = private$nll$fn,
                                                     gradient = private$nll$gr,
                                                     lower = lower.parameter.bound,
                                                     upper = upper.parameter.bound,
                                                     control=private$optim.settings$control.nlminb))
      }

    }, gcFirst = FALSE)

  # add timer to estimation
  private$timers$estimation <- comptime

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
    private$results$opt = NULL
    return(invisible(self))
  }

  # store optimization object
  names(opt$par) <- names(private$model$free.pars)
  private$results$opt = opt

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
  if (any(private$algo.settings$method== c("laplace","laplace.thygesen"))) {
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
  if (is.null(private$results$opt)) {
    return(NULL)
  }

  # clear fit
  private$results$fit = NULL

  # get convergence
  private$results$fit$convergence = private$results$opt$convergence

  # Fit Info -----------------------------------
  compute_mle_gradient_and_hessian(self, private)

  # Parameters and Uncertainties -----------------------------------
  compute_mle_parameters_and_std_errors(self, private)

  if(report){

    if(private$algo.settings$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")) {

      # Call filter with silenced settings
      silent.setting <- private$silent
      on.exit(private$silent <- silent.setting, add=TRUE)

      # perform filtering
      self$filter(data=private$data,
                  pars = NULL,
                  method=private$algo.settings$method,
                  ode.solver=private$algo.settings$ode.solver,
                  ode.timestep=private$algo.settings$ode.timestep,
                  loss=private$algo.settings$loss$loss,
                  loss_c=private$algo.settings$loss$loss_c,
                  ukf.hyperpars=private$algo.settings$ukf.hyperpars,
                  initial.state=private$initial.state,
                  laplace.residuals=laplace.residuals,
                  estimate.initial.state=private$algo.settings$estimate.initial,
                  use.cpp = TRUE,
                  silent = TRUE)

      # add filtered results to fit
      private$results$fit = c(private$results$fit, private$results$filtration)

    }

    if(private$algo.settings$method %in% c("laplace","laplace.thygesen")) {

      laplace_report(self, private, laplace.residuals)

    }

  }


  # snapshot only the fields consumed by S3 methods
  private$results$fit$private <- private$make_private_snapshot()

  # set s3 class -----------------------------------
  class(private$results$fit) = "ctsmTMB.fit"

  # return -----------------------------------
  return(invisible(self))
}
