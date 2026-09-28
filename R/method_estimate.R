create_ad_likelihood_fun <- function(self, private){

  # TMB::openmp(n=1, autopar=TRUE, DLL=private$modelname.with.method)

  # Check for rebuild - exit or continue
  # check_or_save_for_ad_rebuild("check", self, private)
  check_for_ad_rebuild(self, private)
  if (!private$rebuild$ad) return(invisible(self))

  # Rebuild needed - save settings for next time and compile nll fun
  # check_or_save_for_ad_rebuild("save", self, private)
  save_ad_rebuild_fields(self, private)

  if (!private$algo.settings$silent)
    message("Compiling automatic-differentiation graph...")

  if (private$algo.settings$method == "ukf" && !private$algo.settings$silent)
    message("This UKF implementation may be unstable - alternativly try 'method = ukf.cpp'.")

  # We call the precompiled list of
  comptime <- system.time({
    .List_of_MakeADFuns[[private$algo.settings$method]](self, private)
  }, gcFirst = FALSE)

  private$timers$construct_adfun <- comptime

  return(invisible(self))
}

# List of callable MakeADFun functions for the different methods
.List_of_MakeADFuns <- list(
  ekf = MakeADFun_EKF,
  ekf.cpp = MakeADFun_EKF_TMB,
  lkf = MakeADFun_LKF,
  lkf.cpp = MakeADFun_LKF_TMB,
  ukf = MakeADFun_UKF,
  ukf.cpp = MakeADFun_UKF_TMB,
  laplace = MakeADFun_Laplace,
  laplace.thygesen = MakeADFun_Laplace_thygesen
)

#######################################################
# OPTIMISE AD FUN
#######################################################

perform_estimation = function(self, private) {

  if(private$dims$free.pars==0)
    stop("Estimation aborted - there are no (fixed effects) parameters in the model.")

  # if(!private$algo.settings$silent)
  #   message("Minimizing the negative log-likelihood...")
  if(!private$algo.settings$silent)
    message("Estimating...")

  # Parameter Bounds
  initial.parameters <- sapply(private$model$parameters[names(private$model$free.pars)], function(par) par$initial)
  lower.parameter.bound <- sapply(private$model$free.pars, function(par) par$lower)
  upper.parameter.bound <- sapply(private$model$free.pars, function(par) par$upper)

  # unconstrained optimization?
  if(private$optim.settings$unconstrained.optim){
    lower.parameter.bound = -Inf
    upper.parameter.bound = Inf
  }

  kalman.methods <- c("lkf","lkf.cpp","ekf","ekf.cpp","ukf","ukf.cpp")
  laplace.methods <- c("laplace","laplace.thygesen")

  comptime <- system.time({

    # IF METHOD IS KALMAN FILTER
    if (private$algo.settings$method %in% kalman.methods) {

      # use function, gradient and hessian
      if (private$optim.settings$use.hessian) {
        opt <- try_with_warning_recovery(stats::nlminb(start = initial.parameters,
                                                       objective = private$nll$fn,
                                                       gradient = private$nll$gr,
                                                       hessian = private$nll$he,
                                                       lower = lower.parameter.bound,
                                                       upper = upper.parameter.bound,
                                                       control=private$optim.settings$control.nlminb))
      } else {
        # or just function and gradient
        opt <- try_with_warning_recovery(stats::nlminb(start = initial.parameters,
                                                       objective = private$nll$fn,
                                                       gradient = private$nll$gr,
                                                       lower = lower.parameter.bound,
                                                       upper = upper.parameter.bound,
                                                       control=private$optim.settings$control.nlminb))
      }

    } else if (private$algo.settings$method %in% laplace.methods) {
      opt <- try_with_warning_recovery(stats::nlminb(start = initial.parameters,
                                                     objective = private$nll$fn,
                                                     gradient = private$nll$gr,
                                                     lower = lower.parameter.bound,
                                                     upper = upper.parameter.bound,
                                                     control=private$optim.settings$control.nlminb))
    } else {

      stop("No estimation method was recognized.")

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
              4. Consider reducing the 'ode.timestep' (this also reduces the SDE timestep in the laplace methods).
              5. The Kalman filters may benefit from optimization with the hessian i.e. 'use.hessian'
              6. Try using other optimizers by extracting the function handlers via the 'likelihood' method.
              7. Try changing the optimizer tolerances with the 'control' argument. See ?stats::nlminb for details.")
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
  if(!private$algo.settings$silent){
    # if(outer_mgc > 1){
    # message("BEWARE: THE MAXIMUM GRADIENT COMPONENT APPEARS TO BE LARGE ( > 1 ) - THE FOUND OPTIMUM MIGHT BE INVALID.")
    # }

    message(
      sprintf("Estimation results:\n\n"),
      sprintf("  %-30s %s seconds\n", "Elapsed time:", comp.time),
      sprintf("  %-30s %s\n", "Objective value:", format(opt$objective, scientific = TRUE)),
      sprintf("  %-30s %s\n", "Max gradient component:", format(outer_mgc, digits = 2, scientific = TRUE)),
      sprintf("  %-30s %s\n", "Convergence message:", opt$message),
      sprintf("  %-30s %s\n", "Iterations:", opt$iterations),
      sprintf("  %-30s Fun: %s  Grad: %s\n", "Evaluations:", opt$evaluations["function"], opt$evaluations[["gradient"]]),
      "\n  See stats::nlminb for available tolerance/control arguments.\n"
    )
  }

  #
  # For TMB method: run sdreport
  if (any(private$algo.settings$method== c("laplace","laplace.thygesen"))) {
    if(!private$algo.settings$silent) message("Calculating standard deviations...")
    # NOTE: The state covariances can be retrived by inverting sdr$jointPrecision
    # but this takes very long time. Should it be an option?
    private$sdr <- TMB::sdreport(private$nll, getJointPrecision=FALSE)
  }

  # return
  return(invisible(self))
}

create_estimation_return_fit <- function(self, private, report, laplace.residuals){

  if(!private$algo.settings$silent) message("Returning results...")

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
  # NOTE: populates private$results$fit$par.fixed, used by set_parameters below
  compute_mle_parameters_and_std_errors(self, private)

  if(report){

    if(private$algo.settings$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")) {

      # NULL finds the MLE parameters from estimate
      set_parameters(NULL, self, private)

      # Temporarily silence filtering printouts...
      silent.setting <- private$algo.settings$silent
      on.exit(private$algo.settings$silent <- silent.setting, add=TRUE)
      private$algo.settings$silent <- TRUE

      # We run the filter function directly bypassing build, set flags etc
      filter_predict_simulate_smooth(self, private, proc="filter")
      create_filter_results(self, private, laplace.residuals)

      # Merge filter results into fit
      private$results$fit = c(private$results$fit, private$results$filtration)

    }

    if (private$algo.settings$method %in% c("laplace","laplace.thygesen"))
      laplace_report(self, private, laplace.residuals)

  }

  # snapshot only the fields consumed by S3 methods
  private$results$fit$private <- private$make_private_snapshot()

  # set s3 class -----------------------------------
  class(private$results$fit) = "ctsmTMB.fit"

  # return -----------------------------------
  return(invisible(self))
}
