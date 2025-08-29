report_return_fit <- function(self, private, laplace.residuals=NULL) {
  
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
  
  # Get Report -----------------------------------
  rep <- get_state_report(self, private)
  
  # States -----------------------------------
  compute_return_states(rep, self, private)
  
  # Residuals -----------------------------------
  report_residuals(rep, laplace.residuals=laplace.residuals, self, private)
  
  # Observations -----------------------------------
  report_observations(self, private)
  
  # clone private into fit -----------------------------------
  private$fit$private <- self$clone()$.__enclos_env__$private
  
  # set s3 class -----------------------------------
  class(private$fit) = "ctsmTMB.fit"
  
  # return -----------------------------------
  return(invisible(self))
  
}
