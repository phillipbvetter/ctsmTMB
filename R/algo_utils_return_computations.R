#######################################################
#######################################################
# likelihood value, gradient and hessian
#######################################################
#######################################################
compute_mle_gradient_and_hessian <- function(self, private){
  
  # MLE
  private$fit$nll = private$opt$objective
  
  # MLE gradient
  private$fit$nll.gradient = try_with_warning_recovery(
    {
      nll.grad = as.vector(private$nll$gr(private$opt$par))
      names(nll.grad) = names(private$free.pars)
      nll.grad
    }
  )
  if (inherits(private$fit$nll.gradient,"try-error")) {
    private$fit$nll.gradient = NULL
  }
  
  # MLE Hessian
  if(private$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")){
    private$fit$nll.hessian = try_with_warning_recovery(
      {
        nll.hess = private$nll$he(private$opt$par)
        rownames(nll.hess) = names(private$free.pars)
        colnames(nll.hess) = names(private$free.pars)
        nll.hess
      }
    )
    if (inherits(private$fit$nll.hessian, "try-error")) {
      private$fit$nll.hessian = NULL
    }
  }
  
  return(invisible(self))
}

# function that writes parameter estimates and computes std. errors
# from either inverting the hessian (kalman filters) or from 
# using TMB::sdreport (laplace)
compute_mle_parameters_and_std_errors <- function(self, private){
  
  private$fit$par.fixed = private$opt$par
  
  # allocate 
  n.fixed.pars <- private$number.of.fixed.pars
  private$fit$sd.fixed = rep(NA, n.fixed.pars)
  private$fit$cov.fixed = array(NA, dim=rep(n.fixed.pars,2))
  private$fit$tvalue = rep(NA, n.fixed.pars)
  private$fit$Pr.tvalue = rep(NA, n.fixed.pars)
  
  ############ LAPLACE ############
  if(private$method %in% c("laplace","laplace.thygesen")){
    private$fit$cov.fixed <- private$sdr$cov.fixed
    private$fit$sd.fixed <- sqrt(diag(private$fit$cov.fixed))
  }
  
  ############ KALMAN ############
  if(private$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")){
    calculate_covariance_from_hessian(self, private)
  }
  
  # t-values and Pr( t > t_test ) -----------------------------------
  private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
  # The degrees of fredom are number of (non NA) observations minus number of (free) parameters
  freedom.degrees <- sum(!is.na(private$data[private$obs.names])) - private$number.of.free.pars
  private$fit$Pr.tvalue <- rep(NA, length(private$fit$par.fixed))
  if(freedom.degrees > 0.1){
    private$fit$Pr.tvalue = 2 * pt(q=abs(private$fit$tvalue), df=freedom.degrees, lower.tail=FALSE)
  }
  
  # return
  return(invisible(self))
}

#######################################################
#######################################################
# helper for hessian 
#######################################################
#######################################################

calculate_covariance_from_hessian <- function(self, private){
  
  # invert hessian to get covariances
  if(!is.null(private$fit$nll.hessian)){
    
    out <- iterative_hessian_inversion(self, private)
    covariance <- out$covariance
    std.dev <- out$std.dev
    
  } else {
    
    # if the hessian could not be computed, then we may try sdreport instead
    # (sometimes more stable???)
    private$sdr <- TMB::sdreport(private$nll, getJointPrecision=FALSE)
    private$fit$cov.fixed <- private$sdr$cov.fixed
    private$fit$sd.fixed <- sqrt(diag(private$fit$cov.fixed))
    
  }
  
  # write covariance
  if(!inherits(covariance,"try-error")){
    private$fit$cov.fixed <- covariance
    colnames(private$fit$cov.fixed) <- names(private$fit$par.fixed)
    rownames(private$fit$cov.fixed) <- names(private$fit$par.fixed)
    std.dev <- try_with_warning_recovery(sqrt(diag(covariance)))
  }
  
  # write std.dev
  if(!inherits(std.dev,"try-error")){
    private$fit$sd.fixed <- std.dev
    names(private$fit$sd.fixed) <- names(private$fit$par.fixed)
  }
  
  # return
  return(invisible(self))
  
}

#######################################################
#######################################################
# helper for hessian inversion
#######################################################
#######################################################
iterative_hessian_inversion <- function(self, private){
  
  # OPTION 0 -----------------------------------
  # Invert full hessian
  temp.hessian = private$fit$nll.hessian
  covariance = try_with_warning_recovery(solve(temp.hessian))
  std.dev = try_with_warning_recovery(sqrt(diag((temp.hessian))))
  
  # OPTION 1 -----------------------------------
  # Remove all row/cols where the diagonal elements are smalller than threshold
  min.diag = 1e-8
  remove.ids <- which(diag(temp.hessian) < min.diag)
  if(inherits(covariance,"try-error") && any(remove.ids)){
    
    # try to invert reduced hessian
    covariance = try_with_warning_recovery(solve(temp.hessian[-remove.ids, -remove.ids]))
    std.dev = try_with_warning_recovery(sqrt(diag(covariance)))
    
    if(!inherits(covariance,"try-error")){
      covtemp = array(NA, dim=dim(temp.hessian))
      covtemp[-remove.ids, -remove.ids] = covariance
      covariance <- covtemp
    }
    if(!inherits(std.dev,"try-error")){
      stdtemp = rep(NA, length(private$fit$par.fixed))
      stdtemp[-remove.ids] <- std.dev
      std.dev <- stdtemp
    }
  }
  
  # OPTION 2 -----------------------------------
  # Remove small diagonal element one by one until solve is succesful
  failed.to.invert.hessian = TRUE
  id.diag.hess <- order(diag(temp.hessian))
  i = 1
  if(inherits(covariance,"try-error")) {
    while(failed.to.invert.hessian) {
      remove.ids <- id.diag.hess[1:i]
      covariance <- try_with_warning_recovery(solve(temp.hessian[-remove.ids,-remove.ids]))
      std.dev = try_with_warning_recovery(sqrt(diag(covariance)))
      
      # if succesful update results
      if(!inherits(covariance,"try-error")) {
        failed.to.invert.hessian <- FALSE
        # 
        covtemp = array(NA, dim=dim(temp.hessian))
        covtemp[-remove.ids, -remove.ids] = covariance
        covariance <- covtemp
        # 
        if(!inherits(std.dev,"try-error")){
          stdtemp = rep(NA, length(private$fit$par.fixed))
          stdtemp[-remove.ids] <- std.dev
          std.dev <- stdtemp
        }
      }
      
      i = i + 1
      # if unsuccesful break while loop
      if(i == nrow(temp.hessian)) {
        break
      }
    }
  }
  
  # return
  return(list(covariance=covariance, std.dev=std.dev))
}

#######################################################
#######################################################
# laplace report
#######################################################
#######################################################
laplace_report <- function(self, private, laplace.residuals){
  
  # lengths
  n.states <- private$number.of.states
  n.diff <- private$number.of.diffusions
  nobs <- private$number.of.observations
  .colnames <- c("t", private$obs.names)
  
  set.col.names <- function(object, colnames){
    colnames(object) <- colnames
    object
  }
  cbind_t <- function(object, names){
    set.col.names(cbind(private$data$t, object), names)
  }
  
  ############ states ############
  random.ids <- private$ode.timesteps.cumsum+1
  if(private$method %in% c("laplace")){
    
    # Smoothed States -----------------------------------
    temp.states <- data.frame(private$data$t, matrix(private$sdr$par.random, ncol=n.states)[random.ids, ])
    temp.sd <- data.frame(private$data$t, matrix(sqrt(private$sdr$diag.cov.random), ncol=n.states)[random.ids, ])
    names(temp.states) = c("t", private$state.names)
    names(temp.sd) = c("t", private$state.names)
    private$fit$states$mean$smoothed = temp.states
    private$fit$states$sd$smoothed = temp.sd
  }
  
  if(private$method %in% c("laplace.thygesen")){
    
    # Smoothed States -----------------------------------
    # Fill with zeros to complete perfect columns (dB's have 1 missing elements compared to random effect states)
    # for easier column extraction
    par.random <- c(private$sdr$par.random, numeric(n.diff))
    sd.random <- c(sqrt(private$sdr$diag.cov.random), numeric(n.diff))
    temp.states <- data.frame(private$data$t, matrix(par.random, ncol=n.states+n.diff)[random.ids, 1:n.states])
    temp.sd <- data.frame(private$data$t, matrix(sd.random, ncol=n.states+n.diff)[random.ids, 1:n.states])
    names(temp.states) = c("t", private$state.names)
    names(temp.sd) = c("t", private$state.names)
    private$fit$states$mean$smoothed = temp.states
    private$fit$states$sd$smoothed = temp.sd
  }
  
  ############ observations and residuals ############
  if(private$method %in% c("laplace","laplace.thygesen")) {
    
    # compute one-step residuals
    if(laplace.residuals){
      
      message("Calculating one-step ahead residuals...")
      res <- RTMB::oneStepPredict(private$nll,
                                  observation.name="obsMat",
                                  method="oneStepGaussian",
                                  trace=FALSE,
                                  parallel=TRUE)
      private$fit$residuals$residuals = cbind_t(matrix(res[["residual"]], ncol=nobs), .colnames)
      private$fit$observations$mean$smoothed <- cbind_t(matrix(res[["mean"]], ncol=nobs), .colnames)
      private$fit$observations$sd$smoothed <- cbind_t(matrix(res[["sd"]], ncol=nobs), .colnames)
      
    }
  }
  
}

#######################################################
#######################################################
# report states
#######################################################
#######################################################
laplace_states <- function(self, private){
  
  # lengths
  n.states <- private$number.of.states
  n.diff <- private$number.of.diffusions
  
  ############ LAPLACE ############
  
  random.ids <- private$ode.timesteps.cumsum+1
  
  if(private$method %in% c("laplace")){
    
    # Smoothed States -----------------------------------
    temp.states <- data.frame(private$data$t, matrix(private$sdr$par.random, ncol=n.states)[random.ids, ])
    temp.sd <- data.frame(private$data$t, matrix(sqrt(private$sdr$diag.cov.random), ncol=n.states)[random.ids, ])
    names(temp.states) = c("t", private$state.names)
    names(temp.sd) = c("t", private$state.names)
    private$fit$states$mean$smoothed = temp.states
    private$fit$states$sd$smoothed = temp.sd
  }
  
  if(private$method %in% c("laplace.thygesen")){
    
    # Smoothed States -----------------------------------
    # Fill with zeros to complete perfect columns (dB's have 1 missing elements compared to random effect states)
    # for easier column extraction
    par.random <- c(private$sdr$par.random, numeric(n.diff))
    sd.random <- c(sqrt(private$sdr$diag.cov.random), numeric(n.diff))
    temp.states <- data.frame(private$data$t, matrix(par.random, ncol=n.states+n.diff)[random.ids, 1:n.states])
    temp.sd <- data.frame(private$data$t, matrix(sd.random, ncol=n.states+n.diff)[random.ids, 1:n.states])
    names(temp.states) = c("t", private$state.names)
    names(temp.sd) = c("t", private$state.names)
    private$fit$states$mean$smoothed = temp.states
    private$fit$states$sd$smoothed = temp.sd
  }
  
  
  # return
  return(invisible(self))
}

#######################################################
#######################################################
# report residuals
#######################################################
#######################################################
report_residuals_and_observations <- function(laplace.residuals, self , private) {
  
  nobs <- private$number.of.observations
  .colnames <- c("t", private$obs.names)
  
  cbind_t <- function(object, names){
    set.colnames(cbind(private$data$t, object), names)
  }
  
  ############ LAPLACE ############
  if(private$method %in% c("laplace","laplace.thygesen")) {
    
    # compute one-step residuals
    if(laplace.residuals){
      
      message("Calculating one-step ahead residuals...")
      res <- RTMB::oneStepPredict(private$nll,
                                  observation.name="obsMat",
                                  method="oneStepGaussian",
                                  trace=FALSE,
                                  parallel=TRUE)
      private$fit$residuals$residuals = cbind_t(matrix(res[["residual"]], ncol=nobs), .colnames)
      private$fit$observations$mean$smoothed <- cbind_t(matrix(res[["mean"]], ncol=nobs), .colnames)
      private$fit$observations$sd$smoothed <- cbind_t(matrix(res[["sd"]], ncol=nobs), .colnames)
      
    }
  }
  
  # return
  return(invisible(self))
  
}
