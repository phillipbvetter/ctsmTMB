#######################################################
#######################################################
# likelihood value, gradient and hessian
#######################################################
#######################################################
compute_mle_gradient_and_hessian <- function(self, private){
  
  # MLE
  private$results$fit$nll = private$results$opt$objective
  
  # MLE gradient
  private$results$fit$nll.gradient = try_with_warning_recovery(
    {
      nll.grad = as.vector(private$nll$gr(private$results$opt$par))
      names(nll.grad) = names(private$model$free.pars)
      nll.grad
    }
  )
  if (inherits(private$results$fit$nll.gradient,"try-error")) {
    private$results$fit$nll.gradient = NULL
  }
  
  # MLE Hessian
  if(private$algo.settings$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")){
    private$results$fit$nll.hessian = try_with_warning_recovery(
      {
        nll.hess = private$nll$he(private$results$opt$par)
        rownames(nll.hess) = names(private$model$free.pars)
        colnames(nll.hess) = names(private$model$free.pars)
        nll.hess
      }
    )
    if (inherits(private$results$fit$nll.hessian, "try-error")) {
      private$results$fit$nll.hessian = NULL
    }
  }
  
  return(invisible(self))
}

# function that writes parameter estimates and computes std. errors
# from either inverting the hessian (kalman filters) or from 
# using TMB::sdreport (laplace)
compute_mle_parameters_and_std_errors <- function(self, private){
  
  private$results$fit$par.fixed = private$results$opt$par
  
  # allocate 
  n.fixed.pars <- private$dims$fixed.pars
  private$results$fit$sd.fixed = rep(NA, n.fixed.pars)
  private$results$fit$cov.fixed = array(NA, dim=rep(n.fixed.pars,2))
  private$results$fit$tvalue = rep(NA, n.fixed.pars)
  private$results$fit$Pr.tvalue = rep(NA, n.fixed.pars)
  
  ############ LAPLACE ############
  if(private$algo.settings$method %in% c("laplace","laplace.thygesen")){
    private$results$fit$cov.fixed <- private$sdr$cov.fixed
    private$results$fit$sd.fixed <- sqrt(diag(private$results$fit$cov.fixed))
  }
  
  ############ KALMAN ############
  if(private$algo.settings$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")){
    calculate_covariance_from_hessian(self, private)
  }
  
  # t-values and Pr( t > t_test ) -----------------------------------
  private$results$fit$tvalue = private$results$fit$par.fixed / private$results$fit$sd.fixed
  # The degrees of fredom are number of (non NA) observations minus number of (free) parameters
  freedom.degrees <- sum(!is.na(private$data[private$names$obs])) - private$dims$free.pars
  private$results$fit$Pr.tvalue <- rep(NA, length(private$results$fit$par.fixed))
  if(freedom.degrees > 0.1){
    private$results$fit$Pr.tvalue = 2 * pt(q=abs(private$results$fit$tvalue), df=freedom.degrees, lower.tail=FALSE)
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
  if(!is.null(private$results$fit$nll.hessian)){
    
    out <- iterative_hessian_inversion(self, private)
    covariance <- out$covariance
    std.dev <- out$std.dev
    
  } else {
    
    # if the hessian could not be computed, then we may try sdreport instead
    # (sometimes more stable???)
    private$sdr <- TMB::sdreport(private$nll, getJointPrecision=FALSE)
    private$results$fit$cov.fixed <- private$sdr$cov.fixed
    private$results$fit$sd.fixed <- sqrt(diag(private$results$fit$cov.fixed))
    
  }
  
  # write covariance
  if(!inherits(covariance,"try-error")){
    private$results$fit$cov.fixed <- covariance
    colnames(private$results$fit$cov.fixed) <- names(private$results$fit$par.fixed)
    rownames(private$results$fit$cov.fixed) <- names(private$results$fit$par.fixed)
    std.dev <- try_with_warning_recovery(sqrt(diag(covariance)))
  }
  
  # write std.dev
  if(!inherits(std.dev,"try-error")){
    private$results$fit$sd.fixed <- std.dev
    names(private$results$fit$sd.fixed) <- names(private$results$fit$par.fixed)
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
  temp.hessian = private$results$fit$nll.hessian
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
      stdtemp = rep(NA, length(private$results$fit$par.fixed))
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
          stdtemp = rep(NA, length(private$results$fit$par.fixed))
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
  n.states <- private$dims$states
  n.diff <- private$dims$diffusions
  nobs <- private$dims$observations
  .colnames <- c("t", private$names$obs)
  
  set.col.names <- function(object, colnames){
    colnames(object) <- colnames
    object
  }
  cbind_t <- function(object, names){
    set.col.names(cbind(private$data$t, object), names)
  }
  
  ############ states ############
  random.ids <- private$algo.settings$ode.timesteps.cumsum+1
  if(private$algo.settings$method %in% c("laplace")){
    
    # Smoothed States -----------------------------------
    temp.states <- data.frame(private$data$t, matrix(private$sdr$par.random, ncol=n.states)[random.ids, ])
    temp.sd <- data.frame(private$data$t, matrix(sqrt(private$sdr$diag.cov.random), ncol=n.states)[random.ids, ])
    names(temp.states) = c("t", private$names$states)
    names(temp.sd) = c("t", private$names$states)
    private$results$fit$states$mean$smoothed = temp.states
    private$results$fit$states$sd$smoothed = temp.sd
  }
  
  if(private$algo.settings$method %in% c("laplace.thygesen")){
    
    # Smoothed States -----------------------------------
    # Fill with zeros to complete perfect columns (dB's have 1 missing elements compared to random effect states)
    # for easier column extraction
    par.random <- c(private$sdr$par.random, numeric(n.diff))
    sd.random <- c(sqrt(private$sdr$diag.cov.random), numeric(n.diff))
    temp.states <- data.frame(private$data$t, matrix(par.random, ncol=n.states+n.diff)[random.ids, 1:n.states])
    temp.sd <- data.frame(private$data$t, matrix(sd.random, ncol=n.states+n.diff)[random.ids, 1:n.states])
    names(temp.states) = c("t", private$names$states)
    names(temp.sd) = c("t", private$names$states)
    private$results$fit$states$mean$smoothed = temp.states
    private$results$fit$states$sd$smoothed = temp.sd
  }
  
  ############ observations and residuals ############
  if(private$algo.settings$method %in% c("laplace","laplace.thygesen")) {
    
    # compute one-step residuals
    if(laplace.residuals){
      
      message("Calculating one-step ahead residuals...")
      res <- RTMB::oneStepPredict(private$nll,
                                  observation.name="obsMat",
                                  method="oneStepGaussian",
                                  trace=FALSE,
                                  parallel=TRUE)
      private$results$fit$residuals$residuals = cbind_t(matrix(res[["residual"]], ncol=nobs), .colnames)
      private$results$fit$observations$mean$smoothed <- cbind_t(matrix(res[["mean"]], ncol=nobs), .colnames)
      private$results$fit$observations$sd$smoothed <- cbind_t(matrix(res[["sd"]], ncol=nobs), .colnames)
      
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
  n.states <- private$dims$states
  n.diff <- private$dims$diffusions
  
  ############ LAPLACE ############
  
  random.ids <- private$algo.settings$ode.timesteps.cumsum+1
  
  if(private$algo.settings$method %in% c("laplace")){
    
    # Smoothed States -----------------------------------
    temp.states <- data.frame(private$data$t, matrix(private$sdr$par.random, ncol=n.states)[random.ids, ])
    temp.sd <- data.frame(private$data$t, matrix(sqrt(private$sdr$diag.cov.random), ncol=n.states)[random.ids, ])
    names(temp.states) = c("t", private$names$states)
    names(temp.sd) = c("t", private$names$states)
    private$results$fit$states$mean$smoothed = temp.states
    private$results$fit$states$sd$smoothed = temp.sd
  }
  
  if(private$algo.settings$method %in% c("laplace.thygesen")){
    
    # Smoothed States -----------------------------------
    # Fill with zeros to complete perfect columns (dB's have 1 missing elements compared to random effect states)
    # for easier column extraction
    par.random <- c(private$sdr$par.random, numeric(n.diff))
    sd.random <- c(sqrt(private$sdr$diag.cov.random), numeric(n.diff))
    temp.states <- data.frame(private$data$t, matrix(par.random, ncol=n.states+n.diff)[random.ids, 1:n.states])
    temp.sd <- data.frame(private$data$t, matrix(sd.random, ncol=n.states+n.diff)[random.ids, 1:n.states])
    names(temp.states) = c("t", private$names$states)
    names(temp.sd) = c("t", private$names$states)
    private$results$fit$states$mean$smoothed = temp.states
    private$results$fit$states$sd$smoothed = temp.sd
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
  
  nobs <- private$dims$observations
  .colnames <- c("t", private$names$obs)
  
  cbind_t <- function(object, names){
    set.colnames(cbind(private$data$t, object), names)
  }
  
  ############ LAPLACE ############
  if(private$algo.settings$method %in% c("laplace","laplace.thygesen")) {
    
    # compute one-step residuals
    if(laplace.residuals){
      
      message("Calculating one-step ahead residuals...")
      res <- RTMB::oneStepPredict(private$nll,
                                  observation.name="obsMat",
                                  method="oneStepGaussian",
                                  trace=FALSE,
                                  parallel=TRUE)
      private$results$fit$residuals$residuals = cbind_t(matrix(res[["residual"]], ncol=nobs), .colnames)
      private$results$fit$observations$mean$smoothed <- cbind_t(matrix(res[["mean"]], ncol=nobs), .colnames)
      private$results$fit$observations$sd$smoothed <- cbind_t(matrix(res[["sd"]], ncol=nobs), .colnames)
      
    }
  }
  
  # return
  return(invisible(self))
  
}
