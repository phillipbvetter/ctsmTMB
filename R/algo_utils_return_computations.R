#######################################################
#######################################################
# likelihood value, gradient and hessian
#######################################################
#######################################################
compute_mle_gradient_and_hessian <- function(self, private){
  
  # MLE
  private$fit$nll = private$opt$objective
  
  # MLE gradient
  private$fit$nll.gradient = try_withWarningRecovery(
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
    private$fit$nll.hessian = try_withWarningRecovery(
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
    std.dev <- try_withWarningRecovery(sqrt(diag(covariance)))
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
  covariance = try_withWarningRecovery(solve(temp.hessian))
  std.dev = try_withWarningRecovery(sqrt(diag((temp.hessian))))
  
  # OPTION 1 -----------------------------------
  # Remove all row/cols where the diagonal elements are smalller than threshold
  min.diag = 1e-8
  remove.ids <- which(diag(temp.hessian) < min.diag)
  if(inherits(covariance,"try-error") && any(remove.ids)){
    
    # try to invert reduced hessian
    covariance = try_withWarningRecovery(solve(temp.hessian[-remove.ids, -remove.ids]))
    std.dev = try_withWarningRecovery(sqrt(diag(covariance)))
    
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
      covariance <- try_withWarningRecovery(solve(temp.hessian[-remove.ids,-remove.ids]))
      std.dev = try_withWarningRecovery(sqrt(diag(covariance)))
      
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
# compute filters
#######################################################
#######################################################
get_state_report <- function(self, private){
  
  estimated.pars <- self$getParameters()[,"estimate"]
  switch(private$method,
         # kalman filters
         ekf = {rep <- ekf_filter_r(estimated.pars, self, private)},
         # ekf = {rep <- lkf_ekf_ukf_filter_rcpp(estimated.pars, self, private)},
         ekf.cpp = {rep <- private$nll$report(estimated.pars)},
         lkf = {rep <- lkf_filter_r(estimated.pars, self, private)},
         lkf.cpp = {rep <- private$nll$report(estimated.pars)},
         ukf = {rep <- ukf_filter_r(estimated.pars, self, private)},
         ukf.cpp = {rep <- private$nll$report(estimated.pars)},
         # laplace methods
         laplace = {rep <- NULL},
         laplace.thygesen = {rep <- NULL}
  )
  
  return(rep)
}

#######################################################
#######################################################
# report states
#######################################################
#######################################################
compute_return_states <- function(rep, self, private){
  
  # lengths
  n.states <- private$number.of.states
  n.diff <- private$number.of.diffusions
  
  ############# KALMAN ############
  if(private$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")) {
    
    # Prior States
    temp.states = try_withWarningRecovery(data.frame(private$data$t, t(do.call(cbind,rep$xPrior))))
    temp.sd = try_withWarningRecovery(data.frame(private$data$t, sqrt(do.call(rbind,lapply(rep$pPrior,diag)))))
    temp.cov <- rep$pPrior
    names(temp.states) = c("t", private$state.names)
    names(temp.sd) = c("t",private$state.names)
    names(temp.cov) = paste0("t = ", private$data$t)
    private$fit$states$mean$prior = temp.states
    private$fit$states$sd$prior = temp.sd
    private$fit$states$cov$prior = temp.cov
    
    # Posterior States
    temp.states = try_withWarningRecovery(data.frame(private$data$t, t(do.call(cbind,rep$xPost))))
    temp.sd = try_withWarningRecovery(data.frame(private$data$t, sqrt(do.call(rbind,lapply(rep$pPost,diag)))))
    temp.cov <- rep$pPost
    names(temp.states) = c("t",private$state.names)
    names(temp.sd) = c("t",private$state.names)
    names(temp.cov) = paste("t = ",private$data$t,sep="")
    private$fit$states$mean$posterior = temp.states
    private$fit$states$sd$posterior = temp.sd
    private$fit$states$cov$posterior = temp.cov
    
  }
  
  ############ LAPLACE ############
  if(private$method %in% c("laplace","laplace.thygesen")){
    
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
    
  }
  
  # return
  return(invisible(self))
}

#######################################################
#######################################################
# report residuals
#######################################################
#######################################################
report_residuals <- function(rep, laplace.residuals, self , private) {
  
  ############ KALMAN ############
  if(private$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")) {
    rowNAs = as.matrix(!is.na(private$data[private$obs.names]))
    sumrowNAs = rowSums(rowNAs)
    
    temp.res = matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
    temp.var =  matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
    
    innovation = rep$Innovation
    innovation.cov = rep$InnovationCovariance
    names(innovation.cov) = paste("t = ", private$data$t, sep="")
    for (i in 1:nrow(private$data)) {
      if (sumrowNAs[i] > 0) {
        temp.res[i,rowNAs[i,]] = innovation[[i]]
        temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
      }
    }
    temp.res <- data.frame(private$data$t, temp.res)
    temp.sd = data.frame(private$data$t, try_withWarningRecovery(sqrt(temp.var)))
    names(temp.res) = c("t", private$obs.names)
    names(temp.sd) = c("t", private$obs.names)
    # should we remove the empty matrices?
    # innovation.cov = innovation.cov[sumrowNAs!=0]
    
    private$fit$residuals$residuals = temp.res
    private$fit$residuals$sd = temp.sd
    private$fit$residuals$normalized = temp.res
    # normalize residuals (removing time column [,-1])
    private$fit$residuals$normalized[,-1] = temp.res[,-1]/temp.sd[,-1]
    private$fit$residuals$cov = innovation.cov
    
  }
  
  ############ LAPLACE ############
  if(private$method %in% c("laplace","laplace.thygesen")) {
    
    # compute one-step residuals
    if(laplace.residuals){
      
      rowNAs = as.matrix(!is.na(do.call(cbind,private$data[private$obs.names]))[-1,])
      sumrowNAs = rowSums(rowNAs)
      
      message("Calculating one-step ahead residuals...")
      res <- RTMB::oneStepPredict(private$nll,
                                  observation.name="obsMat",
                                  method="oneStepGaussian",
                                  trace=TRUE)
      nr <- nrow(private$data)
      temp.res <- data.frame(private$data$t, matrix(res[["residual"]],ncol=private$number.of.observations))
      temp.sd <- data.frame(private$data$t, matrix(res[["sd"]],ncol=private$number.of.observations))
      names(temp.res) = c("t", private$obs.names)
      names(temp.sd) = c("t", private$obs.names)
      
      private$fit$residuals$residuals <- temp.res
      private$fit$residuals$sd <- temp.sd
      private$fit$residuals$normalized <- temp.res
      # The first column is the time column why it is removed
      private$fit$residuals$normalized[,-1] <- temp.res[,-1]/temp.sd[,-1]
      
    }
  }
  
  # return
  return(invisible(self))
  
}

#######################################################
#######################################################
# report observations
#######################################################
#######################################################
report_observations <- function(self , private) {
  
  ############ KALMAN ############
  if(private$method %in% c("ekf","ekf.cpp","lkf","lkf.cpp","ukf","ukf.cpp")) {
    
    # Observations -----------------------------------
    free.and.fixed.parameters <- self$getParameters()[,"estimate"]
    names(free.and.fixed.parameters) <- private$parameter.names
    
    listofvariables.prior = c(
      as.list(private$fit$states$mean$prior[-1]),
      as.list(free.and.fixed.parameters),
      as.list(private$data)
    )
    listofvariables.posterior = c(
      as.list(private$fit$states$mean$posterior[-1]),
      as.list(free.and.fixed.parameters),
      as.list(private$data)
    )
    obs.df.prior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.prior)})
    )
    obs.df.posterior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.posterior)})
    )
    
    private$fit$observations$mean$prior = data.frame(t=private$data$t, obs.df.prior)
    private$fit$observations$mean$posterior = data.frame(t=private$data$t, obs.df.posterior)
    
    # Observation variances -----------------------------------
    
    # The observation variance (to first order) is:
    # y = h(x) + e -> var(y) = dhdx var(x) dhdx^T + var(e)
    
    n <- private$number.of.states
    m <- private$number.of.observations
    
    # create expression for observation covariance to be 'eval'-uated
    jac.h = c()
    for(i in seq_along(private$obs.names)){
      for(j in seq_along(private$state.names)){
        jac.h = c(jac.h, private$diff.terms.obs[[i]][[j]])
      }
    }
    dhdx <- parse(text=sprintf("matrix(c(%s),nrow=%s,ncol=%s)", paste(jac.h,collapse=","), m, n))[[1]]
    obs.var <- c()
    for(i in seq_along(private$obs.var.trans)){
      obs.var = c(obs.var, private$obs.var.trans[[i]]$rhs)
    }
    obsvar <- parse(text=sprintf("diag(%s)*c(%s)", m, paste(obs.var,collapse=",")))[[1]]
    yCov <- substitute(dhdx %*% xCov %*% t(dhdx) + eCov, list(dhdx=dhdx, eCov=obsvar))
    
    # evaluate prior and posterior variance
    list.of.parameters <- c(
      as.list(private$fit$par.fixed),
      lapply(private$fixed.pars, function(x) x$initial)
    )
    
    # prior
    obsvar.prior <- list()
    for(i in seq_along(private$data$t)){
      obsvar.prior[[i]] <- eval(expr = yCov,
                                envir = c(list.of.parameters,
                                          as.list(private$fit$states$mean$prior[i,-1]),
                                          list(xCov = private$fit$states$cov$prior[[i]]),
                                          as.list(private$data[i,-1])
                                ))
    }
    names(obsvar.prior) <- names(private$fit$states$cov$prior)
    private$fit$observations$cov$prior <- obsvar.prior
    obs.sd.prior <- cbind(private$fit$states$mean$prior["t"], do.call(rbind, lapply(obsvar.prior, diag)))
    rownames(obs.sd.prior) <- NULL
    names(obs.sd.prior) <- c("t",private$obs.names)
    private$fit$observations$sd$prior <- obs.sd.prior
    
    # posterior
    obsvar.post <- list()
    for(i in seq_along(private$data$t)){
      obsvar.post[[i]] <- eval(expr = yCov,
                               envir = c(list.of.parameters,
                                         as.list(private$fit$states$mean$posterior[i,-1]),
                                         list(xCov = private$fit$states$cov$posterior[[i]]),
                                         as.list(private$data[i,-1])
                               ))
    }
    names(obsvar.post) <- names(private$fit$states$cov$posterior)
    private$fit$observations$cov$posterior <- obsvar.post
    obs.sd.post <- cbind(private$fit$states$mean$posterior["t"], do.call(rbind, lapply(obsvar.post, diag)))
    rownames(obs.sd.post) <- NULL
    names(obs.sd.post) <- c("t",private$obs.names)
    private$fit$observations$sd$posterior <- obs.sd.post
    
  }
  
  ############ LAPLACE ############
  if(private$method %in% c("laplace","laplace.thygesen")) {
    
    # Observations -----------------------------------
    n <- private$number.of.states
    m <- private$number.of.observations
    
    listofvariables.smoothed = c(
      # states
      as.list(private$fit$states$mean$smoothed[-1]),
      # estimated free parameters 
      as.list(private$fit$par.fixed),
      # fixed parameters
      lapply(private$fixed.pars, function(x) x$initial),
      # inputs
      as.list(private$data)
    )
    obs.df.smoothed = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.smoothed)})
    )
    private$fit$observations$mean$smoothed = data.frame(t=private$data$t, obs.df.smoothed)
    
    # Observation variances -----------------------------------
    # The observation variance (to first order) is: 
    #
    # y = h(x) + e -> var(y) = dhdx var(x) dhdx^T + var(e)
    jac.h = c()
    for(i in seq_along(private$obs.names)){
      for(j in seq_along(private$state.names)){
        jac.h = c(jac.h, private$diff.terms.obs[[i]][[j]])
      }
    }
    dhdx <- parse(text=sprintf("matrix(c(%s),nrow=%s,ncol=%s)", paste(jac.h,collapse=","), m, n))[[1]]
    obs.var <- c()
    for(i in seq_along(private$obs.var.trans)){
      obs.var = c(obs.var, private$obs.var.trans[[i]]$rhs)
    }
    obsvar <- parse(text=sprintf("diag(%s)*c(%s)", m, paste(obs.var,collapse=",")))[[1]]
    yCov <- substitute(dhdx %*% xCov %*% t(dhdx) + eCov, list(dhdx=dhdx, eCov=obsvar))
    
    # Evaluate prior and posterior variance
    list.of.parameters <- c(
      as.list(private$fit$par.fixed),
      lapply(private$fixed.pars, function(x) x$initial)
    )
    # prior
    # obsvar.smooth <- list()
    # for(i in seq_along(private$data$t)){
    #   obsvar.smooth[[i]] <- eval(expr = yCov,
    #                             envir = c(list.of.parameters,
    #                                       as.list(private$fit$states$mean$smoothed[i,-1]),
    #                                       list(xCov = private$fit$states$cov$smoothed[[i]]),
    #                                       as.list(private$data[i,-1])
    #                             ))
    # }
    # names(obsvar.smooth) <- names(private$fit$states$cov$prior)
    # private$fit$observations$cov$prior <- obsvar.prior
    # obs.sd.prior <- cbind(private$fit$states$mean$prior["t"], do.call(rbind, lapply(obsvar.prior, diag)))
    # rownames(obs.sd.prior) <- NULL
    # names(obs.sd.prior) <- c("t",private$obs.names)
    # private$fit$observations$sd$prior <- obs.sd.prior

  }
  
  # return
  return(invisible(self))
  
}
