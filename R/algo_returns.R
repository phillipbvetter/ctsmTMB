#######################################################
#######################################################
# RETURN FIT FOR EKF ESTIMATION
#######################################################
#######################################################

calculate_fit_statistics_ekf <- function(self, private){
  
  # Initialization and Clearing -----------------------------------
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  # clear fit
  private$fit = NULL
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  
  # Parameters and Uncertainties -----------------------------------
  
  # objective value
  private$fit$nll = private$opt$objective
  
  # gradient
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
  
  # # hessian
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
  
  # parameter estimates and standard deviation
  npars <- length(private$opt$par)
  private$fit$par.fixed = private$opt$par
  private$fit$sd.fixed = rep(NA,npars)
  private$fit$cov.fixed = array(NA,dim=c(npars,npars))
  
  # parameter std. error and full covariance by hessian inversion
  if(!is.null(private$fit$nll.hessian)){
    
    # OPTION 0 -----------------------------------
    # Invert full hessian
    temp.hessian = private$fit$nll.hessian
    covariance = try_withWarningRecovery(solve(temp.hessian))
    
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
    if(inherits(covariance,"try-error")){
      while(failed.to.invert.hessian){
        remove.ids <- id.diag.hess[1:i]
        covariance <- try_withWarningRecovery(solve(temp.hessian[-remove.ids,-remove.ids]))
        std.dev = try_withWarningRecovery(sqrt(diag(covariance)))
        i = i + 1
        
        # if unsuccesful break while loop
        if(i == nrow(temp.hessian)){
          break
        }
        
        # if succesful update results
        if(!inherits(covariance,"try-error")){
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
      }
    }
    
    if(!inherits(covariance,"try-error")){
      private$fit$cov.fixed <- covariance
      colnames(private$fit$cov.fixed) <- names(private$fit$par.fixed)
      rownames(private$fit$cov.fixed) <- names(private$fit$par.fixed)
      std.dev <- try_withWarningRecovery(sqrt(diag(covariance)))
    }
    if(!inherits(std.dev,"try-error")){
      private$fit$sd.fixed <- std.dev
      names(private$fit$sd.fixed) <- names(private$fit$par.fixed)
    }
    
  }
  
  # States -----------------------------------
  
  # Extract reported items from nll
  estimated_pars <- self$getParameters()[,"estimate"]
  rep <- ekf_filter_r(estimated_pars, self, private)
  
  # Prior States
  temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPrior))))
  temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPrior,diag)))))
  
  colnames(temp.states) = c("t", private$state.names)
  colnames(temp.sd) = c("t",private$state.names)
  private$fit$states$mean$prior = as.data.frame(temp.states)
  private$fit$states$sd$prior = as.data.frame(temp.sd)
  private$fit$states$cov$prior = rep$pPrior
  names(private$fit$states$cov$prior) = paste("t = ",private$data$t,sep="")
  
  # Posterior States
  temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPost))))
  temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPost,diag)))))
  colnames(temp.states) = c("t",private$state.names)
  colnames(temp.sd) = c("t",private$state.names)
  private$fit$states$mean$posterior = as.data.frame(temp.states)
  private$fit$states$sd$posterior = as.data.frame(temp.sd)
  private$fit$states$cov$posterior = rep$pPost
  names(private$fit$states$cov$posterior) = paste("t = ",private$data$t,sep="")
  
  # Residuals -----------------------------------
  rowNAs = as.matrix(!is.na(private$data[private$obs.names]))
  sumrowNAs = rowSums(rowNAs)
  
  nr <- nrow(private$data)
  temp.res = matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
  temp.var =  matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
  
  nr <- nrow(private$data)
  innovation = rep$Innovation
  innovation.cov = rep$InnovationCovariance
  names(innovation.cov) = paste("t = ", private$data$t, sep="")
  for (i in 1:nr) {
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
  private$fit$residuals$normalized[,-1] = temp.res[,-1]/temp.sd[,-1]
  private$fit$residuals$cov = innovation.cov
  
  
  # Observations -----------------------------------
  
  # We need all states, inputs and parameter values to evaluate the observation
  # put them in a list
  listofvariables.prior = c(
    # states
    as.list(private$fit$states$mean$prior[-1]),
    # estimated free parameters
    as.list(private$fit$par.fixed),
    # fixed parameters
    lapply(private$fixed.pars, function(x) x$initial),
    # inputs
    as.list(private$data)
  )
  
  listofvariables.posterior = c(
    # states
    as.list(private$fit$states$mean$posterior[-1]),
    # estimated free parameters
    as.list(private$fit$par.fixed),
    # fixed parameters
    lapply(private$fixed.pars, function(x) x$initial),
    # inputs
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
  #
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
  
  # t-values and Pr( t > t_test ) -----------------------------------
  
  private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
  private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
  
  # clone private into fit -----------------------------------
  private$fit$private <- self$clone()$.__enclos_env__$private
  
  # set s3 class -----------------------------------
  class(private$fit) = "ctsmTMB.fit"
  
  # return -----------------------------------
  return(invisible(self))
  
}


#######################################################
#######################################################
# RETURN FIT FOR LKF RTMB
#######################################################
#######################################################

calculate_fit_statistics_lkf <- function(self, private){
  
  # Initialization and Clearing -----------------------------------
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  # clear fit
  private$fit = NULL
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  
  # Parameters and Uncertainties -----------------------------------
  
  # objective value
  private$fit$nll = private$opt$objective
  
  # gradient
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
  
  # # hessian
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
  
  # parameter estimates and standard deviation
  npars <- length(private$opt$par)
  private$fit$par.fixed = private$opt$par
  private$fit$sd.fixed = rep(NA,npars)
  private$fit$cov.fixed = array(NA, dim=c(npars,npars))
  private$fit$tvalue = rep(NA,npars)
  private$fit$Pr.tvalue = rep(NA,npars)
  
  # parameter std. error and full covariance by hessian inversion
  if(!is.null(private$fit$nll.hessian)){
    
    # OPTION 0 -----------------------------------
    # Invert full hessian
    temp.hessian = private$fit$nll.hessian
    covariance = try_withWarningRecovery(solve(temp.hessian))
    
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
    if(inherits(covariance,"try-error")){
      while(failed.to.invert.hessian){
        remove.ids <- id.diag.hess[1:i]
        covariance <- try_withWarningRecovery(solve(temp.hessian[-remove.ids,-remove.ids]))
        std.dev = try_withWarningRecovery(sqrt(diag(covariance)))
        i = i + 1
        
        # if unsuccesful break while loop
        if(i == nrow(temp.hessian)){
          break
        }
        
        # if succesful update results
        if(!inherits(covariance,"try-error")){
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
      }
    }
    
    # save results
    if(!inherits(covariance,"try-error")){
      private$fit$cov.fixed <- covariance
      colnames(private$fit$cov.fixed) <- names(private$fit$par.fixed)
      rownames(private$fit$cov.fixed) <- names(private$fit$par.fixed)
      std.dev <- try_withWarningRecovery(sqrt(diag(covariance)))
    }
    if(!inherits(std.dev,"try-error")){
      private$fit$sd.fixed <- std.dev
      names(private$fit$sd.fixed) <- names(private$fit$par.fixed)
    }
    
  }
  
  # States -----------------------------------
  # Extract reported items from nll
  estimated_pars <- self$getParameters()[,"estimate"]
  rep <- try_withWarningRecovery((lkf_filter_r(estimated_pars, self, private)))
  
  if(!inherits(rep,"try-error")){
    
    # Prior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPrior))))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPrior,diag)))))
    
    colnames(temp.states) = c("t", private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$prior = as.data.frame(temp.states)
    private$fit$states$sd$prior = as.data.frame(temp.sd)
    private$fit$states$cov$prior = rep$pPrior
    names(private$fit$states$cov$prior) = paste("t = ",private$data$t,sep="")
    
    # Posterior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPost))))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPost,diag)))))
    colnames(temp.states) = c("t",private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$posterior = as.data.frame(temp.states)
    private$fit$states$sd$posterior = as.data.frame(temp.sd)
    private$fit$states$cov$posterior = rep$pPost
    names(private$fit$states$cov$posterior) = paste("t = ",private$data$t,sep="")
    
    # Residuals -----------------------------------
    rowNAs = as.matrix(!is.na(private$data[private$obs.names]))
    sumrowNAs = rowSums(rowNAs)
    
    nr <- nrow(private$data)
    temp.res = matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
    temp.var =  matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
    
    nr <- nrow(private$data)
    innovation = rep$Innovation
    innovation.cov = rep$InnovationCovariance
    names(innovation.cov) = paste("t = ", private$data$t, sep="")
    for (i in 1:nr) {
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
    private$fit$residuals$normalized[,-1] = temp.res[,-1]/temp.sd[,-1]
    private$fit$residuals$cov = innovation.cov
    
    
    # Observations -----------------------------------
    
    # We need all states, inputs and parameter values to evaluate the observation
    # put them in a list
    listofvariables.prior = c(
      # states
      as.list(private$fit$states$mean$prior[-1]),
      # estimated free parameters 
      as.list(private$fit$par.fixed),
      # fixed parameters
      lapply(private$fixed.pars, function(x) x$initial),
      # inputs
      as.list(private$data)
    )
    
    listofvariables.posterior = c(
      # states
      as.list(private$fit$states$mean$posterior[-1]),
      # estimated free parameters 
      as.list(private$fit$par.fixed),
      # fixed parameters
      lapply(private$fixed.pars, function(x) x$initial),
      # inputs
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
    #
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
    
    # t-values and Pr( t > t_test ) -----------------------------------
    
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
    
  } else {
    
    message("Unable to perform filtering with the following warning:\n\t", rep)
    
  }
  
  # clone and return -----------------------------------
  
  # clone private and return fit
  private$fit$private <- self$clone()$.__enclos_env__$private
  
  # set s3 class
  class(private$fit) = "ctsmTMB.fit"
  
  return(invisible(self))
}


#######################################################
#######################################################
# RETURN FIT FOR UKF RTMB
#######################################################
#######################################################

calculate_fit_statistics_ukf <- function(self, private){
  
  # Initialization and Clearing -----------------------------------
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  estimated.parameters <- private$opt$par
  
  # clear fit
  private$fit = NULL
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  
  # Parameters and Uncertainties -----------------------------------
  
  # objective value
  private$fit$nll = private$opt$objective
  
  # gradient
  private$fit$nll.gradient = try_withWarningRecovery(
    {
      nll.grad = as.vector(private$nll$gr(estimated.parameters))
      names(nll.grad) = names(private$free.pars)
      nll.grad
    }
  )
  if (inherits(private$fit$nll.gradient,"try-error")) {
    private$fit$nll.gradient = NULL
  }
  
  # # hessian
  private$fit$nll.hessian = try_withWarningRecovery(
    {
      nll.hess = private$nll$he(estimated.parameters)
      rownames(nll.hess) = names(private$free.pars)
      colnames(nll.hess) = names(private$free.pars)
      nll.hess
    }
  )
  if (inherits(private$fit$nll.hessian, "try-error")) {
    private$fit$nll.hessian = NULL
  }
  
  # parameter estimates and standard deviation
  npars <- length(estimated.parameters)
  private$fit$par.fixed = estimated.parameters
  private$fit$sd.fixed = rep(NA,npars)
  private$fit$cov.fixed = array(NA,dim=c(npars,npars))
  
  # parameter std. error and full covariance by hessian inversion
  if(!is.null(private$fit$nll.hessian)){
    
    # OPTION 0 -----------------------------------
    # Invert full hessian
    temp.hessian = private$fit$nll.hessian
    covariance = try_withWarningRecovery(solve(temp.hessian))
    
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
    if(inherits(covariance,"try-error")){
      while(failed.to.invert.hessian){
        remove.ids <- id.diag.hess[1:i]
        covariance <- try_withWarningRecovery(solve(temp.hessian[-remove.ids,-remove.ids]))
        std.dev = try_withWarningRecovery(sqrt(diag(covariance)))
        i = i + 1
        
        # if unsuccesful break while loop
        if(i == nrow(temp.hessian)){
          break
        }
        
        # if succesful update results
        if(!inherits(covariance,"try-error")){
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
      }
    }
    
    # save results
    if(!inherits(covariance,"try-error")){
      private$fit$cov.fixed <- covariance
      colnames(private$fit$cov.fixed) <- names(private$fit$par.fixed)
      rownames(private$fit$cov.fixed) <- names(private$fit$par.fixed)
      std.dev <- try_withWarningRecovery(sqrt(diag(covariance)))
    }
    if(!inherits(std.dev,"try-error")){
      private$fit$sd.fixed <- std.dev
      names(private$fit$sd.fixed) <- names(private$fit$par.fixed)
    }
    
  }
  
  # States -----------------------------------
  
  # Extract reported items from nll
  estimated.and.fixed.pars <- self$getParameters()[,"estimate"]
  rep <- ukf_filter_r(estimated.and.fixed.pars, self, private)
  
  private$fit$rep <- rep
  
  # Prior States
  temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPrior))))
  temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPrior,diag)))))
  
  colnames(temp.states) = c("t", private$state.names)
  colnames(temp.sd) = c("t",private$state.names)
  private$fit$states$mean$prior = as.data.frame(temp.states)
  private$fit$states$sd$prior = as.data.frame(temp.sd)
  private$fit$states$cov$prior = rep$pPrior
  names(private$fit$states$cov$prior) = paste("t = ",private$data$t,sep="")
  
  # Posterior States
  temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPost))))
  temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPost,diag)))))
  colnames(temp.states) = c("t",private$state.names)
  colnames(temp.sd) = c("t",private$state.names)
  private$fit$states$mean$posterior = as.data.frame(temp.states)
  private$fit$states$sd$posterior = as.data.frame(temp.sd)
  private$fit$states$cov$posterior = rep$pPost
  names(private$fit$states$cov$posterior) = paste("t = ",private$data$t,sep="")
  
  # Residuals -----------------------------------
  
  # rowNAs = as.matrix(!is.na(private$data[private$obs.names])[-1,])
  rowNAs = as.matrix(!is.na(private$data[private$obs.names]))
  sumrowNAs = rowSums(rowNAs)
  
  nr <- nrow(private$data)
  # temp.res = matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
  # temp.var =  matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
  temp.res = matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
  temp.var =  matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
  
  nr <- nrow(private$data)
  innovation = rep$Innovation
  innovation.cov = rep$InnovationCovariance
  names(innovation.cov) = paste("t = ", private$data$t, sep="")
  for (i in 1:nr) {
    if (sumrowNAs[i] > 0) {
      temp.res[i,rowNAs[i,]] = innovation[[i]]
      temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
    }
  }
  # temp.res = cbind(private$data$t, temp.res)
  # temp.sd = cbind(private$data$t, try_withWarningRecovery(sqrt(temp.var)))
  temp.res <- data.frame(private$data$t, temp.res)
  temp.sd = data.frame(private$data$t, try_withWarningRecovery(sqrt(temp.var)))
  names(temp.res) = c("t", private$obs.names)
  names(temp.sd) = c("t", private$obs.names)
  # should we remove the empty matrices?
  # innovation.cov = innovation.cov[sumrowNAs!=0]
  
  private$fit$residuals$residuals = temp.res
  private$fit$residuals$sd = temp.sd
  private$fit$residuals$normalized = temp.res
  private$fit$residuals$normalized[,-1] = temp.res[,-1]/temp.sd[,-1]
  private$fit$residuals$cov = innovation.cov
  
  
  # Observations -----------------------------------
  
  # We need all states, inputs and parameter values to evaluate the observation
  # put them in a list
  listofvariables.prior = c(
    # states
    as.list(private$fit$states$mean$prior[-1]),
    # estimated free parameters 
    as.list(private$fit$par.fixed),
    # fixed parameters
    lapply(private$fixed.pars, function(x) x$initial),
    # inputs
    as.list(private$data)
  )
  
  listofvariables.posterior = c(
    # states
    as.list(private$fit$states$mean$posterior[-1]),
    # estimated free parameters 
    as.list(private$fit$par.fixed),
    # fixed parameters
    lapply(private$fixed.pars, function(x) x$initial),
    # inputs
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
  #
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
  
  # t-values and Pr( t > t_test ) -----------------------------------
  
  private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
  private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
  
  # clone and return -----------------------------------
  
  # clone private and return fit
  # self.clone <- self$clone()$.__enclos_env__$private
  # private$fit$private = self.clone$.__enclos_env__$private
  private$fit$private <- self$clone()$.__enclos_env__$private
  
  # set s3 class
  class(private$fit) = "ctsmTMB.fit"
  
  return(invisible(self))
  
}

#######################################################
#######################################################
# RETURN FIT FOR LAPLACE
#######################################################
#######################################################
calculate_fit_statistics_laplace <- function(self, private, laplace.residuals){
  
  # Initialization and clearing -----------------------------------
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  # clear fit
  private$fit = NULL
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  n <- private$number.of.states
  m <- private$number.of.observations
  n.random <- nrow(private$data) - 1 #number of random effects
  
  # Parameters and Uncertainties -----------------------------------
  
  # objective
  private$fit$nll = private$opt$objective
  
  # gradient
  private$fit$nll.gradient = private$sdr$gradient.fixed
  names(private$fit$nll.gradient) = names(private$free.pars)
  
  # parameter estimates and standard deviation
  private$fit$par.fixed = private$opt$par
  private$fit$sd.fixed = diag(private$sdr$cov.fixed)
  
  # parameter covariance
  private$fit$cov.fixed = private$sdr$cov.fixed
  
  # random.ids <- head(private$ode.timesteps.cumsum+1,-1)
  random.ids <- private$ode.timesteps.cumsum+1
  # States (Smoothed) -----------------------------------
  temp.states <- matrix(private$sdr$par.random, ncol=n)[random.ids, ]
  temp.sd <- matrix(sqrt(private$sdr$diag.cov.random), ncol=n)[random.ids, ]
  
  # initialState <- unlist(private$tmb.initial.state[1,])
  # temp.states <- rbind(initialState, temp.states)
  # temp.sd <- rbind(initialState, temp.sd)
  
  temp.states <- cbind(private$data$t, temp.states)
  temp.sd <- cbind(private$data$t, temp.sd)
  
  private$fit$states$mean$smoothed = as.data.frame(temp.states)
  private$fit$states$sd$smoothed = as.data.frame(temp.sd)
  colnames(private$fit$states$sd$smoothed) = c("t",private$state.names)
  colnames(private$fit$states$mean$smoothed) = c("t",private$state.names)
  
  # Residuals -----------------------------------
  rowNAs = as.matrix(!is.na(do.call(cbind,private$data[private$obs.names]))[-1,])
  sumrowNAs = rowSums(rowNAs)
  
  # compute one-step residuals
  if(laplace.residuals){
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
    private$fit$residuals$normalized[,-1] <- temp.res[,-1]/temp.sd[,-1]
  }
  
  # t-values and Pr( t > t_test ) -----------------------------------
  private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
  private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
  
  # Observations -----------------------------------
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
  
  # clone and return -----------------------------------
  # clone private and return fit
  private$fit$private <- self$clone()$.__enclos_env__$private
  
  # set s3 class
  class(private$fit) = "ctsmTMB.fit"
}

#######################################################
#######################################################
# RETURN FIT FOR LAPLACE 2
#######################################################
#######################################################
calculate_fit_statistics_laplace2 <- function(self, private, laplace.residuals){
  
  # Initialization and clearing -----------------------------------
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  # clear fit
  private$fit = NULL
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.diff <- private$number.of.diffusions
  
  n <- private$number.of.states
  m <- private$number.of.observations
  n.random <- nrow(private$data) - 1 #number of random effects
  
  # Parameters and Uncertainties -----------------------------------
  
  # objective
  private$fit$nll = private$opt$objective
  
  # gradient
  private$fit$nll.gradient = private$sdr$gradient.fixed
  names(private$fit$nll.gradient) = names(private$free.pars)
  
  # parameter estimates and standard deviation
  private$fit$par.fixed = private$opt$par
  private$fit$sd.fixed = diag(private$sdr$cov.fixed)
  
  # parameter covariance
  private$fit$cov.fixed = private$sdr$cov.fixed
  
  # States (Smoothed) -----------------------------------
  par.random <- c(private$sdr$par.random, numeric(n.diff))
  sd.random <- c(sqrt(private$sdr$diag.cov.random), numeric(n.diff))
  private$fit$states$mean$smoothed <- data.frame(
    private$data$t,
    matrix(par.random, ncol=n.states+n.diff)[private$ode.timesteps.cumsum+1, 1:n.states]
  )
  private$fit$states$sd$smoothed <- data.frame(
    private$data$t,
    matrix(sd.random, ncol=n.states+n.diff)[private$ode.timesteps.cumsum+1, 1:n.states]
  )
  names(private$fit$states$sd$smoothed) = c("t",private$state.names)
  names(private$fit$states$mean$smoothed) = c("t",private$state.names)
  
  # Residuals -----------------------------------
  rowNAs = as.matrix(!is.na(do.call(cbind,private$data[private$obs.names]))[-1,])
  sumrowNAs = rowSums(rowNAs)
  # compute one-step residuals
  if(laplace.residuals){
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
    private$fit$residuals$normalized[,-1] <- temp.res[,-1]/temp.sd[,-1]
  }
  
  # t-values and Pr( t > t_test ) -----------------------------------
  private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
  private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
  
  # clone and return -----------------------------------
  # clone private and return fit
  private$fit$private <- self$clone()$.__enclos_env__$private
  
  # set s3 class
  class(private$fit) = "ctsmTMB.fit"
}
