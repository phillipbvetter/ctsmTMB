#######################################################
# MAKE RETURN DATA NICE AFTER OPTIMIZATION
#######################################################

create_return_fit = function(self, private, calculate.laplace.onestep.residuals) {
  
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  # clear fit
  private$fit = NULL
  
  # store the provided data in the fit
  private$fit$data = private$data
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  # store the object in fit - these gives access to e.g.
  private$fit$.__object__ = self$clone()
  
  ################################################
  # FOR KALMAN FILTERS
  ################################################
  
  if (any(private$method == c("ekf","ukf"))) {
    
    
    ################################################
    # BASICS
    ################################################
    
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
    
    # compute hessian
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
    
    # parameter estimates
    private$fit$par.fixed = private$opt$par
    
    # parameter std. error and full covariance by hessian inversion
    if(!is.null(private$fit$nll.hessian)){
      
      # Step 1 - invert full hessian
      temp.hessian = private$fit$nll.hessian
      covariance = try(solve(temp.hessian), silent=T)
      
      private$fit$cov.fixed = covariance
      private$fit$sd.fixed = try_withWarningRecovery(sqrt(diag(covariance)))
      
      if(inherits(private$fit$sd.fixed,"try-error")){
        private$fit$sd.fixed = rep(NA,length(private$fit$par.fixed))
      }
      if(inherits(private$fit$cov.fixed,"try-error")){
        private$fit$cov.fixed = matrix(NA,nrow=length(private$fit$par.fixed),ncol=length(private$fit$par.fixed))
      }
      
      # Options 1 - If the above fails, remove all row/cols where the diagonal
      # elements in small than min.diag
      min.diag = 1e-8
      keep.ids = !(diag(temp.hessian) < min.diag)
      if(inherits(covariance,"try-error") && any(keep.ids)){
        
        covariance = temp.hessian[keep.ids, keep.ids]
        covariance = try(solve(covariance, silent=T))
        
        sd.fixed = rep(NA,length(private$fit$par.fixed))
        sd.fixed[keep.ids] = try_withWarningRecovery(sqrt(diag(covariance)))
        private$fit$sd.fixed = sd.fixed
        
        cov.fixed = temp.hessian * NA
        cov.fixed[keep.ids, keep.ids] = covariance
        private$fit$cov.fixed = cov.fixed
      }
      
      # Option 2 - Recursive remove the smallest parameter
      # ids = sort(diag(temp.hessian), index.return=T)$ix
      # i = 0
      # while(inherits(hess,"try-error") && i < length(private$fit$par.fixed)){
      #   i = i + 1
      #   covariance = try(solve(temp.hessian[-ids[1:i],-ids[1:i]]), silent=T)
      # }
      
    }
    
    ################################################
    # STATES, RESIDUALS, OBSERVATIONS ETC.
    ################################################
    
    # Extract reported items from nll
    rep = private$nll$report()
    
    # Prior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, do.call(rbind,rep$xPrior)))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPrior,diag)))))
    
    colnames(temp.states) = c("t", private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$prior = as.data.frame(temp.states)
    private$fit$states$sd$prior = as.data.frame(temp.sd)
    private$fit$states$cov$prior = rep$pPrior
    names(private$fit$states$cov$prior) = paste("t = ",private$data$t,sep="")
    
    # Posterior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, do.call(rbind,rep$xPost)))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPost,diag)))))
    colnames(temp.states) = c("t",private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$posterior = as.data.frame(temp.states)
    private$fit$states$sd$posterior = as.data.frame(temp.sd)
    private$fit$states$cov$posterior = rep$pPost
    names(private$fit$states$cov$posterior) = paste("t = ",private$data$t,sep="")
    
    # Residual
    # rowNAs = as.matrix(!is.na(do.call(cbind, private$data[private$obs.names]))[-1,])
    rowNAs = as.matrix(!is.na(private$data[private$obs.names])[-1,])
    sumrowNAs = rowSums(rowNAs)
    
    innovation = rep$Innovation
    innovation.cov = rep$InnovationCovariance
    innovation[[1]] = NULL
    innovation.cov[[1]] = NULL
    
    temp.res = matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
    temp.var =  matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
    
    # do.call(rbind, lapply(rep$Innovation, "length<-", private$m))
    for (i in seq_along(private$data$t[-1])) {
      if (sumrowNAs[i] > 0) {
        temp.res[i,rowNAs[i,]] = innovation[[i]]
        temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
      }
    }
    temp.res = cbind(private$data$t[-1], temp.res)
    temp.sd = cbind(private$data$t[-1], sqrt(temp.var))
    
    names(innovation.cov) = paste("t = ",private$data$t[-1],sep="")
    
    # should we remove the empty matrices?
    # innovation.cov = innovation.cov[sumrowNAs!=0]
    
    colnames(temp.res) = c("t",private$obs.names)
    colnames(temp.sd) = c("t",private$obs.names)
    private$fit$residuals$mean = as.data.frame(temp.res)
    private$fit$residuals$sd = as.data.frame(temp.sd)
    private$fit$residuals$normalized = as.data.frame(temp.res)
    private$fit$residuals$normalized[,-1] = private$fit$residuals$normalized[,-1]/temp.sd[,-1]
    private$fit$residuals$cov = innovation.cov
    
    
    # Observations
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
      as.list(private$fit$data)
    )
    
    listofvariables.posterior = c(
      # states
      as.list(private$fit$states$mean$posterior[-1]),
      # estimated free parameters 
      as.list(private$fit$par.fixed),
      # fixed parameters
      lapply(private$fixed.pars, function(x) x$initial),
      # inputs
      as.list(private$fit$data)
    )
    obs.df.prior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.prior)})
    )
    obs.df.posterior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.posterior)})
    )
    private$fit$observations$mean$prior = data.frame(t=private$data$t, obs.df.prior)
    private$fit$observations$mean$posterior = data.frame(t=private$data$t, obs.df.posterior)
    
    # t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
    
  }
  
  if (any(private$method == c("ekf_rtmb"))) {
    
    
    ################################################
    # BASICS
    ################################################
    
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
    
    # hessian
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
    private$fit$par.fixed = private$opt$par
    
    # parameter std. error and full covariance by hessian inversion
    if(!is.null(private$fit$nll.hessian)){
      
      # Step 1 - invert full hessian
      temp.hessian = private$fit$nll.hessian
      covariance = try(solve(temp.hessian), silent=T)
      
      private$fit$cov.fixed = covariance
      private$fit$sd.fixed = try_withWarningRecovery(sqrt(diag(covariance)))
      
      if(inherits(private$fit$sd.fixed,"try-error")){
        private$fit$sd.fixed = rep(NA,length(private$fit$par.fixed))
      }
      if(inherits(private$fit$cov.fixed,"try-error")){
        private$fit$cov.fixed = matrix(NA,nrow=length(private$fit$par.fixed),ncol=length(private$fit$par.fixed))
      }
      
      # Options 1 - If the above fails, remove all row/cols where the diagonal
      # elements in small than min.diag
      min.diag = 1e-8
      keep.ids = !(diag(temp.hessian) < min.diag)
      if(inherits(covariance,"try-error") && any(keep.ids)){
        
        covariance = temp.hessian[keep.ids, keep.ids]
        covariance = try(solve(covariance, silent=T))
        
        sd.fixed = rep(NA,length(private$fit$par.fixed))
        sd.fixed[keep.ids] = try_withWarningRecovery(sqrt(diag(covariance)))
        private$fit$sd.fixed = sd.fixed
        
        cov.fixed = temp.hessian * NA
        cov.fixed[keep.ids, keep.ids] = covariance
        private$fit$cov.fixed = cov.fixed
      }
      
      # Option 2 - Recursive remove the smallest parameter
      # ids = sort(diag(temp.hessian), index.return=T)$ix
      # i = 0
      # while(inherits(hess,"try-error") && i < length(private$fit$par.fixed)){
      #   i = i + 1
      #   covariance = try(solve(temp.hessian[-ids[1:i],-ids[1:i]]), silent=T)
      # }
      
    }
    
    ################################################
    # STATES, RESIDUALS, OBSERVATIONS ETC.
    ################################################
    
    # Extract reported items from nll
    rep = private$nll$report()
    
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
    
    # Residual
    # rowNAs = as.matrix(!is.na(do.call(cbind, private$data[private$obs.names]))[-1,])
    rowNAs = as.matrix(!is.na(private$data[private$obs.names])[-1,])
    sumrowNAs = rowSums(rowNAs)
    
    innovation = rep$Innovation
    innovation.cov = rep$InnovationCovariance
    innovation[[1]] = NULL
    innovation.cov[[1]] = NULL
    
    temp.res = matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
    temp.var =  matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
    
    # do.call(rbind, lapply(rep$Innovation, "length<-", private$m))
    for (i in seq_along(private$data$t[-1])) {
      if (sumrowNAs[i] > 0) {
        temp.res[i,rowNAs[i,]] = innovation[[i]]
        temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
      }
    }
    temp.res = cbind(private$data$t[-1], temp.res)
    temp.sd = cbind(private$data$t[-1], sqrt(temp.var))
    
    names(innovation.cov) = paste("t = ",private$data$t[-1],sep="")
    
    # should we remove the empty matrices?
    # innovation.cov = innovation.cov[sumrowNAs!=0]
    
    colnames(temp.res) = c("t",private$obs.names)
    colnames(temp.sd) = c("t",private$obs.names)
    private$fit$residuals$mean = as.data.frame(temp.res)
    private$fit$residuals$sd = as.data.frame(temp.sd)
    private$fit$residuals$normalized = as.data.frame(temp.res)
    private$fit$residuals$normalized[,-1] = private$fit$residuals$normalized[,-1]/temp.sd[,-1]
    private$fit$residuals$cov = innovation.cov
    
    
    # Observations
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
      as.list(private$fit$data)
    )
    
    listofvariables.posterior = c(
      # states
      as.list(private$fit$states$mean$posterior[-1]),
      # estimated free parameters 
      as.list(private$fit$par.fixed),
      # fixed parameters
      lapply(private$fixed.pars, function(x) x$initial),
      # inputs
      as.list(private$fit$data)
    )
    obs.df.prior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.prior)})
    )
    obs.df.posterior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.posterior)})
    )
    private$fit$observations$mean$prior = data.frame(t=private$data$t, obs.df.prior)
    private$fit$observations$mean$posterior = data.frame(t=private$data$t, obs.df.posterior)
    
    # t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
    
  }
  
  ################################################
  # FOR TMB
  ################################################
  
  if (private$method == "laplace") {
    
    n = private$number.of.states
    
    # Objective and Gradient
    private$fit$nll = private$opt$objective
    private$fit$nll.gradient = private$sdr$gradient.fixed
    names(private$fit$nll.gradient) = names(private$free.pars)
    
    # Parameter Estimate
    private$fit$par.fixed = private$opt$par
    private$fit$sd.fixed = diag(private$sdr$cov.fixed)
    
    # Parameter Covariance
    private$fit$cov.fixed = private$sdr$cov.fixed
    
    # Posterior States (Smoothed)
    temp = cbind(matrix(private$sdr$par.random, ncol=n),
                 matrix(sqrt(private$sdr$diag.cov.random), ncol=n))[private$ode.timesteps.cumsum+1, ]
    temp.states = cbind(private$data$t, matrix(temp[,1:n],nrow=length(private$data$t)))
    temp.sd = cbind(private$data$t, matrix(temp[,(n+1):(2*n)],nrow=length(private$data$t)))
    #
    private$fit$states$mean$smoothed = as.data.frame(temp.states)
    private$fit$states$sd$smoothed = as.data.frame(temp.sd)
    colnames(private$fit$states$sd$smoothed) = c("t",private$state.names)
    colnames(private$fit$states$mean$smoothed) = c("t",private$state.names)
    
    # Residuals
    rowNAs = as.matrix(!is.na(do.call(cbind,private$data[private$obs.names]))[-1,])
    sumrowNAs = rowSums(rowNAs)
    
    # compute one-step residuals
    if(calculate.laplace.onestep.residuals){
      message("Calculating one-step ahead residuls...")
      private$fit$residuals = RTMB::oneStepPredict(private$nll,
                                 observation.name="obsMat",
                                 method="oneStepGaussian",
                                 trace=TRUE)
    }
    
    # t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
    
  }
  
  # Set S3 class
  class(private$fit) = "sdeTMB.fit"
  
  return(invisible(self))
}

construct_predict_rcpp_dataframe = function(pars, predict.list, data, return.covariance, return.k.ahead, self, private){
  
  # Simlify variable names
  n = private$number.of.states
  n.ahead = private$n.ahead
  state.names = private$state.names
  last.pred.index = private$last.pred.index
  
  # Create return data.frame
  df.out = data.frame(matrix(nrow=last.pred.index*(n.ahead+1), ncol=5+n+n^2))
  disp_names = sprintf(rep("cor.%s.%s",n^2),rep(state.names, each=n),rep(state.names,n))
  disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  var_bool = !stringr::str_detect(disp_names,"cor")
  if(return.covariance){
    disp_names = sprintf(rep("cov.%s.%s",n^2),rep(state.names,each=n),rep(state.names,n))
    disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  }
  names(df.out) = c("i.","j.","t.i","t.j","k.ahead",state.names,disp_names)
  
  # Fill out data.frame
  ran = 0:(last.pred.index-1)
  df.out["i."] = rep(ran,each=n.ahead+1)
  df.out["j."] = df.out["i."] + rep(0:n.ahead,last.pred.index)
  df.out["t.i"] = rep(data$t[ran+1],each=n.ahead+1)
  df.out["t.j"] = data$t[df.out[,"i."]+1+rep(0:n.ahead,last.pred.index)]
  df.out["k.ahead"] = rep(0:n.ahead,last.pred.index)
  df.obs = df.out[c("i.","j.","t.i","t.j","k.ahead")]
  
  ##### STATES PREDICTIONS ######
  df.out[,state.names] = do.call(rbind,lapply(predict.list$Xpred, function(cur.list) do.call(rbind, cur.list)))
  if(return.covariance){ #covariance
    df.out[,disp_names] = do.call(rbind, lapply(predict.list$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(x))) } ))
  } else { #correlation
    df.out[,disp_names] = do.call(rbind, lapply(predict.list$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(cov2cor(x)))) } )) 
    diag.ids = seq(from=1,to=n^2,by=n+1)
    df.out[,disp_names[diag.ids]] = do.call(rbind, lapply(predict.list$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(diag(x)))) } ))
  }
  
  ##### OBSERVATION PREDICTIONS #####
  # calculate observations at every time-step in predict
  inputs.df = private$data[df.out[,"j."]+1,private$input.names]
  
  named.pars.list = as.list(pars)
  names(named.pars.list) = names(private$free.pars)
  # create environment
  env.list = c(
    # states
    as.list(df.out[state.names]),
    # inputs
    as.list(inputs.df),
    # free parameters
    named.pars.list,
    # fixed parameters
    lapply(private$fixed.pars, function(x) x$initial)
  )
  
  # calculate observations
  obs.df.predict = as.data.frame(
    lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = env.list)})
  )
  # names(obs.df.predict) = paste(private$obs.names,".predict",sep="")
  names(obs.df.predict) = paste(private$obs.names)
  # df.out = cbind(df.out, obs.df.predict)
  
  # add data observation to output data.frame 
  obs.df.data = private$data[df.out[,"j."]+1, private$obs.names, drop=F]
  names(obs.df.data) = paste(private$obs.names,".data",sep="")
  # df.out = cbind(df.out, obs.df)
  
  df.obs = cbind(df.obs, obs.df.predict, obs.df.data)
  
  # return only specific n.ahead
  df.out = df.out[df.out[,"k.ahead"] %in% return.k.ahead,]
  df.obs = df.obs[df.obs[,"k.ahead"] %in% return.k.ahead,]
  
  list.out = list(states = df.out, observations = df.obs)
  class(list.out) = c(class(list.out), "sdeTMB.pred")
  
  return(list.out)
}

construct_simulate_rcpp_dataframe = function(pars, predict.list, data, return.k.ahead, n.sims, self, private){
  
  
  list.out = vector("list",length=private$number.of.states)
  names(list.out) = private$state.names
  
  setRownames = function(obj, nm){rownames(obj) = nm; return(obj)}
  setColnames = function(obj, nm){colnames(obj) = nm; return(obj)}
  

  # for(i in seq_along(list.out)){
  #   list.out[[i]] = stats::setNames(
  #     lapply(predict.list, function(ls.outer){
  #       # setRownames(
  #       t(do.call(cbind, lapply(ls.outer, function(ls.inner) ls.inner[,i])))
  #       # paste0("k.ahead", 0:private$n.ahead)
  #       # )
  #     }),
  #     paste0("t", head(data$t, private$last.pred.index))
  #   )
  # }
  
  # Compute the prediction times for each horizon
  ran = 0:(private$last.pred.index-1)
  t.j = data$t[rep(ran,each=private$n.ahead+1)+1+rep(0:private$n.ahead,private$last.pred.index)]
  t.j.splitlist = split(t.j, ceiling(seq_along(t.j)/(private$n.ahead+1)))
  list.of.time.vectors = lapply(t.j.splitlist, function(x) data.frame(t.j=x))
  # names(list.of.time.vectors) = names(list.out[[1]])
  # list.out2 = c(list.out, list(prediction_times = list.of.time.vectors))
  
  for(i in seq_along(list.out)){
    list.out[[i]] = stats::setNames(
      lapply(predict.list, function(ls.outer){
        # setRownames(
        t(do.call(cbind, lapply(ls.outer, function(ls.inner) ls.inner[,i])))
        # paste0("k.ahead", 0:private$n.ahead)
        # )
      }),
      # paste0("t", head(data$t, private$last.pred.index))
      paste0("i",ran)
    )
  }
  
  for(i in seq_along(list.out)){
    for(j in seq_along(list.out[[i]])){
      list.out[[i]][[j]] = data.frame(i=j-1, 
                                      j=(j-1):(j+private$n.ahead-1), 
                                      t.i=rep(data$t[i],private$n.ahead+1),
                                      t.j=list.of.time.vectors[[j]][,"t.j"], 
                                      k.ahead = 0:private$n.ahead,
                                      list.out[[i]][[j]]
                                      )
      nams = paste0(private$state.names,1:n.sims)
      names(list.out[[i]][[j]]) = c("i","j","t.i","t.j","k.ahead",nams)
    }
  }
  
  list.out2 = list(
    states = list.out,
    observations = list()
  )
  
  return(list.out2)
}


