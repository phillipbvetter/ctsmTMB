#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

construct_makeADFun = function(self, private){
  
  # TMB::openmp(max=TRUE, autopar=TRUE, DLL=private$modelname.with.method)
  
  if(any(private$method == c("ekf","ukf"))){
    comptime <- system.time(construct_kalman_cpp_makeADFun(self, private))
  }
  
  if(private$method == "ekf_rtmb"){
    comptime <- system.time(construct_rtmb_ekf_makeADFun(self, private))
  }
  
  if(private$method=="laplace"){
    comptime <- system.time(construct_rtmb_laplace_makeADFun(self, private))
  }
  
  comptime = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
  if(!private$silent) message("...took: ", comptime, " seconds.")
  
  return(invisible(self))
}

#######################################################
# CONSTRUCT KALMAN MAKEADFUN
#######################################################
construct_kalman_cpp_makeADFun = function(self, private){
  
  ################################################
  # Data
  ################################################
  
  # add mandatory entries to data
  tmb.data = list(
    
    # methods and purpose
    estimation_method = switch(private$method, ekf = 1, ukf = 2),
    ode_solver = private$ode.solver,
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # time-steps
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    ode_cumsum_timesteps = private$ode.timesteps.cumsum,
    
    # loss function
    loss_function = private$loss$loss,
    loss_threshold_value = private$loss$c,
    tukey_loss_parameters = private$tukey.pars,
    
    # system size
    number_of_state_eqs = private$number.of.states,
    number_of_obs_eqs = private$number.of.observations,
    number_of_diffusions = private$number.of.diffusions,
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names])
  )
  
  # unscented parameters
  ukf_hyperpars_list = list()
  if(private$method=="ukf")
  {
    ukf_hyperpars_list = list(
      ukf_alpha = private$ukf_alpha,
      ukf_beta = private$ukf_beta,
      ukf_kappa = private$ukf_kappa
    )
  }
  
  # MAP Estimation?
  tmb.map.data = list(
    MAP_bool = 0L
  )
  if (!is.null(private$map)) {
    bool = self$getParameters()[,"type"] == "free"
    tmb.map.data = list(
      MAP_bool = 1L,
      map_mean__ = private$map$mean[bool],
      map_cov__ = private$map$cov[bool,bool],
      map_ints__ = as.numeric(bool),
      sum_map_ints__ = sum(as.numeric(bool))
    )
  }
  
  # construct final data list
  data = c(tmb.data, private$iobs, tmb.map.data, ukf_hyperpars_list)
  
  ################################################
  # Parameters
  ################################################
  
  parameters = lapply(private$parameters, function(x) x[["initial"]]) # Initial parameter values
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = TMB::MakeADFun(data = data,
                       parameters = parameters,
                       map = lapply(private$fixed.pars, function(x) x$factor),
                       DLL = private$modelname.with.method,
                       silent = TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}

#######################################################
# CONSTRUCT KALMAN MAKEADFUN WITH RTMB
#######################################################

construct_rtmb_ekf_makeADFun = function(self, private)
{
  
  ################################################
  # Data
  ################################################
  
  # methods and purpose
  ode_solver = private$ode.solver
  
  # initial
  stateVec = cbind(private$initial.state$x0)
  covMat = private$initial.state$p0
  
  # loss function
  loss_function = private$loss$loss
  loss_threshold_value = private$loss$c
  tukey_loss_parameters = private$tukey.pars
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # system size
  number_of_state_eqs = private$number.of.states
  number_of_obs_eqs = private$number.of.observations
  number_of_diffusions = private$number.of.diffusions
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # MAP Estimation?
  MAP_bool = 0L
  if (!is.null(private$map)) {
    bool = self$getParameters()[,"type"] == "free"
    MAP_bool = 1L
    map_mean__ = private$map$mean[bool]
    map_cov__ = private$map$cov[bool,bool]
    map_ints__ = as.numeric(bool)
    sum_map_ints__ = sum(as.numeric(bool))
  }
  
  ################################################
  # Initial Parameters
  ################################################
  
  parameters = lapply(private$parameters, function(x) x[["initial"]])
  
  ################################################
  # Define functions
  ################################################
  
  # ODE Solver
  ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt, ode_solver){
    
    # Initials
    X0 = stateVec
    P0 = covMat
    
    # Explicit Forward Euler
    if(ode_solver==1){
      X1 = X0 + f__(stateVec, parVec, inputVec) * dt
      P1 = P0 + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
    }
    
    # Classical 4th Order Runge-Kutta Method
    if(ode_solver==2){
      
      # 1. Approx Slope at Initial Point
      k1 = f__(stateVec, parVec, inputVec)
      c1 = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 2. First Approx Slope at Midpoint
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + 0.5 * dt * k1
      covMat   = P0 + 0.5 * dt * c1
      k2       = f__(stateVec, parVec, inputVec)
      c2       = cov_ode_1step(covMat, stateVec, parVec, inputVec)   
      
      # 3. Second Approx Slope at Midpoint
      stateVec = X0 + 0.5 * dt * k2
      covMat   = P0 + 0.5 * dt * c2
      k3       = f__(stateVec, parVec, inputVec)
      c3       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 4. Approx Slope at End Point
      inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + dt * k3
      covMat   = P0 + dt * c3
      k4       = f__(stateVec, parVec, inputVec)
      c4       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # ODE UPDATE
      X1 = X0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0 * dt
      P1 = P0 + (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0 * dt
    }
    return(invisible(list(X1,P1)))
  }
  
  # Covariance ODE 1-Step
  cov_ode_1step = function(covMat, stateVec, parVec, inputVec){
    A = dfdx__(stateVec, parVec, inputVec)
    G = g__(stateVec, parVec, inputVec)
    AcovMat = A %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  for(i in seq_along(private$rtmb.function.strings)){
    eval(parse(text=private$rtmb.function.strings[[i]]))
  }
  eval(parse(text=private$rtmb.nll.strings$ekf))
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = RTMB::MakeADFun(func = ekf.nll,
                        parameters=parameters,
                        map = lapply(private$fixed.pars, function(x) x$factor),
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}

#######################################################
# CONSTRUCT LAPLACE MAKEADFUN WITH RTMB
#######################################################

construct_rtmb_laplace_makeADFun = function(self, private)
{
  
  ################################################
  # Data
  ################################################
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # system size
  number_of_state_eqs = private$number.of.states
  number_of_obs_eqs = private$number.of.observations
  number_of_diffusions = private$number.of.diffusions
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # indices in state parameter vectors corresponding to indices in observations / inputs
  # add 1 because too lazy to change private$iobs from 0-index to 1-indexed.
  iobs = lapply(private$iobs,function(x) x+1)
  
  ################################################
  # Parameters
  ################################################
  
  parameters = c(
    lapply(private$parameters, function(x) x[["initial"]]),
    private$tmb.initial.state.for.parameters
  )
  
  ################################################
  # Functions
  ################################################
  
  for(i in seq_along(private$rtmb.function.strings)){
    eval(parse(text=private$rtmb.function.strings[[i]]))
  }
  eval(parse(text=private$rtmb.nll.strings$laplace))
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = RTMB::MakeADFun(func = laplace.nll, 
                        parameters=parameters, 
                        random=private$state.names,
                        map = lapply(private$fixed.pars, function(x) x$factor),
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}


#######################################################
# OPTIMISE AD FUN
#######################################################

perform_estimation = function(self, private) {
  
  # Parameter Bounds
  lower.parameter.bound = unlist(lapply(private$free.pars, function(par) par$lower))
  upper.parameter.bound = unlist(lapply(private$free.pars, function(par) par$upper))
  if(private$unconstrained.optim){
    lower.parameter.bound = -Inf
    upper.parameter.bound = Inf
  }
  
  # IF METHOD IS KALMAN FILTER
  if (any(private$method==c("ekf","ukf","ekf_rtmb"))) {
    
    # use function, gradient and hessian
    if (private$use.hessian) {
      comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                           objective = private$nll$fn,
                                                                           gradient = private$nll$gr,
                                                                           hessian = private$nll$he,
                                                                           lower = lower.parameter.bound,
                                                                           upper = upper.parameter.bound,
                                                                           control=private$control.nlminb))
      )
      # or just function and gradient
    } else {
      comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                           objective = private$nll$fn,
                                                                           gradient = private$nll$gr,
                                                                           lower = lower.parameter.bound,
                                                                           upper = upper.parameter.bound,
                                                                           control=private$control.nlminb))
      )
    }
    
  }
  
  # IF METHOD IS TMB
  if (private$method =="laplace") {
    comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                         objective = private$nll$fn,
                                                                         gradient = private$nll$gr,
                                                                         lower = lower.parameter.bound,
                                                                         upper = upper.parameter.bound,
                                                                         control=private$control.nlminb))
    )
  }
  
  # DID THE OPTIMIZATION FAIL?
  if (inherits(opt,"try-error")) {
    
    message("The optimisation failed due to the following error: \n\n\t",opt)
    
    if(stringr::str_detect(opt,"NA/NaN")){
      
      message("You should consider the following to circumvent the error:
              1. Reduce the ODE step-size via argument 'ode.timestep'.
              2. Run the optimization with / without the hessian via argument 'use.hessian'.
              3. Change / explore parameter initial values.
              4. Extract the function handlers with the 'construct_nll' method and investigate outputs, or try other optimizers.
              5. Change the optimization tolerances for 'nlminb' with the 'control' argument.")
      
    }
    
    private$opt = NULL
    
    return(invisible(self))
  }
  
  # store optimization object
  private$opt = opt
  
  # extract maxmimum gradient component, and format computation time to 5 digits
  outer_mgc = max(abs(private$nll$gr(opt$par)))
  comp.time = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
  
  # print convergence and timing result
  if(!private$silent){
    if(outer_mgc > 1){
      message("BEWARE: THE MAXIMUM GRADIENT COMPONENT APPEARS TO BE LARGE ( > 1 ) - THE FOUND OPTIMUM MIGHT BE INVALID.")
    }
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
  if (private$method=="laplace") {
    if(!private$silent) message("Calculating random effects standard deviation...")
    comptime = system.time(
      private$sdr <- RTMB::sdreport(private$nll)
    )
    comptime = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
    if(!private$silent) message("...took: ", comptime, " seconds.")
  }
  
  # return
  return(invisible(self))
}




#######################################################
# Create return fit after calling estimate 
#######################################################

create_fit = function(self, private, calculate.laplace.onestep.residuals) {
  
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
  # private$fit$.__object__ = self$clone()
  
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
  class(private$fit) = "ctsmTMB.fit"
  
  return(invisible(self))
}


