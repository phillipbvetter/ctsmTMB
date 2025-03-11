#######################################################
# CONSTRUCT KALMAN MAKEADFUN WITH RTMB
#######################################################

construct_rtmb_ekf_makeADFun = function(self, private)
{
  
  # Data ----------------------------------------
  
  # methods and purpose
  ode_solver = private$ode.solver
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  # loss function
  # loss_function = private$loss$loss
  # loss_threshold_value = private$loss$c
  # tukey_loss_parameters = private$tukey.pars
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # # MAP Estimation?
  # MAP_bool = 0L
  # if (!is.null(private$map)) {
  #   bool = self$getParameters()[,"type"] == "free"
  #   MAP_bool = 1L
  #   map_mean__ = private$map$mean[bool]
  #   map_cov__ = private$map$cov[bool,bool]
  #   map_ints__ = as.numeric(bool)
  #   sum_map_ints__ = sum(as.numeric(bool))
  # }
  
  # parameters ----------------------------------------
  parameters = lapply(private$parameters, function(x) x[["initial"]])
  
  # adjoints ----------------------------------------
  logdet <- RTMB::ADjoint(
    function(x) {
      dim(x) <- rep(sqrt(length(x)), 2)
      log(det(x))
    },
    function(x, y, dy) {
      dim(x) <- rep(sqrt(length(x)), 2)
      t(RTMB::solve(x)) * dy
    },
    name = "logdet")
  
  kron_left <- RTMB::ADjoint(
    function(x) {
      dim(x) <- rep(sqrt(length(x)), 2)
      i <- diag(sqrt(length(x)))
      kronecker(x,i)
    },
    function(x, y, dy) {
      n <- sqrt(length(x))
      out <- matrix(0,nrow=n,ncol=n)
      for(i in 1:n){
        for(j in 1:n){
          id.seq <- 1+(j-1)*n + (n^2+1) * (1:n-1)
          id.seq <- id.seq + (i-1) * n^3
          out[[i,j]] <- sum(dy[id.seq])
        }
      }
      return(out)
    },
    name = "kron_left")
  
  kron_right <- RTMB::ADjoint(
    function(x) {
      dim(x) <- rep(sqrt(length(x)), 2)
      i <- diag(sqrt(length(x)))
      kronecker(i,x)
    },
    function(x, y, dy) {
      n <- sqrt(length(x))
      out <- matrix(0,nrow=n,ncol=n)
      for(i in 1:n){
        for(j in 1:n){
          id.seq <- j + (n^3+n) * (1:n-1)
          id.seq <- id.seq + (i-1) * n^2
          out[[i,j]] <- sum(dy[id.seq])
        }
      }
      return(out)
    },
    name = "kron_right")
  
  # 1-step covariance ODE ----------------------------------------
  cov_ode_1step = function(covMat, stateVec, parVec, inputVec){
    G <- g__(stateVec, parVec, inputVec)
    AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  
  # forward euler ----------------------------------------
  if(ode_solver==1){
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      X1 = stateVec + f__(stateVec, parVec, inputVec) * dt
      P1 = covMat + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
      
      return(list(X1,P1))
    }
  }
  
  # rk4 ----------------------------------------
  if(ode_solver==2){
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      # Initials
      X0 = stateVec
      P0 = covMat
      
      # Classical 4th Order Runge-Kutta Method
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
      
      return(list(X1,P1))
    }
  }
  
  # user-defined functions ---------------------------
  for(i in seq_along(private$rtmb.function.strings)){
    eval(parse(text=private$rtmb.function.strings[[i]]))
  }
  
  # new functions ----------------------------------------
  # error function ----------------------------------------
  erf = function(x){
    y <- sqrt(2) * x
    2*RTMB::pnorm(y)-1
  }
  
  # global values in likelihood function ------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  estimate.initial <- private$estimate.initial
  
  ekf.nll = function(p){
    
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    
    ####### Parameters into vector #######
    parVec <- do.call(c, p[1:n.pars])
    
    ####### STORAGE #######
    xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
    
    ####### Neg. LogLikelihood #######
    nll <- 0
    
    ####### Pre-Allocated Object #######
    I0 <- RTMB::diag(n.states)
    E0 <- RTMB::diag(n.obs)
    
    ####### INITIAL STATE / COVARIANCE #######
    # The state/covariance is either given by user or obtained from solving the
    # stationary mean, and then solving for the covariance.
    # In principle these are coupled equations, but believe that root-finding
    # both simultaneously can lead to likelihood blow-up.
    inputVec = inputMat[1,]
    if(estimate.initial){
      # 1. Root-find stationary mean
      .F <- RTMB::MakeTape(function(x) sum(f__(x, parVec, inputVec)^2), numeric(n.states))
      stateVec <- .F$newton(1:n.states)(numeric(0))
      # 2. Use stationary mean to solve lyapunov eq. for associated covariance
      A <- dfdx__(stateVec, parVec, inputVec)
      G <- g__(stateVec, parVec, inputVec)
      Q <- G %*% t(G)
      P <- kron_left(A) + kron_right(A)
      X <- -RTMB::solve(P, as.numeric(Q))
      covMat <- RTMB::matrix(X, nrow=n.states)
    }
    xPrior[[1]] <- stateVec
    pPrior[[1]] <- covMat
    
    ######## (PRE) DATA UPDATE ########
    # This is done to include the first measurements in the provided data
    # We update the state and covariance based on the "new" measurement
    obsVec = obsMat[1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      y = obsVec[obsVec_bool]
      E = E0[obsVec_bool,, drop=FALSE]
      C = E %*% dhdx__(stateVec, parVec, inputVec)
      e = y - E %*% h__(stateVec, parVec, inputVec)
      V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
      R = C %*% covMat %*% t(C) + V
      K = covMat %*% t(C) %*% RTMB::solve(R)
      # Likelihood Contribution
      nll = nll - RTMB::dmvnorm(e, Sigma=R, log=TRUE)
      # Update State/Cov
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
      # Store innovation and covariance
      Innovation[[1]] = e
      InnovationCovariance[[1]] = R
    }
    xPost[[1]] <- stateVec
    pPost[[1]] <- covMat
    
    ###### MAIN LOOP START #######
    for(i in 1:(nrow(obsMat)-1)){
      
      inputVec = inputMat[i,]
      dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
      
      ###### TIME UPDATE #######
      # We solve the first two moments forward in time
      for(j in 1:ode_timesteps[i]){
        sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i])
        stateVec = sol[[1]]
        covMat = sol[[2]]
        inputVec = inputVec + dinputVec
      }
      xPrior[[i+1]] = stateVec
      pPrior[[i+1]] = covMat
      
      ######## DATA UPDATE ########
      # We update the state and covariance based on the "new" measurement
      inputVec = inputMat[i+1,]
      obsVec = obsMat[i+1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        y = obsVec[obsVec_bool]
        E = E0[obsVec_bool,, drop=FALSE] #permutation matrix with rows removed
        C = E %*% dhdx__(stateVec, parVec, inputVec)
        e = y - E %*% h__(stateVec, parVec, inputVec)
        V = E %*% hvar__matrix(stateVec, parVec, inputVec) %*% t(E)
        R = C %*% covMat %*% t(C) + V
        K = covMat %*% t(C) %*% RTMB::solve(R)
        # Likelihood Contribution
        nll = nll - RTMB::dmvnorm(e, Sigma=R, log=TRUE)
        # Update State/Cov
        stateVec = stateVec + K %*% e
        covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
        # Store innovation and covariance
        Innovation[[i+1]] <- e
        InnovationCovariance[[i+1]] <- R
      }
      xPost[[i+1]] = stateVec
      pPost[[i+1]] = covMat
    }
    
    ###### MAIN LOOP END #######
    
    # ###### MAXIMUM A POSTERIOR #######
    
    ##### REPORT #######
    RTMB::REPORT(Innovation)
    RTMB::REPORT(InnovationCovariance)
    RTMB::REPORT(xPrior)
    RTMB::REPORT(xPost)
    RTMB::REPORT(pPrior)
    RTMB::REPORT(pPost)
    
    # ###### RETURN #######
    return(nll)
  }
  
  # Create AD-likelihood function ---------------------------------------
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
#######################################################
# RETURN FIT FOR EKF RTMB
#######################################################
#######################################################
create_fit_ekf_rtmb <- function(self, private){
  
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
      covariance = try(solve(covariance), silent=T)
      
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
  
  # States -----------------------------------
  
  # Extract reported items from nll
  comptime <- system.time(rep <- private$nll$report())
  comptime = format(round(as.numeric(comptime["elapsed"])*1e2)/1e2,digits=5,scientific=F)
  if(!private$silent) message("...took ", comptime, " seconds")
  
  # Fix annoying errors when reported objects are advectors
  # which occurs when estimating initial state / cov
  ad2double <- function(x) RTMB:::getValues(RTMB::AD(x))
  bool <- any(sapply(rep$xPrior,function(x) inherits(x,"advector")))
  if(bool){
    rep$xPrior <- lapply(rep$xPrior, ad2double)
    rep$pPrior <- lapply(rep$pPrior, ad2double)
    rep$xPost <- lapply(rep$xPost, ad2double)
    rep$pPost <- lapply(rep$pPost, ad2double)
  }
  
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
  temp.sd = cbind(private$data$t[-1], try_withWarningRecovery(sqrt(temp.var)))
  
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
  
  # Evaluate prior and posterior variance
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
# RETURN FIT FOR EKF CPP
#######################################################
#######################################################
calculate_fit_statistics_ekfukf_cpp <- function(self, private){
  
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
  
  # objective gradient
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
  if (is.null(private$fit$nll.hessian)){
    private$fit$cov.fixed <- NaN * diag(length(private$fit$par.fixed))
    private$fit$sd.fixed <- rep(NaN,length(private$fit$par.fixed))
  } else {
    
    ####### OPTION 0 #######
    # This option tries to invert the full hessain
    full.hess = private$fit$nll.hessian
    covar <- try_withWarningRecovery(solve(full.hess))
    if (inherits(covar,"try-error")) {
      private$fit$cov.fixed = NaN * full.hess
      private$fit$sd.fixed = rep(NaN,length(private$fit$par.fixed))
    } else {
      private$fit$cov.fixed <- covar
      private$fit$sd.fixed <- try_withWarningRecovery(sqrt(diag(private$fit$cov.fixed)))
    }
    
    ####### OPTION 1 #######
    # This option filters out all row/cols of the hessian where the diagonal
    # element is smaller than some set threshold
    min.diag = 1e-8
    keep.ids = !(diag(full.hess) < min.diag)
    if(inherits(covar,"try-error") && any(keep.ids)){
      reduced.hess <- full.hess[keep.ids, keep.ids]
      cov <- try_withWarningRecovery(solve(reduced.hess))
      private$fit$cov.fixed <- NaN * full.hess
      private$fit$sd.fixed <- rep(NaN,length(private$fit$par.fixed))
      if (!inherits(cov,"try-error")) {
        private$fit$cov.fixed[keep.ids, keep.ids] <- cov
        private$fit$sd.fixed[keep.ids] = try_withWarningRecovery(sqrt(diag(cov)))
      }
    }
    
    ####### OPTION 2 #######
    # This option tries to invert the hessian by recursively removing row/cols
    # of the hessian where the diagonal are smallest
    if (inherits(covar,"try-error")) {
      private$fit$cov.fixed <- NaN * full.hess
      private$fit$sd.fixed <- rep(NaN,length(private$fit$par.fixed))
      ids = sort(diag(full.hess), index.return=T)$ix
      for(i in seq_along(ids)){
        id <- ids[1:i]
        hess <- full.hess[-id,-id]
        cov <- try_withWarningRecovery(solve(hess))
        if (!inherits(cov,"try-error")){
          private$fit$cov.fixed[-id,-id] <- cov
          private$fit$sd.fixed[-id] <- try_withWarningRecovery(sqrt(diag(cov)))
          break
        }
      }
    }
  }
  
  # States -----------------------------------
  
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
  
  # Residuals -----------------------------------
  
  rowNAs = as.matrix(!is.na(private$data[private$obs.names]))
  sumrowNAs = rowSums(rowNAs)
  
  innovation = rep$Innovation
  innovation.cov = rep$InnovationCovariance
  
  temp.res = matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
  temp.var =  matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
  
  for (i in seq_along(private$data$t)) {
    if (sumrowNAs[i] > 0) {
      temp.res[i,rowNAs[i,]] = innovation[[i]]
      temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
    }
  }
  temp.res = cbind(private$data$t, temp.res)
  temp.sd = cbind(private$data$t, try_withWarningRecovery(sqrt(temp.var)))
  
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
  
  obs.df.prior = as.data.frame(
    lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.prior)})
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
  
  # Evaluate prior and posterior variance
  list.of.parameters <- c(
    as.list(private$fit$par.fixed),
    lapply(private$fixed.pars, function(x) x$initial)
  )
  
  # prior
  obsvar.prior <- list()
  for(i in seq_along(private$data$t)){
    list.of.states <- as.list(private$fit$states$mean$prior[i,-1,drop=F])
    state.cov <- list(xCov = private$fit$states$cov$prior[[i]])
    list.of.inputs <- as.list(private$data[i,-1,drop=F])
    obsvar.prior[[i]] <- eval(expr = yCov,
                              envir = c(list.of.parameters, 
                                        list.of.states,
                                        state.cov,
                                        list.of.inputs
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
    list.of.states <- as.list(private$fit$states$mean$posterior[i,-1,drop=F])
    state.cov <- list(xCov = private$fit$states$cov$posterior[[i]])
    list.of.inputs <- as.list(private$data[i,-1,drop=F])
    obsvar.post[[i]] <- eval(expr = yCov,
                             envir = c(list.of.parameters, 
                                       list.of.states,
                                       state.cov,
                                       list.of.inputs
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
}
