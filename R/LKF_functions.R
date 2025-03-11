#######################################################
#######################################################
# EKF RTMB-IMPLEMENTATION (FOR OPTIMIZATION)
#######################################################
#######################################################

makeADFun_lkf_rtmb = function(self, private)
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
  
  # user-defined functions ---------------------------
  for(i in seq_along(private$rtmb.function.strings.indexed)){
    eval(parse(text=private$rtmb.function.strings.indexed[[i]]))
  }
  
  # new functions ----------------------------------------
  
  #error function
  erf = function(x){
    y <- sqrt(2) * x
    2*RTMB::pnorm(y)-1
  }
  
  # global values in likelihood function ------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # AD overwrites ----------------------------------------
  f_vec <- RTMB::AD(numeric(n.obs),force=TRUE)
  dfdx_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.states, ncol=n.states),force=TRUE)
  dfdu_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.states, ncol=n.inputs+1),force=TRUE)
  g_mat <- RTMB::AD(RTMB::matrix(0,nrow=n.states, ncol=n.diffusions),force=TRUE)
  h_vec <- RTMB::AD(numeric(n.obs),force=TRUE)
  dhdx_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.obs, ncol=n.states),force=TRUE)
  hvar_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.obs, ncol=n.obs),force=TRUE)
  
  stateVec <- RTMB::AD(stateVec,force=TRUE)
  covMat <- RTMB::AD(covMat,force=TRUE)
  # obsMat <- RTMB::AD(obsMat,force=F)
  # inputMat <- RTMB::AD(inputMat,force=F)
  
  # detect time-variations etc ----------------------
  constant.time.diff <- FALSE
  if(var(diff(ode_timestep_size)) < 1e-15){
    constant.time.diff <- TRUE
    fixed.timestep.size <- ode_timestep_size[1]
  }
  
  # likelihood function --------------------------------
  lkf.nll = function(p){
    
    "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    
    ####### Parameters into vector #######
    parVec <- do.call(c, p[1:n.pars])
    
    ####### Neg. LogLikelihood #######
    nll <- 0
    
    ####### Pre-Allocated Object #######
    I0 <- RTMB::diag(n.states)
    E0 <- RTMB::diag(n.obs)
    
    ####### Compute Matrix Exponentials for 1-Step Mean and Variance #######
    # dX = A*X + B*U + G*dB
    A <- dfdx__(stateVec, parVec, inputVec) #states
    B <- dfdu__(stateVec, parVec, inputVec) #inputs
    G <- g__(stateVec, parVec, inputVec) #diffusions
    H <- dhdx__(stateVec,parVec, inputVec) #observation
    V0 <- hvar__matrix(stateVec, parVec, inputVec) #observation variance
    
    # [A B \\ 0 0]
    Phi1 <- 0.0 * diag(n.states+n.inputs+1)
    Phi1[1:n.states,1:n.states] <- A
    Phi1[1:n.states,(n.states+1):ncol(Phi1)] <- B
    # [-A GG^T \\ 0 A^t]
    Phi2 <- rbind(cbind(-A,G %*% t(G)),cbind(0*A,t(A)))
    ePhi1 <- Matrix::expm(Phi1 * fixed.timestep.size)
    ePhi2 <- Matrix::expm(Phi2 * fixed.timestep.size)
    
    # A and B (for mean calculations)
    Ahat <- ePhi1[1:n.states,1:n.states]
    Ahat_T <- t(Ahat)
    Bhat <- ePhi1[1:n.states, (n.states+1):ncol(ePhi1)]
    # V (for variance)
    Q22 <- ePhi2[(n.states+1):ncol(ePhi2), (n.states+1):ncol(ePhi2)]
    Q12 <- ePhi2[1:n.states, (n.states+1):ncol(ePhi2)]
    Vhat <- t(Q22) %*% Q12
    
    # if(constant.time.diff){
    #   .dt <- ode.timestep.size[1]
    # } else {
    # .dt <- diff(ode.timesteps)
    # RTMB::eigen(A)
    
    ####### INITIAL STATE / COVARIANCE #######
    # The state/covariance is either given by user or obtained from solving the
    # stationary mean, and then solving for the covariance.
    # In principle these are coupled equations, but believe that root-finding
    # both simultaneously can lead to likelihood blow-up.
    # inputVec = inputMat[1,]
    # if(estimate.initial){
    #   # 1. Root-find stationary mean
    #   .F <- RTMB::MakeTape(function(x) sum(f__(x, parVec, inputVec)^2), numeric(n.states))
    #   stateVec <- .F$newton(1:n.states)(numeric(0))
    #   # 2. Use stationary mean to solve lyapunov eq. for associated covariance
    #   A <- dfdx__(stateVec, parVec, inputVec)
    #   G <- g__(stateVec, parVec, inputVec)
    #   Q <- G %*% t(G)
    #   P <- kron_left(A) + kron_right(A)
    #   X <- -RTMB::solve(P, as.numeric(Q))
    #   covMat <- RTMB::matrix(X, nrow=n.states)
    # }
    
    ######## (PRE) DATA UPDATE ########
    # This is done to include the first measurements in the provided data
    # We update the state and covariance based on the "new" measurement
    inputVec = inputMat[1,]
    obsVec = obsMat[1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      y = obsVec[obsVec_bool]
      E = E0[obsVec_bool,, drop=FALSE]
      C = E %*% H
      e = y - C %*% stateVec
      V = E %*% V0 %*% t(E)
      R = C %*% covMat %*% t(C) + V
      K = covMat %*% t(C) %*% RTMB::solve(R)
      # Likelihood Contribution
      nll = nll - RTMB::dmvnorm(e, Sigma=R, log=TRUE)
      # Update State/Cov
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    }
    
    
    ###### MAIN LOOP START #######
    for(i in 1:(nrow(obsMat)-1)){
      
      ###### TIME UPDATE #######
      # augment input vector to account for constant terms in B
      inputVec = c(1,inputMat[i,])
      # Perform one-step prediction of mean and covariance
      stateVec <- Ahat %*% stateVec + Bhat %*% inputVec
      covMat <- Ahat %*% covMat %*% Ahat_T + Vhat
      
      ######## DATA UPDATE ########
      # We update the state and covariance based on the "new" measurement
      inputVec = inputMat[i+1,]
      obsVec = obsMat[i+1,]
      obsVec_bool = !is.na(obsVec)
      if(any(obsVec_bool)){
        y = obsVec[obsVec_bool]
        E = E0[obsVec_bool,, drop=FALSE]
        C = E %*% H
        e = y - C %*% stateVec
        V = E %*% V0 %*% t(E)
        R = C %*% covMat %*% t(C) + V
        K = covMat %*% t(C) %*% RTMB::solve(R)
        # Likelihood Contribution
        nll = nll - RTMB::dmvnorm(e, Sigma=R, log=TRUE)
        # Update State/Cov
        stateVec = stateVec + K %*% e
        covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
      }
    }
    ###### MAIN LOOP END #######
    
    # ###### RETURN #######
    return(nll)
  }
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = RTMB::MakeADFun(func = lkf.nll,
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
# EKF R-IMPLEMENTATION (FOR REPORTING)
#######################################################
#######################################################

lkf_r = function(parVec, self, private)
{
  
  # parameters ----------------------------------------
  if(missing(parVec)){
    parVec = sapply(private$parameters, function(x) x[["initial"]])
  }
  
  # Data ----------------------------------------
  
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
  
  # user-defined functions ---------------------------
  for(i in seq_along(private$rekf.function.strings)){
    eval(parse(text=private$rekf.function.strings[[i]]))
  }
  
  # new functions ----------------------------------------
  #error function
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
  
  # detect time-variations etc ----------------------
  constant.time.diff <- FALSE
  if(var(diff(ode.timestep.size)) < 1e-15){
    constant.time.diff <- TRUE
  }
  
  ####### STORAGE #######
  xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
  
  ####### Neg. LogLikelihood #######
  nll <- 0
  
  ####### Pre-Allocated Object #######
  I0 <- diag(n.states)
  E0 <- diag(n.obs)
  
  ####### Compute Matrix Exponentials for 1-Step Mean and Variance #######
  # dX = A*X + B*U + G*dB
  A <- dfdx__(stateVec, parVec, inputVec) #states
  B <- dfdu__(stateVec, parVec, inputVec) #inputs
  G <- g__(stateVec, parVec, inputVec) #diffusions
  H <- dhdx__(stateVec,parVec, inputVec) #observation
  V0 <- hvar__matrix(stateVec, parVec, inputVec) #observation variance
  
  # [A B \\ 0 0]
  Phi1 <- rbind(cbind(A,B),cbind(0*A, 0*B))
  # [-A GG^T \\ 0 A^t]
  Phi2 <- rbind(cbind(-A,G %*% t(G)),cbind(0*A,t(A)))
  ePhi1 <- Matrix::expm(Phi1*dt.fixed)
  ePhi2 <- Matrix::expm(Phi2*dt.fixed)
  
  Ahat <- ePhi1[1:n.states,1:n.states]
  Ahat_T <- t(Ahat)
  Bhat <- ePhi1[1:n.states, (n.states+1):ncol(ePhi1)]
  Q22 <- ePhi2[(n.states+1):ncol(ePhi1), (n.states+1):ncol(ePhi1)]
  Q12 <- ePhi2[1:n.states, (n.states+1):ncol(ePhi1)]
  Vhat <- t(Q22) %*% Q12
  
  # if(constant.time.diff){
  #   .dt <- ode.timestep.size[1]
  # } else {
  # .dt <- diff(ode.timesteps)
  # RTMB::eigen(A)
  
  ####### INITIAL STATE / COVARIANCE #######
  # The state/covariance is either given by user or obtained from solving the
  # stationary mean, and then solving for the covariance.
  # In principle these are coupled equations, but believe that root-finding
  # both simultaneously can lead to likelihood blow-up.
  inputVec = inputMat[1,]
  if(estimate.initial){
    # 1. Root-find stationary mean
    opt <- stats::nlminb(numeric(n.states), function(x) sum(f__(x, parVec, inputVec)^2))
    stateVec <- opt$par
    # 2. Use stationary mean to solve lyapunov eq. for associated covariance
    A <- dfdx__(stateVec, parVec, inputVec)
    G <- g__(stateVec, parVec, inputVec)
    Q <- G %*% t(G)
    P <- kronecker(A,I0) + kronecker(I0,A)
    X <- -solve(P, as.numeric(Q))
    covMat <- matrix(X, nrow=n.states)
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
    C = E %*% H
    e = y - C %*% stateVec
    V = E %*% V0 %*% t(E)
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
    
    ###### TIME UPDATE #######
    # augment input vector to account for constant terms in B
    inputVec = c(1,inputMat[i,])
    # Perform one-step prediction of mean and covariance
    stateVec <- Ahat %*% stateVec + Bhat %*% inputVec
    covMat <- Ahat %*% covMat %*% Ahat_T + Vhat
    # Store prior predictions
    xPrior[[i+1]] = stateVec
    pPrior[[i+1]] = covMat
    
    ######## DATA UPDATE ########
    # We update the state and covariance based on the "new" measurement
    inputVec = inputMat[i+1,]
    obsVec = obsMat[i+1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      y = obsVec[obsVec_bool]
      E = E0[obsVec_bool,, drop=FALSE]
      C = E %*% H
      e = y - C %*% stateVec
      V = E %*% V0 %*% t(E)
      R = C %*% covMat %*% t(C) + V
      K = covMat %*% t(C) %*% RTMB::solve(R)
      # Likelihood Contribution
      nll = nll - RTMB::dmvnorm(e, Sigma=R, log=TRUE)
      # Update State/Cov
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
      Innovation[[i+1]] <- e
      InnovationCovariance[[i+1]] <- R
    }
    xPost[[i+1]] = stateVec
    pPost[[i+1]] = covMat
  }
  ###### MAIN LOOP END #######
  
  ####### RETURN #######
  returnlist <- list(nll=nll,
                     xPost = xPost,
                     pPost = pPost,
                     xPrior = xPrior,
                     pPrior = pPrior,
                     Innovation = Innovation,
                     InnovationCovariance = InnovationCovariance)
  
  return(invisible(returnlist))
}

#######################################################
#######################################################
# EKF CPP IMPLEMENTATION (FOR OPTIMIZATION)
#######################################################
#######################################################
makeADFun_lkf_cpp = function(self, private){
  
  # Data ----------------------------------------
  
  # add mandatory entries to data
  tmb.data = list(
    
    # methods and purpose
    ode_solver = private$ode.solver,
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # time-steps
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    
    # loss function
    loss_function = private$loss$loss,
    loss_threshold_value = private$loss$c,
    tukey_loss_parameters = private$tukey.pars,
    
    # system size
    number_of_state_eqs = private$number.of.states,
    number_of_obs_eqs = private$number.of.observations,
    number_of_diffusions = private$number.of.diffusions,
    
    # estimate stationary levels
    estimate_stationary_initials = as.numeric(private$estimate.initial),
    initial_variance_scaling = private$initial.variance.scaling,
    
    # parameter bounds
    par_lb = sapply(private$parameters, function(x) x$lower),
    par_ub = sapply(private$parameters, function(x) x$lower),
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names])
  )
  
  # unscented parameters
  ukf_hyperpars = c()
  if(private$method=="ukf"){
    ukf_hyperpars = list(private$ukf_hyperpars)
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
  data = c(tmb.data, tmb.map.data, ukf_pars=ukf_hyperpars)
  
  
  # Parameters ----------------------------------------
  parameters = lapply(private$parameters, function(x) x[["initial"]]) # Initial parameter values
  
  
  # Create AD-likelihood function ---------------------------------------
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
#######################################################
# RETURN FIT FOR EKF RTMB
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
  npars <- length(private$fit$par.fixed)
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
        sdttemp[-remove.ids] <- std.dev
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
      std.dev <- try_withWarningRecovery(sqrt(diag(covariance)))
    }
    if(!inherits(std.dev,"try-error")){
      private$fit$sd.fixed <- std.dev
    }
    
  }
  
  # States -----------------------------------
  
  # Extract reported items from nll
  comptime <- system.time(rep <- ekf_r(model$getParameters()[,"estimate"], self, private))
  comptime = format(round(as.numeric(comptime["elapsed"])*1e2)/1e2,digits=5,scientific=F)
  if(!private$silent) message("...took ", comptime, " seconds")
  
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
  rowNAs = as.matrix(!is.na(private$data[private$obs.names])[-1,])
  sumrowNAs = rowSums(rowNAs)
  
  innovation = rep$Innovation
  innovation.cov = rep$InnovationCovariance
  innovation[[1]] = NULL
  innovation.cov[[1]] = NULL
  
  temp.res = matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
  temp.var =  matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
  
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
