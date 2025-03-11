#######################################################
# CONSTRUCT LAPLACE MAKEADFUN WITH RTMB
#######################################################

makeADFun_laplace_rtmb = function(self, private)
{
  
  # Data ----------------------------------------
  
  # initial states and covariance
  stateVec <- private$initial.state$x0
  covMat <- private$initial.state$p0
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # indices with non-na observations
  iobs <- private$iobs
  
  # parameters ----------------------------------------
  parameters = c(
    lapply(private$parameters, function(x) x[["initial"]]),
    private$tmb.initial.state
  )
  
  # adjoints ----------------------------------------
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
  
  ################################################
  # Functions
  ################################################
  
  for(i in seq_along(private$rtmb.function.strings)){
  eval(parse(text=private$rtmb.function.strings[[i]]))
  }
  
  # error function in terms of pnorm
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
  
  # likelihood function --------------------------------------
  
  laplace.nll = function(p){
    
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    
    # set negative log-likelihood
    nll = 0
    
    # small identity matrix
    small_identity = diag(1e-8, nrow=n.states, ncol=n.states)
    
    # fixed effects parameter vector
    parVec <- do.call(c, p[1:n.pars])
    
    inputVec = inputMat[1,]
    if(estimate.initial){
      # 1. Root-find stationary mean
      .F <- RTMB::MakeTape(function(x) sum(f__(x, parVec, inputVec)^2), numeric(n.states))
      stateVec <- .F$newton(1:n.states)(numeric(0))
      # # 2. Use stationary mean to solve lyapunov eq. for associated covariance
      A <- dfdx__(stateVec, parVec, inputVec)
      G <- g__(stateVec, parVec, inputVec)
      Q <- G %*% t(G)
      P <- kron_left(A) + kron_right(A)
      X <- -RTMB::solve(P, as.numeric(Q))
      covMat <- RTMB::matrix(X, nrow=n.states)
    }
    
    # extract state random effects and fix initial condition
    stateMat <- do.call(cbind, p[(n.pars+1):length(p)])
    
    # prior contribution
    z0 <- stateMat[1,] - stateVec
    nll <- nll -  RTMB::dmvnorm(z0, Sigma=covMat, log=TRUE)
    
    ###### TIME LOOP START #######
    for(i in 1:(nrow(obsMat)-1)) {
      
      # Define inputs and use first order input interpolation
      inputVec = inputMat[i,]
      dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
      
      ###### BETWEEN TIME POINTS LOOP START #######
      for(j in 1:ode_timesteps[i]){
        
        # grab current and next state
        x_now = stateMat[ode_cumsum_timesteps[i]+j,]
        x_next = stateMat[ode_cumsum_timesteps[i]+j+1,]
        
        # compute drift (vector) and diffusion (matrix)
        f = f__(x_now, parVec, inputVec)
        g = g__(x_now, parVec, inputVec)
        inputVec = inputVec + dinputVec
        
        # assume multivariate gauss distribution according to euler-step
        # and calculate the likelihood
        z = x_next - (x_now + f * ode_timestep_size[i])
        v = (g %*% t(g) + small_identity) * ode_timestep_size[i]
        nll = nll - RTMB::dmvnorm(z, Sigma=v, log=TRUE)
      }
      ###### BETWEEN TIME POINTS LOOP END #######
    }
    ###### TIME LOOP END #######
    
    obsMat = RTMB::OBS(obsMat)
    ###### DATA UPDATE START #######
    for(i in 1:n.obs){
      iobs.vec <- iobs[[i]]
      for(j in 1:length(iobs[[i]])){
        
        # Get index where observation is available
        # k = iobs[[i]][[j]]
        k <- iobs.vec[[j]]
        
        # Get corresponding input, state and observation
        inputVec = inputMat[k,]
        stateVec = stateMat[ode_cumsum_timesteps[k]+1,]
        obsScalar = obsMat[k,i]
        
        # Observation equation and variance
        h_x = h__(stateVec, parVec, inputVec)[i]
        hvar_x = hvar__(stateVec, parVec, inputVec)[i]
        
        # likelihood contribution
        nll = nll - RTMB::dnorm(obsScalar, mean=h_x, sd=sqrt(hvar_x), log=TRUE)
      }
    }
    ###### DATA UPDATE END #######
    
    # return
    return(nll)
  }
  
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
#######################################################
# RETURN FIT FOR LAPLACE
#######################################################
#######################################################
calculcate_fit_statistics_laplace <- function(self, private, laplace.residuals){
  
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
    message("Calculating one-step ahead residuls...")
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
