#######################################################
# CONSTRUCT LAPLACE MAKEADFUN WITH RTMB
#######################################################

makeADFun_laplace2_rtmb = function(self, private)
{
  
  # Tape Configration ----------------------
  # The best option was not to change defaults
  RTMB::TapeConfig(atomic="disable")
  RTMB::TapeConfig(vectorize="disable")
  
  # global values in likelihood function ------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  estimate.initial <- private$estimate.initial
  n.dbs <- nrow(private$tmb.initial.state) - 1
  
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
  db.len <- n.diffusions * n.dbs
  dB = matrix(numeric(n.diffusions * n.dbs), nrow=n.dbs, ncol=n.diffusions)
  parameters = c(
    lapply(private$parameters, function(x) x[["initial"]]),
    private$tmb.initial.state,
    dB = list(dB)
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
  
  # likelihood function --------------------------------------
  laplace2.nll = function(p){
    
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    
    # set negative log-likelihood
    nll = 0
    
    # small identity matrix
    # tiny = diag(1e-8, nrow=n.states, ncol=n.states)
    tiny <- 1e-5 * diag(n.states)
    
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
    stateMat <- do.call(cbind, p[(n.pars+1):(length(p)-1)])
    dstateMat <- RTMB::apply(stateMat, 2, diff)
    I0 <- diag(n.diffusions)
    
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
        
        # current index in state matrix
        cur.id <- ode_cumsum_timesteps[i]+j
        
        # compute drift (vector) and diffusion (matrix)
        stateVec <- stateMat[cur.id,]
        f = f__(stateVec, parVec, inputVec)
        g = g__(stateVec, parVec, inputVec)
        inputVec = inputVec + dinputVec
        
        # Compute expected dX from Euler Maruyama
        dstateVecPred <- f * ode_timestep_size[i] + g %*% p$dB[cur.id,]
        
        # Likelihood contribution from state difference (diagonal covariance)
        z <- dstateMat[cur.id,] - dstateVecPred
        nll = nll - RTMB::dmvnorm(z, Sigma=tiny, log=TRUE)
        
        # Likelihood contribution from dBs
        nll = nll - RTMB::dmvnorm(p$dB[cur.id,], Sigma=ode_timestep_size[i]*I0, log=TRUE)
      }
      ###### BETWEEN TIME POINTS LOOP END #######
    }
    ###### TIME LOOP END #######
    
    obsMat = RTMB::OBS(obsMat)
    ##### DATA UPDATE START #######
    for(i in 1:n.obs){
      iobs.vec <- iobs[[i]]
      for(j in 1:length(iobs[[i]])){

        # Get index where observation is available
        k <- iobs.vec[[j]]

        # Get corresponding input, state and observation
        inputVec = inputMat[k,]
        stateVec = stateMat[ode_cumsum_timesteps[k]+1,]
        obsScalar = obsMat[[k,i]]

        # Observation equation and variance
        h_x = h__(stateVec, parVec, inputVec)[[i]]
        hvar_x = hvar__(stateVec, parVec, inputVec)[[i]]

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
  
  nll = RTMB::MakeADFun(func = laplace2.nll, 
                        parameters=parameters, 
                        random=c(private$state.names,"dB"),
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
calculate_fit_statistics_laplace2 <- function(self, private, laplace.residuals){
  
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
  
  # clone and return -----------------------------------
  # clone private and return fit
  private$fit$private <- self$clone()$.__enclos_env__$private
  
  # set s3 class
  class(private$fit) = "ctsmTMB.fit"
}
