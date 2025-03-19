#######################################################
#######################################################
# EKF R-IMPLEMENTATION (FOR REPORTING)
#######################################################
#######################################################

ekf_r_simulation = function(self, private, nsims)
{
  
  
  # parameters ----------------------------------------
  parVec <- private$pars
  
  # Data ----------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # methods and purpose
  ode.solver = private$ode.solver
  k.ahead <- private$n.ahead
  last.pred.index <- private$last.pred.index
  
  # initial
  stateVec = private$pred.initial.state$x0
  covMat = private$pred.initial.state$p0
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  sde_timestep_size <- private$simulation.timestep.size
  sde_timesteps <- private$simulation.timesteps
  
  # 1-step covariance ODE ----------------------------------------
  cov_ode_1step = function(covMat, stateVec, parVec, inputVec){
    G <- g__(stateVec, parVec, inputVec)
    AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  # forward euler ----------------------------------------
  if(ode.solver==1){
    ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      
      X1 = stateVec + f__(stateVec, parVec, inputVec) * dt
      P1 = covMat + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
      
      return(list(X1,P1))
    }
  } else if (ode.solver==2) {
    # rk4 ----------------------------------------
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
  } else {
    # Construct input interpolators
    #---------------------------------
    input.interp.funs <- vector("list",length=n.inputs)
    for(i in 1:length(input.interp.funs)){
      input.interp.funs[[i]] <- stats::approxfun(x=inputMat[,1], y=inputMat[,i], rule=2)
    }
    
    # Construct ode fun for DeSolve
    #---------------------------------
    ode.fun <- function(time, stateVec_and_covMat, parVec){
      inputVec <- sapply(input.interp.funs, function(f) f(time))
      # 
      stateVec <- head(stateVec_and_covMat, n.states)
      covMat <- tail(stateVec_and_covMat, -n.states)
      # 
      G <- g__(stateVec, parVec, inputVec)
      AcovMat = dfdx__(stateVec, parVec, inputVec) %*% covMat
      #
      dX <- f__(stateVec, parVec, inputVec)
      dP <- AcovMat + t(AcovMat) + G %*% t(G)
      return(list(c(dX,dP)))
    }
    
    # Construct function to call in likelihood functon
    #---------------------------------
    ode_integrator <- function(covMat, stateVec, parVec, inputVec, dinputVec, dt){
      out <- deSolve::ode(y = c(stateVec, covMat),
                          times = c(inputVec[1], inputVec[1]+dt),
                          func = ode.fun,
                          parms = parVec,
                          method=ode.solver)[2,-1]
      
      return(
        list(head(out,n.states),
             matrix(tail(out,-n.states),nrow=n.states)
        )
      )
    }
  }
  
  euler_maruyama_simulation <- function(stateMat, parVec, inputVec, sim.dt){
    # We need to perform an euler maruyama step for each column in the stateMat,
    
    # Each column is f(stateVec)
    .F <- apply(stateMat, 2, function(x) f__(x, parVec, inputVec))
    # Each column is g(stateVec) as a vector
    .G <- lapply(1:nsims, function(i) g__(stateMat[,i],parVec, inputVec))
    # Draw nsims random variables for each db
    dB <- matrix(stats::rnorm(n.diffusions*nsims), ncol=nsims)
    # Perform multiplication
    GdB <- sapply(1:nsims, function(i) .G[[i]] %*% dB[,i])
    
    # Produce euler-maruyama step for all states simultaneously
    return(stateMat + .F * sim.dt + GdB * sqrt(sim.dt))
  }
  
  # user-defined functions ---------------------------
  for(i in seq_along(private$rekf.function.strings)){
    eval(parse(text=private$rekf.function.strings[[i]]))
  }
  
  # new functions ----------------------------------------
  # error function ----------------------------------------
  erf = function(x){
    y <- sqrt(2) * x
    2*RTMB::pnorm(y)-1
  }
  
  ####### STORAGE #######
  xSim <- lapply(1:last.pred.index, function(x) vector("list", length=k.ahead+1))
  
  ####### Pre-Allocated Object #######
  I0 <- diag(n.states)
  E0 <- diag(n.obs)
  
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
    K = covMat %*% t(C) %*% solve(R)
    # Update State/Cov
    stateVec = stateVec + K %*% e
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
  }
  
  ###### TIME LOOP ####### 
  for(i in 1:last.pred.index){
    
    # We draw nsims samples of stateVecs from a multivariate normal; z = u + A * dB
    # where u is the mean (stateVec) and A is cholesky factor of covariance matrix (sqrt(covMat))
    # and dB is i.d.d normal vector
    # Each column in stateMat is a sample of a stateVec
    stateMat <- matrix(rep(stateVec, times=nsims), ncol=nsims) + Matrix::chol(covMat) %*% 
      matrix(stats::rnorm(nsims*n.states),ncol=nsims)
    xSim[[i]][[1]] <- t(stateMat)
    
    ###### K-STEP AHEAD LOOP ####### 
    for(k in 1:k.ahead){
      inputVec = inputMat[i+k-1,]
      dinputVec = (inputMat[i+k,] - inputVec)/sde_timesteps[i+k-1]
      
      # Solve moment ODEs
      for(j in 1:ode_timesteps[i+k-1]){
        stateMat <- euler_maruyama_simulation(stateMat, parVec, inputVec, sde_timestep_size[i+k-1])
        inputVec = inputVec + dinputVec
      }
      
      xSim[[i]][[k+1]] <- t(stateMat)
    }
    
    inputVec = inputMat[i,]
    dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
    
    # Solve ODE 1-Step Forward to get next posterior for simulations
    for(j in 1:ode_timesteps[i]){
      sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i])
      stateVec = sol[[1]]
      covMat = sol[[2]]
      inputVec = inputVec + dinputVec
    }
    
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
      K = covMat %*% t(C) %*% solve(R)
      # Update State/Cov
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    }
    
    # End Loop
  }
  ###### MAIN LOOP END #######
  
  ####### RETURN #######
  private$simulation = xSim
  return(invisible(self))
}

# Stochastic Euler-Maruyama simulation function that calls the underlying Rcpp simulation function
rcpp_simulation = function(self, private, n.sims){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # Call C++ function to perform simulation
  mylist = execute_ekf_simulation2(private$rcpp_function_ptr$f,
                                   private$rcpp_function_ptr$g,
                                   private$rcpp_function_ptr$dfdx,
                                   private$rcpp_function_ptr$h,
                                   private$rcpp_function_ptr$dhdx,
                                   private$rcpp_function_ptr$hvar,
                                   obsMat,
                                   inputMat,
                                   private$pars,
                                   private$pred.initial.state$p0,
                                   private$pred.initial.state$x0,
                                   private$ode.timestep.size,
                                   private$ode.timesteps,
                                   private$simulation.timestep.size,
                                   private$simulation.timesteps,
                                   numeric_is_not_na_obsMat,
                                   number_of_available_obs,
                                   private$number.of.states,
                                   private$number.of.observations,
                                   private$number.of.diffusions,
                                   private$last.pred.index,
                                   private$n.ahead,
                                   private$ode.solver,
                                   n.sims)
  
  private$simulation = mylist
  
  return(invisible(NULL))
}

# Generates a user-friendly data.frame of prediction results from private$prediction
create_return_simulation = function(return.k.ahead, n.sims, self, private){
  
  
  list.out = vector("list",length=private$number.of.states)
  names(list.out) = private$state.names
  
  # setRownames = function(obj, nm){rownames(obj) = nm; return(obj)}
  # setColnames = function(obj, nm){colnames(obj) = nm; return(obj)}
  
  # Compute the prediction times for each horizon
  ran = 0:(private$last.pred.index-1)
  t.j = private$data$t[rep(ran,each=private$n.ahead+1)+1+rep(0:private$n.ahead,private$last.pred.index)]
  t.j.splitlist = split(t.j, ceiling(seq_along(t.j)/(private$n.ahead+1)))
  list.of.time.vectors = lapply(t.j.splitlist, function(x) data.frame(t.j=x))
  
  for(i in seq_along(list.out)){
    list.out[[i]] = stats::setNames(
      lapply(private$simulation, function(ls.outer){
        t(do.call(cbind, lapply(ls.outer, function(ls.inner) ls.inner[,i])))
      }),
      paste0("i",ran)
    )
  }
  
  for(i in seq_along(list.out)){
    for(j in seq_along(list.out[[i]])){
      list.out[[i]][[j]] = data.frame(i = j-1, 
                                      j = (j-1):(j+private$n.ahead-1), 
                                      t.i = rep(private$data$t[i],private$n.ahead+1),
                                      t.j = list.of.time.vectors[[j]][,"t.j"], 
                                      k.ahead = 0:private$n.ahead,
                                      list.out[[i]][[j]]
      )
      nams = paste0(private$state.names,1:n.sims)
      names(list.out[[i]][[j]]) = c("i","j","t.i","t.j","k.ahead",nams)
    }
  }
  
  # Observations
  eval(parse(text=private$rekf.function.strings$h))
  
  # We should take the state vector, and pass it through the observation vector, right?
  # The inputs should come from the inputMat for the j'th row...
  # There are so many though...
  
  
  
  
  
  private$simulation = list( states = list.out, observations = list() )
  
  return(invisible(NULL))
}
