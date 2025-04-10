#######################################################
#######################################################
# EKF R-IMPLEMENTATION (FOR REPORTING)
#######################################################
#######################################################

ekf_r_prediction = function(self, private)
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
      covMat <- matrix(tail(stateVec_and_covMat, -n.states),nrow=n.states)
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
  predMats <- lapply(1:last.pred.index, function(x) matrix(NA,nrow=k.ahead+1,ncol=n.states+n.states^2))
  
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
    # Start Loop
    
    predMats[[i]][1,] <- c(stateVec, covMat)
    
    ###### K-STEP AHEAD LOOP ####### 
    for(k in 1:k.ahead){
      
      inputVec = inputMat[i+k-1,]
      dinputVec = (inputMat[i+k,] - inputVec)/ode_timesteps[i+k-1]
      
      # Solve moment ODEs
      for(j in 1:ode_timesteps[i+k-1]){
        sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i+k-1])
        stateVec = sol[[1]]
        covMat = sol[[2]]
        inputVec = inputVec + dinputVec
      }
      
      predMats[[i]][k+1,] <- c(stateVec, covMat)
      
    }
    
    stateVec <- head(predMats[[i]][2,],n.states)
    covMat <- matrix(tail(predMats[[i]][2,],n.states^2),nrow=n.states)
    
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
  
  ####### STORE PREDICTION #######
  # private$prediction = list(Xpred=xPred, Ppred=pPred, predMats=predMats)
  private$prediction = list(predMats=predMats)
  
  ####### RETURN #######
  return(invisible(self))
}

ekf_rcpp_prediction = function(self, private){
  
  if(!any(private$ode.solver==c(1,2))){
    stop("Predictions using C++ currently only support 'euler' or 'rk4' ODE solvers.") 
  }
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # predict using c++ function
  # mylist <- execute_ekf_prediction(private$rcpp_function_ptr$f,
  #                                  private$rcpp_function_ptr$g,
  #                                  private$rcpp_function_ptr$dfdx,
  #                                  private$rcpp_function_ptr$h,
  #                                  private$rcpp_function_ptr$dhdx,
  #                                  private$rcpp_function_ptr$hvar,
  #                                  obsMat,
  #                                  inputMat,
  #                                  private$pars,
  #                                  private$pred.initial.state$p0,
  #                                  private$pred.initial.state$x0,
  #                                  private$ode.timestep.size,
  #                                  private$ode.timesteps,
  #                                  numeric_is_not_na_obsMat,
  #                                  number_of_available_obs,
  #                                  private$number.of.states,
  #                                  private$number.of.observations,
  #                                  private$last.pred.index,
  #                                  private$n.ahead,
  #                                  private$ode.solver)
  mylist <- execute_ekf_prediction2(private$rcpp_function_ptr$f,
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
                                    numeric_is_not_na_obsMat,
                                    number_of_available_obs,
                                    private$number.of.states,
                                    private$number.of.observations,
                                    private$last.pred.index,
                                    private$n.ahead,
                                    private$ode.solver)
  
  ####### STORE PREDICTION #######
  private$prediction = mylist
  
  ####### RETURN #######
  return(invisible(self))
}

create_return_prediction = function(return.covariance, return.k.ahead, use.cpp, self, private){
  
  # Simlify variable names
  n = private$number.of.states
  n.ahead = private$n.ahead
  state.names = private$state.names
  last.pred.index = private$last.pred.index
  
  # Create return data.frame
  df.out = data.frame(matrix(nrow=last.pred.index*(n.ahead+1), ncol=5+n+n^2))
  disp_names = sprintf(rep("cor.%s.%s",n^2),rep(state.names, each=n),rep(state.names,n))
  disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  if(return.covariance){
    disp_names = sprintf(rep("cov.%s.%s",n^2),rep(state.names,each=n),rep(state.names,n))
    disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  }
  names(df.out) = c("i.","j.","t.i","t.j","k.ahead",state.names,disp_names)
  
  # Fill out data.frame
  ran = 0:(last.pred.index-1)
  df.out["i."] = rep(ran,each=n.ahead+1)
  df.out["j."] = df.out["i."] + rep(0:n.ahead,last.pred.index)
  df.out["t.i"] = rep(private$data$t[ran+1],each=n.ahead+1)
  df.out["t.j"] = private$data$t[df.out[,"i."]+1+rep(0:n.ahead,last.pred.index)]
  df.out["k.ahead"] = rep(0:n.ahead,last.pred.index)
  df.obs = df.out[c("i.","j.","t.i","t.j","k.ahead")]
  
  df.out[, c(state.names, disp_names)] <- do.call(rbind, private$prediction$predMats)
  if(!return.covariance){
    diag.ids <- seq(from=1,to=n^2,by=n+1)
    .seq <- seq(from=1,to=n^2,by=1)
    non.diag.ids <- .seq[!(.seq %in% diag.ids)]
    df.out[,disp_names[non.diag.ids]] <- t(apply(df.out, 1, function(x) as.vector(stats::cov2cor(matrix(tail(x, n^2), nrow=n)))))[,non.diag.ids]
  }
  
  ##### OBSERVATION PREDICTIONS #####
  # calculate observations at every time-step in predict
  inputs.df = private$data[df.out[,"j."]+1,private$input.names]
  
  named.pars.list = as.list(private$pars)
  # names(named.pars.list) = names(private$free.pars)
  names(named.pars.list) = private$parameter.names
  # create environment
  env.list = c(
    # states
    as.list(df.out[state.names]),
    # inputs
    as.list(inputs.df),
    # free parameters
    named.pars.list
    # fixed parameters
    # lapply(private$fixed.pars, function(x) x$initial)
  )
  
  # calculate observations
  obs.df.predict = as.data.frame(
    lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = env.list)})
  )
  names(obs.df.predict) = paste(private$obs.names)
  
  # add data observation to output data.frame 
  obs.df.data = private$data[df.out[,"j."]+1, private$obs.names, drop=F]
  names(obs.df.data) = paste(private$obs.names,".data",sep="")
  
  df.obs = cbind(df.obs, obs.df.predict, obs.df.data)
  
  # return only specific n.ahead
  df.out = df.out[df.out[,"k.ahead"] %in% return.k.ahead,]
  df.obs = df.obs[df.obs[,"k.ahead"] %in% return.k.ahead,]
  
  list.out = list(states = df.out, observations = df.obs)
  class(list.out) = c(class(list.out), "ctsmTMB.pred")
  
  private$prediction = list.out
  
  return(list.out)
}
