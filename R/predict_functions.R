#######################################################
#######################################################
# EKF R-IMPLEMENTATION (FOR REPORTING)
#######################################################
#######################################################

ekf_r_prediction = function(self, private)
{
  
  
  # parameters ----------------------------------------
  # parVec = sapply(private$parameters, function(x) x[["initial"]])
  parVec <- private$pars
  
  # Data ----------------------------------------
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # methods and purpose
  ode_solver = private$ode.solver
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
  for(i in seq_along(private$rekf.function.strings)){
    eval(parse(text=private$rekf.function.strings[[i]]))
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
  
  ####### STORAGE #######
  # print(k.ahead)
  # print(last.pred.index)
  xPred <- lapply(1:last.pred.index, function(x) matrix(NA, nrow=k.ahead+1,ncol=n.states))
  pPred <- lapply(1:last.pred.index, function(x) vector("list", length=k.ahead+1))
  
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
    
    xPred[[i]][1,] <- stateVec
    pPred[[i]][[1]] <- covMat
    
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
      
      xPred[[i]][k+1,] <- stateVec
      pPred[[i]][[k+1]] <- covMat
      
    }
    
    stateVec <- xPred[[i]][2,]
    covMat <- pPred[[i]][[2]]
    
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
  return(invisible(list(xPred=xPred, pPred=pPred)))
}

rcpp_prediction = function(self, private){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # predict using c++ function
  mylist <- execute_ekf_prediction(private$rcpp_function_ptr$f,
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
  
  
  private$prediction = mylist
  
  return(invisible(NULL))
}

create_return_prediction = function(return.covariance, return.k.ahead, self, private){
  
  # Simlify variable names
  n = private$number.of.states
  n.ahead = private$n.ahead
  state.names = private$state.names
  last.pred.index = private$last.pred.index
  if(is.null(return.k.ahead)) return.k.ahead <- 0:n.ahead
  
  # Create return data.frame
  df.out = data.frame(matrix(nrow=last.pred.index*(n.ahead+1), ncol=5+n+n^2))
  disp_names = sprintf(rep("cor.%s.%s",n^2),rep(state.names, each=n),rep(state.names,n))
  disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  # var_bool = !stringr::str_detect(disp_names,"cor")
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
  
  ##### STATES PREDICTIONS ######
  df.out[,state.names] = do.call(rbind,lapply(private$prediction$Xpred, function(cur.list) do.call(rbind, cur.list)))
  if(return.covariance){ #covariance
    df.out[,disp_names] = do.call(rbind, lapply(private$prediction$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(x))) } ))
  } else { #correlation
    df.out[,disp_names] = do.call(rbind, lapply(private$prediction$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(cov2cor(x)))) } )) 
    diag.ids = seq(from=1,to=n^2,by=n+1)
    df.out[,disp_names[diag.ids]] = do.call(rbind, lapply(private$prediction$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(diag(x)))) } ))
  }
  
  ##### OBSERVATION PREDICTIONS #####
  # calculate observations at every time-step in predict
  inputs.df = private$data[df.out[,"j."]+1,private$input.names]
  
  named.pars.list = as.list(private$pars)
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
  names(obs.df.predict) = paste(private$obs.names)
  
  # add data observation to output data.frame 
  obs.df.data = private$data[df.out[,"j."]+1, private$obs.names, drop=F]
  names(obs.df.data) = paste(private$obs.names,".data",sep="")
  
  df.obs = cbind(df.obs, obs.df.predict, obs.df.data)
  
  # return only specific n.ahead
  df.out = df.out[df.out[,"k.ahead"] %in% return.k.ahead,]
  df.obs = df.obs[df.obs[,"k.ahead"] %in% return.k.ahead,]
  
  list.out = list(states = df.out, observations = df.obs)
  # class(list.out) = c(class(list.out), "ctsmTMB.pred")
  
  private$prediction = list.out
  
  return(list.out)
}
