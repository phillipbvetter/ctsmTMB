#######################################################
#######################################################
# PREDICTION EKF R IMPLEMENTATION
#######################################################
#######################################################

ekf_r_prediction = function(parVec, self, private)
{
  
  # Data ----------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # initial
  stateVec <- private$initial.state$x0
  covMat = private$initial.state$p0
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  create.state.space.functions.for.filtering()
  
  # create.function.from.string.body("f__", "ans", private$r.function.strings$f)
  # create.function.from.string.body("dfdx__", "ans", private$r.function.strings$dfdx)
  # create.function.from.string.body("g__",  "ans", private$r.function.strings$g)
  # create.function.from.string.body("h__", "ans", private$r.function.strings$h)
  # create.function.from.string.body("dhdx__", "ans", private$r.function.strings$dhdx)
  # create.function.from.string.body("hvar__matrix", "ans", private$r.function.strings$hvar__matrix)
  
  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  getOdeSolvers()
  if(estimate.initial) {
    getInitialStateEstimator()
  }
  getKalmanFunctions()
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  
  # prediction settings
  k.ahead <- private$n.ahead
  last.pred.index <- private$last.pred.index
  
  ####### STORAGE #######
  predMats <- lapply(1:last.pred.index, function(x) matrix(NA,nrow=k.ahead+1,ncol=n.states+n.states^2))
  
  ####### Pre-Allocated Object #######
  I0 <- diag(n.states)
  E0 <- diag(n.obs)
  
  ####### INITIAL STATE / COVARIANCE #######
  inputVec = inputMat[1,]
  if(estimate.initial){
    stateVec <- f.initial.state.newton(c(parVec, inputVec))
    covMat <- f.initial.covar.solve(stateVec, parVec, inputVec)
  }
  
  ######## (PRE) DATA UPDATE ########
  obsVec = obsMat[1,]
  obsVec_bool = !is.na(obsVec)
  if(any(obsVec_bool)){
    data.update <- kalman.data.update(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
    stateVec <- data.update[[1]]
    covMat <- data.update[[2]]
  }
  
  for(i in 1:last.pred.index){ ###### TIME LOOP ####### 
    
    predMats[[i]][1,] <- c(stateVec, covMat)
    
    for(k in 1:k.ahead){ ###### K-STEP AHEAD LOOP ####### 
      inputVec = inputMat[i+k-1,]
      dinputVec = (inputMat[i+k,] - inputVec)/ode_timesteps[i+k-1]
      
      # Solve moment ODEs
      for(j in 1:ode_timesteps[i+k-1]){
        sol = ode.integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i+k-1])
        stateVec = sol[[1]]
        covMat = sol[[2]]
        inputVec = inputVec + dinputVec
      }
      predMats[[i]][k+1,] <- c(stateVec, covMat)
      
    }
    stateVec <- head(predMats[[i]][2,], n.states)
    covMat <- matrix(tail(predMats[[i]][2,], n.states^2), nrow=n.states)
    
    ######## DATA UPDATE ########
    inputVec = inputMat[i+1,]
    obsVec = obsMat[i+1,]
    obsVec_bool = !is.na(obsVec)
    if(any(obsVec_bool)){
      data.update <- kalman.data.update(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
      stateVec <- data.update[[1]]
      covMat <- data.update[[2]]
    }
    
  } ###### MAIN LOOP END #######
  
  ####### STORE PREDICTION #######
  private$prediction = list(predMats=predMats)
  
  ####### RETURN #######
  return(invisible(self))
}

#######################################################
#######################################################
# PREDICTION EKF C++ IMPLEMENTATION
#######################################################
#######################################################

ekf_rcpp_prediction = function(parVec, self, private){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # predict using c++ function
  mylist <- execute_ekf_prediction2(private$rcpp_function_ptr$f,
                                    private$rcpp_function_ptr$g,
                                    private$rcpp_function_ptr$dfdx,
                                    private$rcpp_function_ptr$h,
                                    private$rcpp_function_ptr$dhdx,
                                    private$rcpp_function_ptr$hvar,
                                    obsMat,
                                    inputMat,
                                    parVec,
                                    private$initial.state$p0,
                                    private$initial.state$x0,
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

create_ekf_predict_return = function(return.covariance, return.k.ahead, self, private){
  
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
  names(named.pars.list) = private$parameter.names
  # create environment
  env.list = c(
    # states
    as.list(df.out[state.names]),
    # inputs
    as.list(inputs.df),
    # free and fixed parameters
    named.pars.list
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
