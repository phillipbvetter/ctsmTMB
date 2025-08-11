#######################################################
#######################################################
# EKF SIMULATION R
#######################################################
#######################################################

ekf_r_simulation = function(parVec, self, private, nsims)
{
  
  # Data ----------------------------------------
  n.states <- private$number.of.states
  n.obs <- private$number.of.observations
  n.pars <- private$number.of.pars
  n.diffusions <- private$number.of.diffusions
  n.inputs <- private$number.of.inputs
  estimate.initial <- private$estimate.initial
  
  # initial
  stateVec = private$initial.state$x0
  covMat = private$initial.state$p0
  
  create.function.from.string.body("f__", "ans", private$r.function.strings$f)
  create.function.from.string.body("dfdx__", "ans", private$r.function.strings$dfdx)
  create.function.from.string.body("g__",  "ans", private$r.function.strings$g)
  create.function.from.string.body("h__", "ans", private$r.function.strings$h)
  create.function.from.string.body("dhdx__", "ans", private$r.function.strings$dhdx)
  create.function.from.string.body("hvar__matrix", "ans", private$r.function.strings$hvar__matrix)
  
  # various utility functions for likelihood calculations ---------------------
  # Note - order can be important here
  getOdeSolvers()
  if(estimate.initial) {
    getInitialStateEstimator()
  }
  getKalmanFunctions()
  getSimulationFunctions()
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  # prediction settings
  k.ahead <- private$n.ahead
  last.pred.index <- private$last.pred.index
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  sde_timestep_size <- private$simulation.timestep.size
  sde_timesteps <- private$simulation.timesteps
  
  ####### STORAGE #######
  xSim <- lapply(1:last.pred.index, function(x) vector("list", length=k.ahead+1))
  
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
  
  ###### TIME LOOP ####### 
  for(i in 1:last.pred.index){
    
    # We draw nsims samples of stateVecs from a multivariate normal; z = u + A * dB
    # where u is the mean (stateVec) and A is cholesky factor of covariance matrix (sqrt(covMat))
    # and dB is i.d.d normal vector
    # Each column in stateMat is a sample of a stateVec
    stateMat <- matrix(rep(stateVec, times=nsims), ncol=nsims) + 
      Matrix::chol(covMat) %*% matrix(stats::rnorm(nsims*n.states), ncol=nsims)
    xSim[[i]][[1]] <- t(stateMat)
    
    ###### K-STEP AHEAD LOOP ####### 
    for(k in 1:k.ahead){
      inputVec = inputMat[i+k-1,]
      dinputVec = (inputMat[i+k,] - inputVec)/sde_timesteps[i+k-1]
      
      # Solve moment ODEs
      for(j in 1:ode_timesteps[i+k-1]){
        stateMat <- euler.maruyama.simulation(stateMat, parVec, inputVec, sde_timestep_size[i+k-1])
        inputVec = inputVec + dinputVec
      }
      
      xSim[[i]][[k+1]] <- t(stateMat)
    }
    
    inputVec = inputMat[i,]
    dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
    
    # Solve ODE 1-Step Forward to get next posterior for simulations
    for(j in 1:ode_timesteps[i]){
      sol = ode.integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i])
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
      data.update <- kalman.data.update(stateVec, covMat, parVec, inputVec, obsVec, obsVec_bool, E0, I0)
      stateVec <- data.update[[1]]
      covMat <- data.update[[2]]
    }
    
    # End Loop
  }
  ###### MAIN LOOP END #######
  
  ####### RETURN #######
  private$simulation = xSim
  return(invisible(self))
}

#######################################################
#######################################################
# EKF SIMULATION C++
#######################################################
#######################################################

# Stochastic Euler-Maruyama simulation function that calls the underlying Rcpp simulation function
ekf_rcpp_simulation = function(parVec, self, private, n.sims){
  
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
                                   parVec,
                                   private$initial.state$p0,
                                   private$initial.state$x0,
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

create_ekf_simulation_return = function(return.k.ahead, n.sims, self, private){
  
  list.out = vector("list",length=private$number.of.states)
  names(list.out) = private$state.names
  
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
  # eval(parse(text=private$rekf.function.strings$h))
  
  # We should take the state vector, and pass it through the observation vector, right?
  # The inputs should come from the inputMat for the j'th row...
  # There are so many though...
  
  private$simulation = list( states = list.out, observations = list() )
  
  return(invisible(NULL))
}
