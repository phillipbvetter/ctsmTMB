#######################################################
#######################################################
# EKF LKF UKF SIMULATION R
#######################################################
#######################################################

ekf_lkf_ukf_simulate_r = function(parVec, self, private, nsims)
{
  
  # Data ----------------------------------------
  getSystemDimensions()
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  # prediction settings
  k.ahead <- private$n.ahead
  last.pred.index <- private$last.pred.index
  # time-steps
  sde_timestep_size <- private$simulation.timestep.size
  sde_timesteps <- private$simulation.timesteps
  
  # various utility functions for likelihood calculations ---------------------
  create.state.space.functions.for.filtering()
  getSimulationFunctions()
  
  # Get filtered states
  filt <- switch(private$method,
                 lkf = lkf_filter_r(parVec, self, private),
                 ekf = ekf_filter_r(parVec, self, private),
                 ukf = ukf_filter_r(parVec, self, private)
  )
  xPost <- filt$xPost
  pPost <- filt$pPost
  
  ####### STORAGE #######
  xSim <- lapply(1:last.pred.index, function(x) vector("list", length=k.ahead+1))
  
  ###### TIME LOOP ####### 
  for(i in 1:last.pred.index){
    
    stateVec <- xPost[[i]]
    covMat <- matrix(pPost[[i]],nrow=n.states)
    
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
      for(j in 1:sde_timesteps[i+k-1]){
        stateMat <- euler.maruyama.simulation(stateMat, parVec, inputVec, sde_timestep_size[i+k-1])
        inputVec = inputVec + dinputVec
      }
      xSim[[i]][[k+1]] <- t(stateMat)
    }
  }
  ###### MAIN LOOP END #######
  
  ####### STORE SIMULATION #######
  private$simulation = xSim
  
  ####### RETURN #######
  return(invisible(self))
}

#######################################################
#######################################################
# EKF SIMULATION C++
#######################################################
#######################################################

lkf_ekf_ukf_simulate_rcpp <- function(pars, self, private, n.sims){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  output <- NULL
  if(private$method=="lkf"){
    output = lkf_simulate_rcpp(private$rcpp_function_ptr$f,
                               private$rcpp_function_ptr$g,
                               private$rcpp_function_ptr$dfdx,
                               private$rcpp_function_ptr$h,
                               private$rcpp_function_ptr$dhdx,
                               private$rcpp_function_ptr$hvar,
                               private$rcpp_function_ptr$dfdu,
                               obsMat,
                               inputMat,
                               pars,
                               private$initial.state$p0,
                               private$initial.state$x0,
                               private$ode.timestep.size,
                               private$simulation.timestep.size,
                               private$simulation.timesteps,
                               numeric_is_not_na_obsMat,
                               number_of_available_obs,
                               private$number.of.diffusions,
                               private$last.pred.index,
                               private$n.ahead,
                               n.sims)
  } 
  if(private$method=="ekf"){
    output = ekf_simulate_rcpp(private$rcpp_function_ptr$f,
                               private$rcpp_function_ptr$g,
                               private$rcpp_function_ptr$dfdx,
                               private$rcpp_function_ptr$h,
                               private$rcpp_function_ptr$dhdx,
                               private$rcpp_function_ptr$hvar,
                               obsMat,
                               inputMat,
                               pars,
                               private$initial.state$p0,
                               private$initial.state$x0,
                               private$ode.timestep.size,
                               private$ode.timesteps,
                               private$simulation.timestep.size,
                               private$simulation.timesteps,
                               numeric_is_not_na_obsMat,
                               number_of_available_obs,
                               private$number.of.diffusions,
                               private$last.pred.index,
                               private$n.ahead,
                               private$ode.solver,
                               n.sims)
  }
  if(private$method == "ukf"){
    output <- ukf_simulate_rcpp(private$rcpp_function_ptr$f,
                                private$rcpp_function_ptr$g,
                                private$rcpp_function_ptr$dfdx,
                                private$rcpp_function_ptr$h,
                                private$rcpp_function_ptr$dhdx,
                                private$rcpp_function_ptr$hvar,
                                obsMat,
                                inputMat,
                                pars,
                                private$initial.state$p0,
                                private$initial.state$x0,
                                private$ode.timestep.size,
                                private$ode.timesteps,
                                private$simulation.timestep.size,
                                private$simulation.timesteps,
                                numeric_is_not_na_obsMat,
                                number_of_available_obs,
                                private$ukf_hyperpars,
                                private$number.of.diffusions,
                                private$last.pred.index,
                                private$n.ahead,
                                private$ode.solver,
                                n.sims)
  }
  
  private$simulation = output
  
  return(invisible(NULL))
}
