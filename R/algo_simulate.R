#######################################################
#######################################################
# EKF LKF UKF SIMULATION R
#######################################################
#######################################################

ekf_lkf_ukf_simulate_r = function(parVec, self, private, nsims)
{
  
  # Data ----------------------------------------
  get_sys_dims()
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  # prediction settings
  k.ahead <- private$k.ahead
  last.pred.index <- private$last.pred.index
  # time-steps
  sde_timestep_size <- private$simulation.timestep.size
  sde_timesteps <- private$simulation.timesteps
  
  # various utility functions for likelihood calculations ---------------------
  create_state_space_functions_for_filtering()
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
  private$simulation.raw = xSim
  
  ####### RETURN #######
  return(invisible(self))
}
