# This function is the main caller for carrying out simulations
perform_simulation <- function(self, private, use.cpp, n.sims){
  
  
  if(private$method=="laplace"){
    stop("Simulations are not available for the Laplace method")
  }
  
  if(private$method=="laplace.thygesen"){
    stop("Simulations are not available for the Laplace-Thygesen method")
  }
  
  # time prediction
  comptime <- system.time(
    {
      # Predict with Rcpp implementation
      if(use.cpp){
        
        if(!private$silent) message("Simulating with C++...")
        
        lkf_ekf_ukf_simulate_rcpp(private$pars, self, private, n.sims)
        
        # Predict with R implementation
      } else {
        
        if(!private$silent) message("Simulating with R...")
        
        ekf_lkf_ukf_simulate_r(private$pars, self, private, n.sims)
        
      }
      
    }, gcFirst = FALSE)
  
  private$timer_simulation <- comptime
  
  return(invisible(self))
}

lkf_ekf_ukf_simulate_rcpp <- function(pars, self, private, n.sims){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  ids <- 1:private$number.of.observations - 1 #minus 1 for 0 indexing
  non_na_ids <- apply(obsMat, 1, function(x) ids[!is.na(x)], simplify = FALSE)
  any_available_obs <- sapply(non_na_ids, function(x) length(x) != 0)
  
  output <- NULL
  if(private$method=="lkf"){
    output <- lkf_simulate_rcpp(private$rcpp_function_ptr,
                                obsMat,
                                inputMat,
                                pars,
                                private$initial.state$p0,
                                private$initial.state$x0,
                                private$ode.timestep.size,
                                private$simulation.timestep.size,
                                private$simulation.timesteps,
                                any_available_obs,
                                non_na_ids,
                                private$number.of.diffusions,
                                private$last.pred.index,
                                private$k.ahead,
                                n.sims,
                                private$seed$state.seed)
  } 
  if(private$method=="ekf"){
    output <- ekf_simulate_rcpp(private$rcpp_function_ptr,
                                obsMat,
                                inputMat,
                                pars,
                                private$initial.state$p0,
                                private$initial.state$x0,
                                private$ode.timestep.size,
                                private$ode.timesteps,
                                private$simulation.timestep.size,
                                private$simulation.timesteps,
                                any_available_obs,
                                non_na_ids,
                                private$ode.solver,
                                private$last.pred.index,
                                private$k.ahead,
                                private$number.of.diffusions,
                                n.sims,
                                private$seed$state.seed)
  }
  if(private$method == "ukf"){
    output <- ukf_simulate_rcpp(private$rcpp_function_ptr,
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
                                private$ukf.hyperpars,
                                private$number.of.diffusions,
                                private$last.pred.index,
                                private$k.ahead,
                                private$ode.solver,
                                n.sims,
                                private$seed$state.seed)
  }
  
  private$simulation.raw <- output
  
  return(invisible(NULL))
}


# This function returns simulation results to the user
create_return_simulation <- function(return.k.ahead, n.sims, self, private){
  
  if(!private$silent) message("Returning results...")
  
  # create names for inner list
  inner.names <- paste0("i", 0:(private$last.pred.index-1))
  
  # Build returnlist for states
  state.list <- build_simulation_returnlist(private$simulation.raw,
                                            private$data$t,
                                            private$number.of.states,
                                            private$k.ahead,
                                            n.sims)
  # set state names
  names(state.list) <- private$state.names
  for(i in seq_along(state.list)){
    names(state.list[[i]]) <- inner.names
  }
  
  # create index/time list
  time.list <- build_simulation_timelists(private$data$t,
                                          private$last.pred.index,
                                          private$k.ahead)
  names(time.list) <- inner.names
  
  # Build returnlist for observations
  # First we must calculate the simulated observation trajectories
  simulation.raw.obs <- calculate_simulation_observations(
    private$simulation.raw,
    private$rcpp_function_ptr,
    t(as.matrix(private$data[private$input.names])),
    private$pars,
    private$number.of.states,
    private$number.of.observations,
    private$k.ahead,
    n.sims,
    private$seed$obs.seed
  )
  
  # # Now we can build the returnlist
  obs.list <- build_simulation_returnlist(simulation.raw.obs,
                                          private$data$t,
                                          private$number.of.observations,
                                          private$k.ahead,
                                          n.sims)
  names(obs.list) <- private$obs.names
  for(i in seq_along(obs.list)){
    names(obs.list[[i]]) <- inner.names
  }
  
  private$simulation <- list(states=state.list, observations=obs.list, times=time.list)
  
  return(invisible(self))
}
