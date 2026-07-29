# This function is the main caller for carrying out simulations
perform_simulation <- function(self, private, use.cpp, n.sims){


  if(private$algo.settings$method=="laplace"){
    stop("Simulations are not available for the Laplace method")
  }

  if(private$algo.settings$method=="laplace.thygesen"){
    stop("Simulations are not available for the Laplace-Thygesen method")
  }

  # time prediction
  comptime <- system.time(
    {
      # Predict with Rcpp implementation
      if(use.cpp){

        if(!private$algo.settings$silent) message("Simulating with C++...")
        lkf_ekf_ukf_simulate_rcpp(private$algo.settings$argument.parameters, self, private, n.sims)

        # Predict with R implementation
      } else {

        if(!private$algo.settings$silent) message("Simulating with R...")
        ekf_lkf_ukf_simulate_r(private$algo.settings$argument.parameters, self, private, n.sims)

      }

    }, gcFirst = FALSE)

  private$timers$simulation <- comptime

  return(invisible(self))
}

lkf_ekf_ukf_simulate_rcpp <- function(pars, self, private, n.sims){

  # check for null pointers
  check_if_rcpp_pointers_are_valid(self, private)

  # observation/input matrix
  obsMat = as.matrix(private$data[private$names$obs])
  inputMat = as.matrix(private$data[private$names$inputs])

  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)

  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)

  ids <- 1:private$dims$observations - 1 #minus 1 for 0 indexing
  non_na_ids <- apply(obsMat, 1, function(x) ids[!is.na(x)], simplify = FALSE)
  any_available_obs <- sapply(non_na_ids, function(x) length(x) != 0)

  output <- NULL
  if(private$algo.settings$method=="lkf"){
    output <- lkf_simulate_rcpp(private$model$rcpp_function_ptr,
                                obsMat,
                                inputMat,
                                pars,
                                private$algo.settings$initial.state$p0,
                                private$algo.settings$initial.state$x0,
                                private$algo.settings$ode.timestep.size,
                                private$algo.settings$ode.timesteps,
                                private$algo.settings$simulation.timestep.size,
                                private$algo.settings$simulation.timesteps,
                                any_available_obs,
                                non_na_ids,
                                private$dims$diffusions,
                                private$algo.settings$last.pred.index,
                                private$algo.settings$k.ahead,
                                n.sims,
                                private$algo.settings$seed$state.seed,
                                private$algo.settings$first.order.input.hold)
  }
  if(private$algo.settings$method=="ekf"){
    output <- ekf_simulate_rcpp(private$model$rcpp_function_ptr,
                                obsMat,
                                inputMat,
                                pars,
                                private$algo.settings$initial.state$p0,
                                private$algo.settings$initial.state$x0,
                                private$algo.settings$ode.timestep.size,
                                private$algo.settings$ode.timesteps,
                                private$algo.settings$simulation.timestep.size,
                                private$algo.settings$simulation.timesteps,
                                any_available_obs,
                                non_na_ids,
                                private$algo.settings$ode.solver,
                                private$algo.settings$last.pred.index,
                                private$algo.settings$k.ahead,
                                private$dims$diffusions,
                                n.sims,
                                private$algo.settings$seed$state.seed,
                                private$algo.settings$first.order.input.hold)
  }
  if(private$algo.settings$method == "ukf"){
    output <- ukf_simulate_rcpp(private$model$rcpp_function_ptr,
                                obsMat,
                                inputMat,
                                pars,
                                private$algo.settings$initial.state$p0,
                                private$algo.settings$initial.state$x0,
                                private$algo.settings$ode.timestep.size,
                                private$algo.settings$ode.timesteps,
                                private$algo.settings$simulation.timestep.size,
                                private$algo.settings$simulation.timesteps,
                                numeric_is_not_na_obsMat,
                                number_of_available_obs,
                                private$algo.settings$ukf.hyperpars,
                                private$dims$diffusions,
                                private$algo.settings$last.pred.index,
                                private$algo.settings$k.ahead,
                                private$algo.settings$ode.solver,
                                n.sims,
                                private$algo.settings$seed$state.seed,
                                private$algo.settings$first.order.input.hold)
  }

  private$results$simulation.raw <- output

  return(invisible(NULL))
}


# This function returns simulation results to the user
create_return_simulation <- function(return.k.ahead, n.sims, self, private){

  if(!private$algo.settings$silent) message("Returning results...")

  # create names for inner list
  inner.names <- paste0("i", 0:(private$algo.settings$last.pred.index-1))

  # Build returnlist for states
  state.list <- build_simulation_returnlist(private$results$simulation.raw,
                                            private$data$t,
                                            private$dims$states,
                                            private$algo.settings$k.ahead,
                                            n.sims)
  # set state names
  names(state.list) <- private$names$states
  for(i in seq_along(state.list)){
    names(state.list[[i]]) <- inner.names
  }

  # create index/time list
  time.list <- build_simulation_timelists(private$data$t,
                                          private$algo.settings$last.pred.index,
                                          private$algo.settings$k.ahead)
  names(time.list) <- inner.names

  # Build returnlist for observations
  # First we must calculate the simulated observation trajectories
  simulation.raw.obs <- calculate_simulation_observations(
    private$results$simulation.raw,
    private$model$rcpp_function_ptr,
    t(as.matrix(private$data[private$names$inputs])),
    private$algo.settings$argument.parameters,
    private$dims$states,
    private$dims$observations,
    private$algo.settings$k.ahead,
    n.sims,
    private$algo.settings$seed$obs.seed
  )

  # # Now we can build the returnlist
  obs.list <- build_simulation_returnlist(simulation.raw.obs,
                                          private$data$t,
                                          private$dims$observations,
                                          private$algo.settings$k.ahead,
                                          n.sims)
  names(obs.list) <- private$names$obs
  for(i in seq_along(obs.list)){
    names(obs.list[[i]]) <- inner.names
  }

  private$results$simulation <- list(states=state.list, observations=obs.list, times=time.list)

  return(invisible(self))
}
