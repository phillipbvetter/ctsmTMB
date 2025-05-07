
set_flags = function(proc, args, self, private){
  
  # estimation, prediction, simulation, construction
  private$set_procedure(proc)
  
  if(private$procedure == "estimation"){

    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_control(args$control)
    private$use_hessian(args$use.hessian)
    private$set_unconstrained_optim(args$unconstrained.optim)
    private$set_silence(args$silent)
    # private$set_loss(args$loss, args$loss_c)
    # private$set_ukf_hyperpars(args$unscented_hyperpars)
    private$set_initial_state_estimation(args$estimate.initial.state)
    
  }
  
  if(private$procedure == "construction"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)  
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    # private$set_loss(args$loss, args$loss_c)
    # private$set_ukf_hyperpars(args$unscented_hyperpars)
    private$set_initial_state_estimation(args$estimate.initial.state)
    
  }
  
  if(private$procedure == "prediction"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_pred_initial_state(args$initial.state)
    # private$set_ukf_hyperpars(args$unscented_hyperpars)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    # private$set_k_ahead(args$k.ahead)
    
  }
  
  if(private$procedure == "simulation"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_pred_initial_state(args$initial.state)
    # private$set_ukf_hyperpars(args$unscented_hyperpars)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    # private$set_k_ahead(args$k.ahead)
    private$set_cpp_seed(args$cpp.seed)
    
  }
  
  if(private$procedure == "filter"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    
  }
  
  if(private$procedure == "smoother"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    
  }
  
}

save_settings_for_comparison_next_time <- function(self, private){
  
  cloned.self <- self$clone(deep=TRUE)
  cloned.private <- cloned.self$.private()
  
  # private$old.data = list()
  # private$old.data$data = private$data
  private$old.data$method = cloned.private$method
  private$old.data$ode.solver = cloned.private$ode.solver
  private$old.data$ode.timestep = cloned.private$ode.timestep
  private$old.data$loss = cloned.private$loss
  private$old.data$estimate.initial = cloned.private$estimate.initial
  
  # private$simulation.timestep.size
  # private$simulation.timesteps
  
  return(invisible(self))
}

check_for_ADfun_rebuild <- function(self, private){
  
  # print(str(private$old.data))
  
  # We perform checks against the old data on the entries that would
  # require a new call to RTMB::MakeADFun (i.e. entries that affect the
  # calculations in the likelihood function)
  bool <- c(
    private$rebuild.ad,
    !identical(private$old.data$method, private$method),
    !identical(private$old.data$ode.solver, private$ode.solver),
    !identical(private$old.data$ode.timestep, private$ode.timestep),
    !identical(private$old.data$loss, private$loss),
    !identical(private$old.data$estimate.initial, private$estimate.initial)
  )
  
  private$rebuild.ad <- any(bool)
  
  return(invisible(self))
}

check_for_data_rebuild <- function(data, self, private){
  
  # We check if the data needs to reset
  bool <- c(
    private$rebuild.data,
    !identical(private$old.data$entry.data, data)
  )
  
  private$rebuild.data <- any(bool)
  
  return(invisible(self))
}
