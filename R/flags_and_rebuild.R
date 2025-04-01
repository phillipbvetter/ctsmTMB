
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
    
  }
  
}

temporary_save_old_data <- function(self, private){
  
  newobj <- self$clone(deep=TRUE)
  private <- self$.private()
  
  private$old.data = list()
  private$old.data$data = private$data
  private$old.data$method = private$method
  private$old.data$ode.solver = private$ode.solver
  private$old.data$ode.timestep = private$ode.timestep
  private$old.data$loss$loss = private$loss$loss
  private$old.data$estimate.initial = private$estimate.initial
  
  # private$simulation.timestep.size
  # private$simulation.timesteps
  
  return(invisible(self))
}

check_for_data_rebuild <- function(self, private){
  
  bool <- c(
    private$rebuild.data,
    !identical(private$old.data$data, private$data),
    !identical(private$old.data$method, private$method),
    !identical(private$old.data$ode.solver, private$ode.solver),
    !identical(private$old.data$ode.timestep, private$ode.timestep),
    !identical(private$old.data$loss$loss, private$loss$loss),
    !identical(private$old.data$estimate.initial, private$estimate.initial)
  )
  
  private$rebuild.data = any(bool)
  
  return(invisible(self))
}
