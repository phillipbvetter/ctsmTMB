
set_flags = function(proc, args, self, private){
  
  # estimation, prediction, simulation, construction
  private$set_procedure(proc)
  
  if(private$procedure == "estimation"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_loss(args$loss, args$loss_c)
    private$set_ukf_hyperpars(args$unscented_hyperpars)
    private$set_control(args$control)
    private$use_hessian(args$use.hessian)
    private$set_unconstrained_optim(args$unconstrained.optim)
    private$set_compile(args$compile)
    private$set_silence(args$silent)
    # private$estimate.initial = args$estimate.initial.state
    
  }
  
  if(private$procedure == "construction"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)  
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_loss(args$loss, args$loss_c)
    private$set_ukf_hyperpars(args$unscented_hyperpars)
    private$set_compile(args$compile)
    private$set_silence(args$silent)
    # private$estimate.initial = estimate.initial.state
    
  }
  
  if(private$procedure == "prediction"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_pred_initial_state(args$initial.state)
    private$set_ukf_hyperpars(args$unscented_hyperpars)
    private$set_silence(args$silent)
    
  }
  
  if(private$procedure == "simulation"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_pred_initial_state(args$initial.state)
    private$set_ukf_hyperpars(args$unscented_hyperpars)
    private$set_silence(args$silent)
    
  }
  
}

temporary_save_old_data <- function(self, private){
  
  newobj <- self$clone(deep=TRUE)
  priv <- newobj$.__enclos_env__$private
  
  private$old.data = list()
  private$old.data$data = priv$data
  private$old.data$method = priv$method
  private$old.data$ode.solver = priv$ode.solver
  private$old.data$ode.timestep = priv$ode.timestep
  
  # private$simulation.timestep.size
  # private$simulation.timesteps
  
  return(invisible(self))
}

check_for_data_rebuild <- function(self, private){
  
  bool <- any(c(
    private$rebuild.data,
    !identical(private$old.data$data, private$data),
    !identical(private$old.data$method, private$method),
    !identical(private$old.data$ode.solver, private$ode.solver),
    !identical(private$old.data$ode.timestep, private$ode.timestep)
  ))
  private$rebuild.data = bool
  
  return(invisible(self))
}
