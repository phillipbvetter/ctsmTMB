
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
