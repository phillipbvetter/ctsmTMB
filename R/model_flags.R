
set_flags = function(proc, args, self, private){
  
  private$set_procedure(proc)
  
  if(private$procedure == "construction"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)  
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_compile(args$compile)
    private$set_initial_state(args$initial.state)
    
  }
  
  if(private$procedure == "estimation"){

    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_control(args$control)
    private$use_hessian(args$use.hessian)
    private$set_unconstrained_optim(args$unconstrained.optim)
    private$set_silence(args$silent)
    private$set_compile(args$compile)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_initial_state(args$initial.state)
    
  }
  
  if(private$procedure == "filtration"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_initial_state(args$initial.state)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    
  }
  
  if(private$procedure == "smoother"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_initial_state(args$initial.state)
    
  }
  
  if(private$procedure == "prediction"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_initial_state(args$initial.state)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    
  }
  
  if(private$procedure == "simulation"){
    
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_cpp_seed(args$cpp.seed)
    private$set_initial_state(args$initial.state)
    
  }
  
}
