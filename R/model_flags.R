
set_flags = function(proc, args, self, private){

  private$set_procedure(proc)

  if(private$procedure == "likelihood"){
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_ode_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_compile(args$compile)
    private$set_initial_state(args$initial.state)
    private$set_first_order_hold(args$first.order.input.hold)
  }

  if(private$procedure == "estimation"){
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_ode_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_control(args$control)
    private$use_hessian(args$use.hessian)
    private$set_unconstrained_optim(args$unconstrained.optim)
    private$set_silence(args$silent)
    private$set_compile(args$compile)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_initial_state(args$initial.state)
    private$set_first_order_hold(args$first.order.input.hold)
  }

  if(private$procedure == "filtration"){
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_ode_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_initial_state(args$initial.state)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    private$set_first_order_hold(args$first.order.input.hold)
  }

  if(private$procedure == "prediction"){
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_ode_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_initial_state(args$initial.state)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    private$set_first_order_hold(args$first.order.input.hold)
  }

  if(private$procedure == "simulation"){
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_ode_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$simulation.timestep)
    private$set_silence(args$silent)
    private$set_ukf_hyperpars(args$ukf.hyperpars)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_cpp_seed(args$cpp.seed)
    private$set_initial_state(args$initial.state)
    private$set_first_order_hold(args$first.order.input.hold)
  }

  if(private$procedure == "smoothing"){
    private$set_method(args$method)
    private$set_ode_solver(args$ode.solver)
    private$set_ode_timestep(args$ode.timestep)
    private$set_simulation_timestep(args$ode.timestep)
    private$set_silence(args$silent)
    private$set_initial_state_estimation(args$estimate.initial.state)
    private$set_initial_state(args$initial.state)
    # private$set_first_order_hold(args$first.order.input.hold)
  }

}
