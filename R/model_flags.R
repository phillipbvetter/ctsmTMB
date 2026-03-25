
set_flags = function(proc, args, self, private){

  private$set_procedure(proc)

  # Common core — applies to every procedure
  private$set_method(args$method)
  private$set_ode_solver(args$ode.solver)
  private$set_timestep("ode", args$ode.timestep)
  private$set_silence(args$silent)
  private$set_initial_state_estimation(args$estimate.initial.state)
  private$set_initial_state(args$initial.state)

  # All procedures except smoother support UKF hyperparameters
  if (proc != "smoother") {
    private$set_ukf_hyperpars(args$ukf.hyperpars)
  }

  # Simulation timestep: simulation uses its own arg; all others mirror ode.timestep
  sim_dt <- if (proc == "simulation") args$simulation.timestep else args$ode.timestep
  private$set_timestep("simulation", sim_dt)

  # Procedure-specific extras
  if (proc %in% c("likelihood", "estimation")) {
    private$set_compile(args$compile)
  }
  if (proc == "estimation") {
    private$set_control(args$control)
    private$use_hessian(args$use.hessian)
    private$set_unconstrained_optim(args$unconstrained.optim)
  }
  if (proc == "simulation") {
    private$set_cpp_seed(args$cpp.seed)
  }

}
