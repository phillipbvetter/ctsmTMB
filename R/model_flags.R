
set_flags = function(proc, args, self, private){

  # Procedure is "estimation", "filtration", "simulation" etc...
  private$set_procedure(proc)

  # EKF, LKF, UKF, Laplace, Laplace Thygesen
  private$set_method(args$method)

  # The inital state value is the prior state for Kalman filtering, or the initial guess for all states in the Laplace methods
  private$set_initial_state(args$initial.state)
  private$set_initial_state_estimation(args$estimate.initial.state)

  # ODE / SDE settings
  private$set_ode_solver(args$ode.solver)
  private$set_first_order_hold(args$first.order.input.hold)
  private$set_timestep("ode", args$ode.timestep)
  # TODO
  # We only need this for simulate, right?
  sim.dt <- if (proc == "simulation") args$simulation.timestep else args$ode.timestep
  private$set_timestep("simulation", sim.dt)

  # Utilities
  private$set_silence(args$silent)

  # likelihood
  # ----------------------------------- #
  if (proc == "likelihood") {
    private$set_compile(args$compile)
  }

  # estimate
  # ----------------------------------- #
  if (proc == "estimation") {
    private$set_compile(args$compile)
    private$set_control(args$control)
    private$use_hessian(args$use.hessian)
    private$set_unconstrained_optim(args$unconstrained.optim)
  }

  # simulate
  # ----------------------------------- #
  if (proc == "simulation") {
    private$set_cpp_seed(args$cpp.seed)
  }

  # not smoother
  # ----------------------------------- #
  if (proc != "smoother") {
    private$set_ukf_hyperpars(args$ukf.hyperpars)
  }


}
