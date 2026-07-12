
# Hacks to solve issues with R CMD check
# Prevent 'no visible binding for global variables'

utils::globalVariables(
  c(
    "f__",
    "dfdx__",
    "g__",
    "dhdx__",
    "dfdu__",
    "hvar__matrix",
    "hvar__",
    "h__",
    "h.sigma"
  )
)

utils::globalVariables(
  c(
    "f.initial.state.newton",
    "f.initial.covar.solve",
    "kalman.data.update",
    "kalman.data.update.with.nll",
    "kalman.no.update.with.nll",
    "ode.integrator",
    "euler.maruyama.simulation",
    "kron.left",
    "kron.right",
    "logdet",
    "loss.function",
    "sigma2chol",
    "create.sigmaPoints"
  )
)

utils::globalVariables(
  c(
    "private",
    "n.states",
    "n.pars",
    "n.inputs",
    "n.obs",
    "n.diffusions",
    "n.zeros",
    "n.sigmapoints",
    "nn",
    "nsims",
    "sqrt_c",
    "W.m",
    "W"
  )
)

utils::globalVariables(
  c(
    "x",
    "level",
    "sd"
  )
)

# C++ function from Rcpp::sourceCpp()
utils::globalVariables(
  "get_sysfun_cpp_function_ptrs"
)
