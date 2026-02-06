
# Hacks to solve issues with R CMD check

# 1. Prevent 'no visile binding for global variables' for functions defined
# by eval(parse(text=function_as_string) type constructions used in RTMB

# f__ <- function(stateVec, parVec, inputVec) NULL
# dfdx__ <- function(stateVec, parVec, inputVec) NULL
# g__ <- function(stateVec, parVec, inputVec) NULL
# dhdx__ <- function(stateVec, parVec, inputVec) NULL
# dfdu__ <- function(stateVec, parVec, inputVec) NULL
# h__ <- function(stateVec, parVec, inputVec) NULL
# hvar__ <- function(stateVec, parVec, inputVec) NULL
# hvar__matrix <- function(stateVec, parVec, inputVec) NULL
# x <- NULL
# level <- NULL

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
