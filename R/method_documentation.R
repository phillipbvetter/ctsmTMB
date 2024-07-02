# This document contains empty functions with names similar to the methods in sdeTMB
# This allow users to easily find method documentation using ?add_systems, instead of having
# to manually retrieve them from ?sdeTMB and then finding the method among the large list of other methods


########################################################################
# ADD SYSTEMS
########################################################################
#' @title Add state equations to model object
#' 
#' @description 
#' Add stochastic differential equation to the \code{sdeTMB} model-object that 
#' governs the differential evolution of states in the specified model.
#' 
#' @usage
#' add_systems(form,
#'             ...)
#' 
#' @examples 
#' # Example 1 - Linear System
#' add_systems(dx ~ theta * (mu - x + u) * dt + sigma * dw)
#' 
#' # Example 2 - Non-Linear System
#' add_systems(dx ~ theta * (mu - exp(x)^2 + u) * dt + sigma * x * (1-x) * dw)
#' 
#' @details
#' ## Usable functions
#' The formulas can contain most elementary functions such as \code{log}, 
#' \code{exp}, \code{logit} and \code{invlogit}. In general the supported
#' functions are only those that are both 1) defined in the derivative table 
#' of the \code{Deriv} package, and 2) undestood by *TMB* on the *C++* side.
#' 
#'              
#' @param form formula specifying the stochastic differential equation to be 
#' added to the system.
#' @param ... formulas similar to \code{form}, used to allow specifying
#' multiple formulas by comma-seperation rather than providing a list.
#' 
add_systems = function(form,...) {
  NULL
}

########################################################################
# ADD OBSERVATIONS
########################################################################
#' @title Add observation equations to model object
#' 
#' @description
#' Add an observation equation to the \code{sdeTMB} model-object that links 
#' states to an observed variable.
#'  
#' @examples
#' # Example 1
#' add_observations(y ~ x)
#' 
#' # Example 2
#' Specify observation as the sum of exponentials of two latent states
#' add_observations(y ~ exp(x1) + exp(x2))
#' 
#' # Example 3
#' Specify observation as the sum of exponentials of two latent states
#' add_observations(log(y) ~ x, obsnames = "logy")
#' 
#' @usage
#' add_observations(form,
#'                  ...,
#'                  obsnames=NULL)
#' 
#' @param form formula specifying the observation equation to be added to the 
#' system.
#' 
#' @param ... formulas similar to \code{form}, used to allow specifying
#' multiple formulas by comma-seperation rather than providing a list.
#' 
#' @param obsnames character vector specifying observation names only used when 
#' when the observation left-hand side is a function call. See details.
#' 
#' @details
#' ## \code{obsnames}
#' The \code{obsnames} argument is used when the left-hand side of \code{form}
#' is a function of a variable i.e. \code{log(y)} (when its of 
#' class 'call' instead of 'name'), as in Example 3. The user should then only 
#' provide data for \code{y}, and the log-transformation
#' is then handled internally.
#' 
#' The supported functions are those discussed in the \code{\link{add_systems}}.
#' 
add_observations = function(form,...,obsnames=NULL) {
  NULL
}

########################################################################
# ADD OBSERVATION VARIANCES
########################################################################
#' @title Add observation variances to the model object.
#' 
#' @description
#' Specify the observation variance of an existing observation equation.
#' 
#' @usage
#' add_observation_variances(form,
#'                           ...)
#' 
#' @examples
#' # Example 1
#' add_observation_variances(y ~ sigma_y^2). 
#'
#' # Example 2 
#' add_observation_variances(y ~ 0.1 + x * exp(logsigma_y)^2). 
#' 
#' @param form formula class specifying the obsevation equation to be added 
#' to the system.
#' @param ... formulas similar to \code{form}, used to allow specifying
#' multiple formulas by comma-seperation rather than providing a list.
#' 
add_observation_variances = function(form,...) {
  NULL
}

########################################################################
# ADD INPUTS
########################################################################
#' @title Specify input variables in the model object.
#'
#' @description Declare whether a variable contained in system, observation or observation 
#' variance equations is an input variable.
#' 
#' 
#' @usage 
#' add_inputs(...)
#' 
#' @examples
#' # Example 1
#' add_inputs(u)
#' 
#' # Example 2
#' add_inputs(u1, u2, u3)
#' 
#' @param ... a series of variable names (unquouted) that match variable names
#' in the defined system which should be treated as input variables.
#' 
add_inputs =  function(...) {
  NULL
}

########################################################################
# ADD PARAMETERS
########################################################################
#' @title Specify parameters in the model object
#' 
#' @description Declare which variables that are (fixed effects) parameters in
#' the specified model, and specify the initial optimizer values, as well as
#' lower / upper bounds. Parameters can be declared either as vectors or as
#' matrices. The first entry is the initial value, the second entry is the lower
#' bound and the third entry is the upper bound. Providing only a first entry
#' fixes the particular parameter to that value.
#' 
#' @examples
#' # Example 1
#' add_parameters(
#' alpha = c(initial=0, lower=-10, upper=10),
#' beta = c(2, 0, 5),
#' gamma = 5)
#' )
#' 
#' # Example 2
#' parmat = matrix(rep(c(0,-10,10),times=3),ncol=3,byrow=T)
#' rownames(parmat) = c("a","b","c")
#' colnames(parmat) = c("initial","lower","upper")
#' add_parameters(parmat)
#'
#' @param ... a comma-seperated series of vectors/matrix entries
add_parameters = function(...) {
  NULL
}

########################################################################
# ADD ALGEBRAICS
########################################################################
#' @title Add algebraic relationships to the model object.
#'
#' @description
#' Algebraic relations is a convenient way to transform parameters in your equations,
#' to reduce clutter when specying the various equations, for instance to ensure
#' positivity (log-transform). 
#' 
#' @examples
#' # Example 1
#' add_algebraics(
#' sigma ~ exp(logsigma),
#' theta ~ invlogit(alpha + beta)
#' )
#' 
#' @details
#' The left-hand side of the provided formula specifies which parameter that
#' should be overwritten with the expression on the right-hand side. This also
#' means that the left-hand side parameter will vanish from the model formulation
#' and \code{link{add_parameters}} should therefore specify values for the 
#' new parameters.
#' 
#' @usage
#' add_algebraics(form,
#'                ...)
#' 
#' @param form formula specifying algebraic relation.
#' @param ... formulas similar to \code{form}, used to allow specifying
#' multiple formulas by comma-seperation rather than providing a list.
#' 
add_algebraics = function(form,...) {
  NULL
}

########################################################################
# SET INITIAL STATE
########################################################################
#' Set initial state mean and covariance 
#'
#' @description Declare the initial state values i.e. mean and covariance for the system states.
#' 
#' @param initial.state a named list of two entries 'x0' and 'p0' containing the initial state and covariance of the state
#' @param estimate boolean value which indicates whether or not the initial conditions
#' shall be estimated as fixed effects parameters. The provided mean and covariance are then
#' used as initial guesses
#' 
set_initial_state = function(initial.state, estimate=FALSE) {
  NULL
}

