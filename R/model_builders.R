
build_model = function(self, private) {
  
  # Check if model is already built
  if(!private$rebuild.model) return(invisible(self))
  private$rebuild.model <- FALSE
  private$rebuild.data <- TRUE
  private$rebuild.ad <- TRUE
  
  # Print
  if(!private$silent) message("Checking model components...")
  
  # check_model
  basic_model_check(self, private)
  
  # last check
  final_build_check(self, private)
  
  # return
  return(invisible(self))
}

#######################################################
# FIRST FUNCTION TO RUN WHEN BUILDING
#######################################################

basic_model_check = function(self, private) {
  
  # system eqs
  if(private$number.of.states == 0) {
    stop("There were no specified system equations - use 'addSystem'.")
  }
  
  # obs eqs
  if(private$number.of.observations == 0) {
    stop("There were no specified observation equations - use 'addObs'.")
  }
  
  # obs var
  missing.var = sapply(private$obs.var, length) == 0
  if (any(missing.var)) {
    missing.names = paste(private$obs.names[missing.var], collapse=", ")
    stop("There are no observation variances specified for the observation(s): \n\t", missing.names)
  }
  
  # parameters
  if(private$number.of.pars == 0) {
    stop("There are no parameters in the model.")
  }
  
  return(invisible(self))
}

#######################################################
# LAST CHECK BEFORE COMPILING
#######################################################

# We perform a series of basics checks to see if the model satisfies assumings e.g.
# 1. The observation equation must have at least one state on the rhs
# 2. There are no observations on the rhs of an observation equation
# 3. We check that all variables in the model are declared as states, inputs or parameters


final_build_check = function(self, private) {
  
  # Verify that all observations relate to a state (have a state on their rhs)
  bool = unlist(lapply(private$obs.eqs.trans, function(x) any(private$state.names %in% x$allvars)))
  if (any(!bool)) {
    stop("Error: There are no states on the right-hand side of the following observation equation(s) : \n\t ",paste(private$obs.names[!bool],collapse=", "))
  }
  
  # Verify that observations dont have other observations on their rhs
  bool = unlist(lapply(private$obs.eqs.trans, function(x) any(private$obs.names %in% x$allvars)))
  if (any(bool)) {
    stop("The following observation(s) attempt to observe other observations! \n\t ",paste(private$obs.names[!bool],collapse=", "))
  }
  
  # Verify that all variables on the rhs of system, obs and obs.variance equations
  # have been provided. They can either be inputs, parameters and states
  vars = list()
  vars[[1]] = unlist(lapply(private$sys.eqs.trans, function(x) x$allvars))
  vars[[2]] = unlist(lapply(private$obs.eqs.trans, function(x) x$allvars))
  vars[[3]] = unlist(lapply(private$obs.var.trans, function(x) x$allvars))
  rhs.vars = unique(unlist(vars))
  given.vars = c("pi", private$parameter.names, private$input.names, private$state.names)
  bool = rhs.vars %in% given.vars
  if (any(!bool)) {
    stop("Error: The following variables(s) in the model have not been declared as parameters, inputs or states: \n\t ",
         paste(rhs.vars[!bool],collapse=", "))
  }
  
  ##### NOTE::: Do we need to remove this? doesnt really matter right? ####
  # Verify that all input and parameters are used in the model
  # given.vars = c( private$parameter.names, private$input.names[-1])
  # bool = given.vars %in% rhs.vars
  # if (any(!bool)) {
  #   stop("The following variables(s) are unused: \n\t ", paste(given.vars[!bool],collapse=", "))
  # }
  
  # return
  return(invisible(self))
}
