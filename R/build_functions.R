# These functions are helper functions used when calling the ctsmrTMB method
# 'build_model'.

#######################################################
# MAIN BUILDING FUNCTION THAT CALLS ALL OTHER FUNCTIONS
#######################################################

build_model = function(self, private) {
  
  # check_model
  basic_model_check(self, private)
  
  # set dimensions, diff processes, etc...
  set_model_settings(self, private)
  
  # 1) apply algebraics, and define new trans_system
  # 2) calculate new diff terms
  apply_algebraics_and_define_trans_equations(self, private)
  calculate_diff_terms(self, private)
  
  # apply lamperti and update diff terms
  apply_lamperti(self, private) 
  calculate_diff_terms(self, private)
  
  # create rtmb functions
  create_rtmb_function_strings(self, private)
  create_rtmb_laplace_string(self, private)
  create_rtmb_ekf_string(self, private)
  
  # rcpp
  create_rcpp_function_strings(self, private)
  
  # last check
  final_build_check(self, private)
  
  # set initial state
  if(any(private$procedure == c("estimation","construction"))){
    if(is.null(private$estimate.initial)){
      self$setInitialState(list(rep(0,private$number.of.states), 1*diag(private$number.of.states)), TRUE)
    }
  }
  
  # return
  return(invisible(self))
}

#######################################################
# FIRST FUNCTION TO RUN WHEN BUILDING
#######################################################

basic_model_check = function(self, private) {
  
  # system eqs
  if (length(private$sys.eqs) == 0) {
    stop("There were no specified system equations - use 'addSystem'")
  }
  
  # obs eqs
  if (length(private$obs.eqs) == 0) {
    stop("There were no specified observation equations - use 'addObs'")
  }
  
  # obs var
  missing.var = sapply(private$obs.var, length) == 0
  if (any(missing.var)) {
    missing.names = paste(private$obs.names[missing.var],collapse=", ")
    stop("There are no observation variances specified for the observation(s): \n\t", missing.names)
  }
  
  # parameters
  if(is.null(self$getParameters())){
    stop("No parameters were specified.")
  }
  
  # initial state for estimation
  # The same is handled in set_pred_initial_state for prediction/simulation
  # if(any(private$procedure == c("estimation","construction"))){
  #   if (is.null(private$initial.state)) {
  #     stop("You must set an initial state estimate and covariance")
  #   }
  # }
  
  return(invisible(self))
}

#######################################################
# SET SYSTEM VARIABLES
#######################################################
set_model_settings = function(self, private){
  
  # set system size variables
  private$diff.processes = unique(unlist(lapply(private$sys.eqs, function(x) x$diff)))
  private$number.of.states = length(private$sys.eqs)
  private$number.of.observations = length(private$obs.eqs)
  private$number.of.diffusions =  length(private$diff.processes) - 1 # minus 1 to remove 'dt'
  private$number.of.pars = length(private$parameters)
  private$number.of.inputs = length(private$inputs)
  
  
  # update free parameter list contains: initial, lower, upper for each parameter
  private$free.pars = private$parameters[!(private$parameter.names %in% names(private$fixed.pars))]
  
  return(invisible(self))
}

#######################################################
# APPLY ALGEBRAIC TRANSFORMATIONS
#######################################################

apply_algebraics_and_define_trans_equations = function(self, private) {
  
  # extract rhs's
  sys.rhs = lapply(private$sys.eqs,function(x) x$rhs)
  obs.rhs = lapply(private$obs.eqs,function(x) x$rhs)
  obs.var.rhs = lapply(private$obs.var,function(x) x$rhs)
  alg.rhs = lapply(private$alg.eqs, function(x) x$rhs)
  
  # apply algebraics with substitute
  sys.rhs = lapply(sys.rhs, function(x) do.call(substitute, list(x,alg.rhs)))
  obs.rhs = lapply(obs.rhs, function(x) do.call(substitute, list(x,alg.rhs)))
  obs.var.rhs = lapply(obs.var.rhs, function(x) do.call(substitute, list(x,alg.rhs)))
  
  # replace rhs in the already defined system, obs, obs.var equations and
  # add these to the transformed systems/obs/obs.var private fields
  
  # system
  for(i in seq_along(private$sys.eqs)) {
    temp.form = private$sys.eqs[[i]]$form
    temp.form[[3]] = sys.rhs[[i]] #[[3]] is rhs of a formula
    temp.list = list(form=temp.form, name=private$state.names[i])
    private$add_trans_systems(temp.list)
  }
  
  # observations
  for(i in seq_along(private$obs.eqs)) {
    temp.form = private$obs.eqs[[i]]$form
    temp.form[[3]] = obs.rhs[[i]] #[[3]] is rhs of a formula
    temp.list = list(form=temp.form, name=private$obs.names[i])
    private$add_trans_observations(temp.list)
  }
  
  # observation variances
  for(i in seq_along(private$obs.var)) {
    temp.form = private$obs.var[[i]]$form
    temp.form[[3]] = obs.var.rhs[[i]] #[[3]] is rhs of a formula
    temp.list = list(form=temp.form, name=private$obs.names[i])
    private$add_trans_observation_variances(temp.list)
  }
  
  return(invisible(self))
}

#######################################################
# UPDATE DIFF TERMS TO FIT ALGEBRAIC EQUATIONS
#######################################################

calculate_diff_terms = function(self, private) {
  
  # SYSTEM EQUATIONS
  # Calculate the differential terms in front of 'dt' and diffusion 'dw...'
  for (i in seq_along(private$sys.eqs.trans)) {
    
    # must use cache.exp=FALSE to prevent renaming variables
    private$diff.terms[[i]] = lapply(private$diff.processes, 
                                     function(x) Deriv::Deriv(private$sys.eqs.trans[[i]]$rhs, x, cache.exp=FALSE))
    names(private$diff.terms[[i]]) = private$diff.processes
  }
  names(private$diff.terms) = private$state.names
  
  # OBSERVATION EQUATIONS WRT STATES
  for(i in seq_along(private$obs.eqs.trans)){
    private$diff.terms.obs[[i]] = lapply(private$state.names,
                                         function(x) Deriv::Deriv(private$obs.eqs.trans[[i]]$rhs, x, cache.exp=F))
    names(private$diff.terms.obs[[i]]) = private$state.names
  }
  names(private$diff.terms.obs) = private$obs.names
}

#######################################################
# APPLY LAMPERTI TRANSFORM IF SPECIFIED
#######################################################

# This function applies the set lampeti transformation to each state of the
# system, but only states that have 1 diffusion process.
#
# Note: The observation equation is not transformed - should we also do that?
# i.e. the states that are present in the observations should be transformed
apply_lamperti = function(self, private) {
  
  ##### Extract ####
  # get the states to transform, and the transformations to use on those states
  transforms =  private$lamperti$transforms
  states = private$lamperti$states
  
  if(all(transforms == "identity")){
    return(invisible(self))
  }
  
  ##### Initial Filtering ####
  
  # Remove states with multiple diffusion terms (lamperti only works in 1D)
  # check how many diff.terms is non-zero. Each non-zero is a diffusion process
  nonzero.diffterms = lapply(private$diff.terms, 
                             function(x) unlist(lapply(x, function(y) y==0))) #state names are retained
  number.of.diffterms = lapply(nonzero.diffterms, sum)
  diffusion.id = numeric(private$number.of.states)
  bool = transforms != "identity" #this removes identity transforms altogether
  
  for(i in seq_along(states)){
    #remove states with more than one dw process (2 because dt and 1 dw process)
    if(number.of.diffterms[[states[i]]] > 2.5){
      warning("The lamperti transformation on ", states[i], " was aborted, because there are multiple diffusion terms")
      bool[i] = FALSE
    } else {
      # which dw diffusion process was active?
      diffusion.id[i] = which(!nonzero.diffterms[[states[i]]])[2]
      
    }
  }
  states = states[bool]
  transforms = transforms[bool]
  diffusion.id = diffusion.id[bool]
  
  ##### List of Available Transformations ####
  # Define and select lamperti transform and 1st and 2nd derivative: list( x(z) , dzdz(x) , dz2/dx2(x) )
  # Note that the first entry is in terms of z, not x. We must therefore
  # substitute x(z) into dzdx(x) and d2z/dx2(x) to get the expression in terms of the new state variable z.
  psi.all = list(
    "log" = list( quote(exp(x)), 
                  quote(1/x), 
                  quote(-1/x^2)
    ),
    # 
    "logit" = list( 
      quote(exp(x/(1+x))), 
      quote(1/(x*(1-x))), 
      quote((2*x-1)/(x^2*(x-1)^2))
    ),
    # 
    "sqrt-logit" = list( 
      quote(0.5*(sin(x)+1)), 
      quote(1/(sqrt(x*(1-x)))), 
      quote(0.5*(2*x-1)/(x*(1-x)*sqrt(x*(1-x))))
    )
  )
  # name each list in psi.all
  for(i in seq_along(psi.all)){ 
    names(psi.all[[i]]) = c("..psi..","..dpsidx..","..d2psidx2..") 
  }
  
  
  ############### MAIN LOOP FOR STATES ###############
  # Perform lamperti transformation with substitutions
  for(i in seq_along(states)){
    
    # select current state and transformation
    state = states[i]
    transform = transforms[i]
    
    # get the current transformation
    psi = psi.all[[transform]]
    
    # substitute the state name into the psi_transformation list instead of the placeholder 'x'
    templist = list(x=parse(text=state)[[1]])
    psi.correct.state = lapply(psi, function(x) do.call(substitute, list(x, templist)))
    
    # get the single diffusion process used in the current system equation
    dw.list = list(dw=parse(text=private$diff.processes[diffusion.id[i]])[[1]])
    
    # get drift and diffusion
    fg.list = list(
      f = private$diff.terms[[state]][["dt"]],
      g = private$diff.terms[[state]][[dw.list$dw]]
    )
    
    # apply the lamperti substitution
    lamperti.formula = quote((f * ..dpsidx.. + 0.5 * g^2 * ..d2psidx2.. ) * dt + g * ..dpsidx.. * dw)
    substitute.list = c(psi.correct.state[-1], dw.list, fg.list)
    lamperti.equation = do.call(substitute, list(lamperti.formula, substitute.list))
    lamperti.simplified = Deriv::Simplify(lamperti.equation)
    
    # The state equation is now in terms of the original state variable 'x', but we want it
    # in terms of the new state variable 'z'
    templist = psi.correct.state[1]
    names(templist) = state
    lamperti.complete = do.call(substitute, list(lamperti.simplified, templist))
    
    # replace the original SDE RHS with the new transformed one
    form = as.formula(paste(
      parse(text=paste("d",state,sep=""))[[1]],
      paste(deparse1(lamperti.complete),collapse=""),
      sep ="~"
    ))
    
    # add the new transformed equation to trans_system
    private$add_trans_systems(list(form=form, name=state))
  }
  
  ############### MAIN LOOP FOR OBSERVATIONS AND VARIANCES ###############
  # Transform the state entries in observation equations from e.g. x to exp(x)
  for(j in seq_along(private$obs.eqs)){
    
    # Get observation and variance rhs
    obs.rhs = private$obs.eqs.trans[[j]]$rhs
    obs.var.rhs = private$obs.var.trans[[j]]$rhs
    
    # copy rhs for repeated substitutions
    transformed.obs.rhs = obs.rhs
    transformed.obs.var.rhs = obs.var.rhs
    
    for(i in seq_along(states)){
      
      # select current state and transformation
      state = states[i]
      transform = transforms[i]
      
      # get the current transformation
      psi = psi.all[[transform]][1]
      names(psi) = state
      
      # substitute the state name into the psi_transformation list instead of the placeholder 'x'
      templist = list(x=parse(text=state)[[1]])
      psi.correct.state = lapply(psi, function(x) do.call(substitute, list(x, templist)))
      
      # substitute into new obs rhs
      transformed.obs.rhs = do.call(substitute, list(transformed.obs.rhs, psi.correct.state))
      transformed.obs.var.rhs = do.call(substitute, list(transformed.obs.var.rhs, psi.correct.state))
    }
    
    # obsname
    obsname = private$obs.names[j]
    # get lhs of observation (not the same as obsname for "complex" obs e.g. log(y))
    obslhs = private$obs.eqs[[j]]$lhs
    
    # Create formula and add the new observation rhs
    transformed.form.obs = as.formula(paste(c(obslhs, "~", transformed.obs.rhs),collapse=" "))
    private$add_trans_observations(list(form=transformed.form.obs, name=obsname))
    
    # Create formula and add the new observation variance rhs
    transformed.form.obs.var = as.formula(paste(c(obslhs, "~", transformed.obs.var.rhs),collapse=" "))
    private$add_trans_observation_variances(list(form=transformed.form.obs.var, name=obsname))
    
  }
  
  # return
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
  given.vars = c( private$parameter.names, private$input.names, private$state.names)
  bool = rhs.vars %in% given.vars
  if (any(!bool)) {
    stop("Error: The following variables(s) in the model have not been declared as parameters, inputs or states: \n\t ",
         paste(rhs.vars[!bool],collapse=", "))
  }
  
  ##### NOTE::: Do we need to remove this? doesnt really matter right? ####
  # Verify that all input and parameters are used in the model
  given.vars = c( private$parameter.names, private$input.names[-1])
  bool = given.vars %in% rhs.vars
  if (any(!bool)) {
    stop("The following variables(s) are unused: \n\t ", paste(given.vars[!bool],collapse=", "))
  }
  
  # return
  return(invisible(self))
}

