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
  if(any(private$procedure == c("estimation","construction"))){
    if (is.null(private$initial.state)) {
      stop("You must set an initial state estimate and covariance")
    }
  }
  
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

#######################################################
# CONSTRUCT DRIFT, DIFF, OBS FUNCTIONS FOR RTMB
#######################################################

create_rtmb_function_strings = function(self, private)
{
  
  # Create substitution translation list
  obsList = lapply(seq_along(private$obs.names), function(id) substitute(obsVec[i],list(i=as.numeric(id))))
  parList = lapply(seq_along(private$parameter.names), function(id) substitute(parVec[i],list(i=as.numeric(id))))
  stateList = lapply(seq_along(private$state.names), function(id) substitute(stateVec[i],list(i=as.numeric(id))))
  inputList = lapply(seq_along(private$input.names), function(id) substitute(inputVec[i],list(i=as.numeric(id))))
  names(obsList) = private$obs.names
  names(parList) = private$parameter.names
  names(stateList) = private$state.names
  names(inputList) = private$input.names
  subsList = c(obsList, parList, stateList, inputList)
  
  ##################################################
  # drift
  ##################################################
  f.elements = sapply( seq_along(private$state.names),
                       function(i){
                         drift.term = private$diff.terms[[i]]$dt
                         deparse1(do.call(substitute, list(drift.term, subsList)))
                       })
  
  f.function.text = paste('
  f__ = function(stateVec, parVec, inputVec){
    ans = c(F_ELEMENTS)
    return(ans)
  }')
  
  private$rtmb.function.strings$f = stringr::str_replace_all(f.function.text, 
                                                             pattern="F_ELEMENTS", 
                                                             replacement=paste(f.elements,collapse=","))
  
  ##################################################
  # drift jacobian
  ##################################################
  dfdx.elements = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term = Deriv::Deriv(private$diff.terms[[i]]$dt,x=private$state.names[j], cache.exp=F)
      dfdx.elements = c(dfdx.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  dfdx.function.text = paste('
  dfdx__ = function(stateVec, parVec, inputVec){
    ans = RTMB::matrix(c(DFDX_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_STATES, byrow=T)
    return(ans)
  }')
  
  dfdx.function.text = stringr::str_replace_all(dfdx.function.text, 
                                                pattern="NUMBER_OF_STATES", 
                                                replacement=deparse(as.numeric(private$number.of.states)))
  
  dfdx.function.text = stringr::str_replace_all(dfdx.function.text, 
                                                pattern="DFDX_ELEMENTS", 
                                                replacement=paste(dfdx.elements,collapse=","))
  
  private$rtmb.function.strings$dfdx = dfdx.function.text
  
  ##################################################
  # diffusion
  ##################################################
  g.elements = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term = private$diff.terms[[i]][[j+1]]
      g.elements = c(g.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  g.function.text = paste('
  g__ = function(stateVec, parVec, inputVec){
    ans = RTMB::matrix(c(G_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_DIFFUSIONS, byrow=T)
    return(ans)
  }')
  
  g.function.text = stringr::str_replace_all(g.function.text, 
                                             pattern="NUMBER_OF_STATES", 
                                             replacement=deparse(as.numeric(private$number.of.states)))
  
  g.function.text = stringr::str_replace_all(g.function.text, 
                                             pattern="NUMBER_OF_DIFFUSIONS", 
                                             replacement=deparse(as.numeric(private$number.of.diffusions)))
  
  g.function.text = stringr::str_replace_all(g.function.text, 
                                             pattern="G_ELEMENTS", 
                                             replacement=paste(g.elements,collapse=","))
  
  private$rtmb.function.strings$g = g.function.text
  
  
  ##################################################
  # observation
  ##################################################
  h.elements = sapply(seq_along(private$obs.names), 
                      function(i){
                        term = private$obs.eqs.trans[[i]]$rhs
                        deparse1(do.call(substitute, list(term, subsList)))
                      })
  
  h.function.text = paste('
  h__ = function(stateVec, parVec, inputVec){
    ans = c(H_ELEMENTS)
    return(ans)
  }')
  
  private$rtmb.function.strings$h = stringr::str_replace_all(h.function.text, 
                                                             pattern="H_ELEMENTS", 
                                                             replacement=paste(h.elements,collapse=","))
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dhdx.elements = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = private$diff.terms.obs[[i]][[j]]
      dhdx.elements = c(dhdx.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  dhdx.function.text = paste('
  dhdx__ = function(stateVec, parVec, inputVec){
    ans = RTMB::matrix(c(DHDX_ELEMENTS),nrow=NUMBER_OF_OBSERVATIONS, ncol=NUMBER_OF_STATES, byrow=T)
    return(ans)
  }')
  
  dhdx.function.text = stringr::str_replace_all(dhdx.function.text, 
                                                pattern="NUMBER_OF_STATES", 
                                                replacement=deparse(as.numeric(private$number.of.states)))
  
  dhdx.function.text= stringr::str_replace_all(dhdx.function.text, 
                                               pattern="NUMBER_OF_OBSERVATIONS", 
                                               replacement=deparse(as.numeric(private$number.of.observations)))
  
  dhdx.function.text = stringr::str_replace_all(dhdx.function.text, 
                                                pattern="DHDX_ELEMENTS", 
                                                replacement=paste(dhdx.elements,collapse=","))
  
  private$rtmb.function.strings$dhdx = dhdx.function.text
  
  ##################################################
  # observation variance
  ##################################################
  
  ### VECTOR FORM ### (for laplace)
  
  hvar.elements = sapply(seq_along(private$obs.var.trans), 
                         function(i) {
                           term = private$obs.var.trans[[i]]$rhs
                           deparse1(do.call(substitute, list(term, subsList)))
                         })
  
  hvar.function.text = paste('
  hvar__ = function(stateVec, parVec, inputVec){
    ans = c(HVAR_ELEMENTS)
    return(ans)
  }')
  
  private$rtmb.function.strings$hvar = stringr::str_replace_all(hvar.function.text, 
                                                                pattern="HVAR_ELEMENTS", 
                                                                replacement=paste(hvar.elements,collapse=","))
  
  ### MATRIX FORM ### (for kalman)
  
  hvar.function.text = paste('
  hvar__matrix = function(stateVec, parVec, inputVec){
    ans = RTMB::matrix(c(HVAR_ELEMENTS),nrow=NUMBER_OF_OBSERVATIONS, ncol=NUMBER_OF_OBSERVATIONS, byrow=T)
    return(ans)
  }')
  
  hvar.function.text = stringr::str_replace_all(hvar.function.text, 
                                                pattern="NUMBER_OF_OBSERVATIONS", 
                                                replacement=deparse(as.numeric(private$number.of.observations)))
  
  hvar.function.text = stringr::str_replace_all(hvar.function.text, 
                                                pattern="HVAR_ELEMENTS", 
                                                replacement=paste(hvar.elements,collapse=","))
  
  private$rtmb.function.strings$hvar__matrix = hvar.function.text
  
  return(invisible(NULL))
}

create_rtmb_laplace_string = function(self, private){
  
  
  # define the negative loglikelihood 
  laplace.text = paste('
  laplace.nll = function(p){
    # set negative log-likelihood
    nll = 0
    
    # small identity matrix
    small_identity = diag(1e-8, nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_STATES)
    
    # extract state random effects
    stateMat = cbind(STATE_NAMES)
    
    # fixed effects parameter vector
    parVec = cbind(FIXED_PARAMETERS)
    
    ###### TIME LOOP START #######
    for(i in 1:(nrow(obsMat)-1)) {
    
    # Define inputs and use first order input interpolation
    inputVec = inputMat[i,]
    dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
    
    ###### BETWEEN TIME POINTS LOOP START #######
    for(j in 1:ode_timesteps[i]){
    
    # grab current and next state
    x_now = stateMat[ode_cumsum_timesteps[i]+j,]
    x_next = stateMat[ode_cumsum_timesteps[i]+j+1,]
    
    # compute drift (vector) and diffusion (matrix)
    f = f__(x_now, parVec, inputVec)
    g = g__(x_now, parVec, inputVec)
    inputVec = inputVec + dinputVec
    
    # assume multivariate gauss distribution according to euler-step
    # and calculate the likelihood
    z = x_next - (x_now + f * ode_timestep_size[i])
    v = (g %*% t(g) + small_identity) * ode_timestep_size[i]
    nll = nll - RTMB::dmvnorm(z, Sigma=v, log=TRUE)
    }
    ###### BETWEEN TIME POINTS LOOP END #######
    }
    ###### TIME LOOP END #######
    
    obsMat = RTMB::OBS(obsMat)
    ###### DATA UPDATE START #######
    for(i in 1:NUMBER_OF_OBSERVATIONS){
    for(j in 1:length(iobs[[i]])){
    
    # Get index where observation is available
    k = iobs[[i]][j]
    
    # Get corresponding input, state and observation
    inputVec = inputMat[k,]
    stateVec = stateMat[ode_cumsum_timesteps[k]+1,]
    obsScalar = obsMat[k,i]
    
    # Observation equation and varianace
    h_x = h__(stateVec, parVec, inputVec)[i]
    hvar_x = hvar__(stateVec, parVec, inputVec)[i]
    
    # likelihood contribution
    nll = nll - RTMB::dnorm(obsScalar, mean=h_x, sd=sqrt(hvar_x), log=TRUE)
    }}
    ###### DATA UPDATE END #######
    
    # return
    return(invisible(nll))
}')
  
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="NUMBER_OF_STATES", 
                                          replacement=deparse(as.numeric(private$number.of.states)))
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="NUMBER_OF_DIFFUSIONS", 
                                          replacement=deparse(as.numeric(private$number.of.diffusions)))
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="NUMBER_OF_OBSERVATIONS", 
                                          replacement=deparse(as.numeric(private$number.of.observations)))
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="STATE_NAMES", 
                                          replacement=paste("p$",private$state.names,sep="",collapse=","))
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="FIXED_PARAMETERS", 
                                          replacement=paste("p$",private$parameter.names,sep="",collapse=","))
  
  # Store likelihood function
  private$rtmb.nll.strings$laplace = laplace.text
  
  return(invisible(NULL))
}

create_rtmb_ekf_string = function(self, private){
  
  ################################################
  ekf.text = paste('
  ekf.nll = function(p){
    ####### INITIALIZATION #######
    # storage variables
    xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
    xPrior[[1]] <- xPost[[1]] <- stateVec
    pPrior[[1]] <- pPost[[1]] <- covMat
    
    # constants
    halflog2pi = log(2*pi)/2
    I0 <- diag(NUMBER_OF_STATES)
    E0 <- diag(NUMBER_OF_OBSERVATIONS)
    
    # set negative log-likelihood
    nll = 0
    
    # fixed effects parameter vector
    parVec = c(FIXED_PARAMETERS)
    
    # ###### TIME LOOP START #######
    for(i in 1:(nrow(obsMat)-1)){
    
    # # Define inputs and use first order input interpolation
    inputVec = inputMat[i,]
    dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
    
    ###### TIME UPDATE - ODE SOLVER #######
    for(j in 1:ode_timesteps[i]){
    sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i], ode_solver)
    stateVec = sol[[1]]
    covMat = sol[[2]]
    }
    xPrior[[i+1]] = stateVec
    pPrior[[i+1]] = covMat
    
    ######## DATA UPDATE - KALMAN ALGORITHM ########
    obsVec = obsMat[i+1,]
    inputVec = inputMat[i+1,]
    obsVec_bool = !is.na(obsVec)
    
    ######## DATA UPDATE - IF ANY DATA AVAILABLE ######## 
    if(any(obsVec_bool)){
    s = sum(obsVec_bool)
    y = obsVec[obsVec_bool]
    E = E0[obsVec_bool,obsVec_bool,drop=FALSE]
    H = h__(stateVec, parVec, inputVec)
    C = E %*% dhdx__(stateVec, parVec, inputVec)
    e = y - E %*% H
    V0 = hvar__matrix(stateVec, parVec, inputVec)
    # V0 = RTMB::diag(hvar__(stateVec, parVec, inputVec))
    V = E %*% V0 %*% t(E)
    R = C %*% covMat %*% t(C) + V
    Ri = RTMB::solve(R)
    K = covMat %*% t(C) %*% Ri
    
    # Update State
    stateVec = stateVec + K %*% e
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    nll = nll + 0.5 * logdet(R) + 0.5 * t(e) %*% Ri %*% e + halflog2pi * s
    
    # Store innovation and covariance
    Innovation[[i+1]] = e
    InnovationCovariance[[i+1]] = R
    }
    xPost[[i+1]] = stateVec
    pPost[[i+1]] = covMat
    }
    
    # ###### MAXIMUM A POSTERIOR #######
    
    # ###### REPORT #######
    RTMB::REPORT(Innovation)
    RTMB::REPORT(InnovationCovariance)
    RTMB::REPORT(xPrior)
    RTMB::REPORT(xPost)
    RTMB::REPORT(pPrior)
    RTMB::REPORT(pPost)
    
    # ###### RETURN #######
    return(nll)
}')
  
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="NUMBER_OF_STATES", 
                                      replacement=deparse(as.numeric(private$number.of.states)))
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="NUMBER_OF_DIFFUSIONS", 
                                      replacement=deparse(as.numeric(private$number.of.diffusions)))
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="NUMBER_OF_OBSERVATIONS", 
                                      replacement=deparse(as.numeric(private$number.of.observations)))
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="STATE_NAMES", 
                                      replacement=paste("p$",private$state.names,sep="",collapse=","))
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="FIXED_PARAMETERS", 
                                      replacement=paste("p$",private$parameter.names,sep="",collapse=","))
  
  private$rtmb.nll.strings$ekf = ekf.text
  
  
}

create_rcpp_function_strings = function(self, private){
  
  # Create substitution translation list
  obsList = lapply(seq_along(private$obs.names), function(id) substitute(obsVec(i),list(i=as.numeric(id-1))))
  parList = lapply(seq_along(private$parameter.names), function(id) substitute(parVec(i),list(i=as.numeric(id-1))))
  stateList = lapply(seq_along(private$state.names), function(id) substitute(stateVec(i),list(i=as.numeric(id-1))))
  inputList = lapply(seq_along(private$input.names), function(id) substitute(inputVec(i),list(i=as.numeric(id-1))))
  names(obsList) = private$obs.names
  names(parList) = private$parameter.names
  names(stateList) = private$state.names
  names(inputList) = private$input.names
  subsList = c(obsList, parList, stateList, inputList)
  
  
  ##################################################
  # drift
  ##################################################
  
  # Perform substitution of parameters, inputs and states
  f = sapply( seq_along(private$state.names),
              function(i){
                drift.term = hat2pow(private$diff.terms[[i]]$dt)
                new.drift.term = do.call(substitute, list(drift.term, subsList))
                sprintf("f(%i) = %s;",i-1,deparse1(new.drift.term))
              })
  code = sprintf("SEXP f(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd f(%s); 
                 %s
                 return Rcpp::wrap(f);
                 }",private$number.of.states, paste(f,collapse=""))
  
  private$Rcppfunction_f = code
  
  ##################################################
  # drift jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dfdx = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term = hat2pow(Deriv::Deriv(private$diff.terms[[i]]$dt,x=private$state.names[j], cache.exp=F))
      new.term = do.call(substitute, list(term, subsList))
      dfdx = c(dfdx, sprintf("dfdx(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  code = sprintf("SEXP dfdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dfdx(%s,%s);
                 %s
                 return Rcpp::wrap(dfdx);
                 }",private$number.of.states, private$number.of.states, paste(dfdx,collapse=""))
  
  private$Rcppfunction_dfdx = code
  
  ##################################################
  # diffusion
  ##################################################
  
  # calculate all the terms and substitute variables
  g = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term = hat2pow(private$diff.terms[[i]][[j+1]])
      new.term = do.call(substitute, list(term, subsList))
      g = c(g, sprintf("g(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  code = sprintf("SEXP g(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd g(%s,%s);
                 %s
                 return Rcpp::wrap(g);
                 }",private$number.of.states, private$number.of.diffusions, paste(g,collapse=""))
  
  private$Rcppfunction_g <- code
  
  ##################################################
  # observation
  ##################################################
  
  # calculate all the terms and substitute variables
  h = sapply(seq_along(private$obs.names), 
             function(i){
               term = hat2pow(private$obs.eqs.trans[[i]]$rhs)
               new.term = do.call(substitute, list(term, subsList))
               sprintf("h(%s) = %s;",i-1, deparse1(new.term))
             }) 
  
  code = sprintf("SEXP h(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd h(%s);
                 %s
                 return Rcpp::wrap(h);
                 }",private$number.of.observations, paste(h,collapse=""))
  
  private$Rcppfunction_h <- code
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dhdx = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = hat2pow(private$diff.terms.obs[[i]][[j]])
      new.term = do.call(substitute, list(term, subsList))
      dhdx = c(dhdx, sprintf("dhdx(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  code = sprintf("SEXP dhdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dhdx(%s,%s);
                 %s
                 return Rcpp::wrap(dhdx);
                 }", private$number.of.observations, private$number.of.states, paste(dhdx,collapse=""))
  
  private$Rcppfunction_dhdx <- code
  
  
  ##################################################
  # observation variance
  ##################################################
  
  hvar = lapply(seq_along(private$obs.var.trans), 
                function(i) {
                  term = hat2pow(private$obs.var.trans[[i]]$rhs)
                  new.term = do.call(substitute, list(term, subsList))
                  sprintf("hvar(%s,%s) = %s;", i-1, i-1, deparse1(new.term))
                })
  
  code = sprintf("SEXP hvar(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd hvar(%s,%s);
                 hvar.setZero();
                 %s
                 return Rcpp::wrap(hvar);
                 }", private$number.of.observations, private$number.of.observations, paste(hvar,collapse=""))
  
  private$Rcppfunction_hvar <- code
  
  
}
