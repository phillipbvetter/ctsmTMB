# These functions are helper functions used when calling the ctsmrTMB method


# Applies algebraics, lamperti transformation and diff terms
apply_algebraics_and_lamperti <- function(self, private){
  apply_algebraics_and_define_trans_equations(self, private)
  calculate_diff_terms(self, private)
  apply_lamperti(self, private) 
  calculate_diff_terms(self, private)
}

create_state_space_function_strings <- function(self, private){
  create.stace.space.function.strings(self, private)
  create.rcpp.stace.space.function.strings(self, private)
}

#######################################################
# APPLY ALGEBRAIC TRANSFORMATIONS
#######################################################

apply_algebraics_and_define_trans_equations = function(self, private) {
  
  # extract rhs's
  alg.rhs = lapply(private$alg.eqs, function(x) x$rhs)
  sys.rhs = lapply(private$sys.eqs,function(x) x$rhs)
  obs.rhs = lapply(private$obs.eqs,function(x) x$rhs)
  obs.var.rhs = lapply(private$obs.var,function(x) x$rhs)
  
  # substitute algebraics into rhs's
  for(i in seq_along(alg.rhs)){
    this.alg <- alg.rhs[i]
    alg.rhs <- lapply(alg.rhs, function(x) do.call(substitute, list(x, this.alg)))
    sys.rhs = lapply(sys.rhs, function(x) do.call(substitute, list(x, this.alg)))
    obs.rhs = lapply(obs.rhs, function(x) do.call(substitute, list(x, this.alg)))
    obs.var.rhs = lapply(obs.var.rhs, function(x) do.call(substitute, list(x, this.alg)))
  }

  # system
  # for(i in seq_along(private$sys.eqs)) {
  for(i in seq_along(sys.rhs)) {
    temp.form = private$sys.eqs[[i]]$form
    temp.form[[3]] = sys.rhs[[i]] #[[3]] is rhs of a formula
    temp.list = list(form=temp.form, name=private$state.names[i])
    private$add_trans_systems(temp.list)
  }

  # observations
  # for(i in seq_along(private$obs.eqs)) {
  for(i in seq_along(obs.rhs)) {
    temp.form = private$obs.eqs[[i]]$form
    temp.form[[3]] = obs.rhs[[i]] #[[3]] is rhs of a formula
    temp.list = list(form=temp.form, name=private$obs.names[i])
    private$add_trans_observations(temp.list)
  }
  
  # observation variances
  # for(i in seq_along(private$obs.var)) {
  for(i in seq_along(obs.var.rhs)) {
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
  
  # Calculate drift and diffusion terms (differentiate w.r.t dt and dw)
  for (i in seq_along(private$sys.eqs.trans)) {
    private$diff.terms[[i]] = lapply(private$diff.processes, function(x) ctsmTMB.Deriv(f=private$sys.eqs.trans[[i]]$rhs, x=x))
    names(private$diff.terms[[i]]) = private$diff.processes
  }
  names(private$diff.terms) = private$state.names
  
  
  # dfdx
  for(i in seq_along(private$sys.eqs.trans)){
    private$diff.terms.drift[[i]] = lapply(private$state.names, function(x) ctsmTMB.Deriv(f=private$sys.eqs.trans[[i]]$diff.dt, x=x))
    names(private$diff.terms.drift[[i]]) = private$state.names
  }
  names(private$diff.terms.drift) = private$state.names
  
  # dhdx
  for(i in seq_along(private$obs.eqs.trans)){
    private$diff.terms.obs[[i]] = lapply(private$state.names, function(x) ctsmTMB.Deriv(f=private$obs.eqs.trans[[i]]$rhs, x=x))
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
  nonzero.diffterms = lapply(private$diff.terms, function(x) unlist(lapply(x, function(y) y==0))) #state names are retained
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
  
  # EXPLANATION / NOTE NOTE NOTE NOTE:::
  # The exp(x) is actually exp(z),  but we keep same names. 
  # This ensures that e.g if the original equation had state dep. diffusion:
  
  # dx ~ ... sigma * x * dw 
  # y ~ x + e
  
  # then the lamperti transform is z = log(x), so x = exp(z) so the transformed 
  # system is:
  
  # dz ~ ... + sigma * dw
  # but we just write x instead of z:
  # dx ~ ... + sigma * dw

  # This ensures that the transformed observation equation rhs is untransformed i.e.
  # y ~ exp(z) + e = x + e
  # where x is strictly positive so it is now the users responsibility to log
  # transform y i.e. the user should know to input the appropriate observation
  # equation from the start, namely:
  
  # log(y) ~ x + e
  # EXPLANATION / NOTE NOTE NOTE NOTE:::
  
  
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

