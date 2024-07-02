#######################################################
# CHECK EXISTENCE OF NAME FUNCTION
#######################################################

check_if_name_is_overwritable = function(newvar, type, self, private) {
  # A state (or observation or input) can only overwrite another states (or 
  # observation or input), but not an input or an observation.
  
  if (newvar %in% names(private$sys.eqs) & type!="state") {
    stop("The variable ", newvar, " is already in use (state).")
  }
  
  if (newvar %in% names(private$obs.eqs) & type!="obs" & type!="obsvar") {
    stop("The variable ", newvar, " is already in use (observation).")
  }
  
  if (newvar %in% names(private$inputs) & type!="input") {
    stop("The variable ", newvar, " is already in use (input).")
  }
  
  return(invisible(self))
}

#######################################################
# CHECK SYSTEM EQUATION FUNCTION
#######################################################

check_system_eqs = function(form, self, private) {
  
  # CHECK "FORM"
  #######################################################
  
  if(!inherits(form,"formula")){
    stop("The system equation should be a formula e.g. dx ~ ... * dt + ... * dw1 + ... * dw2")
  }
  
  lhs = form[[2]]
  rhs = form[[3]]
  
  # If  LHS has length one (single term) and class "name"
  if (!(length(lhs) == 1)) {
    stop("You have multiple terms on the left-hand side")
  }
  
  # Is the state name valid?
  state = stringr::str_match(deparse1(lhs),"^d([a-zA-Z]+[a-zA-Z0-9]*)$")[2]
  if (is.na(state)) {
    stop("That state name is not allowed - use d followed by any number of letters, followed by any number of digits")
  }
  
  # dont use dt or dw in the state name
  match = stringr::str_match(deparse1(lhs),"^(?!d[tw])")
  if (is.na(match)) {
    stop("The state name can't begin with dt or dw")
  }
  
  # CHECK VARIABLES IN FORM SDE
  #######################################################
  
  # Find all diffusion processes
  diff.proc = unique(unlist(stringr::str_extract_all(paste(deparse1(rhs),collapse=""),"dw([a-zA-Z0-9]*)")))
  diff.processes = c("dt",diff.proc)
  diff.terms = lapply(diff.processes, function(x) { Deriv::Deriv(rhs, x, cache.exp=FALSE) })
  names(diff.terms) = diff.processes
  
  # There must be diffusion processes
  if(length(diff.proc)==0){
    stop("You are not allowed to specify processes without any diffusion dw(...). You can use 0 * dw if no diffusion is desired.")
  }
  
  # Check for dt/dw cross terms
  valid = all(unlist(lapply(diff.terms, function(x) all(is.na(match(diff.processes, all.vars(x)))))))
  if (!valid) { stop("There are illegal dt and dw cross terms") }
  
  # Check if any variables are outside scope like c in: f(.) * dt + g(.) * dw + c
  pars = unique(all.vars(rhs))
  pars.after = as.vector(c(diff.processes, unlist(lapply(diff.terms, all.vars))))
  if (any(is.na(match(pars, pars.after)))) {
    stop("You have illegally specified terms other than drifts (... * dt) and diffusions (... * dw).")
  }
  # The above does not capture constants like the 5 in ' ... * dw + 5' so check that.
  # This is a bit annoying because you can't see numerics directly. We work around by
  # multiplying a variable onto, and the checking if there are any variables left after
  # setting the diffusions equal to zero
  rhs.checker1 = parse(text=paste("(",deparse1(rhs),")","*K__",sep=""))[[1]]
  zero.list = as.list(numeric(length(diff.processes)))
  names(zero.list) = diff.processes
  rhs.checker2 = do.call(substitute,list(rhs.checker1,zero.list))
  rhs.checker3 = all.vars(Deriv::Simplify(rhs.checker2))
  if(length(rhs.checker3)>0){
    stop("There are illegal terms outside of the drifts dt or diffusions dw(s).")
  }
  
  # Check if any variables are called dt(something)
  pars = unique(all.vars(rhs))
  pars = pars[!(pars %in% "dt")]
  if(any(stringr::str_detect(pars,"dt.*"))){
    stop("There are illegal variable names apart from 'dt' that begins with dt.")
  }
  
  # extract all variables
  bool = unique(all.vars(rhs)) %in% diff.processes
  variables = unique(all.vars(rhs))[!bool]
  
  # return
  result = list(name=state, form=form, rhs=rhs, allvars=variables, diff=diff.processes)
  return(result)
}

#######################################################
# CHECK OBSERVATION EQUATION FUNCTION
#######################################################

check_observation_eqs = function(forms, self, private) {
  
  form = forms$form
  obsname = forms$name
  
  if(!inherits(form,"formula")){
    stop("The observation equation should be a formula e.g. 'y ~ ...")
  }
  
  lhs = form[[2]]
  rhs = form[[3]]
  
  # if the observation is complex (of class 'call') then we must have a name provided
  if(inherits(lhs,"call")){
    if(is.null(obsname)){
      stop("You must provide argument 'obsnames' for observations with complex left-hand sides.")
    }
  } 
  # if the observation is of class 'name' then just grab that variable name with deparse
  else {
    obsname = deparse1(lhs) 
  }
  
  # Check if the observation name is OK
  bool = stringr::str_detect(obsname,"^[a-zA-Z]+")
  if(!bool){
    stop("The observation name must begin with a letter")
  }
  
  # you cannot observe a differential process
  bool = stringr::str_detect(obsname,"^(?!d[tw])[[:alnum:]]*")
  if (!bool) {
    stop("You can't observe a differential process.")
  }
  
  # extract all variables
  variables = unique(all.vars(rhs))
  
  # return
  result = list(name=obsname,form=form,rhs=rhs,lhs=lhs,allvars=variables)
  return(result)
}

#######################################################
# CHECK OBSERVATION EQUATION FUNCTION
#######################################################

check_observation_variance_eqs <- function(form, self, private) {
  
  if(!inherits(form,"formula")){
    stop("The observation variance equation should be a formula whose left-hand side is the name of a previously specified observation e.g. y ~ ...")
  }
  
  lhs = form[[2]]
  rhs = form[[3]]
  obsname = deparse1(lhs)
  
  if(inherits(lhs, "call")){
    stop("The left-hand side of an observation variance equation can only be a single variable, not a function expression")
  }
  
  # Is there an observation with that name?
  if (!(obsname %in% names(private$obs.eqs))) {
    stop("Please add an observation equation for ", deparse1(lhs), " before specifying its variance")
  }
  
  # extract all variables
  variables = unique(all.vars(rhs))
  
  # overwrite the lhs side of form when its complex e.g. log(y) ~ ...
  form[[2]] = private$obs.eqs[[obsname]]$form[[2]]
  
  # return
  result = list(name=obsname, form=form, rhs=rhs, allvars=variables)
  return(result)
}

#######################################################
# CHECK DATA INPUT FUNCTION
#######################################################

check_inputs <- function(input, self, private) {
  
  # Check for correct input class
  if (!is.name(input)) {
    stop("The inputs should be of class 'name' i.e. use $add_inputs(a)")
  }
  
  name = deparse1(input)
  
  # Does the input name start with dt or dw?
  valid = !is.na(stringr::str_match(name,"^(?!d[tw])[[:alnum:]]*"))
  if (!valid) {
    stop("Input names are not allowed to start with dt or dw")
  }
  
  # Reserved input names
  valid = !(name == "t")
  if (!valid) {
    stop("The name 't' is already reserved for the time vector")
  }
  
  result = list(name=name, input=input)
  return(result)
}

check_parameter_vector = function(par, parname, self, private) {
  
  # check if numeric
  if(!is.numeric(par)){
    stop(sprintf("The parameter %s gave an error:
                 Please provide a numeric vector",parname))
  }
  
  # must be length 1 or 3
  if(!any(length(par) == c(1,3))){
    stop("The parameter vector must have length 1 or 3")
  }

  # the parameter name strings must start with a character
  bool = stringr::str_detect(parname,"^[[:alpha:]][[:alnum:]_-]*$")
  if(!bool){
    stop("The parameter name ",parname, " is not valid. The name must begin with a letter, 
         and can only contain numerals, letters, underscore (_) and dash (-).")
  }
  
  # parameter name can't begin with dw or dt
  bool = stringr::str_detect(parname,"^(?!d[tw])[[:alnum:]]*")
  if(!bool){
    stop("The parameter names are not allowed to start with dt or dw, since these are reserved for differentials")
  }
  
  # set expected names names
  expected.names = c("initial","lower","upper")
  
  # if the vector has length 1, then set to length 3 and set names
  if(length(par)==1){
    length(par) = 3
    names(par) = expected.names
  }
  
  # Is the 3-vector named?, otherwise name it
  if (is.null(names(par))){
    if(length(par)==3) names(par) = c("initial","lower","upper")
  }
  
  # if the 3-vector is already named, are all names present?
  if(!all(expected.names %in% names(par))){
    stop("The parameter ", parname, " gave an error because it was named but did not contain all three required names  'initial', 'lower' and 'upper'")
  }
  
  # the initial value can't be NA
  if(is.na(par["initial"])){
    stop("The parameter ", parname, " gave an error because the initial value was NA") 
  }
  
  # if either of lower or upper are NA, then set both as NA
  if(any(is.na(par[c("lower","upper")]))){
    par[c("lower","upper")] = NA
  }
  
  # check if the values are ascending lower <= initial <= upper
  if(!all(is.na(par[c("lower","upper")]))){
    if(any(diff(par[c("lower","initial","upper")]) < 0)){
      stop("The parameter ", parname, " does not have ascending bounds i.e. lower bound <= initial value <= upper bound.")
    }
  }
  
  # IS THE PARAMETER IN THE MODEL?
  # the parameter name must be present in the object already - check all entries
  # but disregard parameter names on LHS of algebraics that will be replaced by
  # the algebraic RHS
  all.names = unique(unlist(c(
    lapply(private$sys.eqs, function(x) x$allvars),
    lapply(private$obs.eqs, function(x) x$allvars),
    lapply(private$obs.var, function(x) x$allvars),
    lapply(private$alg.eqs, function(x) all.vars(x$rhs))
  ))) 
  bool = all.names %in% names(private$alg.eqs)
  all.names = all.names[!bool]
  check.bool = parname %in% all.names
  if(!check.bool){
    stop("The following parameter is missing from the defined model (after applying the algebraic substitutions): ", parname)
  }
  
  return(invisible(par))
  
  }

check_parameter_matrix <- function(parmat, self, private) {
  
  # set column names if 3 columns and no column names
  if(is.null(colnames(parmat)) & ncol(parmat)==3){
    colnames(parmat) = c("initial","lower","upper")
    message("Note: No colnames were provided in parameter matrix - assuming order 'initial', 'lower', 'upper'")
  }
  
  # are column names initial, lower and upper present?
  col.names = colnames(parmat)
  expected.names = c("initial","lower","upper")
  bool = expected.names %in% col.names
  if(!all(bool)){
    stop(sprintf("Missing column names: %s", paste(expected.names[!bool],collapse=", ")))
  }
  
  # extract relevant columns
  parmat = as.matrix(parmat[,c("initial","lower","upper"),drop=FALSE])
  
  # is numerics?
  if(!is.numeric(parmat)){
    stop("The parameter matrix values must be numerics")
  }
  
  # has 3 columns?
  if (nrow(parmat)==0) {
    stop("The parameter matrix must have at least one row")
  }
  
  # are parameter names supplied?
  parnames = rownames(parmat)
  if (is.null(parnames)) {
    stop("You have not supplied any parameter names. Use rownames")
  }
  
  # the parameter name strings must start with a character
  bool = stringr::str_detect(parnames,"^[[:alpha:]][[:alnum:]_-]*$")
  if (sum(bool) != length(bool)) {
    stop("The parameter names ",paste(parnames[!bool],collapse = ", "), " are not valid")
  }
  
  # parameter name can't begin with dw or dt
  bool = stringr::str_detect(parnames,"^(?!d[tw])[[:alnum:]]*")
  if (sum(bool) != length(bool)) {
    stop("Parameter names are not allowed to start with dt or dw")
  }
  
  # the parameter name must be present in the object already - check all entries
  all.names = unique(unlist(c(
    lapply(private$sys.eqs, function(x) x$allvars),
    lapply(private$obs.eqs, function(x) x$allvars),
    lapply(private$obs.var, function(x) x$allvars),
    lapply(private$alg.eqs, function(x) all.vars(x$rhs))
  ))) 
  
  bool = all.names %in% names(private$alg.eqs)
  all.names = all.names[!bool]
  check.bool = parnames %in% all.names
  if(!all(check.bool)){
    stop("The following parameter is not a part of the current model, after applying the algebraic substitutions: ", paste(parnames[!check.bool],collapse=", "))
  }
  
  # result = list(parnames)
  # return(invisible(result))
  return(parmat)
}

#######################################################
# CHECK ALGEBRAIC RELATIONS
#######################################################

check_algebraics = function(form, self, private) {
  
  if(!inherits(form,"formula")){
    stop("The algebraic relation should be a formula e.g. 'theta ~ exp(log_theta) or x ~ logit(z)")
  }
  
  lhs = form[[2]]
  rhs = form[[3]]
  
  # Only single terms on LHS
  if (!(length(lhs) == 1)) {
    stop("You have multiple terms on the left-hand side")
  }
  
  name = deparse1(lhs)
  deparse_rhs = deparse1(rhs)
  
  # You can't redefine differentials
  bool = stringr::str_match(name,"^(?!d[tw])[[:alnum:]]*")
  if (is.na(bool)) {
    stop("You are not allowed to redefine differential processes.")
  }
  
  # You can't have differentials on the RHS
  bool = stringr::str_match(deparse_rhs,"^(?!d[tw])[[:alnum:]]*")
  if (is.na(bool)) {
    stop("You are not allowed to have differential processes on the right-hand side of an algebraic relation.")
  }
  
  # You can not apply algebraics to a state
  if (name %in% private$state.names) {
    stop("Redefining a state is not allowed: ", deparse1(form))
  }
  
  # You can't apply algebraics to an input
  if (name %in% private$input.names) {
    stop("Redefining an input is not allowed: ", deparse1(form))
  }
  
  # You can't apply algebraics to an observation
  if (name %in% private$obs.names) {
    stop("Redefining an observation is not allowed: ", deparse1(form))
  }
  
  result = list(name=name, form=form,rhs=rhs)
  return(result)
}

#######################################################
# REMOVE PARAMETER
#######################################################

remove_parameter = function(parname, self, private) {
  
  # remove parameter from parameter list
  bool = !(private$parameter.names %in% parname)
  private$parameters = private$parameters[bool]
  
  # update parameter names
  private$parameter.names = names(private$parameters)
  
  # remove parameter from fixed parameter list
  bool = !(names(private$fixed.pars) %in% parname)
  private$fixed.pars = private$fixed.pars[bool]
  
  return(invisible(self))
}
