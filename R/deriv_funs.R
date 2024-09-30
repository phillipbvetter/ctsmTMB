#' @title Create Deriv::Deriv enviroment
#' @description
#' This functions returns an environment for Deriv, used in the myDeriv function.
#' @details
#' Deriv requires e.g. that the "erf" function is specified, even though it does
#' not need to do numeric calculations (in our case) - just symbolics. This is 
#' why we only need to specify some arbitrary function return (NULL).
get.Deriv.env = function(){
  e <- new.env()
  
  # add custom functions
  
  # return
  return(e)
}

#' @title Create Deriv::Deriv custom drule
#' @description
#' This functions returns a table of derivatives for symbolic differentiation
#' with Deriv, used in the myDeriv function.
get.Deriv.drule = function(){
  
  # get standard library
  drule <- Deriv::drule
  
  # add custom entries
  drule$erf <- alist(x = 2/sqrt(pi)*exp(-x^2))
  
  # return
  return(drule)
}


#' @title Modify Deriv::Deriv
#' @description
#' We create our own Deriv function based on Deriv::Deriv with custom environment and
#' drules to support new functions e.g erf
ctsmTMB.Deriv = function(
    f,
    x = if (is.function(f)) NULL else all.vars(if (is.character(f)) parse(text = f) else
      f),
    env = get.Deriv.env(),
    use.D = FALSE,
    cache.exp = FALSE,
    nderiv = NULL,
    combine = "c",
    drule = get.Deriv.drule()
){
  Deriv::Deriv(
    f = f,
    x = x,
    env = env,
    use.D = use.D,
    cache.exp = cache.exp,
    nderiv = nderiv,
    combine = combine,
    drule = drule
  )
}

#' @title Error Function
#' @description
#' Error function used internally in the package when returning results etc...matches error function on C++ side.
erf = function(x) 2 * pnorm(sqrt(2)*x) - 1

#' @title Logit
logit = function(x) log(x/(1-x))

#' @title Inverse Logit
invlogit = function(x) 1/(1+exp(-x))

# NOTE:
# WE COULD ADD THE FUNCTIONS TO THE ENVIROMENT OF DERIV BUT
# BECAUSE THEY ARE NEEDED ALSO FOR RTMB AND FOR RETURNING FIT RESULTS
# WE CAN DEFINE THEM GLOBALLY IN THE PACKAGE AS ABOVE AND THEREBY
# WE DO NOT NEED TO SPECIFY THEM IN THE DERIV ENVIROMENT...
