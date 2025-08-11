# Deriv requires e.g. that the "erf" function is specified, even though it does
# not need to do numeric calculations (in our case) - just symbolics. This is 
# why we only need to specify some arbitrary function return (NULL).
get_Deriv_environment = function(){
  e <- new.env()
  
  # add custom functions
  
  # return
  return(e)
}

# This functions returns a table of derivatives for symbolic differentiation
# with Deriv, used in the ctsmTMB.Deriv function.
get_Deriv_drules = function(){
  
  # get standard library
  drule <- Deriv::drule
  
  # add custom entries
  drule$erf <- alist(x = 2/sqrt(pi)*exp(-x^2))
  
  # return
  return(drule)
}


# We create our own Deriv function based on Deriv::Deriv with custom environment and
# drules to support new functions e.g erf
ctsmTMB.Deriv = function(
    f,
    x = if (is.function(f)) NULL else all.vars(if (is.character(f)) parse(text = f) else
      f),
    env = get_Deriv_environment(),
    use.D = FALSE,
    cache.exp = FALSE,
    nderiv = NULL,
    combine = "c",
    drule = get_Deriv_drules()
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



#################################################################################
#################################################################################

# These are additional functions needed.

# NOTE:
# Rather than specifying these functions only in the environment of Deriv,
# they are specified globally here. That way Deriv sees them, but they can also
# be seen by RTMB, and also when computing various statistics of the return fit

# erf <- function(x) 2 * pnorm(sqrt(2)*x) - 1
erf <- function(x) 2 * RTMB::pnorm(sqrt(2)*x) - 1

logit <- function(x) log(x/(1-x))

invlogit <- function(x) 1/(1+exp(-x))

invlogit2 <- function(x,a,b) 1/(1+exp(-a*(x-b))) # not useable on the C++ side

