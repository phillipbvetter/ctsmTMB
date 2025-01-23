#' @title Create Deriv::Deriv enviroment
#' @description
#' This functions returns an environment for Deriv, used in the ctsmTMB.Deriv function.
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
#' with Deriv, used in the ctsmTMB.Deriv function.
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
#' @param f see \link[Deriv]{Deriv} documentation
#' @param x see \link[Deriv]{Deriv} documentation
#' @param env see \link[Deriv]{Deriv} documentation
#' @param use.D see \link[Deriv]{Deriv} documentation
#' @param cache.exp see \link[Deriv]{Deriv} documentation
#' @param nderiv see \link[Deriv]{Deriv} documentation
#' @param combine see \link[Deriv]{Deriv} documentation
#' @param drule see \link[Deriv]{Deriv} documentation
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
