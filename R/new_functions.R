# NOTE:
# Rather than specifying these functions only in the environment of Deriv,
# they are specified globally here. That way Deriv sees them, but they can also
# be seen by RTMB, and also when computing various statistics of the return fit

#' @title Error Function
#' @param x Input
erf = function(x) 2 * pnorm(sqrt(2)*x) - 1

#' @title Logit
#' @param x Input
logit = function(x) log(x/(1-x))

#' @title Inverse Logit
#' @param x Input
invlogit = function(x) 1/(1+exp(-x))

#' @title Inverse Logit with Shift and Scale
#' @param x Input
#' @param a Scaling
#' @param b Shifting
invlogit2 = function(x,a,b) 1/(1+exp(-a*(x-b)))
