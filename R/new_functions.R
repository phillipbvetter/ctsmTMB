# NOTE:
# Rather than specifying these functions only in the environment of Deriv,
# they are specified globally here. That way Deriv sees them, but they can also
# be seen by RTMB, and also when computing various statistics of the return fit

erf <- function(x) 2 * pnorm(sqrt(2)*x) - 1

logit <- function(x) log(x/(1-x))

invlogit <- function(x) 1/(1+exp(-x))

invlogit2 <- function(x,a,b) 1/(1+exp(-a*(x-b)))
