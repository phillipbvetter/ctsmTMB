#######################################################
# PASTE WITH COLLAPSE
#######################################################
paste00 = function(str_vec){
  paste(str_vec,collapse=", ")
}

#######################################################
# SIMPLIFY FORMULA
#######################################################

simplify_formula = function(form) {
  
  form = stats::as.formula(paste(
    form[[2]],
    paste(deparse(Deriv::Simplify(form[[3]])),collapse=""),
    sep="~"
  ))
  
}

#######################################################
# EXTRA FUNCTIONS FOR DERIV THAT ALSO WORKS IN TMB
#######################################################

# this only works because logit and invlogit is already defined in TMB
# if not they must be defined in the c++ file

logit = function(x) log(x/(1-x))
invlogit = function(x) 1/(1+exp(-x))
invlogit2 = function(x,a,b) 1/(1+exp(-a*(x-b)))

#######################################################
# CHANGE FROM R POWER NOTATION TO C++
#######################################################

# Changes powers in expressions from R-style (e.g. x^2) to pow(x,2) which is interpretable by C++
# from https://stackoverflow.com/questions/40606723/substitute-the-power-symbol-with-cs-pow-syntax-in-mathematical-expression

hat2pow <- function(e) {
  #check if you are at the end of the tree's branch
  if (is.name(e) || is.atomic(e)) {
    #replace ^
    if (e == quote(`^`)) return(quote(pow))
    return(e)
  }
  #follow the tree with recursion
  for (i in seq_along(e)) e[[i]] <- hat2pow(e[[i]])
  return(e)
}

###########################################################
# TRY THAT SUPPRESSES AND RECOVERS ANSWER DESPITE WARNINGS
###########################################################

try_withWarningRecovery = function(expr, silent=TRUE){
  output = try(withCallingHandlers(
    {
      expr
    },
    warning = function(w) {
      invokeRestart("muffleWarning")
    }
  ),silent=silent)
  
  return(output)
}

###########################################################
# ATOMIC LOG-DETERMINANT FOR RTMB
###########################################################

logdet <- RTMB::ADjoint(
  function(x) {
    dim(x) <- rep(sqrt(length(x)), 2)
    determinant(x, log=TRUE)$modulus
  },
  function(x, y, dy) {
    dim(x) <- rep(sqrt(length(x)), 2)
    t(RTMB::solve(x)) * dy
  },
  name = "logdet")

###########################################################
# SOLVE LYAPUNOV EQUATION
###########################################################

# This code is grabbed from Uffe Thygesens SDEtools package
# https://github.com/Uffe-H-Thygesen/SDEtools/tree/main

lyap_solver = function (A, G) {
  # Solves the Lyapunov Equation AP + PA' + GG' = 0 for P.
  
  Q <- G %*% t(G)
  A <- as.matrix(A)
  I <- diag(rep(1, nrow(A)))
  P <- kronecker(I, A) + kronecker(A, I)
  X <- -solve(P, as.numeric(Q))
  
  return(matrix(X, nrow = nrow(A)))
}

mean_solver = function(A, B, u) {
  # Solves the Mean Stationarity Equation Ax + Bu = 0 for x
  # Matrices: A, B
  # Vectors: x, u
  
  x <- solve(A, - B %*% u)
  
  return(x)
}
