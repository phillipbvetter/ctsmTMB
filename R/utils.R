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
