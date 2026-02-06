#######################################################
# CHANGE FROM R POWER NOTATION TO C++
#######################################################

# Changes powers in expressions from R-style (e.g. x^2) to pow(x,2) which is interpretable by C++
# from https://stackoverflow.com/questions/40606723/substitute-the-power-symbol-with-cs-pow-syntax-in-mathematical-expression

# This function gave an error for very few people which the
# other function did not...
# hat2pow <- function(e) {
#   #check if you are at the end of the tree's branch
#   if (is.name(e) || is.atomic(e)) {
#     #replace ^
#     if (e == quote(`^`)) return(quote(pow))
#     return(e)
#   }
#   #follow the tree with recursion
#   for (i in seq_along(e)) e[[i]] <- hat2pow(e[[i]])
#   return(e)
# }

hat2pow <- function(x) {
  if(is.call(x)) {
    # Check if the operator is "^" and replace with "pow"
    if(identical(x[[1L]], as.name("^"))) {
      x[[1L]] <- as.name("pow")
    }
    # Recursively travel through the expression
    if(length(x) > 1L) {
      x[2L:length(x)] <- lapply(x[2L:length(x)], hat2pow)
    }
  }
  
  return(x)
}
