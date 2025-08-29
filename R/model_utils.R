#######################################################
# CHANGE FROM R POWER NOTATION TO C++
#######################################################

# Changes powers in expressions from R-style (e.g. x^2) to pow(x,2) which is interpretable by C++
# from https://stackoverflow.com/questions/40606723/substitute-the-power-symbol-with-cs-pow-syntax-in-mathematical-expression

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
  ), silent=silent)
  
  return(output)
}

###########################################################
# Upper-case first letter of string
###########################################################
capitalize_first <- function(s) {
  if (nchar(s) == 0) return(s)  # Handle empty string case
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

#######################################################
# GGPLOT2 FUNCTIONS FOR USE IN PLOTS
#######################################################

getggplot2theme = function() {
  mytheme =
    ggplot2::theme_minimal() +
    ggplot2::theme(
      # text = ggplot2::element_text("Avenir Next Condensed",size=12),
      text = ggplot2::element_text(size=12),
      legend.text = ggplot2::element_text(size=12),
      axis.text = ggplot2::element_text(size=12),
      strip.text = ggplot2::element_text(face="bold",size=12),
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.position = "top",
      plot.title = ggplot2::element_text(hjust=0.5),
      plot.subtitle = ggplot2::element_text(hjust=0.5)
    )
  return(mytheme)
}

# getggplot2colors = function(n) {
#   hues = seq(15, 375, length = n + 1)
#   ggcolors = grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
#   return(ggcolors)
# }

###########################################################
# Shortcut for model creation
###########################################################
#' Create a ctsmTMB model faster avoiding $...
#' @returns Print of ctsmTMB model object
#' @export
new.model <- function(){
  ctsmTMB$new()
}
