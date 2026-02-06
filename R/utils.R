###########################################################
# Shortcut for model creation
###########################################################
#' Create a ctsmTMB model faster avoiding $...
#' @returns Print of ctsmTMB model object
#' @export
newModel <- function(){
  ctsmTMB$new()
}

###########################################################
# TRY THAT SUPPRESSES AND RECOVERS ANSWER DESPITE WARNINGS
###########################################################

try_with_warning_recovery = function(expr, silent=TRUE){
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

###########################################################
# set colum names of matrix
###########################################################
set.colnames <- function(object, colnames){
  colnames(object) <- colnames
  object
}
