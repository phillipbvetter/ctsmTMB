#######################################################
# Print - S3 Method
#######################################################

#' Basic print of ctsmTMB objects
#' @param x  an object of class 'ctsmTMB'
#' @param ... additional arguments (not in use)
#'
#' @returns Print of ctsmTMB model object
#' @export
print.ctsmTMB = function(x,...) {
  obj <- x
  output <- obj$print()

  return(invisible(output))
}

#' Basic print of objects ctsmTMB fit objects
#' @param x a ctsmTMB fit object
#' @param ... additional arguments
#'
#' @returns Print of ctsmTMB fit object
#' @export
print.ctsmTMB.fit = function(x,...) {
  fit <- x
  if(is.null(fit$sd.fixed)) {
    fit$sd.fixed <- rep(NA,length(fit$par.fixed))
    fit$tvalue <- rep(NA,length(fit$par.fixed))
    fit$Pr.tvalue <- rep(NA,length(fit$par.fixed))
  }
  mat = cbind(fit$par.fixed, fit$sd.fixed, fit$tvalue, fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)

  return(invisible(mat))
}

