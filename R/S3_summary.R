#' Basic summary of ctsmTMB fit object
#' @param object a ctsmTMB fit object
#' @param correlation boolean indicating whether or not to display the
#' parameter correlation structure
#' @param ... additional arguments
#' @returns a summary of the estimated ctsmTMB model fit
#' @export
summary.ctsmTMB.fit = function(object,
                               correlation = FALSE,
                               ...) {

  fit <- object #method consistency (argument must be called 'object')

  if (!is.logical(correlation)) {
    stop("correlation must be logical")
  }

  mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)
  if (correlation){
    cat("\nCorrelation Matrix\n")
    S = diag(1/fit$sd.fixed)
    cor = S %*% fit$cov.fixed %*% S
    rownames(cor) = names(fit$par.fixed)
    colnames(cor) = names(fit$par.fixed)
    cor <- format(round(cor, 2), nsmall = 2, digits = 6)
    cor[!lower.tri(cor,diag=T)] <- ""
    print(cor, drop = FALSE, quote = FALSE)
    # cor[!lower.tri(cor)] <- ""
    # print(cor[-1, -(dim(cor)[1]), drop = FALSE], quote = FALSE)
  }

  return(invisible(list(parameters=mat)))
}
