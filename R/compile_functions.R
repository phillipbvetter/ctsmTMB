compile_rcpp_functions = function(self, private){
  
  .depends <- c("RcppEigen", "ctsmTMB")
  
  private$rcpp_function_ptr$f <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$f, depends=.depends)
  
  private$rcpp_function_ptr$dfdx <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$dfdx, depends=.depends)
  
  private$rcpp_function_ptr$g <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$g, depends=.depends)
  
  private$rcpp_function_ptr$h <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$h, depends=.depends)
  
  private$rcpp_function_ptr$dhdx <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$dhdx, depends=.depends)
  
  private$rcpp_function_ptr$hvar <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$hvar, depends=.depends)
  
}
