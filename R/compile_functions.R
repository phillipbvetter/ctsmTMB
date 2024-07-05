perform_compilation = function(self, private, type="estimation"){
  
  # estimation
  if(type=="estimation"){
    compile_cppfile(self, private)
  }
  
  
  # prediction or simulation
  if(any(type==c("prediction","simulation"))){
    compile_rcpp_functions(self, private)
  }
  
  return(invisible(NULL))
}

#######################################################
# COMPILE CPP FUNCTION
#######################################################

# This function decides whether to compile the C++ function, and performs
# reloading of the dynamic library.
# NOTE::: The manual dyn.unload dyn.load must be performed, otherwise TMB
# will reuse the old model. This is a TMB bug. 

compile_cppfile = function(self, private) {
  
  # Start Check:
  # - Exit if the method uses RTMB and does not need C++ compilation.
  bool = any(private$method == c("laplace", "ekf_rtmb"))
  if(bool){
    return(invisible(self))
  }
  
  ############################################################################
  ############################################################################
  
  # If the user requested a compilaton
  if(private$compile){
    
    # Create folder if it doesnt exist
    if(!dir.exists(private$cppfile.directory)){
      dir.create(private$cppfile.directory)
    } 
    
    # Create the C++ file
    write_method_cppfile(self, private)
    
    # Compile the C++ file with TMB
    was.there.an.error = tryCatch(
      TMB::compile(file = paste(private$cppfile.path.with.method,".cpp",sep=""), 
                   framework = "TMBad",
                   openmp = TRUE
      ),
      error = function(e){
        message("----------------------")
        message("A compilation error occured. Here is the original error message: \n\t", 
                conditionMessage(e))
        message("Do you perhaps use spaces in the file-path? Please remove these.")
        return(TRUE)
      })
    if(was.there.an.error){
      stop("Compilation stopped.")
    }
    
    # reload shared libraries
    # Suppress TMB output 'removing X pointers' with capture.output
    utils::capture.output(try(dyn.unload(TMB::dynlib(private$cppfile.path.with.method)),silent=T))
    utils::capture.output(try(dyn.load(TMB::dynlib(private$cppfile.path.with.method)),silent=T))
  }
  
  # If compilation not requested
  if (!private$compile) {
    
    # Unix/MAC-OS platforms: check that the C++ file exists
    if (.Platform$OS.type=="unix") {
      
      # get the compiled dynamic library file
      modelpath.so = paste(private$cppfile.path.with.method,".so",sep="")
      
      # If the model does not exist, then set compile flag and call function again
      if (!file.exists(modelpath.so)) {
        private$set_compile(TRUE)
        compile_cppfile(self, private)
      }
    }
    
    # Windows platforms : check that the C++ file exists
    if (.Platform$OS.type=="windows") {
      
      # get the compiled dynamic library file
      modelpath.dll = paste(private$cppfile.path.with.method,".dll",sep="")
      
      # If the model does not exist, then set compile flag and call function again
      if (!file.exists(modelpath.dll)) {
        private$set_compile(TRUE)
        compile_cppfile(self, private)
      }
    }
    
    # reload shared libraries
    # Suppress TMB output 'removing X pointers' with capture.output
    utils::capture.output(try(dyn.unload(TMB::dynlib(private$cppfile.path.with.method)),silent=T))
    utils::capture.output(try(dyn.load(TMB::dynlib(private$cppfile.path.with.method)),silent=T))
  }
  return(invisible(self))
}

compile_rcpp_functions = function(self, private){
  
  private$Rcppfunction_f <- RcppXPtrUtils::cppXPtr(private$Rcppfunction_f, depends=c("RcppEigen","ctsmTMB"))
  
  private$Rcppfunction_dfdx <- RcppXPtrUtils::cppXPtr(private$Rcppfunction_dfdx, depends="RcppEigen")
  
  private$Rcppfunction_g <- RcppXPtrUtils::cppXPtr(private$Rcppfunction_g, depends=c("RcppEigen"))
  
  private$Rcppfunction_h <- RcppXPtrUtils::cppXPtr(private$Rcppfunction_h, depends="RcppEigen")
  
  private$Rcppfunction_dhdx <- RcppXPtrUtils::cppXPtr(private$Rcppfunction_dhdx, depends="RcppEigen")
  
  private$Rcppfunction_hvar <- RcppXPtrUtils::cppXPtr(private$Rcppfunction_hvar, depends="RcppEigen")
  
}
