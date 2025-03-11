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
  bool <- !any(private$method == c("ekf_cpp", "ukf"))
  if(bool){
    return(invisible(self))
  }
  
  ############################################################################
  ############################################################################
  
  # If the user requested a compilaton
  if(private$compile){
    
    # Create folder if it doesnt exist
    # This is necessary because no folder was created if model$setCppfilesDirectory 
    # was not envoked by the user 
    if(!dir.exists(private$cppfile.directory)){
      dir.create(private$cppfile.directory, recursive=TRUE)
    }
    
    # Write the C++ file
    write_method_cppfile(self, private)
    
    # Compile C++ File
    # we need optimization level 1 to avoid compilation errors on windows(?)
    if (.Platform$OS.type=="windows") {
      comptime = tryCatch(
        system.time(
          TMB::compile(file = paste(private$cppfile.path.with.method,".cpp",sep=""),
                       flags = paste("-O1",paste0("-I", system.file("include", package = "ctsmTMB"))),
                       framework = "TMBad",
                       openmp = TRUE
          )
        ),
        error = function(e){
          message("----------------------")
          message("A compilation error occured with the following error message: \n\t", 
                  conditionMessage(e))
          return(e)
        })
    }
    
    # No optimization flag needed on unix
    if (.Platform$OS.type=="unix") {
      comptime = tryCatch(
        system.time(
          TMB::compile(file = paste(private$cppfile.path.with.method,".cpp",sep=""),
                       flags = paste0("-I", system.file("include", package = "ctsmTMB")),
                       framework = "TMBad",
                       openmp = TRUE
          )
        ),
        error = function(e){
          message("----------------------")
          message("A compilation error occured with the following error message: \n\t", 
                  conditionMessage(e))
          return(e)
        })
    }
    
    if(inherits(comptime,"error")){
      stop("Stopping because compilation failed.")
    }
    
    comptime = format(round(as.numeric(comptime["elapsed"])*1e2)/1e2,digits=5,scientific=F)
    if(!private$silent) message("...took ", comptime, " seconds")
    
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
  
  .depends <- c("RcppEigen", "ctsmTMB")
  
  private$rcpp_function_ptr$f <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$f, depends=.depends)
  
  private$rcpp_function_ptr$dfdx <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$dfdx, depends=.depends)
  
  private$rcpp_function_ptr$g <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$g, depends=.depends)
  
  private$rcpp_function_ptr$h <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$h, depends=.depends)
  
  private$rcpp_function_ptr$dhdx <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$dhdx, depends=.depends)
  
  private$rcpp_function_ptr$hvar <- RcppXPtrUtils::cppXPtr(private$rcpp.function.strings$hvar, depends=.depends)
  
}
