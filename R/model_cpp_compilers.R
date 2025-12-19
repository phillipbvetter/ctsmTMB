#######################################################
# COMPILE RCPP DRIFT DIFFUSION ETC FUNCTIONS
#######################################################
compile_rcpp_functions = function(self, private){
  
  if(!private$silent) message("Compiling C++ function pointers...")
  
  # # Settings
  # .depends <- c("Rcpp", "RcppEigen", "ctsmTMB")
  # 
  # # COMMENT:
  # # These .includes should in principle be added via the .depends = 'ctsmTMB'.
  # # This automatically fetches the code in inst/include/ctsmTMB.h
  # # This does not work however when using devtools::load_all.
  # .includes <- c(
  #   "double invlogit(double x){return 1/(1 + exp(-x));}",
  #   "const double pi = 3.14159265358979323846;"
  # )
  # 
  # # Create XPtr's
  # print(system.time(
  #   outlist <- lapply(private$rcpp.function.strings, 
  #                     function(s) RcppXPtrUtils::cppXPtr(s, 
  #                                                        depends=.depends, 
  #                                                        includes = .includes)
  #   )
  # ))
  # 
  # # Add to private fields
  # nams <- c("f","dfdx","g","h","dhdx","hvar","dfdu")
  # names(outlist) <- nams
  # private$rcpp_function_ptr[nams] <- outlist[nams]
  
  # Read lines from inst/include template
  txt <- readLines(system.file("include/template_user_functions.h", package="ctsmTMB"))
  # Insert our created system function strings
  txt[which(txt %in% "// INSERT F")] <- private$rcpp.function.strings$f
  txt[which(txt %in% "// INSERT DFDX")] <- private$rcpp.function.strings$dfdx
  txt[which(txt %in% "// INSERT G")] <- private$rcpp.function.strings$g
  txt[which(txt %in% "// INSERT H")] <- private$rcpp.function.strings$h
  txt[which(txt %in% "// INSERT DHDX")] <- private$rcpp.function.strings$dhdx
  txt[which(txt %in% "// INSERT HVAR")] <- private$rcpp.function.strings$hvar
  txt[which(txt %in% "// INSERT DFDU")] <- private$rcpp.function.strings$dfdu
  
  # We compile using 'code' over 'file' in sourceCpp. 
  # This is better for several reasons:
    # 1. The caching is automatically handled, so second time is faster
    # 2. Rstudio does not enter the "Source Cpp" tab automatically (annoying for users)
    # 3. This solution is faster than using RcppXptrUtils
  Rcpp::sourceCpp(
    code=paste(txt, collapse=" \n "), 
    verbose=FALSE, 
    showOutput = FALSE
    )
  private$rcpp_function_ptr <- get_sysfun_cpp_function_ptrs()
  
  # return
  return(invisible(self))
  
}

#######################################################
# COMPILE CPP FUNCTION
#######################################################

# This function decides whether to compile the C++ function, and performs
# reloading of the dynamic library.
# NOTE::: The manual dyn.unload dyn.load must be performed, otherwise TMB
# will reuse the old model. This is a TMB bug. 

compile_cppfile <- function(self, private) {
  
  # Start Check:
  # - Exit if the method uses RTMB and does not need C++ compilation.
  bool <- any(private$method == c("lkf", "ekf", "ukf", "laplace", "laplace.thygesen"))
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
    write_cppfile(self, private)
    
    # Compile C++ File
    # we need optimization level 1 to avoid compilation errors on windows(?)
    
    comptime <- system.time(
      {
        if (.Platform$OS.type=="windows") {
          out <- tryCatch(
            TMB::compile(file = paste(private$cppfile.path.with.method,".cpp",sep=""),
                         flags = paste("-O1",paste0("-I", system.file("include", package = "ctsmTMB"))),
                         framework = "TMBad",
                         openmp = TRUE),
            error = function(e){
              message("----------------------")
              message("A compilation error occured with the following error message: \n\t", 
                      conditionMessage(e))
              return(e)
            })
        }
        
        # No optimization flag needed on unix
        if (.Platform$OS.type=="unix") {
          out <- tryCatch(
            TMB::compile(file = paste(private$cppfile.path.with.method,".cpp",sep=""),
                         flags = paste0("-I", system.file("include", package = "ctsmTMB")),
                         # for some reason the flag above disabled -02 which is good for
                         # compilation speed (from 40s to 10 s)
                         framework = "TMBad",
                         openmp = TRUE),
            error = function(e){
              message("----------------------")
              message("A compilation error occured with the following error message: \n\t", 
                      conditionMessage(e))
              return(e)
            })
        }
        
      }, gcFirst = FALSE)
    
    private$timer_cppbuild <- comptime
    
    if(inherits(out,"error")){
      stop("Stopping because compilation failed.")
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
      file.end <- ".so"
    }
    if (.Platform$OS.type=="windows") {
      file.end <- ".dll"
    }
    
    # get the compiled dynamic library file
    model.dyn.path <- paste0(private$cppfile.path.with.method, file.end)
    
    # If the model exists, check that dimensions are correct
    if (file.exists(model.dyn.path)) {
      model.cpp.path <- paste0(private$cppfile.path.with.method, ".cpp")
      out <- readLines(model.cpp.path)
      model.dims <- as.numeric(stringr::str_extract(out[1:4], ".*:(\\d+)", group=1))
      bool <- !all(model.dims[1] == private$number.of.states,
                   model.dims[2] == private$number.of.observations,
                   model.dims[3] == private$number.of.inputs,
                   model.dims[4] == private$number.of.pars
      )
      # if the dimensions are wrong recompile the c++ file
      if(bool){
        message("Recompiling C++ file because the current is inconsistent with the specified model.")
        private$set_compile(TRUE)
        compile_cppfile(self, private)
      }
    }
    
    # If the model does not exist, then set compile flag and call function again
    if (!file.exists(model.dyn.path)) {
      private$set_compile(TRUE)
      compile_cppfile(self, private)
    }
    
    # reload shared libraries
    # Suppress TMB output 'removing X pointers' with capture.output
    utils::capture.output(try(dyn.load(TMB::dynlib(private$cppfile.path.with.method)), silent=T))
  }
  
  # return
  return(invisible(self))
}
