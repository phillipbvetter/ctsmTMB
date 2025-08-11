

get_substitution_list = function(self, private){
  
  # Create substitution translation list
  obsList = lapply(seq_along(private$obs.names), function(id) substitute(obsVec(i),list(i=as.numeric(id-1))))
  parList = lapply(seq_along(private$parameter.names), function(id) substitute(parVec(i),list(i=as.numeric(id-1))))
  stateList = lapply(seq_along(private$state.names), function(id) substitute(stateVec(i),list(i=as.numeric(id-1))))
  inputList = lapply(seq_along(private$input.names), function(id) substitute(inputVec(i),list(i=as.numeric(id-1))))
  names(obsList) = private$obs.names
  names(parList) = private$parameter.names
  names(stateList) = private$state.names
  names(inputList) = private$input.names
  subsList = c(obsList, parList, stateList, inputList)
  
  return(subsList)
}

##################################################
# drift
##################################################

write_f = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  f <- c()
  for(i in seq_along(private$state.names)){
    drift.term <- Deriv::Simplify(private$diff.terms[[i]]$dt)
    if(!(drift.term==0)){
      drift.term = hat2pow(private$diff.terms[[i]]$dt)
      new.drift.term = do.call(substitute, list(drift.term, subsList))
      f <- c(f, sprintf("f__(%i) = %s;",i-1, deparse1(new.drift.term)))
    }
  }
  
  newtxt = "\n//////////// drift function //////////
  template<class Type>
  vector<Type> f__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> f__(%s);
    f__.setZero();
    %s
    return f__;
  }"
  
  newtxt = sprintf(newtxt, private$number.of.states, paste(f,collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# drift jacobian
##################################################

write_jac_f = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  jac.f = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term <- Deriv::Simplify(private$diff.terms.drift[[i]][[j]])
      if(!(term==0)){
        term = hat2pow(term)
        new.term = do.call(substitute, list(term, subsList))
        jac.f = c(jac.f, sprintf("dfdx__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
      }
    }
  }
  
  newtxt = "\n//////////// jacobian of drift function ///////////
  template<class Type>
  matrix<Type> dfdx__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> dfdx__(%s, %s);
    dfdx__.setZero();
    %s
    return dfdx__;
  }"
  
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.states, paste(jac.f, collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# diffusion
##################################################

write_g = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  g = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term <- Deriv::Simplify(private$diff.terms[[i]][[j+1]])
      if(!(term==0)){
        term = hat2pow(term)
        new.term = do.call(substitute, list(term, subsList))
        g = c(g, sprintf("g__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
      }
    }
  }
  
  newtxt = "\n//////////// diffusion function ///////////
  template<class Type>
  matrix<Type> g__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> g__(%s, %s);
    g__.setZero();
    %s
    return g__;
  }"
  
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.diffusions, paste(g,collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# observation
##################################################

write_h = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  h <- c()
  for(i in seq_along(private$obs.names)){
    term <- Deriv::Simplify(private$obs.eqs.trans[[i]]$rhs)
    if(term!=0){
      term = hat2pow(term)
      new.term = do.call(substitute, list(term, subsList))
      h <- c(h, sprintf("h__(%s) = %s;",i-1, deparse1(new.term)))
    }
  }
  
  newtxt = "\n//////////// observation function ///////////
  template<class Type>
  vector<Type> h__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> h__(%s);
    h__.setZero();
    %s
    return h__;
  }"
  
  newtxt = sprintf(newtxt, private$number.of.observations, paste(h,collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# observation jacobian
##################################################

write_jac_h = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  jac.h = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = Deriv::Simplify(private$diff.terms.obs[[i]][[j]])
      if(term!=0){
        term = hat2pow(term)
        new.term = do.call(substitute, list(term, subsList))
        jac.h = c(jac.h, sprintf("dhdx__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
      }
    }
  }
  
  newtxt = "\n//////////// jacobian of obs function ///////////
  template<class Type>
  matrix<Type> dhdx__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> dhdx__(%s, %s);
    dhdx__.setZero();
    %s
    return dhdx__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, private$number.of.states, paste(jac.h,collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# derivative w.r.t inputs (for linear kalman filter)
##################################################

write_jac_f_wrt_u <- function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  # 1. We first find constant terms
  # this corresponds to an input that is always 1 (first element of U)
  zero.list <- as.list(numeric(private$number.of.inputs + private$number.of.states))
  names(zero.list) <-c(private$input.names, private$state.names)
  constant.terms <- sapply(private$sys.eqs.trans, function(x) Deriv::Simplify(do.call(substitute, list(x$diff.dt, zero.list))))
  if(inherits(constant.terms, "try-error")){
    constant.terms <- rep(0, private$number.of.states)
  }
  
  jac.f.wrt.u = c()
  # 2. Now we find input-terms by differentiation
  for(i in seq_along(private$state.names)){
    #its inputs + 1 below because the first element is for constants
    for(j in 1:(private$number.of.inputs+1)){
      if(j==1){
        term <- constant.terms[[i]]
      } else {
        term <- ctsmTMB.Deriv(f=private$sys.eqs.trans[[i]]$diff.dt, x=private$input.names[j-1])
      }
      # skip if zero
      if(term=="0") next
      term = hat2pow(term)
      new.term = do.call(substitute, list(term, subsList))
      jac.f.wrt.u = c(jac.f.wrt.u, sprintf("dfdu__(%s, %s) = %s;", i-1, j-1, deparse1(new.term)))
    }
  }
  
  newtxt = "\n//////////// jacobian of drift function wrt inputs ///////////
  template<class Type>
  matrix<Type> dfdu__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> dfdu__(%s, %s);
    dfdu__.setZero();
    %s
    return dfdu__;
  }"
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.inputs+1, paste(jac.f.wrt.u, collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# observation variance
##################################################

write_h_var = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  hvar <- c()
  for(i in seq_along(private$obs.var.trans)){
    term <- Deriv::Simplify(private$obs.var.trans[[i]]$rhs)
    if(term!=0){
      term <- hat2pow(term)
      new.term = do.call(substitute, list(term, subsList))
      hvar <- c(hvar, sprintf("hvar__(%s) = %s;", i-1, deparse1(new.term)))
    }
  }
  
  newtxt = "\n//////////// observation variance matrix function ///////////
  template<class Type>
  vector<Type> hvar__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> hvar__(%s);
    hvar__.setZero();
    %s
    return hvar__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(hvar,collapse="\n\t\t"))
  
  return(newtxt)
}

write_cppfile <- function(self, private){
  
  if(private$method == "lkf.cpp") write_lkf_cppfile(self, private)
  if(private$method == "ekf.cpp") write_ekf_cppfile(self, private)
  if(private$method == "ukf.cpp") write_ukf_cppfile(self, private)
  
}

write_lkf_cppfile <- function(self, private){
  
  # Get template
  txt <- readLines(system.file("include/template_lkf.h", package="ctsmTMB"))
  
  # Embed system info
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_STATES")] <- sprintf("// STATES:%s", private$number.of.states)
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_OBS")] <- sprintf("// OBS:%s", private$number.of.observations)
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_INPUTS")] <- sprintf("// INPUTS:%s", private$number.of.inputs)
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_PARS")] <- sprintf("// PARS:%s", private$number.of.pars)
  
  # Insert user functions
  txt[which(txt %in% "// INSERT F")] <- write_f(self, private)
  txt[which(txt %in% "// INSERT DFDX")] <- write_jac_f(self, private)
  txt[which(txt %in% "// INSERT G")] <- write_g(self, private)
  txt[which(txt %in% "// INSERT H")] <- write_h(self, private)
  txt[which(txt %in% "// INSERT HVAR")] <- write_h_var(self, private)
  txt[which(txt %in% "// INSERT DHDX")] <- write_jac_h(self, private)
  txt[which(txt %in% "// INSERT DFDU")] <- write_jac_f_wrt_u(self, private)
  
  # Open, write and close new file
  fileconn = file(paste0(private$cppfile.path.with.method,".cpp"))
  writeLines(txt, fileconn)
  close(fileconn)
  
}

write_ekf_cppfile <- function(self, private){
  
  # Get template
  txt <- readLines(system.file("include/template_ekf.h", package="ctsmTMB"))
  
  # Embed system info
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_STATES")] <- sprintf("// STATES:%s", private$number.of.states)
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_OBS")] <- sprintf("// OBS:%s", private$number.of.observations)
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_INPUTS")] <- sprintf("// INPUTS:%s", private$number.of.inputs)
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_PARS")] <- sprintf("// PARS:%s", private$number.of.pars)
  
  # Insert user functions
  txt[which(txt %in% "// INSERT F")] <- write_f(self, private)
  txt[which(txt %in% "// INSERT DFDX")] <- write_jac_f(self, private)
  txt[which(txt %in% "// INSERT G")] <- write_g(self, private)
  txt[which(txt %in% "// INSERT H")] <- write_h(self, private)
  txt[which(txt %in% "// INSERT DHDX")] <- write_jac_h(self, private)
  txt[which(txt %in% "// INSERT HVAR")] <- write_h_var(self, private)
  
  # Open, write and close new file
  fileconn = file(paste0(private$cppfile.path.with.method,".cpp"))
  writeLines(txt, fileconn)
  close(fileconn)
  
}

write_ukf_cppfile <- function(self, private){
  
  # Get template
  txt <- readLines(system.file("include/template_ukf.h", package="ctsmTMB"))
  
  # Embed system info
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_STATES")] <- sprintf("// STATES:%s", private$number.of.states)
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_OBS")] <- sprintf("// OBS:%s", private$number.of.observations)
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_INPUTS")] <- sprintf("// INPUTS:%s", private$number.of.inputs)
  txt[which(txt %in% "// SYSINFO: NUMBER_OF_PARS")] <- sprintf("// PARS:%s", private$number.of.pars)
  
  # Insert user functions
  txt[which(txt %in% "// INSERT F")] <- write_f(self, private)
  txt[which(txt %in% "// INSERT DFDX")] <- write_jac_f(self, private)
  txt[which(txt %in% "// INSERT G")] <- write_g(self, private)
  txt[which(txt %in% "// INSERT H")] <- write_h(self, private)
  txt[which(txt %in% "// INSERT DHDX")] <- write_jac_h(self, private)
  txt[which(txt %in% "// INSERT HVAR")] <- write_h_var(self, private)
  
  # Open, write and close new file
  fileconn = file(paste0(private$cppfile.path.with.method,".cpp"))
  writeLines(txt, fileconn)
  close(fileconn)
  
}

#######################################################
# MAIN WRITER FUNCTION FOR WRITING C++ FILE
#######################################################
# This is the main writer function for creating the likelihood C++ file
write_cppfile = function(self, private) {
  
  if(private$method == "ekf.cpp"){
    write_ekf_cppfile(self, private)
  }
  
  if(private$method == "lkf.cpp"){
    stop("lkf.cpp is not implemented yet.")
    write_lkf_cppfile(self, private)
  }
  
  if(private$method == "ukf.cpp"){
    write_ukf_cppfile(self, private)
  }
  
  # Return
  return(invisible(self))
}


# #######################################################
# # MAIN WRITER FUNCTION FOR WRITING C++ FILE
# #######################################################
# # This is the main writer function for creating the likelihood C++ file
# write_cppfile = function(self, private) {
#   
#   #Initialize C++ file
#   fileconn = file(paste0(private$cppfile.path.with.method,".cpp"))
#   
#   # Header etc...
#   txt = c(
#     "#include <TMB.hpp>",
#     "using namespace density;"
#   )
#   
#   # Define constants
#   newtxt = "const double pi = M_PI;"
#   txt = c(txt, newtxt)
#   
#   # Custom extra functions
#   newtxt = write_custom_cppfunctions()
#   txt = c(txt, newtxt)
#   
#   # Various helper functions
#   if(any(private$method == c("ekf_cpp","ukf_cpp"))){
#     newtxt = write_helper_cppfunctions()
#     txt = c(txt, newtxt)
#   }
#   
#   # Specific method functions
#   if(private$method=="ekf_cpp"){
#     newtxt = write_ekf_cppfunctions(self, private)
#     txt = c(txt, newtxt)
#   }
#   if(private$method=="ukf_cpp"){
#     newtxt = write_ukf_cppfunctions(self,private)
#     txt = c(txt,newtxt)
#   }
#   
#   # Initialize TMB Objective Function
#   
#   txt = c(txt,"template<class Type>\nType objective_function<Type>::operator() ()\n{")
#   txt = c(txt, "Type nll__ = 0;")
#   
#   # Specific estimation method
#   if(private$method=="ekf_cpp"){
#     newtxt = write_ekf_cpp_likelihood(self, private)
#     txt = c(txt, newtxt)
#   }
#   if(private$method=="ukf_cpp"){
#     newtxt = write_ukf_cpp_likelihood(self, private)
#     txt = c(txt, newtxt)
#   }
#   
#   # Return statement
#   txt = c(txt, "return nll__;")
#   txt = c(txt, "}")
#   
#   # Write cpp function and close file connection
#   writeLines(txt,fileconn)
#   close(fileconn)
#   
#   # Return
#   return(invisible(self))
# }
