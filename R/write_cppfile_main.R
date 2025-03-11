#
#
#
# This is the main writer function for creating the likelihood C++ file
write_method_cppfile = function(self, private) {
  
  #Initialize C++ file
  fileconn = file(paste0(private$cppfile.path.with.method,".cpp"))
  
  
  # Header etc...
  txt = c(
    "#include <TMB.hpp>",
    '#include <ctsmTMB_customfuns.hpp>',
    "using namespace density;"
  )
  
  # Define constants
  newtxt = "const double pi = M_PI;"
  txt = c(txt, newtxt)
  
  # Custom extra functions
  # newtxt = write_custom_functions()
  # txt = c(txt, newtxt)
  
  # Various helper functions
  if(any(private$method == c("ekf_cpp","ukf"))){
    newtxt = write_helper_cppfunctions()
    txt = c(txt, newtxt)
  }
  
  # Specific method functions
  if(private$method=="ekf_cpp"){
    newtxt = write_ekf_functions(self, private)
    txt = c(txt, newtxt)
  }
  if(private$method=="ukf"){
    newtxt = write_ukf_functions(self,private)
    txt = c(txt,newtxt)
  }
  
  # Initialize TMB Objective Function
  
  txt = c(txt,"template<class Type>\nType objective_function<Type>::operator() ()\n{")
  txt = c(txt, "Type nll__ = 0;")
  
  # Specific estimation method
  if(private$method=="ekf_cpp"){
    newtxt = write_ekf_estimate(self, private)
    txt = c(txt, newtxt)
  }
  if(private$method=="ukf"){
    newtxt = write_ukf_estimate(self, private)
    txt = c(txt, newtxt)
  }
  
  # Return statement
  txt = c(txt, "return nll__;")
  txt = c(txt, "}")
  
  # Write cpp function and close file connection
  writeLines(txt,fileconn)
  close(fileconn)
  
  # Return
  return(invisible(self))
}
