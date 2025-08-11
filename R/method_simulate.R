#######################################################
# MAIN SIMULATE FUNCTION THAT CALL OTHERS
#######################################################
perform_simulation <- function(self, private, use.cpp, n.sims){
  
  
  # Laplace cant produce predictions currently...
  if(private$method=="laplace"){
    stop("Simulations arent available for laplace")
  }
  
  if(private$method=="laplace.thygesen"){
    stop("Simulations arent available for laplace")
  }
  
  # time prediction
  comptime <- system.time(
    {
      
      # Predict with Rcpp implementation
      if(use.cpp){
        
        if(!private$silent) message("Simulating with C++...")
        
        
        if(any(private$method == c("lkf","lkf.cpp"))){
          stop("Simulations arent available for LKF")
        }
        
        if(any(private$method == c("ekf","ekf.cpp"))){
          ekf_rcpp_simulation(private$pars, self, private, n.sims)
        }
        
        if(any(private$method == c("ukf","ukf.cpp"))){
          stop("Simulations arent available for UKF")
        }
        
        # Predict with R implementation
      } else {
        
        if(!private$silent) message("Predicting with R...")
        
        
        if(any(private$method == c("lkf","lkf.cpp"))){
          stop("Simulations arent available for LKF")
        }
        
        if(any(private$method == c("ekf","ekf.cpp"))){
          ekf_r_simulation(private$pars, self, private, n.sims)
        }
        
        if(any(private$method == c("ukf","ukf.cpp"))){
          stop("Simulations arent available for UKF")
        }
        
      }
      
    }, gcFirst = FALSE)
  
  private$timer_prediction <- comptime
  
  return(invisible(self))
}

#######################################################
# MAIN RETURN PREDICTIONS FUNCTION THAT CALL OTHERS
#######################################################
create_return_simulation <- function(return.k.ahead, n.sims, self, private){
  
  if(!private$silent) message("Returning results...")
  
  if(any(private$method == c("lkf","lkf.cpp"))){
    
  }
  
  if(any(private$method == c("ekf","ekf.cpp"))){
    create_ekf_simulation_return(return.k.ahead, n.sims, self, private)
  }
  
  if(any(private$method == c("ukf","ukf.cpp"))){
    
  }
  
  return(invisible(self))
  
}
