#######################################################
# MAIN PREDICTION FUNCTION THAT CALL OTHERS
#######################################################
perform_prediction <- function(self, private, use.cpp){
  
  
  # Laplace cant produce predictions currently...
  if(private$method=="laplace"){
    stop("Predictions arent available for laplace")
  }
  
  if(private$method=="laplace2"){
    stop("Predictions arent available for laplace")
  }
  
  # time prediction
  comptime <- system.time(
    {
      
      # Predict with Rcpp implementation
      if(use.cpp){
        
        if(!private$silent) message("Predicting with C++...")
        
        
        if(any(private$method == c("lkf","lkf_cpp"))){
          stop("Predictions arent available for LKF")
        }
        
        if(any(private$method == c("ekf","ekf_cpp"))){
          ekf_rcpp_prediction(self, private)
        }
        
        if(any(private$method == c("ukf","ukf_cpp"))){
          stop("Predictions arent available for UKF")
        }
        
        # Predict with R implementation
      } else {
        
        if(!private$silent) message("Predicting with R...")
        
        
        if(any(private$method == c("lkf","lkf_cpp"))){
          stop("Predictions arent available for LKF")
        }
        
        if(any(private$method == c("ekf","ekf_cpp"))){
          ekf_r_prediction(self, private)
        }
        
        if(any(private$method == c("ukf","ukf_cpp"))){
          stop("Predictions arent available for UKF")
        }
        
      }
      
    }, gcFirst = FALSE)
  
  private$timer_prediction <- comptime
  
  return(invisible(self))
}

#######################################################
# MAIN RETURN PREDICTIONS FUNCTION THAT CALL OTHERS
#######################################################
create_return_prediction <- function(return.covariance, return.k.ahead, self, private){
  
  if(!private$silent) message("Returning results...")
  
  if(any(private$method == c("lkf","lkf_cpp"))){
    
  }
  
  if(any(private$method == c("ekf","ekf_cpp"))){
    create_ekf_predict_return(return.covariance, return.k.ahead, self, private)
  }
  
  if(any(private$method == c("ukf","ukf_cpp"))){
    
  }
  
  return(invisible(self))
  
}
