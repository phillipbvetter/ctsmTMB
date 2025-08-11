#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

perform_filtering = function(self, private, use.cpp){
  
  if(!private$silent) message("Filtering...")
  
  if(private$method == c("laplace")){
    stop("The Laplace method is a smoothing method. Use 'smooth' method instead.")
  }
  
  if(private$method == c("laplace.thygesen")){
    stop("The Laplace method is a smoothing method. Use 'smooth' method instead.")
  }
  
  
  # Linear Kalman Filter
  if(private$method == "lkf"){
    if(use.cpp) {
      stop("C++ filtering is not available for the LKF method")
    } else {
      filt <- lkf_filter_r(private$pars, self, private)
    }
  }
  
  # Extended Kalman Filter
  if(private$method == "ekf"){
    if(use.cpp) {
      filt <- ekf_filter_rcpp(self, private)
    } else {
      filt <- ekf_filter_r(private$pars, self, private)
    }
  }
  
  # Unscented Kalman Filter
  if(private$method == "ukf"){
    if(use.cpp) {
      stop("C++ filtering is not available for the UKF method")
    } else {
      filt <- ukf_filter_r(private$pars, self, private)
    }
  }
  
  private$filt <- filt
  
  return(invisible(self))
}

ekf_filter_rcpp = function(self, private){
  
  if(!any(private$ode.solver==c(1,2))){
    stop("Filtering using C++ currently only support 'euler' or 'rk4' ODE solvers.") 
  }
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # predict using c++ function
  mylist <- ekf_filter_cpp(private$rcpp_function_ptr$f,
                           private$rcpp_function_ptr$g,
                           private$rcpp_function_ptr$dfdx,
                           private$rcpp_function_ptr$h,
                           private$rcpp_function_ptr$dhdx,
                           private$rcpp_function_ptr$hvar,
                           obsMat,
                           inputMat,
                           private$pars,
                           private$initial.state$p0,
                           private$initial.state$x0,
                           private$ode.timestep.size,
                           private$ode.timesteps,
                           numeric_is_not_na_obsMat,
                           number_of_available_obs,
                           private$number.of.states,
                           private$number.of.observations,
                           private$ode.solver)
  
  ####### RETURN #######
  return(invisible(mylist))
}


create_filter_results <- function(self, private, laplace.residuals){
  
  if(!private$silent) message("Returning results...")
  
  filt <- list()
  
  if(any(private$method == c("lkf","ekf","ukf"))){
    
    # extract results from previous
    rep <- private$filt
    
    # States -----------------------------------
    
    # Prior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPrior))))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPrior,diag)))))
    
    colnames(temp.states) = c("t", private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    filt$states$mean$prior = as.data.frame(temp.states)
    filt$states$sd$prior = as.data.frame(temp.sd)
    # covariance
    filt$states$cov$prior = as.data.frame(cbind(private$data$t, do.call(rbind,lapply(rep$pPrior, c))))
    cov.names <- with(
      expand.grid(private$state.names, private$state.names), 
      paste0(paste0("cov_",Var1),"_",Var2)
    )
    names(filt$states$cov$prior) = c("t",cov.names)
    
    # Posterior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, t(do.call(cbind,rep$xPost))))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPost,diag)))))
    colnames(temp.states) = c("t",private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    filt$states$mean$posterior = as.data.frame(temp.states)
    filt$states$sd$posterior = as.data.frame(temp.sd)
    # covariance
    filt$states$cov$posterior = as.data.frame(cbind(private$data$t, do.call(rbind,lapply(rep$pPost, c))))
    cov.names <- with(
      expand.grid(private$state.names, private$state.names), 
      paste0(paste0("cov_",Var1),"_",Var2)
    )
    names(filt$states$cov$posterior) = c("t",cov.names)
    
    # Residuals -----------------------------------
    rowNAs = as.matrix(!is.na(private$data[private$obs.names]))
    sumrowNAs = rowSums(rowNAs)
    
    nr <- nrow(private$data)
    temp.res = matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
    temp.var =  matrix(nrow=length(private$data$t), ncol=private$number.of.observations)
    
    nr <- nrow(private$data)
    innovation = rep$Innovation
    innovation.cov = rep$InnovationCovariance
    names(innovation.cov) = paste("t = ", private$data$t, sep="")
    for (i in 1:nr) {
      if (sumrowNAs[i] > 0) {
        temp.res[i,rowNAs[i,]] = innovation[[i]]
        temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
      }
    }
    temp.res <- data.frame(private$data$t, temp.res)
    temp.sd = data.frame(private$data$t, try_withWarningRecovery(sqrt(temp.var)))
    names(temp.res) = c("t", private$obs.names)
    names(temp.sd) = c("t", private$obs.names)
    
    filt$residuals$residuals = temp.res
    filt$residuals$sd = temp.sd
    filt$residuals$normalized = temp.res
    filt$residuals$normalized[,-1] = temp.res[,-1]/temp.sd[,-1]
    filt$residuals$cov = innovation.cov
    
  }
  
  private$filt <- filt
  
  return(invisible(self))
}

