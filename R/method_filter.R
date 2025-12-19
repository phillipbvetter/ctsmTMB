#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

perform_filtering = function(self, private, use.cpp){
  
  if(private$method == c("laplace")){
    stop("The Laplace method is a smoothing method. Use 'smooth' method instead.")
  }
  
  if(private$method == c("laplace.thygesen")){
    stop("The Laplace method is a smoothing method. Use 'smooth' method instead.")
  }
  
  comptime <- system.time({
      if(use.cpp){
        
        if(!private$silent) message("Filtering with C++...")
        lkf_ekf_ukf_filter_rcpp(private$pars, self, private)
        
      } else {
        
        if(!private$silent) message("Filtering with R...")
        lkf_ekf_ukf_filter_r(private$pars, self, private)
        
      }
    }, gcFirst = FALSE)
  
  private$timer_filtration <- comptime
  
  return(invisible(self))
}

lkf_ekf_ukf_filter_r <- function(pars, self, private){
  
  filt <- switch(private$method,
                 lkf = lkf_filter_r(pars, self, private),
                 ekf = ekf_filter_r(pars, self, private),
                 ukf = ukf_filter_r(pars, self, private),
  )
  
  private$filtration <- filt
  
  return(invisible(self))
}

lkf_ekf_ukf_filter_rcpp <- function(pars, self, private){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  output.list <- list()
  if(private$method == "lkf"){
    output.list <- lkf_filter_rcpp(private$rcpp_function_ptr$f,
                                   private$rcpp_function_ptr$g,
                                   private$rcpp_function_ptr$dfdx,
                                   private$rcpp_function_ptr$h,
                                   private$rcpp_function_ptr$dhdx,
                                   private$rcpp_function_ptr$hvar,
                                   private$rcpp_function_ptr$dfdu,
                                   obsMat,
                                   inputMat,
                                   pars,
                                   private$initial.state$p0,
                                   private$initial.state$x0,
                                   private$ode.timestep.size,
                                   numeric_is_not_na_obsMat,
                                   number_of_available_obs)
  }
  if(private$method=="ekf"){
    output.list <- ekf_filter_rcpp(private$rcpp_function_ptr$f,
                                   private$rcpp_function_ptr$g,
                                   private$rcpp_function_ptr$dfdx,
                                   private$rcpp_function_ptr$h,
                                   private$rcpp_function_ptr$dhdx,
                                   private$rcpp_function_ptr$hvar,
                                   obsMat,
                                   inputMat,
                                   pars,
                                   private$initial.state$p0,
                                   private$initial.state$x0,
                                   private$ode.timestep.size,
                                   private$ode.timesteps,
                                   numeric_is_not_na_obsMat,
                                   number_of_available_obs,
                                   private$ode.solver)
  }
  if(private$method == "ukf"){
    output.list <- ukf_filter_rcpp(private$rcpp_function_ptr$f,
                                   private$rcpp_function_ptr$g,
                                   private$rcpp_function_ptr$dfdx,
                                   private$rcpp_function_ptr$h,
                                   private$rcpp_function_ptr$dhdx,
                                   private$rcpp_function_ptr$hvar,
                                   obsMat,
                                   inputMat,
                                   pars,
                                   private$initial.state$p0,
                                   private$initial.state$x0,
                                   private$ode.timestep.size,
                                   private$ode.timesteps,
                                   numeric_is_not_na_obsMat,
                                   number_of_available_obs,
                                   private$ukf_hyperpars,
                                   private$ode.solver)
  }
  
  private$filtration <- output.list
  
  return(invisible(self))
}

create_filter_results <- function(self, private, laplace.residuals){
  
  if(!private$silent) message("Returning results...")
  
  filt <- list()
  
  if(any(private$method == c("lkf","ekf","ukf"))){
    
    # extract results from previous
    rep <- private$filtration
    
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
  
  private$filtration <- filt
  
  return(invisible(self))
}

