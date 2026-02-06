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
  
  private$filtration.raw <- filt
  
  return(invisible(self))
}

lkf_ekf_ukf_filter_rcpp <- function(pars, self, private){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  ids <- 1:private$number.of.observations - 1 #minus 1 for 0 indexing
  non_na_ids <- apply(obsMat, 1, function(x) ids[!is.na(x)], simplify = FALSE)
  any_available_obs <- sapply(non_na_ids, function(x) !is.null(x))
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  output.list <- list()
  if(private$method == "lkf"){
    output.list <- lkf_filter_rcpp(private$rcpp_function_ptr,
                                   obsMat,
                                   inputMat,
                                   pars,
                                   private$initial.state$p0,
                                   private$initial.state$x0,
                                   private$ode.timestep.size,
                                   any_available_obs,
                                   non_na_ids)
  }
  if(private$method=="ekf"){
    output.list <- ekf_filter_rcpp(private$rcpp_function_ptr,
                                   obsMat,
                                   inputMat,
                                   pars,
                                   private$initial.state$p0,
                                   private$initial.state$x0,
                                   private$ode.timestep.size,
                                   private$ode.timesteps,
                                   any_available_obs,
                                   non_na_ids,
                                   private$ode.solver)
  }
  if(private$method == "ukf"){
    output.list <- ukf_filter_rcpp(private$rcpp_function_ptr,
                                   obsMat,
                                   inputMat,
                                   pars,
                                   private$initial.state$p0,
                                   private$initial.state$x0,
                                   private$ode.timestep.size,
                                   private$ode.timesteps,
                                   numeric_is_not_na_obsMat,
                                   number_of_available_obs,
                                   private$ukf.hyperpars,
                                   private$ode.solver)
  }
  
  private$filtration.raw <- output.list
  
  return(invisible(self))
}

create_filter_results <- function(self, private, laplace.residuals, return.directly=FALSE, silent=private$silent){
  
  if(!silent) message("Returning results...")
  
  
  if(any(private$method == c("lkf","ekf","ukf"))){
    
    rep <- private$filtration.raw
    filt <- list()
    
    ##### helper functions #####
    rbind_vectors <- function(list, colnames=NULL, extra=TRUE){
      if(extra){
        try_with_warning_recovery(
          set.colnames(
            cbind(private$data$t, do.call(rbind, list)),
            colnames
          )
        )
      } else {
        do.call(rbind, list)
      }
    }
    rbind_matrices_flatten <- function(list, colnames=NULL, fn=identity, extra=TRUE){
      if(extra){
        try_with_warning_recovery(
          set.colnames(
            cbind(private$data$t, fn(do.call(rbind, lapply(list,c)))),
            colnames
          )
        )
      } else {
        fn(do.call(rbind, lapply(list,c)))
      }
    }
    rbind_matrices_diag <- function(list, colnames=NULL, fn=identity, extra=TRUE){
      if(extra){
        try_with_warning_recovery(
          set.colnames(
            cbind(private$data$t, fn(do.call(rbind, lapply(list, diag)))),
            colnames
          )
        )
      } else {
        fn(do.call(rbind, lapply(list, diag)))
      }
    }
    
    # column names
    .colnames <- c("t", private$state.names)
    .obs.colnames <- c("t", private$obs.names)
    .covnames <- c("t",as.vector(outer(private$state.names, private$state.names, function(a, b) paste0(a,b))))
    .covnames.list <- paste("t = ", private$data$t, sep="")
    
    ##### priors
    filt$states$mean$prior = rbind_vectors(rep$xPrior, .colnames)
    filt$states$sd$prior = rbind_matrices_diag(rep$pPrior, .colnames, sqrt)
    filt$states$cov$prior = rep$pPrior
    names(filt$states$cov$prior) = .covnames.list
    
    ##### posteriors
    filt$states$mean$posterior = rbind_vectors(rep$xPost, .colnames)
    filt$states$sd$posterior = rbind_matrices_diag(rep$pPost, .colnames, sqrt)
    filt$states$cov$posterior = rep$pPost
    names(filt$states$cov$posterior) = .covnames.list
    
    ##### residuals
    obsMat = as.matrix(private$data[private$obs.names])
    ids <- 1:private$number.of.observations
    non.na.ids <- apply(obsMat, 1, function(x) ids[!is.na(x)], simplify = FALSE)
    length.non.na.ids <- lapply(non.na.ids, length)
    
    
    # pre-allocate
    filt$residuals <- lapply(1:3, function(x) set.colnames(cbind(private$data$t, matrix(NA, nrow=nrow(obsMat), ncol=ncol(obsMat))), .obs.colnames))
    names(filt$residuals) = c("residuals", "sd", "normalized")
    
    # Take care rows where all observations are present
    ids.full.obs <- length.non.na.ids == private$number.of.observations
    filt$residuals$residuals[ids.full.obs,-1] <- rbind_vectors(rep$Innovation[ids.full.obs], extra=FALSE)
    filt$residuals$sd[ids.full.obs,-1] <- rbind_matrices_diag(rep$InnovationCovariance[ids.full.obs], fn=sqrt, extra=FALSE)
    # Take care of all other rows with some missing observations
    for(i in seq_along(private$obs.names[-1])){
      ids <- which(length.non.na.ids == i)
      for(j in ids){
        filt$residuals$residuals[j, non.na.ids[[j]]+1] <- rep$Innovation[[j]]
        filt$residuals$sd[j, non.na.ids[[j]]+1] <- sqrt(diag(rep$InnovationCovariance[[j]]))
      }
    }
    filt$residuals$normalized = filt$residuals$residuals
    filt$residuals$normalized[,-1] <- filt$residuals$normalized[,-1]/filt$residuals$sd[,-1] 
    filt$residuals$cov <- rep$InnovationCovariance
    names(filt$residuals$cov) = .covnames.list
    
    ##### observations
    # call c++ function
    observations <- calculate_filtering_observations(private$filtration.raw,
                                                     private$rcpp_function_ptr,
                                                     as.matrix(private$data[private$input.names]),
                                                     private$pars,
                                                     private$number.of.observations
    )
    colnames(observations$mean$prior) <- .obs.colnames
    colnames(observations$mean$posterior) <- .obs.colnames
    colnames(observations$sd$prior) <- .obs.colnames
    colnames(observations$sd$posterior) <- .obs.colnames
    names(observations$cov$prior) = .covnames.list
    names(observations$cov$posterior) = .covnames.list
    filt$observations <- observations
    
  }
  
  # When we use for estimate we don't want to store in the filtration
  if(return.directly){
    return(filt)
  }
  
  # store
  private$filtration <- filt
  
  # return
  return(invisible(self))
}

