#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

perform_smoothing = function(self, private){
  
  if(!private$silent) message("Smoothing...")
  
  if(private$algo.settings$method=="lkf"){
    stop("Smoothing with the lkf method is not available. Use the 'filter' method instead.")
  }
  if(private$algo.settings$method=="ekf"){
    stop("Smoothing with the ekf method is not available. Use the 'filter' method instead.")
  }
  if(private$algo.settings$method=="ukf"){
    stop("Smoothing with the ukf method is not available. Use the 'filter' method instead.")
  }
  
  # Kalman smoothers
  private$smooth <- switch(private$algo.settings$method,
                           # These must be changed to e.g. lkf_smoother_r when available
                           # lkf = lkf_smoother_r(private$model$argument.parameters, self, private),
                           # ekf = ekf_smoother_r(private$model$argument.parameters, self, private),
                           # ukf = ukf_smoother_r(private$model$argument.parameters, self, private),
                           laplace = laplace_smooth(self, private),
                           laplace.thygesen = laplace_smooth(self, private)
  )
  
  return(invisible(self))
}

laplace_smooth <- function(self, private){
  
  # Integrate out random effect states (perform laplace approx)
  par.type.free <- self$getParameters()[,"type"] == "free"
  free.pars <- private$model$argument.parameters[par.type.free]
  private$nll$fn(free.pars)
  
  # Report result
  sdr <- RTMB::sdreport(private$nll,
                        ignore.parm.uncertainty=T,
                        skip.delta.method=T)
  private$sdr <- sdr
}

create_smooth_results <- function(self, private, laplace.residuals){
  
  smooth <- list()
  
  if(any(private$algo.settings$method == c("laplace","laplace.thygesen"))){
    
    random.effects.ids <- private$ode.timesteps.cumsum + 1
  
    # States (Smoothed) -----------------------------------
    temp.states <- matrix(private$sdr$par.random, ncol=private$dimensions$states)[random.effects.ids, ]
    temp.sd <- matrix(sqrt(private$sdr$diag.cov.random), ncol=private$dimensions$states)[random.effects.ids, ]
    temp.states <- cbind(private$data$t, temp.states)
    temp.sd <- cbind(private$data$t, temp.sd)
    
    smooth$states$mean$smoothed = as.data.frame(temp.states)
    smooth$states$sd$smoothed = as.data.frame(temp.sd)
    colnames(smooth$states$sd$smoothed) = c("t",private$names$states)
    colnames(smooth$states$mean$smoothed) = c("t",private$names$states)
    
    # Residuals -----------------------------------
    rowNAs = as.matrix(!is.na(do.call(cbind, private$data[private$names$obs]))[-1,])
    sumrowNAs = rowSums(rowNAs)
    
    # compute one-step residuals
    if(laplace.residuals){
      message("Calculating one-step ahead residuals...")
      res <- RTMB::oneStepPredict(private$nll,
                                  observation.name="obsMat",
                                  method="oneStepGaussian",
                                  trace=TRUE)
      nr <- nrow(private$data)
      temp.res <- data.frame(private$data$t, matrix(res[["residual"]],ncol=private$dimensions$observations))
      temp.sd <- data.frame(private$data$t, matrix(res[["sd"]],ncol=private$dimensions$observations))
      names(temp.res) = c("t", private$names$obs)
      names(temp.sd) = c("t", private$names$obs)
      
      smooth$residuals$residuals <- temp.res
      smooth$residuals$sd <- temp.sd
      smooth$residuals$normalized <- temp.res
      smooth$residuals$normalized[,-1] <- temp.res[,-1]/temp.sd[,-1]
    }
    
  }
  
  private$smooth <- smooth
  
}
