#######################################################
# MAIN PREDICTION FUNCTION THAT CALL OTHERS
#######################################################
perform_prediction <- function(self, private, use.cpp){
  
  
  # Laplace cant produce predictions currently...
  if(private$method=="laplace"){
    stop("Predictions arent available for laplace yet.")
  }
  
  if(private$method=="laplace.thygesen"){
    stop("Predictions arent available for laplace.thygesen yet.")
  }
  
  # time prediction
  comptime <- system.time(
    {
      
      # Predict with Rcpp implementation
      if(use.cpp){
        
        if(!private$silent) message("Predicting with C++...")
        
        lkf_ekf_ukf_predict_rcpp(private$pars, self, private)
        
        # Predict with R implementation
      } else {
        
        if(!private$silent) message("Predicting with R...")
        
        lkf_ekf_ukf_predict_r(private$pars, self, private)
        
      }
      
    }, gcFirst = FALSE)
  
  private$timer_prediction <- comptime
  
  return(invisible(self))
}

lkf_ekf_ukf_predict_r <- function(pars, self, private){
  
  output <- switch(private$method,
                   lkf = lkf_predict_r(pars, self, private),
                   ekf = ekf_predict_r(pars, self, private),
                   ukf = ukf_predict_r(pars, self, private),
  )
  
  private$prediction <- output
  
  return(invisible(self))
}

lkf_ekf_ukf_predict_rcpp <- function(pars, self, private){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  output <- NULL
  # predict using c++ function
  if(private$method=="lkf"){
    output <- lkf_predict_rcpp(private$rcpp_function_ptr$f,
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
                               number_of_available_obs,
                               private$last.pred.index,
                               private$n.ahead)
  }
  if(private$method == "ekf"){
    output <- ekf_predict_rcpp(private$rcpp_function_ptr$f,
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
                               private$last.pred.index,
                               private$n.ahead,
                               private$ode.solver)
  }
  if(private$method == "ukf"){
    output <- ukf_predict_rcpp(private$rcpp_function_ptr$f,
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
                               private$last.pred.index,
                               private$n.ahead,
                               private$ode.solver)
  }
  
  ####### STORE PREDICTION #######
  private$prediction <- output[[1]]
  
  ####### RETURN #######
  return(invisible(self))
  
}

#######################################################
# MAIN RETURN PREDICTIONS FUNCTION THAT CALL OTHERS
#######################################################
create_return_prediction <- function(return.covariance, return.k.ahead, self, private){
  
  if(!private$silent) message("Returning results...")
  
  # Simlify variable names
  n = private$number.of.states
  n.ahead = private$n.ahead
  state.names = private$state.names
  last.pred.index = private$last.pred.index
  
  # Create return data.frame
  df.out = data.frame(matrix(nrow=last.pred.index*(n.ahead+1), ncol=5+n+n^2))
  disp_names = sprintf(rep("cor.%s.%s",n^2),rep(state.names, each=n),rep(state.names,n))
  disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  if(return.covariance){
    disp_names = sprintf(rep("cov.%s.%s",n^2),rep(state.names,each=n),rep(state.names,n))
    disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  }
  names(df.out) = c("i.","j.","t.i","t.j","k.ahead",state.names,disp_names)
  
  # Fill out data.frame
  ran = 0:(last.pred.index-1)
  df.out["i."] = rep(ran,each=n.ahead+1)
  df.out["j."] = df.out["i."] + rep(0:n.ahead,last.pred.index)
  df.out["t.i"] = rep(private$data$t[ran+1],each=n.ahead+1)
  df.out["t.j"] = private$data$t[df.out[,"i."]+1+rep(0:n.ahead,last.pred.index)]
  df.out["k.ahead"] = rep(0:n.ahead,last.pred.index)
  df.obs = df.out[c("i.","j.","t.i","t.j","k.ahead")]
  
  df.out[, c(state.names, disp_names)] <- do.call(rbind, private$prediction)
  if(!return.covariance){
    diag.ids <- seq(from=1,to=n^2,by=n+1)
    .seq <- seq(from=1,to=n^2,by=1)
    non.diag.ids <- .seq[!(.seq %in% diag.ids)]
    df.out[,disp_names[non.diag.ids]] <- t(apply(df.out, 1, function(x) as.vector(stats::cov2cor(matrix(tail(x, n^2), nrow=n)))))[,non.diag.ids]
  }
  
  ##### OBSERVATION PREDICTIONS #####
  # calculate observations at every time-step in predict
  inputs.df = private$data[df.out[,"j."]+1,private$input.names]
  
  named.pars.list = as.list(private$pars)
  names(named.pars.list) = private$parameter.names
  # create environment
  env.list = c(
    # states
    as.list(df.out[state.names]),
    # inputs
    as.list(inputs.df),
    # free and fixed parameters
    named.pars.list
  )
  
  # calculate observations
  obs.df.predict = as.data.frame(
    lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = env.list)})
  )
  names(obs.df.predict) = paste(private$obs.names)
  
  # add data observation to output data.frame 
  obs.df.data = private$data[df.out[,"j."]+1, private$obs.names, drop=F]
  names(obs.df.data) = paste(private$obs.names,".data",sep="")
  
  df.obs = cbind(df.obs, obs.df.predict, obs.df.data)
  
  # return only specific n.ahead
  df.out = df.out[df.out[,"k.ahead"] %in% return.k.ahead,]
  df.obs = df.obs[df.obs[,"k.ahead"] %in% return.k.ahead,]
  
  list.out = list(states = df.out, observations = df.obs)
  class(list.out) = c(class(list.out), "ctsmTMB.pred")
  
  private$prediction = list.out
  
  return(invisible(self))
  
}
