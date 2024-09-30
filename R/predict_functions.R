#' Prediction function that calls the underlying Rcpp prediction function
#' 
#' @param private model object private fields
#' @param self model object
rcpp_prediction = function(self, private){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # predict using c++ function
  mylist <- execute_ekf_prediction(private$Rcppfunction_f,
                                   private$Rcppfunction_g,
                                   private$Rcppfunction_dfdx,
                                   private$Rcppfunction_h,
                                   private$Rcppfunction_dhdx,
                                   private$Rcppfunction_hvar,
                                   obsMat,
                                   inputMat,
                                   private$pars,
                                   private$pred.initial.state$p0,
                                   private$pred.initial.state$x0,
                                   private$ode.timestep.size,
                                   private$ode.timesteps,
                                   numeric_is_not_na_obsMat,
                                   number_of_available_obs,
                                   private$number.of.states,
                                   private$number.of.observations,
                                   private$last.pred.index,
                                   private$n.ahead,
                                   private$ode.solver)
  
  
  private$prediction = mylist
  
  return(invisible(NULL))
}


#' Generates a user-friendly data.frame of prediction results from private$prediction
#' 
#' @param private model object private fields
#' @param self model object
create_return_prediction = function(return.covariance, return.k.ahead, self, private){
  
  # Simlify variable names
  n = private$number.of.states
  n.ahead = private$n.ahead
  state.names = private$state.names
  last.pred.index = private$last.pred.index
  
  # Create return data.frame
  df.out = data.frame(matrix(nrow=last.pred.index*(n.ahead+1), ncol=5+n+n^2))
  disp_names = sprintf(rep("cor.%s.%s",n^2),rep(state.names, each=n),rep(state.names,n))
  disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  var_bool = !stringr::str_detect(disp_names,"cor")
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
  
  ##### STATES PREDICTIONS ######
  df.out[,state.names] = do.call(rbind,lapply(private$prediction$Xpred, function(cur.list) do.call(rbind, cur.list)))
  if(return.covariance){ #covariance
    df.out[,disp_names] = do.call(rbind, lapply(private$prediction$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(x))) } ))
  } else { #correlation
    df.out[,disp_names] = do.call(rbind, lapply(private$prediction$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(cov2cor(x)))) } )) 
    diag.ids = seq(from=1,to=n^2,by=n+1)
    df.out[,disp_names[diag.ids]] = do.call(rbind, lapply(private$prediction$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(diag(x)))) } ))
  }
  
  ##### OBSERVATION PREDICTIONS #####
  # calculate observations at every time-step in predict
  inputs.df = private$data[df.out[,"j."]+1,private$input.names]
  
  named.pars.list = as.list(private$pars)
  names(named.pars.list) = names(private$free.pars)
  # create environment
  env.list = c(
    # states
    as.list(df.out[state.names]),
    # inputs
    as.list(inputs.df),
    # free parameters
    named.pars.list,
    # fixed parameters
    lapply(private$fixed.pars, function(x) x$initial)
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
  # class(list.out) = c(class(list.out), "ctsmTMB.pred")
  
  private$prediction = list.out
  
  return(list.out)
}
