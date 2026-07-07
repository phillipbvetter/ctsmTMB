#######################################################
# MAIN PREDICTION FUNCTION THAT CALL OTHERS
#######################################################
perform_prediction <- function(self, private, use.cpp){


  # Laplace cant produce predictions currently...
  if(private$algo.settings$method=="laplace"){
    stop("Predictions arent available for laplace yet.")
  }

  if(private$algo.settings$method=="laplace.thygesen"){
    stop("Predictions arent available for laplace.thygesen yet.")
  }

  # time prediction
  comptime <- system.time(
    {

      # Predict with Rcpp implementation
      if(use.cpp){

        if(!private$silent) message("Predicting with C++...")

        lkf_ekf_ukf_predict_rcpp(private$algo.settings$argument.parameters, self, private)

        # Predict with R implementation
      } else {

        if(!private$silent) message("Predicting with R...")

        lkf_ekf_ukf_predict_r(private$algo.settings$argument.parameters, self, private)

      }

    }, gcFirst = FALSE)

  private$timers$prediction <- comptime

  return(invisible(self))
}

lkf_ekf_ukf_predict_r <- function(pars, self, private){

  output <- switch(private$algo.settings$method,
                   lkf = lkf_predict_r(pars, self, private),
                   ekf = ekf_predict_r(pars, self, private),
                   ukf = ukf_predict_r(pars, self, private),
  )

  # private$results$prediction <- output
  private$results$prediction.raw <- output

  return(invisible(self))
}

lkf_ekf_ukf_predict_rcpp <- function(pars, self, private){

  # check for null pointers
  check_if_rcpp_pointers_are_valid(self, private)

  # observation/input matrix
  obsMat = as.matrix(private$data[private$names$obs])
  inputMat = as.matrix(private$data[private$names$inputs])

  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)

  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)

  ids <- 1:private$dims$observations - 1 #minus 1 for 0 indexing
  non_na_ids <- apply(obsMat, 1, function(x) ids[!is.na(x)], simplify = FALSE)
  any_available_obs <- sapply(non_na_ids, function(x) !is.null(x))

  output <- NULL
  # predict using c++ function
  if(private$algo.settings$method=="lkf"){
    output <- lkf_predict_rcpp(private$model$rcpp_function_ptr,
                               obsMat,
                               inputMat,
                               pars,
                               private$initial.state$p0,
                               private$initial.state$x0,
                               private$ode.timestep.size,
                               private$ode.timesteps,
                               any_available_obs,
                               non_na_ids,
                               private$algo.settings$last.pred.index,
                               private$algo.settings$k.ahead,
                               private$algo.settings$first.order.input.hold)
  }
  if(private$algo.settings$method == "ekf"){
    output <- ekf_predict_rcpp(private$model$rcpp_function_ptr,
                               obsMat,
                               inputMat,
                               pars,
                               private$initial.state$p0,
                               private$initial.state$x0,
                               private$ode.timestep.size,
                               private$ode.timesteps,
                               any_available_obs,
                               non_na_ids,
                               private$algo.settings$ode.solver,
                               private$algo.settings$last.pred.index,
                               private$algo.settings$k.ahead,
                               private$algo.settings$first.order.input.hold)
  }
  if(private$algo.settings$method == "ukf"){
    output <- ukf_predict_rcpp(private$model$rcpp_function_ptr,
                               obsMat,
                               inputMat,
                               pars,
                               private$initial.state$p0,
                               private$initial.state$x0,
                               private$ode.timestep.size,
                               private$ode.timesteps,
                               numeric_is_not_na_obsMat,
                               number_of_available_obs,
                               private$algo.settings$ukf.hyperpars,
                               private$algo.settings$last.pred.index,
                               private$algo.settings$k.ahead,
                               private$algo.settings$ode.solver,
                               private$algo.settings$first.order.input.hold)
  }

  ####### STORE PREDICTION #######
  private$results$prediction.raw <- output[[1]]

  ####### RETURN #######
  return(invisible(self))

}

#######################################################
# MAIN RETURN PREDICTIONS FUNCTION THAT CALL OTHERS
#######################################################
create_return_prediction <- function(reported.dispersion.type, return.k.ahead, self, private){

  if(!private$silent) message("Returning results...")

  # Simlify variable names
  n               <- private$dims$states
  n.obs           <- private$dims$observations
  k.ahead         <- private$algo.settings$k.ahead
  state.names     <- private$names$states
  last.pred.index <- private$algo.settings$last.pred.index
  diag.ids        <- seq(from=1, to=n^2, by=n+1)
  diag.ids.obs    <- seq(from=1, to=n.obs^2, by=n.obs+1)
  rbinded.predmat <- do.call(rbind, private$results$prediction.raw)

  # time-related entries
  m.state = matrix(nrow=last.pred.index*(k.ahead+1), ncol=5+n)
  colnames(m.state) = c("i.","j.","t.i","t.j","k.ahead", private$names$states)
  ran = 0:(last.pred.index-1)
  m.state[,"i."] <- rep(ran, each=k.ahead+1)
  m.state[,"k.ahead"] <- rep(0:k.ahead, last.pred.index)
  m.state[,"j."] <- m.state[,"i."] + m.state[,"k.ahead"]
  m.state[,"t.i"] <- private$data$t[m.state[,"i."]+1]
  m.state[,"t.j"] <- private$data$t[m.state[,"j."]+1]
  m.obs <- m.state[,1:5] #for observations further below

  ##### STATE PREDICTIONS #####
  # predicted state means
  m.state[,state.names] = rbinded.predmat[,1:n, drop=FALSE]
  # predicted state dispersions
  if(reported.dispersion.type != "none"){
    if(reported.dispersion.type == "marginal"){
      m.disp <- rbinded.predmat[, n+diag.ids, drop=FALSE]
      colnames(m.disp) <- sprintf(rep("var.%s",n), state.names)
    } else {
      m.disp <- rbinded.predmat[,-(1:n), drop=FALSE]
      colnames(m.disp) <- sprintf(rep("cov.%s.%s",n^2), rep(state.names,each=n), rep(state.names,n))
      if(reported.dispersion.type == "correlation"){
        m.disp <- t(apply(m.disp, 1, function(x) as.vector(cov2cor(matrix(x, nrow=n)))))
        colnames(m.disp) <- sprintf(rep("cor.%s.%s",n^2), rep(state.names,each=n), rep(state.names,n))
      }
      colnames(m.disp)[diag.ids] <- sprintf(rep("var.%s",n), state.names)
    }
    m.state <- cbind(m.state, m.disp)
  }

  ##### OBSERVATION PREDICTIONS #####
  # calculate observations
  m.obs.pred <- calculate_prediction_observations(rbinded.predmat,
                                                  private$model$rcpp_function_ptr,
                                                  as.matrix(private$data[private$names$inputs]),
                                                  private$algo.settings$argument.parameters,
                                                  private$dims$states,
                                                  private$dims$observations,
                                                  private$algo.settings$last.pred.index,
                                                  private$algo.settings$k.ahead,
                                                  reported.dispersion.type != "none")
  # set name etc
  if(reported.dispersion.type != "none"){
    m.obs.pred.mean <- m.obs.pred[,1:n.obs,drop=FALSE]
    colnames(m.obs.pred.mean) <- private$names$obs
    if(reported.dispersion.type == "marginal"){
      m.obs.pred.disp <- m.obs.pred[,n.obs+diag.ids.obs, drop=FALSE]
      colnames(m.obs.pred.disp) <- sprintf(rep("var.%s",n.obs), private$names$obs)
    } else {
      m.obs.pred.disp <- m.obs.pred[,-c(1:n.obs), drop=FALSE]
      colnames(m.obs.pred.disp) <- sprintf(rep("cov.%s.%s",n.obs^2), rep(private$names$obs,each=n.obs), rep(private$names$obs,n.obs))
      if(reported.dispersion.type == "correlation"){
        m.obs.pred.disp <- t(apply(m.obs.pred.disp, 1, function(x) as.vector(cov2cor(matrix(x, nrow=n.obs)))))
        colnames(m.obs.pred.disp) <- sprintf(rep("cor.%s.%s",n.obs^2), rep(private$names$obs,each=n.obs), rep(private$names$obs,n.obs))
      }
      colnames(m.disp)[diag.ids] <- sprintf(rep("var.%s",n.obs), private$names$obs)
    }
    m.obs.pred <- cbind(m.obs.pred.mean, m.obs.pred.disp)
  }

  m.obs.data = as.matrix(private$data[m.state[,"j."]+1, private$names$obs, drop=F])
  colnames(m.obs.data) = paste(private$names$obs,".data",sep="")

  # cbind all obs-related mats
  m.obs <- cbind(m.obs, m.obs.pred, m.obs.data)

  # return only specific k.ahead
  keep.rows <- m.state[,"k.ahead"] %in% return.k.ahead
  m.state = m.state[keep.rows,]
  m.obs = m.obs[keep.rows,]

  list.out = list(states = m.state, observations = m.obs)
  class(list.out) = c(class(list.out), "ctsmTMB.pred")

  private$results$prediction = list.out

  return(invisible(self))

}
