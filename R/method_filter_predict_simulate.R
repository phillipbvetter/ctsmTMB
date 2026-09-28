filter_predict_simulate_smooth <- function(self, private, proc, n.sims=NULL) {

  ########################################################################
  # Pre-computations
  ########################################################################

  parVec <- private$algo.settings$argument.parameters

  # Note:
  # The UKF method still uses the old matrices 'numeric_is_not_na_obsMat' and 'number_of_available_obs'

  # Convert from data.frames to matrices
  obsMat <- as.matrix(private$data[private$names$obs])
  inputMat <- as.matrix(private$data[private$names$inputs])

  # Find the missing observations
  not_na <- !is.na(obsMat)

  # Convert logical entries to numerical without dropping dimensions
  numeric_is_not_na_obsMat <- 1 * not_na

  # For each time-step find dimension of obs vector (how many available observations there are)
  number_of_available_obs <- rowSums(not_na)

  # Also find the index of these observations - we will pick out those indices in the algorithm
  # non_na_ids <- lapply(seq_len(nrow(obsMat)), function(i) which(not_na[i, ]) - 1L)
  # The above line is the understandable code - this one below is just 10 times faster...
  non_na_ids <- unname(split(
    col(obsMat)[not_na] - 1L,
    factor(row(obsMat)[not_na], levels = seq_len(nrow(obsMat)))
  ))

  # We need to know if we can skip the update step if all observations are missing
  any_available_obs <- number_of_available_obs > 0

  # Calculating initial state
  stateVec <- private$algo.settings$initial.state$x0
  covMat <- private$algo.settings$initial.state$p0
  if (private$algo.settings$estimate.initial)
    stateVec <- compute_initial_state_estimate(parVec, inputMat[1,], self, private)


  output <- NULL
  ########################################################################
  # FILTERING
  ########################################################################
  if (proc == "filter") {

    if(!private$algo.settings$silent) message("Filtering...")

    if (private$algo.settings$method %in% c("laplace","laplace.thygesen"))
      stop("The filter method is not available for the Laplace methods")

    # Linear Kalman Filter
    if ( private$algo.settings$method == "lkf" )
      output <- lkf_filter_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        any_available_obs,
        non_na_ids,
        private$algo.settings$first.order.input.hold
      )

    # Extended Kalman Filter
    if ( private$algo.settings$method == "ekf" )
      output <- ekf_filter_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        any_available_obs,
        non_na_ids,
        private$algo.settings$ode.solver,
        private$algo.settings$first.order.input.hold
      )

    # Sarkka Unscented Kalman Filter
    if ( private$algo.settings$method == "ukf" )
      output <- ukf_filter_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        numeric_is_not_na_obsMat,
        number_of_available_obs,
        private$algo.settings$ukf.hyperpars,
        private$algo.settings$ode.solver,
        private$algo.settings$first.order.input.hold
      )

    private$results$filtration.raw <- output

  }


  ########################################################################
  # PREDICTIONS
  ########################################################################
  if (proc == "predict") {

    if(!private$algo.settings$silent) message("Predicting...")

    if (private$algo.settings$method %in% c("laplace","laplace.thygesen"))
      stop("The predict method is not available for the Laplace methods")

    if ( private$algo.settings$method=="lkf" )
      output <- lkf_predict_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        any_available_obs,
        non_na_ids,
        private$algo.settings$last.pred.index,
        private$algo.settings$k.ahead,
        private$algo.settings$first.order.input.hold
      )

    if(private$algo.settings$method == "ekf")
      output <- ekf_predict_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        any_available_obs,
        non_na_ids,
        private$algo.settings$ode.solver,
        private$algo.settings$last.pred.index,
        private$algo.settings$k.ahead,
        private$algo.settings$first.order.input.hold
      )

    if (private$algo.settings$method == "ukf")
      output <- ukf_predict_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        numeric_is_not_na_obsMat,
        number_of_available_obs,
        private$algo.settings$ukf.hyperpars,
        private$algo.settings$last.pred.index,
        private$algo.settings$k.ahead,
        private$algo.settings$ode.solver,
        private$algo.settings$first.order.input.hold
      )

    private$results$prediction.raw <- output[[1]]


  }

  ########################################################################
  # SIMULATIONS
  ########################################################################

  if (proc == "simulate") {

    if(!private$algo.settings$silent) message("Simulating...")

    if (private$algo.settings$method %in% c("laplace","laplace.thygesen"))
      stop("The simulate method is not available for the Laplace methods")

    if (private$algo.settings$method=="lkf")
      output <- lkf_simulate_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        private$algo.settings$sim.stepsizes,
        private$algo.settings$sim.number.of.steps,
        any_available_obs,
        non_na_ids,
        private$dims$diffusions,
        private$algo.settings$last.pred.index,
        private$algo.settings$k.ahead,
        n.sims,
        private$algo.settings$seed$state.seed,
        private$algo.settings$first.order.input.hold
      )

    if (private$algo.settings$method=="ekf")
      output <- ekf_simulate_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        private$algo.settings$sim.stepsizes,
        private$algo.settings$sim.number.of.steps,
        any_available_obs,
        non_na_ids,
        private$algo.settings$ode.solver,
        private$algo.settings$last.pred.index,
        private$algo.settings$k.ahead,
        private$dims$diffusions,
        n.sims,
        private$algo.settings$seed$state.seed,
        private$algo.settings$first.order.input.hold
      )

    if (private$algo.settings$method == "ukf")
      output <- ukf_simulate_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        private$algo.settings$sim.stepsizes,
        private$algo.settings$sim.number.of.steps,
        numeric_is_not_na_obsMat,
        number_of_available_obs,
        private$algo.settings$ukf.hyperpars,
        private$dims$diffusions,
        private$algo.settings$last.pred.index,
        private$algo.settings$k.ahead,
        private$algo.settings$ode.solver,
        n.sims,
        private$algo.settings$seed$state.seed,
        private$algo.settings$first.order.input.hold
      )

    private$results$simulation.raw <- output

  }


  if( proc == "smooth") {

    if(!private$algo.settings$silent) message("Smoothing...")

    if (private$algo.settings$method %in% c("laplace","laplace.thygesen"))
      perform_smoothing(self, private)

    # Linear Kalman Filter
    if ( private$algo.settings$method == "lkf" ) {
      stop("LKF smoothing not implemented")
      output <- lkf_filter_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        any_available_obs,
        non_na_ids,
        private$algo.settings$first.order.input.hold
      )
    }

    # Extended Kalman Filter
    if ( private$algo.settings$method == "ekf" ) {
      output <- ekf_smooth_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        any_available_obs,
        non_na_ids,
        private$algo.settings$ode.solver,
        private$algo.settings$first.order.input.hold
      )
    }

    # Sarkka Unscented Kalman Filter
    if ( private$algo.settings$method == "ukf" ) {
      stop("LKF smoothing not implemented")
      output <- ukf_filter_rcpp(
        private$model$rcpp_function_ptr,
        obsMat,
        inputMat,
        parVec,
        covMat,
        stateVec,
        private$algo.settings$ode.stepsizes,
        private$algo.settings$ode.number.of.steps,
        numeric_is_not_na_obsMat,
        number_of_available_obs,
        private$algo.settings$ukf.hyperpars,
        private$algo.settings$ode.solver,
        private$algo.settings$first.order.input.hold
      )
    }

    private$results$smooth.raw <- output

  }

  return(invisible(self))

}

#######################################################
# RETURN FILTERING
#######################################################
create_filter_results <- function(self, private, laplace.residuals){

  if(!private$algo.settings$silent) message("Returning results...")

  rtmb.kalman.methods <- c("lkf","ekf","ukf")
  all.kalman.methods <- c(rtmb.kalman.methods, paste0(rtmb.kalman.methods,".cpp"))

  filt <- list()

  if(private$algo.settings$method %in% all.kalman.methods){

    rep <- private$results$filtration.raw

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
    .colnames <- c("t", private$names$states)
    .obs.colnames <- c("t", private$names$obs)
    .covnames <- c("t",as.vector(outer(private$names$states, private$names$states, function(a, b) paste0(a,b))))
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
    obsMat = as.matrix(private$data[private$names$obs])
    not.na <- !is.na(obsMat)
    non.na.ids <- unname(split(
      col(obsMat)[not.na],
      factor(row(obsMat)[not.na], levels = seq_len(nrow(obsMat)))
    ))
    length.non.na.ids <- lapply(non.na.ids, length)

    # pre-allocate
    filt$residuals <- lapply(1:3, function(x) set.colnames(cbind(private$data$t, matrix(NA, nrow=nrow(obsMat), ncol=ncol(obsMat))), .obs.colnames))
    names(filt$residuals) = c("residuals", "sd", "normalized")

    # Take care of rows where all observations are present
    ids.full.obs <- unlist(length.non.na.ids) == private$dims$observations
    if(any(ids.full.obs)){
      filt$residuals$residuals[ids.full.obs,-1] <- rbind_vectors(rep$Innovation[ids.full.obs], extra=FALSE)
      filt$residuals$sd[ids.full.obs,-1] <- rbind_matrices_diag(rep$InnovationCovariance[ids.full.obs], fn=sqrt, extra=FALSE)
    }
    # Take care of all other rows with some missing observations
    for(i in seq_along(private$names$obs[-1])){
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
    # Computed in the C++ backend:
    observations <- calculate_filtering_observations(private$results$filtration.raw,
                                                     private$model$rcpp_function_ptr,
                                                     as.matrix(private$data[private$names$inputs]),
                                                     private$algo.settings$argument.parameters,
                                                     private$dims$observations
    )
    colnames(observations$mean$prior) <- .obs.colnames
    colnames(observations$mean$posterior) <- .obs.colnames
    colnames(observations$sd$prior) <- .obs.colnames
    colnames(observations$sd$posterior) <- .obs.colnames
    names(observations$cov$prior) = .covnames.list
    names(observations$cov$posterior) = .covnames.list
    filt$observations <- observations
  }

  # store
  private$results$filtration <- filt

  # return
  return(invisible(self))
}

#######################################################
# RETUNR PREDICTIONS
create_return_prediction <- function(reported.dispersion.type, return.k.ahead, self, private){

  if (!private$algo.settings$silent) message("Returning results...")

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

#######################################################
# RETURN SIMULATIONS
create_return_simulation <- function(return.k.ahead, n.sims, self, private){

  if(!private$algo.settings$silent) message("Returning results...")

  # create names for inner list
  inner.names <- paste0("i", 0:(private$algo.settings$last.pred.index-1))

  # Build returnlist for states
  state.list <- build_simulation_returnlist(
    private$results$simulation.raw,
    private$data$t,
    private$dims$states,
    private$algo.settings$k.ahead,
    n.sims
  )
  names(state.list) <- private$names$states
  for(i in seq_along(state.list)){
    names(state.list[[i]]) <- inner.names
  }

  # create index/time list
  time.list <- build_simulation_timelists(
    private$data$t,
    private$algo.settings$last.pred.index,
    private$algo.settings$k.ahead
  )
  names(time.list) <- inner.names

  # Build returnlist for observations
  # First we must calculate the simulated observation trajectories
  simulation.raw.obs <- calculate_simulation_observations(
    private$results$simulation.raw,
    private$model$rcpp_function_ptr,
    t(as.matrix(private$data[private$names$inputs])),
    private$algo.settings$argument.parameters,
    private$dims$states,
    private$dims$observations,
    private$algo.settings$k.ahead,
    n.sims,
    private$algo.settings$seed$obs.seed
  )

  # # Now we can build the returnlist
  obs.list <- build_simulation_returnlist(
    simulation.raw.obs,
    private$data$t,
    private$dims$observations,
    private$algo.settings$k.ahead,
    n.sims
  )
  names(obs.list) <- private$names$obs
  for(i in seq_along(obs.list)){
    names(obs.list[[i]]) <- inner.names
  }

  private$results$simulation <- list(states=state.list, observations=obs.list, times=time.list)

  return(invisible(self))
}
