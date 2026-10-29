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
  if (private$algo.settings$estimate.initial) {
    stationary.solver <- construct.initial.state.estimators(self, private)
    get.cov <- TRUE
    out <- stationary.solver(parVec, inputMat[1,],get.cov=get.cov)
    stateVec <- out[,1]
    if(get.cov) {
      covMat <- out[,-1]
    }
  }


  laplace.methods <- get.inference.algos("laplace")
  output <- NULL
  ########################################################################
  # FILTERING
  ########################################################################
  if (proc == "filter") {

    if(!private$algo.settings$silent) message("Filtering...")

    if (private$algo.settings$method %in% laplace.methods)
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
        any_available_obs,
        non_na_ids,
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
        any_available_obs,
        non_na_ids,
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
        any_available_obs,
        non_na_ids,
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

    if (private$algo.settings$method %in% c("laplace","laplace.thygesen")) {
      laplace_smooth(self, private)
    }

    # Linear Kalman Filter
    if ( private$algo.settings$method == "lkf" ) {
      output <- lkf_smooth_rcpp(
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
      stop("UKF smoothing not yet implemented.")
      # output <- ukf_filter_rcpp(
      #   private$model$rcpp_function_ptr,
      #   obsMat,
      #   inputMat,
      #   parVec,
      #   covMat,
      #   stateVec,
      #   private$algo.settings$ode.stepsizes,
      #   private$algo.settings$ode.number.of.steps,
      #   numeric_is_not_na_obsMat,
      #   number_of_available_obs,
      #   private$algo.settings$ukf.hyperpars,
      #   private$algo.settings$ode.solver,
      #   private$algo.settings$first.order.input.hold
      # )
    }

    private$results$smooth.raw <- output

  }

  return(invisible(self))

}

#######################################################
# RETURN FILTERING
#######################################################
create_return_filter <- function(self, private, laplace.residuals){

  if(!private$algo.settings$silent) message("Returning results...")

  # Extract and pre-allocate
  filt <- list()
  rep <- private$results$filtration.raw

  # Helper functions
  split.state.matrix <- function(mat, names, time.column=private$data[["t"]]){
    # mat has construction [mean, cov]
    # Extract into mean, standard dev and covariance columns
    n <- length(names)
    diag.ids <- seq(1, n^2, by=n+1)
    mean <- mat[, seq_len(n), drop = FALSE]
    cov  <- mat[, n + seq_len(n^2), drop = FALSE]
    sd   <- sqrt(cov[, diag.ids, drop = FALSE])

    # Construct [time, mean, sd] matrix
    summary <- cbind(time.column, mean, sd)
    colnames(summary) <- c("t", names, paste0(names, ".sd"))

    # Construct [time, cov] matrix
    covnames <- apply(expand.grid(names,names),1,\(a) paste0(a[1],".",a[2],".cov"))
    cov <- cbind(time.column, cov)
    colnames(cov) <- c("t", covnames)

    # Return list
    x <- list(mean_and_sd = summary, cov = cov)
    lapply(x, as.data.frame)
  }

  split.residual.matrix <- function(mat, names, time.column=private$data[["t"]]){
    # mat has construction [mean, cov]
    # Extract into mean, standard dev and covariance columns
    n <- length(names)
    diag.ids <- seq(1, n^2, by=n+1)
    mean <- mat[, seq_len(n), drop = FALSE]
    cov  <- mat[, n + seq_len(n^2), drop = FALSE]
    sd   <- sqrt(cov[, diag.ids, drop = FALSE])


    # Construct [time, mean, normalized] matrix
    summary <- cbind(time.column, mean, sd)
    colnames(summary) <- c("t", names, paste0(names, ".sd"))

    # Construct [time, cov] matrix
    covnames <- apply(expand.grid(names,names),1,\(a) paste0(a[1],".",a[2],".cov"))
    cov <- cbind(time.column, cov)
    colnames(cov) <- c("t", covnames)

    marginal.standardized <- cbind(time.column, mean/sd)
    colnames(marginal.standardized) <- c("t", names)

    # Return list
    x <- list(
      residuals = summary,
      normalized = marginal.standardized,
      cov = cov
    )
    lapply(x, as.data.frame)
  }

  # States
  filt$states$prior <- split.state.matrix(rep$priorMat, private$names$states)
  filt$states$posterior <- split.state.matrix(rep$postMat, private$names$states)

  # Observations
  # We do not include observation noise by default i.e.
  # cov(h(x)) = H %*% P %*% t(H)
  # not
  # cov(y) = H %*% P %*% t(H) + V
  obs.noise <- FALSE
  obs.calculate.helper <- function(mat) from_states_to_observations(private$model$rcpp_function_ptr,
                                                                    mat,
                                                                    as.matrix(private$data[private$names$inputs]),
                                                                    private$algo.settings$argument.parameters,
                                                                    private$dims$observations,
                                                                    private$dims$states,
                                                                    obs.noise)
  obs.prior <- obs.calculate.helper(rep$priorMat)
  obs.post <- obs.calculate.helper(rep$postMat)
  filt$observations$prior <- split.state.matrix(obs.prior$obsMat, private$names$obs)
  filt$observations$posterior <- split.state.matrix(obs.post$obsMat, private$names$obs)

  # Residuals
  filt$residuals <- split.residual.matrix(rep$innoMat, private$names$obs)

  # Store
  private$results$filtration <- filt

  # return
  return(invisible(self))
}

create_return_smoother <- function(self, private, laplace.residuals){

  out <- list()
  method <- private$algo.settings$method
  state.names <- private$names$states
  obs.names <- private$names$obs
  n.obs <- private$dims$observations

  # Helper functions
  split.matrix <- function(mat, names, time.column=private$data[["t"]]){
    # mat has construction [mean, cov]
    # Extract into mean, standard dev and covariance columns
    n <- length(names)
    diag.ids <- seq(1, n^2, by=n+1)
    mean <- mat[, seq_len(n), drop = FALSE]
    cov  <- mat[, n + seq_len(n^2), drop = FALSE]
    sd   <- sqrt(cov[, diag.ids, drop = FALSE])

    # Construct [time, mean, sd] matrix
    summary <- cbind(time.column, mean, sd)
    colnames(summary) <- c("t", names, paste0(names, ".sd"))

    # Construct [time, cov] matrix
    covnames <- apply(expand.grid(names,names),1,\(a) paste0(a[1],".",a[2],".cov"))
    covnames[diag.ids] <- paste0(names,".var")
    cov <- cbind(time.column, cov)
    colnames(cov) <- c("t", covnames)

    # Return list
    x <- list(mean_and_sd = summary, cov = cov)
    lapply(x, as.data.frame)
  }

  if(method %in% get.inference.algos("laplace")) {

    smooth <- list()

    random.effects.ids <- private$algo.settings$ode.number.of.steps.cumsum + 1

    # States (Smoothed) -----------------------------------
    temp.states <- matrix(private$sdr$par.random, ncol=private$dims$states)[random.effects.ids, ]
    temp.sd <- matrix(sqrt(private$sdr$diag.cov.random), ncol=private$dims$states)[random.effects.ids, ]
    temp.df <- data.frame(private$data$t, temp.states, temp.sd)
    names(temp.df) <- c("t", state.names, paste0(state.names,".sd"))
    smooth$states$mean_and_sd <- temp.df

    # Observations
    obs.noise <- FALSE
    observations <- from_states_to_observations_no_state_cross_covariances(
      private$model$rcpp_function_ptr,
      cbind(temp.states, temp.sd^2),
      as.matrix(private$data[private$names$inputs]),
      private$algo.settings$argument.parameters,
      private$dims$observations,
      private$dims$states,
      obs.noise)

    smooth$observations <- split.matrix(observations$obsMat, obs.names)

    # Residuals -----------------------------------
    # compute one-step residuals
    if(laplace.residuals){
      message("Calculating one-step ahead residuals...")
      res <- RTMB::oneStepPredict(private$nll,
                                  observation.name="obsMat",
                                  method="oneStepGaussian",
                                  trace=TRUE)
      nr <- nrow(private$data)
      temp.res <- data.frame(t=private$data$t, matrix(res[["residual"]],ncol=private$dims$observations))
      temp.sd <- as.data.frame(matrix(res[["sd"]],ncol=private$dims$observations))
      names(temp.res) = c("t", private$names$obs)
      names(temp.sd) = private$names$obs
      smooth$residuals$residuals <- cbind(temp.res, temp.sd)
      smooth$residuals$normalized <- temp.res
      smooth$residuals$normalized[,-1] <- temp.res[,-1]/temp.sd

    }


    private$results$smooth <- smooth

  } else if (method %in% get.inference.algos("kalman")) {

    # Get smoothing results from algos
    rep <- private$results$smooth.raw

    # States
    out$states <- split.matrix(rep$smoothMat, private$names$states)

    # Observations
    # We do not include observation noise by default i.e.
    # cov(h(x)) = H %*% P %*% t(H)
    # not
    # cov(y) = H %*% P %*% t(H) + V
    obs.noise <- FALSE
    observations <- from_states_to_observations(private$model$rcpp_function_ptr,
                                                rep$smoothMat,
                                                as.matrix(private$data[private$names$inputs]),
                                                private$algo.settings$argument.parameters,
                                                private$dims$observations,
                                                private$dims$states,
                                                obs.noise)
    out$observations <- split.matrix(observations$obsMat, private$names$obs)

    #store
    private$results$smooth <- out

  }

  # return
  return(invisible(self))
}

#######################################################
# RETUNR PREDICTIONS
create_return_prediction <- function(reported.dispersion.type, return.k.ahead, self, private){

  if (!private$algo.settings$silent) message("Returning results...")

  t               <- private$data$t
  n               <- private$dims$states
  n.obs           <- private$dims$observations
  k.ahead         <- private$algo.settings$k.ahead
  state.names     <- private$names$states
  obs.names       <- private$names$obs
  last.pred.index <- private$algo.settings$last.pred.index
  diag.ids        <- seq(from=1, to=n^2, by=n+1)
  diag.ids.obs    <- seq(from=1, to=n.obs^2, by=n.obs+1)
  predmat         <- do.call(rbind, private$results$prediction.raw)
  pred.means      <- predmat[, seq_len(n), drop = FALSE]
  pred.covs       <- predmat[, n + seq_len(n^2), drop = FALSE]

  ##### TIME / INDICES #####
  i <- rep(seq_len(last.pred.index)-1, each = k.ahead+1)
  k <- rep(0:k.ahead, times = last.pred.index)
  j <- i + k
  time.mat <- cbind(i. = i, j. = j, t.i = t[i + 1], t.j = t[j + 1], k.ahead = k)
  keep.rows <- k %in% return.k.ahead

  # Format the dispersion part of a [mean | flattened cov] matrix
  format_prediction_matrix <- function(mat, names, type = c("covariance", "correlation")) {
    type <- match.arg(type)
    n <- length(names)
    diag.ids <- seq(1, n^2, by = n + 1)
    disp.type <- ifelse(type=="covariance", "cov", "cor")
    disp.names <- apply(expand.grid(names, names),1,\(a) paste0(a[1],".",a[2],".", disp.type))
    means      <- mat[, seq_len(n), drop = FALSE]
    cov       <- mat[, n + seq_len(n^2), drop = FALSE]

    # Means
    m <- means
    colnames(m) <- names

    # Standard deviations
    sd <- sqrt(cov[, diag.ids, drop = FALSE])
    colnames(sd) <- paste0(names, ".sd")

    # Set covariances
    disp <- cov
    colnames(disp) <- disp.names
    if (type == "correlation") {
      # vectorized cov2cor: cor_ij = cov_ij / (sd_i * sd_j), variances kept on the diagonal
      row.id <- rep(seq_len(n), times = n)
      col.id <- rep(seq_len(n), each = n)
      disp <- cov / (sd[, row.id, drop = FALSE] * sd[, col.id, drop = FALSE])
      disp[, diag.ids] <- cov[, diag.ids]
      colnames(disp) <- disp.names
    }
    colnames(disp)[diag.ids] <- paste0(names, ".var")

    # Create return list
    state.list <- list(
      mean_and_sd = cbind(time.mat, m, sd),
      dispersion = cbind(time.mat, disp)
    )
    names(state.list)[2] <- disp.type
    state.list <- lapply(state.list, function(x) x[keep.rows,])
    state.list <- lapply(state.list, as.data.frame)
  }

  ##### STATES #####
  state.list <- format_prediction_matrix(predmat, names=state.names, type=reported.dispersion.type)

  ##### OBS #####
  obs.predmat <- calculate_prediction_observations(predmat,
                                                   private$model$rcpp_function_ptr,
                                                   as.matrix(private$data[private$names$inputs]),
                                                   private$algo.settings$argument.parameters,
                                                   private$dims$states,
                                                   private$dims$observations,
                                                   private$algo.settings$last.pred.index,
                                                   private$algo.settings$k.ahead,
                                                   TRUE)
  obs.list <- format_prediction_matrix(obs.predmat, names=obs.names, type=reported.dispersion.type)

  #### DATA ####
  data.matrix <- cbind(time.mat, private$data[j+1, obs.names,drop=FALSE])


  # Store in private
  list.out <- list(states = state.list, observations=obs.list, data=data.matrix)
  class(list.out) = c(class(list.out), "ctsmTMB.pred")
  attr(list.out, "state.names") <- state.names
  attr(list.out, "obs.names") <- obs.names

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
    private$dims$states,
    private$algo.settings$k.ahead,
    n.sims
  )
  names(state.list) <- private$names$states
  for(i in seq_along(state.list)){
    names(state.list[[i]]) <- inner.names
  }

  build_simulation_timelists <- function(time_vector, n_forecasts, k_ahead) {
    ka <- 0:k_ahead
    lapply(seq_len(n_forecasts) - 1L, function(k) {
      cbind(i=k, j=k+ka, t.i=time_vector[k + 1L],t.j=time_vector[k + ka + 1L],k.ahead = ka)
    })
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
