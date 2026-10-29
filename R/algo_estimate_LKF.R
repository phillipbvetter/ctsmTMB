construct.matrix.exponentials <- function(funs,dims,x,p,u,dt,first.order.hold){

  ##### Sometimes necessary to avoid rtmb errors #######
  "[<-" <- RTMB::ADoverload("[<-")
  # "diag<-" <- RTMB::ADoverload("diag<-")
  # "c" <- RTMB::ADoverload("c")

  n.states <- dims$n.states
  n.inputs <- dims$n.inputs

  ####### Compute Matrix Exponentials for 1-Step Mean and Variance #######
  # dX = A*X*dt + B*U*dt + G*dB
  out <- lapply(funs[c("dfdx","dfdu","g")], function(f) do.call(f,list(x,p,u)))
  A <- out[[1]]; B <- out[[2]]; G <- out[[3]]

  # Compute drift matrices Ahat * X and Bhat * U
  # ------------------------------------
  if(first.order.hold){
    # First Order Hold: expm([A B 0; 0 0 I; 0 0 0]) -> [expm(A), B0, B1; 0 0 0; 0 0 0]
    Phi1 <- 0.0 * RTMB::diag(n.states+2*(n.inputs+1))
    Phi1[1:n.states,1:n.states] <- A
    Phi1[1:n.states,(n.states+1):(n.states+1+n.inputs)] <- B
    Phi1[(n.states+1):(n.states+1+n.inputs),(n.states+n.inputs+2):ncol(Phi1)] <- RTMB::diag(n.inputs+1)
    ePhi1 <- as.matrix(Matrix::expm(Phi1 * dt))
    # A, B0 and B1
    Ahat <- ePhi1[1:n.states, 1:n.states, drop=FALSE]
    Ahat_T <- t(Ahat)
    Bhat <- ePhi1[1:n.states, (n.states+1):(n.states+n.inputs+1), drop=FALSE]
    B1hat <- ePhi1[1:n.states, (n.states+n.inputs+2):ncol(Phi1), drop=FALSE]
  } else {
    # Zero Order Hold: expm([A B; 0 0]) -> [expm(A), B0; 0 0]
    # Mean calculations
    Phi1 <- 0.0 * RTMB::diag(n.states+n.inputs+1)
    Phi1[1:n.states,1:n.states] <- A
    Phi1[1:n.states,(n.states+1):ncol(Phi1)] <- B
    ePhi1 <- as.matrix(Matrix::expm(Phi1 * dt))
    #
    Ahat <- ePhi1[1:n.states, 1:n.states, drop=FALSE]
    Ahat_T <- t(Ahat)
    Bhat <- ePhi1[1:n.states, (n.states+1):ncol(ePhi1), drop=FALSE]
    B1hat <- NULL
  }

  # Compute variance matrix Vhat
  # ------------------------------------
  # expm([-A GG^T \\ 0 A^t]) = [Q11 Q12; Q21 Q22]
  # with variance V = Q22^{T} * Q12
  Phi2 <- rbind(cbind(-A, G %*% t(G)),cbind(0*A,t(A)))
  ePhi2 <- as.matrix(Matrix::expm(Phi2 * dt))
  Q22 <- ePhi2[(n.states+1):ncol(ePhi2), (n.states+1):ncol(ePhi2), drop=FALSE]
  Q12 <- ePhi2[1:n.states, (n.states+1):ncol(ePhi2), drop=FALSE]
  Vhat <- t(Q22) %*% Q12

  return(list(Ahat = Ahat, Vhat = Vhat, Bhat = Bhat,  B1hat = B1hat))
}

# integrator for time update
# ------------------------------------
lkf.construct.ode.integrator <- function(self, private){

  # Unpack dimensions
  n.states <- private$dims$states

  # # tape the ode integrator (3 times faster adfun compilation)
  arglist <- list(
    Ahat=diag(n.states),
    x=numeric(n.states),
    bu=numeric(n.states),
    Vhat=diag(n.states),
    P=diag(n.states)
  )

  taped.integrator <- RTMB::MakeTape(function(x) {
    # x = list(Ahat, stateVec, inputTerm, Vhat, covMat)
    cbind(
      # Ahat %*% stateVec + (B %*% inputVec + B1 %*% dinputVec)
      # Note that the last term above is pre-computed and passed as the vector x[[3]]
      x[[1]] %*% x[[2]] + x[[3]],
      # Ahat %*% covMat %*% t(Ahat) + Vhat
      x[[1]] %*% x[[5]] %*% t(x[[1]]) + x[[4]]
    )
  }, arglist)

  time.update <- function(x,P,bu,Ahat,Vhat){
    sol <- taped.integrator(list(Ahat, x, bu, Vhat, P))
    list(x=sol[,1], P=sol[,-1,drop=FALSE])
  }

  return(time.update)
}


# data update step
# ------------------------------------
lkf.construct.data.updater <- ekf.construct.data.updater


# setup function
# ------------------------------------
lkf.setup <- function(self, private){

  dims <- list(
    n.states = private$dims$states,
    n.pars   = private$dims$pars,
    n.inputs = private$dims$inputs
  )

  settings <- list(
    first.order.hold = private$algo.settings$first.order.input.hold,
    estimate.initial.state = private$algo.settings$estimate.initial
  )

  advanced.settings <- list(
    nllreport = private$algo.settings$advanced.settings$nllreport,
    force.ad.on.input = private$algo.settings$advanced.settings$force.ad.on.input
  )

  data <- list(
    x = private$algo.settings$initial.state$x0,
    P = private$algo.settings$initial.state$p0,
    U = as.matrix(private$data[private$names$inputs]),
    Y = as.matrix(private$data[private$names$obs])
  )
  data$dt <- diff(data$U[1:2, 1])
  data$U_aug <- cbind(1,data$U)
  data$dU_aug <- with(c(data,settings),
                      if(first.order.hold) rbind(diff(U_aug),0) / dt else NULL)

  algorithms <- list(
    time.update = lkf.construct.ode.integrator(self, private),
    data.update = lkf.construct.data.updater(self, private),
    stationary.solver = construct.initial.state.estimators(self, private)
  )

  # final check?
  with(data, if (diff(range(diff(U[,1]))) > 1e-10 * dt) stop("The LKF requires equidistant observation times."))

  # return
  return(named_list(dims, data, settings, advanced.settings, algorithms))
}

# likelihood constructor
# ------------------------------------
MakeADFun_LKF = function(self, private)
{
  # Settings
  configure_ad_tape("RTMB", self, private)
  setup <- lkf.setup(self, private)
  list2env(setup$data, envir=environment())
  list2env(setup$settings, envir=environment())
  list2env(setup$advanced.settings, envir=environment())
  list2env(setup$algorithms, envir=environment())
  funs <- private$model$state.space.funs

  # likelihood function --------------------------------
  lkf.nll = function(p){

    ####### Sometimes necessary to avoid rtmb errors #######
    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")
    "c" <- RTMB::ADoverload("c")

    # Convert parameter list to vector
    p <- do.call(c, p[1:private$dims$pars])

    # Set likelihood to zero
    nll <- 0

    # Cast initial state and inputs to AD
    x <- RTMB::AD(x)
    P <- RTMB::AD(P)
    if (force.ad.on.input) U <- RTMB::AD(U)

    # Construct matrix exponentials for exact linear solution
    mats <- construct.matrix.exponentials(funs, setup$dims,x,p,U[1,],dt,first.order.hold)

    # We now compute input term for all time indices
    BU <- mats$Bhat %*% t(U_aug) # Zero Order Hold
    if (first.order.hold) {
      BU <- BU + mats$B1hat %*% t(dU_aug)
    }
    BU <- t(BU)

    # Estimate initial conditions as stationary solution with first input
    # This estimate is, as with the user provided initial state, assumed to be
    # the prior mean and covariance x_{0 | -1} and P_{0 | -1}
    if (estimate.initial.state) {
      get.cov <- TRUE
      sol <- stationary.solver(p,U[1,],get.cov=get.cov)
      x <- sol[,1]
      if(get.cov) P <- sol[,-1, drop=FALSE]
    }

    # First data update
    post <- data.update(x,P,p,U[1,],Y[1,])
    x <- post$x
    P <- post$P
    nll <- nll + post$nll

    for(i in seq_len(nrow(Y)-1)){
      ###### time update #######
      sol <- time.update(x,P,BU[i,],mats$Ahat,mats$Vhat)
      x <- sol$x; P <- sol$P
      ######## data update ########
      post <- data.update(x,P,p,U[i+1,],Y[i+1,])
      x <- post$x; P <- post$P
      nll <- nll + post$nll
    }

    # return
    return(nll)
  }

  ################################################
  # Construct Neg. Log-Likelihood
  ################################################

  parameters <- lapply(private$model$parameters, function(x) x[["initial"]])
  nll = RTMB::MakeADFun(func = lkf.nll,
                        parameters = parameters,
                        map = lapply(private$model$fixed.pars, function(x) x$factor),
                        silent=TRUE)

  # save objective function
  private$nll = nll

  # return
  return(invisible(self))

}
