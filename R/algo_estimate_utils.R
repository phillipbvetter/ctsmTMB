# loss function
# ------------------------------------
construct.loss.function <- function(dim = c("scalar","joint"), type = c("quadratic", "huber", "tukey"), c = 3) {
  type <- match.arg(type)
  dim <- match.arg(dim)

  # quadratic
  if (type == "quadratic") {
    loss <- switch(
      dim,
      scalar = function(e,R) -RTMB::dnorm(e, sd=sqrt(R), log=TRUE),
      joint = function(e, R) -RTMB::dmvnorm(e, Sigma = R, log = TRUE)
    )
  } else if (type %in% c("huber","tukey")) {
    # tukey and huber
    blend <- function(r2) 1 / (1 + exp(-5 * (sqrt(r2) - c)))
    rho <- switch(
      type,
      huber = function(r2) { s <- blend(r2); (1 - s) * r2 + s * c * (2 * sqrt(r2) - c) },
      tukey = function(r2) { s <- blend(r2); (1 - s) * r2 + s * c^2 }
    )
    loss <- switch(
      dim,
      scalar = function(e,R) 0.5 * (log(R) + length(e) * log(2*pi) + rho(e^2/R)),
      joint = function(e,R) 0.5 * (logdet(R) + length(e) * log(2*pi) + rho(sum(e*RTMB::solve(R,e))))
    )
  }

  return(loss)
}

# stationary solution solver
# ------------------------------------
construct.initial.state.estimators <- function(self, private){

  n.states <- private$dims$states
  n.pars   <- private$dims$pars
  n.inputs <- private$dims$inputs
  funs <- private$model$state.space.funs
  settings <- private$algo.settings

  # STATIONARY SOLVER FOR THE MEAN
  # --------------------------------
  arglist <- list(rep(1,n.states),rep(1,n.pars),rep(1,n.inputs))
  f <- RTMB::MakeTape(function(z) sum(do.call(funs$f,z)^2), arglist)
  stationary.state <- f$newton(seq_len(n.states), grad_tol=1e-20)

  # STATIONARY SOLVER FOR COVARIANCE
  # # --------------------------------
  stationary.covariance <- function(A, G) {
    P <- solve_lyapunov(c(A, G %*% t(G)))
    dim(P) <- dim(A)
    P
  }

  # Return solver
  # get.mean and get.cov are controlled by the entries (state mean, state covariance)
  # of the algo.settings$estimate.initial field set via set_initial_state_estimation.
  # When the state mean is not estimated (get.mean=FALSE) the covariance is
  # computed around the user-provided initial state x0.
  stationary.solve <- function(p,u,
                               get.mean = settings$estimate.initial[1],
                               get.cov = settings$estimate.initial[2],
                               x0 = NULL){
    z <- list(p,u)
    x <- if (get.mean) stationary.state(z) else x0
    if(get.cov){
      z <- list(x,p,u)
      A <- do.call(funs$dfdx, z)
      G <- do.call(funs$g, z)
      P <- stationary.covariance(A,G)
      return(cbind(x,P))
    }
    cbind(x)
  }

  return(stationary.solve)
}

create_model_funs <- function(self, private){

  # We need the parent to be ctsmTMB to find custom functions like
  # e.g. invlogit
  fun.env <- new.env(parent = asNamespace("ctsmTMB"))
  as_function <- function(body) {
    eval(parse(text = sprintf(
      "function(x,p,u) {%s; ans}", body)),
      envir = fun.env)
  }

  fun.strs <- create.state.space.string.bodies(self, private)
  funs <- list(
    f    = as_function(fun.strs$f),    # drift
    dfdx = as_function(fun.strs$dfdx), # drift jacobian (w.r.t statea)
    g    = as_function(fun.strs$g),    # diffusion
    h    = as_function(fun.strs$h),    # observation
    dhdx = as_function(fun.strs$dhdx), # observation jacobian (w.r.t states)
    hvar = as_function(fun.strs$hvar),  # observation variances (diagonal)
    hvar_matrix = as_function(fun.strs$hvar_matrix), #observation variances (matrix)
    dfdu = as_function(fun.strs$dfdu) # drift jacobian w.r.t inputs
  )
  private$model$state.space.funs <- funs
}

create_rcpp_model_funs <- function(self, private){

  fun.args <-"const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec"

  as.vec.fun.string <- function(name,body){
    head <- "Eigen::VectorXd"
    paste(sprintf("%s %s(%s){%s;return(ans);}",head,name,fun.args,body), collapse="")
  }
  as.mat.fun.string <- function(name,body){
    head <- "Eigen::MatrixXd"
    paste(sprintf("%s %s(%s){%s;return(ans);}",head,name,fun.args,body), collapse="")
  }
  as.array.fun.string <- function(name,body){
    head <- "Eigen::ArrayXd"
    paste(sprintf("%s %s(%s){%s;return(ans);}",head,name,fun.args,body), collapse="")
  }

  fun.strs <- create.rcpp.state.space.string.bodies(self, private)
  fun.list <- list(
    # vectors
    f_const = as.vec.fun.string("f_const", fun.strs$f),
    h_const = as.vec.fun.string("h_const", fun.strs$h),
    # matrices
    dfdx_const = as.mat.fun.string("dfdx_const", fun.strs$dfdx),
    g_const = as.mat.fun.string("g_const", fun.strs$g),
    dhdx_const = as.mat.fun.string("dhdx_const", fun.strs$dhdx),
    hvar_const = as.mat.fun.string("hvar_const", fun.strs$hvar),
    dfdu_const = as.mat.fun.string("dfdu_const", fun.strs$dfdu),
    # arrays 
    hvar_array_const = as.array.fun.string("hvar_array_const", fun.strs$hvar_array)
  )
  private$model$rcpp.function.strings <- fun.list

  return(invisible(self))
}

create.rcpp.state.space.string.bodies = function(self, private){

  n.states <- private$dims$states
  n.obs <- private$dims$observations
  n.pars <- private$dims$pars
  n.inputs <- private$dims$inputs
  n.diffusions <- private$dims$diffusions

  # Create named list of replacement names i.e. the state name x1 becomes stateVec(0)
  create.substitution.list <- function(names, arg.vector.name){
    fun <- function(i) parse(text=sprintf("%s(%i)", arg.vector.name, i-1))[[1]]
    calls <- lapply(seq_along(names), fun)
    setNames(calls, names)
  }
  substition.list <- c(
    create.substitution.list(private$names$obs, "obsVec"),
    create.substitution.list(private$names$parameters, "parVec"),
    create.substitution.list(private$names$states, "stateVec"),
    create.substitution.list(private$names$inputs, "inputVec")
  )

  # helpers
  ##################################################

  # substitute variable names with vector indexing and deparse
  to_term <- function(expr) deparse1(do.call(substitute, list(hat2pow(expr), substition.list)))

  vec_body <- function(terms) {
    nz <- which(terms != "0")
    paste(c(
      sprintf("Eigen::VectorXd ans = Eigen::VectorXd::Zero(%i)", length(terms)),
      sprintf("ans(%i) = %s", nz - 1, terms[nz])
    ), collapse = "; ")
  }

  vec_as_array_body <- function(terms) {
    nz <- which(terms != "0")
    paste(c(
      sprintf("Eigen::ArrayXd ans = Eigen::ArrayXd::Zero(%i)", length(terms)),
      sprintf("ans(%i) = %s", nz - 1, terms[nz])
    ), collapse = "; ")
  }

  mat_body <- function(terms) {
    nz <- which(terms != "0", arr.ind = TRUE)
    paste(c(
      sprintf("Eigen::MatrixXd ans = Eigen::MatrixXd::Zero(%i, %i)", nrow(terms), ncol(terms)),
      sprintf("ans(%i, %i) = %s", nz[,1] - 1, nz[,2] - 1, terms[nz])
    ), collapse = "; ")
  }

  # Functions
  ##################################################

  # drift
  f.terms <- vapply(seq_len(n.states), function(i) to_term(private$model$diff.terms[[i]]$dt), "")
  private$model$rcpp.state.space.fun.strings$f <- vec_body(f.terms)


  # drift jacobian
  dfdx.terms <- matrix("0", n.states, n.states)
  for(i in seq_len(n.states)){
    for(j in seq_len(n.states)){
      dfdx.terms[i,j] <- to_term(private$model$diff.terms.drift[[i]][[j]])
    }
  }
  private$model$rcpp.state.space.fun.strings$dfdx <- mat_body(dfdx.terms)

  # diffusion
  g.terms <- matrix("0", n.states, n.diffusions)
  for(i in seq_len(n.states)){
    for(j in seq_len(n.diffusions)){
      g.terms[i,j] <- to_term(private$model$diff.terms[[i]][[j+1]])
    }
  }
  private$model$rcpp.state.space.fun.strings$g <- mat_body(g.terms)

  # observation
  h.terms <- vapply(seq_len(n.obs), function(i) to_term(private$model$obs.eqs.trans[[i]]$rhs), "")
  private$model$rcpp.state.space.fun.strings$h <- vec_body(h.terms)

  # observation jacobian
  dhdx.terms <- matrix("0", n.obs, n.states)
  for(i in seq_len(n.obs)){
    for(j in seq_len(n.states)){
      dhdx.terms[i,j] <- to_term(private$model$diff.terms.obs[[i]][[j]])
    }
  }
  private$model$rcpp.state.space.fun.strings$dhdx <- mat_body(dhdx.terms)

  # observation variance
  hvar.terms <- vapply(seq_len(n.obs), function(i) to_term(private$model$obs.var.trans[[i]]$rhs), "")
  hvar.mat.terms <- matrix("0", n.obs, n.obs)
  diag(hvar.mat.terms) <- hvar.terms
  private$model$rcpp.state.space.fun.strings$hvar <- mat_body(hvar.mat.terms)
  private$model$rcpp.state.space.fun.strings$hvar_array <- vec_as_array_body(hvar.terms)

  # input matrix (Linear Kalman Filter)
  # derivative w.r.t inputs (for linear kalman filter)
  # B in dX = (A*X + B*U) dt + G dW with U = [1, u1(t), u2(t), ...]
  # first column: constant terms (all inputs and states set to zero)
  zero.list <- as.list(numeric(n.inputs + n.states))
  names(zero.list) <- c(private$names$inputs, private$names$states)
  constant.terms <- try_with_warning_recovery(
    unname(sapply(
      sapply(private$model$sys.eqs.trans, function(x) Deriv::Simplify(do.call(substitute, list(x$diff.dt, zero.list)))),
      function(x) deparse1(do.call(substitute, list(x, subsList)))
    ))
  )
  if(inherits(constant.terms, "try-error")){
    constant.terms <- rep("0", n.states)
  }

  # remaining columns: derivatives w.r.t. each input
  dfdu.terms <- matrix("0", n.states, n.inputs + 1)
  dfdu.terms[,1] <- constant.terms
  for(i in seq_len(n.states)){
    for(j in seq_len(n.inputs)){
      tempterm <- ctsmTMB_Deriv(f=private$model$sys.eqs.trans[[i]]$diff.dt, x=private$names$inputs[j])
      dfdu.terms[i,j+1] <- to_term(tempterm)
    }
  }
  private$model$rcpp.state.space.fun.strings$dfdu <- mat_body(dfdu.terms)

  # return
  return(private$model$rcpp.state.space.fun.strings)
}



create.state.space.string.bodies = function(self, private){

  n.states <- private$dims$states
  n.obs <- private$dims$observations
  n.pars <- private$dims$pars
  n.inputs <- private$dims$inputs
  n.diffusions <- private$dims$diffusions

  # Create named list of replacement names i.e. the state name x1 becomes x[[1]] in the
  # state space function body, input "Ta" becomes u[[1]] etc.
  create.substitution.list <- function(names, arg.vector.name){
    calls <- lapply(seq_along(names), function(i) call("[[", as.name(arg.vector.name), as.numeric(i)))
    setNames(calls, names)
  }
  subsList <- c(
    create.substitution.list(private$names$obs, "y"),
    create.substitution.list(private$names$parameters, "p"),
    create.substitution.list(private$names$states, "x"),
    create.substitution.list(private$names$inputs, "u")
  )

  # helpers
  ##################################################

  # substitute variable names with vector indexing and deparse
  to_term <- function(expr) deparse1(do.call(substitute, list(expr, subsList)))

  # terms that depend on no state/parameter/input are plain numerics, and assigning
  # those into an advector throws a lost AD error, so wrap them in RTMB::AD
  wrap_const <- function(term) {
    if (grepl("\\b[xpu]\\[\\[", term)) term else sprintf("RTMB::AD(%s)", term)
  }

  # allocate an AD vector and fill non-zero entries
  vec_body <- function(terms) {
    nz <- which(terms != "0")
    paste(c(
      sprintf("ans <- RTMB::AD(numeric(%i))", length(terms)),
      sprintf("ans[[%i]] <- %s", nz, vapply(terms[nz], wrap_const, "", USE.NAMES = FALSE))
    ), collapse = "; ")
  }

  # allocate an AD matrix and fill non-zero entries
  mat_body <- function(terms) {
    nz <- which(terms != "0", arr.ind = TRUE)
    paste(c(
      sprintf("ans <- RTMB::AD(matrix(0, nrow=%i, ncol=%i))", nrow(terms), ncol(terms)),
      sprintf("ans[[%i, %i]] <- %s", nz[,1], nz[,2], vapply(terms[nz], wrap_const, "", USE.NAMES = FALSE))
    ), collapse = "; ")
  }

# Functions
  ##################################################

  # drift
  f.terms <- vapply(seq_len(n.states), function(i) to_term(private$model$diff.terms[[i]]$dt), "")
  private$model$state.space.fun.strings$f <- vec_body(f.terms)

  # drift jacobian
  dfdx.terms <- matrix("0", n.states, n.states)
  for(i in seq_len(n.states)){
    for(j in seq_len(n.states)){
      dfdx.terms[i,j] <- to_term(private$model$diff.terms.drift[[i]][[j]])
    }
  }
  private$model$state.space.fun.strings$dfdx <- mat_body(dfdx.terms)

  # diffusion
  g.terms <- matrix("0", n.states, n.diffusions)
  for(i in seq_len(n.states)){
    for(j in seq_len(n.diffusions)){
      g.terms[i,j] <- to_term(private$model$diff.terms[[i]][[j+1]])
    }
  }
  private$model$state.space.fun.strings$g <- mat_body(g.terms)

  # observation
  h.terms <- vapply(seq_len(n.obs), function(i) to_term(private$model$obs.eqs.trans[[i]]$rhs), "")
  private$model$state.space.fun.strings$h <- vec_body(h.terms)

  # observation jacobian
  dhdx.terms <- matrix("0", n.obs, n.states)
  for(i in seq_len(n.obs)){
    for(j in seq_len(n.states)){
      dhdx.terms[i,j] <- to_term(private$model$diff.terms.obs[[i]][[j]])
    }
  }
  private$model$state.space.fun.strings$dhdx <- mat_body(dhdx.terms)

  # observation variance (vector and diagonal matrix form)
  hvar.terms <- vapply(seq_len(n.obs), function(i) to_term(private$model$obs.var.trans[[i]]$rhs), "")
  private$model$state.space.fun.strings$hvar <- vec_body(hvar.terms)

  hvar.mat.terms <- matrix("0", n.obs, n.obs)
  diag(hvar.mat.terms) <- hvar.terms
  private$model$state.space.fun.strings$hvar_matrix <- mat_body(hvar.mat.terms)

  # derivative w.r.t inputs (for linear kalman filter)
  # B in dX = (A*X + B*U) dt + G dW with U = [1, u1(t), u2(t), ...]
  # first column: constant terms (all inputs and states set to zero)
  zero.list <- as.list(numeric(n.inputs + n.states))
  names(zero.list) <- c(private$names$inputs, private$names$states)
  constant.terms <- try_with_warning_recovery(
    unname(sapply(
      sapply(private$model$sys.eqs.trans, function(x) Deriv::Simplify(do.call(substitute, list(x$diff.dt, zero.list)))),
      function(x) deparse1(do.call(substitute, list(x, subsList)))
    ))
  )
  if(inherits(constant.terms, "try-error")){
    constant.terms <- rep("0", n.states)
  }

  # remaining columns: derivatives w.r.t. each input
  dfdu.terms <- matrix("0", n.states, n.inputs + 1)
  dfdu.terms[,1] <- constant.terms
  for(i in seq_len(n.states)){
    for(j in seq_len(n.inputs)){
      tempterm <- ctsmTMB_Deriv(f=private$model$sys.eqs.trans[[i]]$diff.dt, x=private$names$inputs[j])
      dfdu.terms[i,j+1] <- to_term(tempterm)
    }
  }
  private$model$state.space.fun.strings$dfdu <- mat_body(dfdu.terms)

  # return
  return(private$model$state.space.fun.strings)
}

configure_ad_tape <- function(context, self, private){

  if (context=="RTMB"){
    if(is.null(private$algo.settings$advanced.settings$rtmb.tapeconfig)){
      # Default settings
      # atomic: disabling yields 7x optimization speed
      x <- list(atomic="disable", comparison="NA", vectorize="NA")
      do.call(RTMB::TapeConfig, x)
    } else {
      # User-provided settings
      x <- as.list(private$algo.settings$advanced.settings$rtmb.tapeconfig)
      nams <- names(formals(RTMB::TapeConfig))
      arg.names <- nams[nams %in% names(x)]
      x <- x[arg.names]
      do.call(RTMB::TapeConfig, x)
    }
  }

  if(context=="TMB"){
    tmb.config.names <- c(
      "trace.atomic",
      "trace.optimize",
      "nthreads",
      "trace.parallel",
      "tmbad.atomic_sparse_log_determinant",
      "tmbad_deterministic_hash",
      "tape.parallel",
      "optimize.parallel",
      "autopar",
      "optimize.instantly",
      "debug.getListElement",
      "tmbad.sparse_hessian_compress"
    )
    if (is.null(private$algo.settings$advanced.settings$tmb.tapeconfig)) {
      # Default settings
      x <- list(
        trace.atomic = 1,
        trace.optimize = 1,
        nthreads = 1,
        trace.parallel = 1,
        tmbad.atomic_sparse_log_determinant = 1,
        tmbad_deterministic_hash = 1,
        tape.parallel = 1,
        optimize.parallel = 0,
        autopar = 0,
        optimize.instantly = 1,
        debug.getListElement = 0,
        tmbad.sparse_hessian_compress = 0
      )
      do.call(TMB::config, c(x, list(DLL=private$modelname.with.method)))
    } else {
      # User-provided settings
      x <- as.list(private$algo.settings$advanced.settings$tmb.tapeconfig)
      arg.names <- tmb.config.names[tmb.config.names %in% names(x)]
      x <- x[arg.names]
      do.call(TMB::config, c(x, list(DLL=private$modelname.with.method)))
    }
  }

  return(invisible(self))
}