# loss function
# ------------------------------------
construct.loss.function <- function(dim = c("scalar","multi"), type = c("quadratic", "huber", "tukey"), c = 3) {
  type <- match.arg(type)
  dim <- match.arg(dim)

  # quadratic
  if (type == "quadratic") {
    loss <- switch(
      dim,
      scalar = function(e,R) -RTMB::dnorm(e, sd=sqrt(R), log=TRUE),
      multi = function(e, R) -RTMB::dmvnorm(e, Sigma = R, log = TRUE)
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
      multi = function(e,R) 0.5 * (logdet(R) + length(e) * log(2*pi) + rho(sum(e*RTMB::solve(R,e))))
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

  # STATIONARY SOLVER FOR THE MEAN
  # --------------------------------
  arglist <- list(rep(1,n.states),rep(1,n.pars),rep(1,n.inputs))
  f <- RTMB::MakeTape(function(z) sum(do.call(funs$f,z)^2), arglist)
  stationary.state <- f$newton(seq_len(n.states),grad_tol=1e-20)

  # STATIONARY SOLVER FOR COVARIANCE
  # # --------------------------------
  stationary.covariance <- function(A, G) {
    P <- solve_lyapunov(c(A, G %*% t(G)))
    dim(P) <- dim(A)
    P
  }

  # Return solver
  # get.cov is controlled by the second entry (state covariance) of the
  # algo.settings$estimate.initial field set via set_initial_state_estimation
  stationary.solve <- function(p,u, get.cov=private$algo.settings$estimate.initial[2]){
    z <- list(p,u)
    x <- stationary.state(z)
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

create.model.funs <- function(self, private){

  # We need the parent to be ctsmTMB to find custom functions like
  # e.g. invlogit
  fun.env <- new.env(parent = asNamespace("ctsmTMB"))
  as_function <- function(body) {
    eval(parse(text = sprintf(
      "function(x,p,u) {%s; ans}", body)),
      envir = fun.env)
  }

  fun.strs <- private$model$state.space.fun.strings
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

  return(funs)
}

create.state.space.funs = function(self, private)
{

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

  ##################################################
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

  ##################################################
  # drift
  ##################################################
  f.terms <- vapply(seq_len(n.states), function(i) to_term(private$model$diff.terms[[i]]$dt), "")
  private$model$state.space.fun.strings$f <- vec_body(f.terms)

  ##################################################
  # drift jacobian
  ##################################################
  dfdx.terms <- matrix("0", n.states, n.states)
  for(i in seq_len(n.states)){
    for(j in seq_len(n.states)){
      dfdx.terms[i,j] <- to_term(private$model$diff.terms.drift[[i]][[j]])
    }
  }
  private$model$state.space.fun.strings$dfdx <- mat_body(dfdx.terms)

  ##################################################
  # diffusion
  ##################################################
  g.terms <- matrix("0", n.states, n.diffusions)
  for(i in seq_len(n.states)){
    for(j in seq_len(n.diffusions)){
      g.terms[i,j] <- to_term(private$model$diff.terms[[i]][[j+1]])
    }
  }
  private$model$state.space.fun.strings$g <- mat_body(g.terms)

  ##################################################
  # observation
  ##################################################
  h.terms <- vapply(seq_len(n.obs), function(i) to_term(private$model$obs.eqs.trans[[i]]$rhs), "")
  private$model$state.space.fun.strings$h <- vec_body(h.terms)

  ##################################################
  # observation jacobian
  ##################################################
  dhdx.terms <- matrix("0", n.obs, n.states)
  for(i in seq_len(n.obs)){
    for(j in seq_len(n.states)){
      dhdx.terms[i,j] <- to_term(private$model$diff.terms.obs[[i]][[j]])
    }
  }
  private$model$state.space.fun.strings$dhdx <- mat_body(dhdx.terms)

  ##################################################
  # observation variance (vector and diagonal matrix form)
  ##################################################
  hvar.terms <- vapply(seq_len(n.obs), function(i) to_term(private$model$obs.var.trans[[i]]$rhs), "")
  private$model$state.space.fun.strings$hvar <- vec_body(hvar.terms)

  hvar.mat.terms <- matrix("0", n.obs, n.obs)
  diag(hvar.mat.terms) <- hvar.terms
  private$model$state.space.fun.strings$hvar_matrix <- mat_body(hvar.mat.terms)

  ##################################################
  # derivative w.r.t inputs (for linear kalman filter)
  ##################################################
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

  return(invisible(NULL))
}
