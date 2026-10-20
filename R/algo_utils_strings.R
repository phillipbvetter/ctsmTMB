apply_forceAD_to_some_variables <- function(force.ad, .envir=parent.frame()) {

  list2env(as.list(.envir), envir = environment())

  method <- private$algo.settings$method

  # UKF related
  ukf.methods <- get.inference.algos("kalman.ukf")

  if (method %in% ukf.methods) {
    Fsigma <- array(list(), c(1,n.sigmapoints))
    Hsigma <- array(list(),c(1,n.sigmapoints))
    assign("Fsigma", Fsigma, envir = .envir)
    assign("Hsigma", Hsigma, envir = .envir)
  }

  # This "hack" where zero-matrices/vectors are created globally rather than inside the state space functions
  # reduces MakeADFun compilation time by roughly 20%. We must use force=TRUE to force the variables into ad
  # context. This works for the estimation, but breaks the report functionality because variables become
  # advectors regardless of calling method.
  if(force.ad) {

    if(method == "ekf"){
      stateVec <- RTMB::AD(stateVec, force=TRUE)
      covMat <- RTMB::AD(covMat, force=TRUE)
      inputMat <- RTMB::AD(inputMat,force=TRUE)

      assign("stateVec", stateVec, envir = .envir)
      assign("covMat", covMat, envir = .envir)
      assign("inputMat", inputMat, envir = .envir)
    }

    if(method == "lkf") {
      stateVec <- RTMB::AD(stateVec, force=TRUE)
      covMat <- RTMB::AD(covMat, force=TRUE)

      assign("stateVec", stateVec, envir = .envir)
      assign("covMat", covMat, envir = .envir)
    }

    # if(private$algo.settings$method == "ukf"){
    if(method %in% ukf.methods){
      stateVec <- RTMB::AD(stateVec, force=TRUE)
      covMat <- RTMB::AD(covMat, force=TRUE)
      inputMat <- RTMB::AD(inputMat,force=TRUE)
      # Fsigma <- array(list(),c(1,n.sigmapoints))
      # Hsigma <- array(list(),c(1,n.sigmapoints))

      assign("stateVec", stateVec, envir = .envir)
      assign("covMat", covMat, envir = .envir)
      assign("inputMat", inputMat, envir = .envir)
      # assign("Fsigma", Fsigma, envir = .envir)
      # assign("Hsigma", Hsigma, envir = .envir)
    }

  }

  return(invisible(NULL))

}

create_function_from_string_body <- function(fun.name, return.name, body.string, .envir=parent.frame()){
  func_string <- sprintf("%s <- function(stateVec, parVec, inputVec) { %s; return(%s) }", fun.name, body.string, return.name)
  func <- eval(parse(text = func_string), envir=.envir)
  return(NULL)
}

create_state_space_functions_for_estimation <- function(force.ad, .envir=parent.frame()){

  list2env(as.list(.envir), envir = environment())

  # This "hack" where zero-matrices/vectors are created globally rather than inside the state space functions
  # reduces MakeADFun compilation time by roughly 20%. We must use force=TRUE to force the variables into ad
  # context. This works for the estimation, but breaks the report functionality.
  if(force.ad) {

    # State space functions
    f_vec <- RTMB::AD(numeric(n.states), force=TRUE)
    dfdx_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.states, ncol=n.states), force=TRUE)
    g_mat <- RTMB::AD(RTMB::matrix(0,nrow=n.states, ncol=n.diffusions), force=TRUE)
    h_vec <- RTMB::AD(numeric(n.obs), force=TRUE)
    dhdx_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.obs, ncol=n.states), force=TRUE)
    hvar_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.obs, ncol=n.obs), force=TRUE)
    hvar_vec <- RTMB::AD(numeric(n.obs), force=TRUE)
    dfdu_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.states, ncol=n.inputs+1), force=TRUE)
    hvar_vec <- RTMB::AD(numeric(n.obs), force=TRUE)

    assign("f_vec", f_vec, envir = .envir)
    assign("dfdx_mat", dfdx_mat, envir = .envir)
    assign("g_mat", g_mat, envir = .envir)
    assign("h_vec", h_vec, envir = .envir)
    assign("dhdx_mat", dhdx_mat, envir = .envir)
    assign("hvar_mat", hvar_mat, envir = .envir)
    assign("hvar_vec", hvar_vec, envir = .envir)
    assign("dfdu_mat", dfdu_mat, envir = .envir)

    create_function_from_string_body("f__", "f_vec", private$model$rtmb.function.strings.indexed2$f, .envir=.envir)
    create_function_from_string_body("dfdx__", "dfdx_mat", private$model$rtmb.function.strings.indexed2$dfdx, .envir=.envir)
    create_function_from_string_body("g__",  "g_mat", private$model$rtmb.function.strings.indexed2$g, .envir=.envir)
    create_function_from_string_body("h__", "h_vec", private$model$rtmb.function.strings.indexed2$h, .envir=.envir)
    create_function_from_string_body("dhdx__", "dhdx_mat", private$model$rtmb.function.strings.indexed2$dhdx, .envir=.envir)
    create_function_from_string_body("hvar__matrix", "hvar_mat", private$model$rtmb.function.strings.indexed2$hvar__matrix, .envir=.envir)
    create_function_from_string_body("dfdu__", "dfdu_mat", private$model$rtmb.function.strings.indexed2$dfdu, .envir=.envir)
    create_function_from_string_body("hvar__", "hvar_vec", private$model$rtmb.function.strings.indexed2$hvar, .envir=.envir)

  } else {

    create_function_from_string_body("f__", "ans", private$model$r.function.strings$f, .envir=.envir)
    create_function_from_string_body("dfdx__", "ans", private$model$r.function.strings$dfdx, .envir=.envir)
    create_function_from_string_body("g__",  "ans", private$model$r.function.strings$g, .envir=.envir)
    create_function_from_string_body("h__", "ans", private$model$r.function.strings$h, .envir=.envir)
    create_function_from_string_body("dhdx__", "ans", private$model$r.function.strings$dhdx, .envir=.envir)
    create_function_from_string_body("hvar__matrix", "ans", private$model$r.function.strings$hvar__matrix, .envir=.envir)
    create_function_from_string_body("dfdu__", "ans", private$model$r.function.strings$dfdu, .envir=.envir)
    create_function_from_string_body("hvar__", "ans", private$model$r.function.strings$hvar, .envir=.envir)

  }

  return(invisible(NULL))
}

# Create and assign functions to caller environment
unpack_state_space_functions <- function(self, private, .envir=parent.frame()){

  fun_from_string <- function(fun.name, return.name, body.string, .envir=parent.frame()){
    func_string <- sprintf("%s <- function(stateVec, parVec, inputVec) { %s; return(%s) }", fun.name, body.string, return.name)
    func <- eval(parse(text = func_string), envir=.envir)
    return(NULL)
  }

  fun.strs <- private$model$r.function.strings
  names(fun.strs) <- paste0(names(fun.strs),"__")
  names(fun.strs)[names(fun.strs) == "hvar__matrix__"] = "hvar__matrix"
  nams <- names(fun.strs)

  e <- new.env()
  local({
    for(i in seq_along(nams)){
      do.call(fun_from_string, list(nams[i], "ans", fun.strs[[nams[i]]]))
    }
  },
  envir=e)

  # Unpack into the caller environment
  list2env(as.list(e, all.names=TRUE), envir=.envir)
  return(invisible(e))
}


create_state_space_functions_for_filtering <- function(.envir=parent.frame()){

  list2env(as.list(.envir), envir = environment())

  # UKF related
  if(private$algo.settings$method=="ukf"){
    Fsigma <- array(list(), c(1,n.sigmapoints))
    Hsigma <- array(list(),c(1,n.sigmapoints))
    assign("Fsigma", Fsigma, envir = .envir)
    assign("Hsigma", Hsigma, envir = .envir)
  }

  create_function_from_string_body("f__", "ans", private$model$r.function.strings$f, .envir=.envir)
  create_function_from_string_body("dfdx__", "ans", private$model$r.function.strings$dfdx, .envir=.envir)
  create_function_from_string_body("g__",  "ans", private$model$r.function.strings$g, .envir=.envir)
  create_function_from_string_body("h__", "ans", private$model$r.function.strings$h, .envir=.envir)
  create_function_from_string_body("dhdx__", "ans", private$model$r.function.strings$dhdx, .envir=.envir)
  create_function_from_string_body("hvar__matrix", "ans", private$model$r.function.strings$hvar__matrix, .envir=.envir)
  create_function_from_string_body("dfdu__", "ans", private$model$r.function.strings$dfdu, .envir=.envir)
  create_function_from_string_body("hvar__", "ans", private$model$r.function.strings$hvar, .envir=.envir)

  return(invisible(NULL))
}

create_state_space_function_strings_new = function(self, private)
{

  obsList = lapply(seq_along(private$names$obs), function(id) substitute(obsVec[[i]],list(i=as.numeric(id))))
  parList = lapply(seq_along(private$names$parameters), function(id) substitute(parVec[[i]],list(i=as.numeric(id))))
  stateList = lapply(seq_along(private$names$states), function(id) substitute(stateVec[[i]],list(i=as.numeric(id))))
  inputList = lapply(seq_along(private$names$inputs), function(id) substitute(inputVec[[i]],list(i=as.numeric(id))))
  names(obsList) = private$names$obs
  names(parList) = private$names$parameters
  names(stateList) = private$names$states
  names(inputList) = private$names$inputs
  subsList = c(obsList, parList, stateList, inputList)

  n.states <- private$dims$states
  n.obs <- private$dims$observations
  n.diffusions <- length(private$model$diff.processes[-1])
  n.inputs <- private$dims$inputs

  ##################################################
  # helpers
  ##################################################

  # substitute variable names with vector indexing and deparse
  to_term <- function(expr) deparse1(do.call(substitute, list(expr, subsList)))

  # terms that depend on no state/parameter/input are plain numerics, and assigning
  # those into an advector throws a lost AD error, so wrap them in RTMB::AD
  wrap_const <- function(term) {
    if (grepl("stateVec|parVec|inputVec", term)) term else sprintf("RTMB::AD(%s)", term)
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
  private$model$r.function.strings$f <- vec_body(f.terms)

  ##################################################
  # drift jacobian
  ##################################################
  dfdx.terms <- matrix("0", n.states, n.states)
  for(i in seq_len(n.states)){
    for(j in seq_len(n.states)){
      dfdx.terms[i,j] <- to_term(private$model$diff.terms.drift[[i]][[j]])
    }
  }
  private$model$r.function.strings$dfdx <- mat_body(dfdx.terms)

  ##################################################
  # diffusion
  ##################################################
  g.terms <- matrix("0", n.states, n.diffusions)
  for(i in seq_len(n.states)){
    for(j in seq_len(n.diffusions)){
      g.terms[i,j] <- to_term(private$model$diff.terms[[i]][[j+1]])
    }
  }
  private$model$r.function.strings$g <- mat_body(g.terms)

  ##################################################
  # observation
  ##################################################
  h.terms <- vapply(seq_len(n.obs), function(i) to_term(private$model$obs.eqs.trans[[i]]$rhs), "")
  private$model$r.function.strings$h <- vec_body(h.terms)

  ##################################################
  # observation jacobian
  ##################################################
  dhdx.terms <- matrix("0", n.obs, n.states)
  for(i in seq_len(n.obs)){
    for(j in seq_len(n.states)){
      dhdx.terms[i,j] <- to_term(private$model$diff.terms.obs[[i]][[j]])
    }
  }
  private$model$r.function.strings$dhdx <- mat_body(dhdx.terms)

  ##################################################
  # observation variance (vector and diagonal matrix form)
  ##################################################
  hvar.terms <- vapply(seq_len(n.obs), function(i) to_term(private$model$obs.var.trans[[i]]$rhs), "")
  private$model$r.function.strings$hvar <- vec_body(hvar.terms)

  hvar.mat.terms <- matrix("0", n.obs, n.obs)
  diag(hvar.mat.terms) <- hvar.terms
  private$model$r.function.strings$hvar__matrix <- mat_body(hvar.mat.terms)

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
  private$model$r.function.strings$dfdu <- mat_body(dfdu.terms)

  return(invisible(NULL))
}


create_state_space_function_strings = function(self, private)
{

  obsList = lapply(seq_along(private$names$obs), function(id) substitute(obsVec[[i]],list(i=as.numeric(id))))
  parList = lapply(seq_along(private$names$parameters), function(id) substitute(parVec[[i]],list(i=as.numeric(id))))
  stateList = lapply(seq_along(private$names$states), function(id) substitute(stateVec[[i]],list(i=as.numeric(id))))
  inputList = lapply(seq_along(private$names$inputs), function(id) substitute(inputVec[[i]],list(i=as.numeric(id))))
  names(obsList) = private$names$obs
  names(parList) = private$names$parameters
  names(stateList) = private$names$states
  names(inputList) = private$names$inputs
  subsList = c(obsList, parList, stateList, inputList)

  ##################################################
  # drift
  ##################################################
  f.elements <- c()
  f.elements2 <- numeric(private$dims$states)
  for(i in seq_along(private$names$states)){
    term <- deparse1(do.call(substitute, list(private$model$diff.terms[[i]]$dt, subsList)))
    # skip if zero
    if(term=="0") next
    f.elements2[i] <- term
    # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error so wrap in RTMB::AD
    out <- suppressWarnings(as.numeric(term))
    if(!is.na(out)){
      f.elements <- c(f.elements, sprintf("f_vec[[%s]] <- RTMB::AD(%s)",i,term))
    } else{
      f.elements <- c(f.elements, sprintf("f_vec[[%s]] <- %s",i,term))
    }
  }
  private$model$rtmb.function.strings.indexed2$f = paste(f.elements,collapse="; ")
  private$model$r.function.strings$f <- paste("ans <- c(",paste(f.elements2,collapse=", "),")")

  ##################################################
  # drift jacobian
  ##################################################
  f.elements = c()
  f.elements2 <- 0*diag(private$dims$states)
  for(i in seq_along(private$names$states)){
    for(j in seq_along(private$names$states)){
      term <- deparse1(do.call(substitute, list(private$model$diff.terms.drift[[i]][[j]], subsList)))
      # skip if zero
      if(term=="0") next
      f.elements2[i,j] <- term
      # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
      # so wrap in RTMB::AD
      out <- suppressWarnings(as.numeric(term))
      if(!is.na(out)){
        f.elements <- c(f.elements, sprintf("dfdx_mat[[%s,%s]] <- RTMB::AD(%s)",i,j,term))
      } else{
        f.elements <- c(f.elements, sprintf("dfdx_mat[[%s,%s]] <- %s",i,j,term))
      }
    }
  }
  private$model$rtmb.function.strings.indexed2$dfdx = paste(f.elements, collapse=";")
  private$model$r.function.strings$dfdx <- paste("ans <- RTMB::matrix(c(",paste(c(f.elements2),collapse=", "),"), nrow=n.states)")

  ##################################################
  # diffusion
  ##################################################
  f.elements = c()
  f.elements2 <- matrix(0, private$dims$states, private$dims$diffusions)
  for(i in seq_along(private$names$states)){
    for(j in seq_along(private$model$diff.processes[-1])){
      term <- deparse1(do.call(substitute, list(private$model$diff.terms[[i]][[j+1]], subsList)))
      # skip if zero
      if(term=="0") next
      f.elements2[i,j] <- term
      # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
      # so wrap in RTMB::AD
      out <- suppressWarnings(as.numeric(term))
      if(!is.na(out)){
        f.elements <- c(f.elements, sprintf("g_mat[[%s,%s]] <- RTMB::AD(%s)",i,j, term))
      } else{
        f.elements <- c(f.elements, sprintf("g_mat[[%s,%s]] <- %s",i,j, term))
      }
    }
  }
  private$model$rtmb.function.strings.indexed2$g = paste(f.elements,collapse=";")
  private$model$r.function.strings$g <- paste("ans <- RTMB::matrix(c(",paste(c(f.elements2),collapse=", "),"), nrow=n.states)")

  ##################################################
  # observation
  ##################################################
  f.elements <- c()
  f.elements2 <- numeric(private$dims$observations)
  for(i in seq_along(private$names$obs)){
    term <- deparse1(do.call(substitute, list(private$model$obs.eqs.trans[[i]]$rhs, subsList)))
    # skip if zero
    if(term=="0") next
    f.elements2[i] <- term
    # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
    # so wrap in RTMB::AD
    out <- suppressWarnings(as.numeric(term))
    if(!is.na(out)){
      # hack
      f.elements <- c(f.elements, sprintf("h_vec[[%s]] <- RTMB::AD(%s)",i, out))
    } else{
      f.elements <- c(f.elements, sprintf("h_vec[[%s]] <- %s",i, term))
    }
  }
  private$model$rtmb.function.strings.indexed2$h = paste(f.elements, collapse=";")
  private$model$r.function.strings$h <- paste("ans <- c(",paste(f.elements2,collapse=", "),")")
  # private$model$r.function.strings$h <- paste("ans <- RTMB::AD(c(",paste(f.elements2,collapse=", "),"))")

  ##################################################
  # observation jacobian
  ##################################################

  # calculate all the terms and substitute variables
  f.elements = c()
  f.elements2 <- matrix(0, private$dims$observations, private$dims$states)
  for(i in seq_along(private$names$obs)){
    for(j in seq_along(private$names$states)){
      term <- deparse1(do.call(substitute, list(private$model$diff.terms.obs[[i]][[j]], subsList)))
      # skip if zero
      if(term=="0") next
      f.elements2[i,j] <- term
      # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
      # so wrap in RTMB::AD
      out <- suppressWarnings(as.numeric(term))
      if(!is.na(out)){
        f.elements <- c(f.elements, sprintf("dhdx_mat[[%s,%s]] <- RTMB::AD(%s)",i,j,term))
      } else{
        f.elements <- c(f.elements, sprintf("dhdx_mat[[%s,%s]] <- %s",i,j,term))
      }
    }
  }
  private$model$rtmb.function.strings.indexed2$dhdx = paste(f.elements,collapse=";")
  private$model$r.function.strings$dhdx <- paste("ans <- RTMB::matrix(c(",paste(c(f.elements2),collapse=", "),"), ncol=n.states)")

  ##################################################
  # observation variance
  ##################################################

  f.elements <- c()
  f.elements2 <- numeric(private$dims$observations)
  for(i in seq_along(private$names$obs)){
    term <- deparse1(do.call(substitute, list(private$model$obs.var.trans[[i]]$rhs, subsList)))
    # skip if zero
    if(term=="0") next
    f.elements2[i] <- term
    # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
    # so wrap in RTMB::AD
    out <- suppressWarnings(as.numeric(term))
    if(!is.na(out)){
      # hack
      f.elements <- c(f.elements, sprintf("hvar_vec[[%s]] <- RTMB::AD(%s)",i, out))
    } else{
      f.elements <- c(f.elements, sprintf("hvar_vec[[%s]] <- %s",i, term))
    }
  }
  private$model$rtmb.function.strings.indexed2$hvar = paste(f.elements,collapse=";")
  private$model$r.function.strings$hvar <- paste("ans <- c(",paste(c(f.elements2),collapse=", "),")")
  # private$model$r.function.strings$hvar <- paste("ans <- RTMB::AD(c(",paste(f.elements2,collapse=", "),"))")

  ### MATRIX FORM ### (for kalman)
  f.elements <- c()
  f.elements2 <- numeric(private$dims$observations)
  for(i in seq_along(private$names$obs)){
    term <- deparse1(do.call(substitute, list(private$model$obs.var.trans[[i]]$rhs, subsList)))
    # skip if zero
    if(term=="0") next
    f.elements2[i] <- term
    # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
    # so wrap in RTMB::AD
    out <- suppressWarnings(as.numeric(term))
    if(!is.na(out)){
      f.elements <- c(f.elements, sprintf("hvar_mat[[%s,%s]] <- RTMB::AD(%s)",i,i,term))
    } else{
      f.elements <- c(f.elements, sprintf("hvar_mat[[%s,%s]] <- %s",i,i,term))
    }
  }
  private$model$rtmb.function.strings.indexed2$hvar__matrix = paste(f.elements,collapse=";")
  private$model$r.function.strings$hvar__matrix <- paste("ans <- RTMB::diag(c(",paste(c(f.elements2),collapse=", "),"), nrow=n.obs)")

  ##################################################
  # derivative w.r.t inputs (for linear kalman filter)
  ##################################################

  # This derivative is concerned with the linear system
  # dX = (A*X + B*U) dt + G dW
  # We need to find B. The vector U = [1, u1(t), u2(t),...]
  # We first scalar covers "constants and parameters e.g. in the Ornstein Uhlenbeck
  # dx = theta * (mu - x + u(t))*dt + sigma * dw
  # then B = [theta*mu, 0, theta]
  # and U = Transpose([1, t, u(t)])
  # Therefore: For each system equation - find all constants, and derivatives w.r.t inputs

  # 1. We first find constant terms
  # this corresponds to an input that is always 1 (first element of U)
  zero.list <- as.list(numeric(private$dims$inputs + private$dims$states))
  names(zero.list) <-c(private$names$inputs, private$names$states)
  constant.terms <- try_with_warning_recovery(
    unname(sapply(
      sapply(private$model$sys.eqs.trans, function(x) Deriv::Simplify(do.call(substitute, list(x$diff.dt, zero.list)))),
      function(x) deparse1(do.call(substitute, list(x, subsList)))
    ))
  )
  if(inherits(constant.terms, "try-error")){
    constant.terms <- rep(0, private$dims$states)
  }

  f.elements <- c()
  # f.elements2 <- matrix(0, private$dims$states, private$dims$inputs+1)
  f.elements2 <- matrix("RTMB::AD(0)", private$dims$states, private$dims$inputs+1)
  # 2. Now we find input-terms by differentiation
  for(i in seq_along(private$names$states)){
    #its inputs + 1 below because the first element is for constants
    for(j in 1:(private$dims$inputs+1)){

      if(j==1){
        term <- constant.terms[i]
      } else {
        tempterm <- ctsmTMB_Deriv(f=private$model$sys.eqs.trans[[i]]$diff.dt, x=private$names$inputs[j-1])
        term <- deparse1(do.call(substitute, list(tempterm, subsList)))
      }

      # skip if zero
      if(term=="0") next
      f.elements2[i,j] <- term

      # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
      # so wrap in RTMB::AD
      out <- suppressWarnings(as.numeric(term))
      if(!is.na(out)){
        f.elements <- c(f.elements, sprintf("dfdu_mat[[%s,%s]] <- RTMB::AD(%s)",i,j,term))
      } else{
        f.elements <- c(f.elements, sprintf("dfdu_mat[[%s,%s]] <- %s",i,j,term))
      }
    }
  }
  private$model$rtmb.function.strings.indexed2$dfdu = paste(f.elements,collapse=";")
  private$model$r.function.strings$dfdu <- paste("ans <- RTMB::matrix(c(",paste(c(f.elements2),collapse=", "),"), nrow=n.states)")

  return(invisible(NULL))
}

##########################################################
##########################################################
##########################################################
# Construct system functions for Rcpp implementations
##########################################################
##########################################################
##########################################################

create_rcpp_state_space_function_strings = function(self, private){

  # Create substitution translation list
  obsList = lapply(seq_along(private$names$obs), function(id) substitute(obsVec(i),list(i=as.numeric(id-1))))
  parList = lapply(seq_along(private$names$parameters), function(id) substitute(parVec(i),list(i=as.numeric(id-1))))
  stateList = lapply(seq_along(private$names$states), function(id) substitute(stateVec(i),list(i=as.numeric(id-1))))
  inputList = lapply(seq_along(private$names$inputs), function(id) substitute(inputVec(i),list(i=as.numeric(id-1))))
  names(obsList) = private$names$obs
  names(parList) = private$names$parameters
  names(stateList) = private$names$states
  names(inputList) = private$names$inputs
  subsList = c(obsList, parList, stateList, inputList)
  # subsList <- get_substitution_list(self, private)


  #################################################
  # drift
  ##################################################

  # Perform substitution of parameters, inputs and states
  f = sapply( seq_along(private$names$states),
              function(i){
                drift.term = hat2pow(private$model$diff.terms[[i]]$dt)
                new.drift.term = do.call(substitute, list(drift.term, subsList))
                sprintf("f(%i) = %s;",i-1,deparse1(new.drift.term))
              })

  code = sprintf("Eigen::VectorXd f(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd f(%s);
                 %s
                 return f;
                 }",private$dims$states, paste(f,collapse=""))

  private$model$rcpp.function.strings$f = code

  code = sprintf("Eigen::VectorXd f_const(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec){
                 Eigen::VectorXd f(%s);
                 %s
                 return f;
                 }",private$dims$states, paste(f,collapse="\n "))
  private$model$rcpp.function.strings$f_const = code

  f = sapply( seq_along(private$names$states),
              function(i){
                drift.term = hat2pow(private$model$diff.terms[[i]]$dt)
                new.drift.term = do.call(substitute, list(drift.term, subsList))
                sprintf("out(%i) = %s;",i-1,deparse1(new.drift.term))
              })
  code = sprintf("void f_prealoc(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec, Eigen::VectorXd& out){
                 %s
                 }", paste(f, collapse="\n "))
  private$model$rcpp.function.strings$f_prealoc = code

  ##################################################
  # drift jacobian
  ##################################################

  # calculate all the terms and substitute variables
  dfdx = c()
  for(i in seq_along(private$names$states)){
    for(j in seq_along(private$names$states)){
      term = hat2pow(ctsmTMB_Deriv(f=private$model$diff.terms[[i]]$dt, x=private$names$states[j]))
      new.term = do.call(substitute, list(term, subsList))
      dfdx = c(dfdx, sprintf("dfdx(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  code = sprintf("Eigen::MatrixXd dfdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dfdx(%s,%s);
                 %s
                 return dfdx;
                 }",private$dims$states, private$dims$states, paste(dfdx,collapse=""))

  private$model$rcpp.function.strings$dfdx = code

  code = sprintf("Eigen::MatrixXd dfdx_const(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec){
                 Eigen::MatrixXd dfdx(%s,%s);
                 %s
                 return dfdx;
                 }",private$dims$states, private$dims$states, paste(dfdx,collapse=""))
  private$model$rcpp.function.strings$dfdx_const = code

  dfdx = c()
  for(i in seq_along(private$names$states)){
    for(j in seq_along(private$names$states)){
      term = hat2pow(ctsmTMB_Deriv(f=private$model$diff.terms[[i]]$dt, x=private$names$states[j]))
      new.term = do.call(substitute, list(term, subsList))
      dfdx = c(dfdx, sprintf("out(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  code = sprintf("void dfdx_prealoc(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec, Eigen::MatrixXd& out){
                 %s
                 }", paste(dfdx,collapse=""))
  private$model$rcpp.function.strings$dfdx_prealoc = code

  ##################################################
  # diffusion
  ##################################################

  # calculate all the terms and substitute variables
  g = c()
  for(i in seq_along(private$names$states)){
    for(j in seq_along(private$model$diff.processes[-1])){
      term = hat2pow(private$model$diff.terms[[i]][[j+1]])
      new.term = do.call(substitute, list(term, subsList))
      g = c(g, sprintf("g(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  code = sprintf("Eigen::MatrixXd g(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd g(%s,%s);
                 %s
                 return g;
                 }",private$dims$states, private$dims$diffusions, paste(g,collapse=""))

  private$model$rcpp.function.strings$g = code

  code = sprintf("Eigen::MatrixXd g_const(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec){
                 Eigen::MatrixXd g(%s,%s);
                 %s
                 return g;
                 }",private$dims$states, private$dims$diffusions, paste(g,collapse=""))
  private$model$rcpp.function.strings$g_const = code

  g = c()
  for(i in seq_along(private$names$states)){
    for(j in seq_along(private$model$diff.processes[-1])){
      term = hat2pow(private$model$diff.terms[[i]][[j+1]])
      new.term = do.call(substitute, list(term, subsList))
      g = c(g, sprintf("out(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  code = sprintf("void g_prealoc(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec, Eigen::MatrixXd& out) {
                 %s
                 }", paste(g,collapse="\n "))

  private$model$rcpp.function.strings$g_prealoc = code

  ##################################################
  # observation
  ##################################################

  # calculate all the terms and substitute variables
  h = sapply(seq_along(private$names$obs),
             function(i){
               term = hat2pow(private$model$obs.eqs.trans[[i]]$rhs)
               new.term = do.call(substitute, list(term, subsList))
               sprintf("h(%s) = %s;",i-1, deparse1(new.term))
             })
  code = sprintf("Eigen::VectorXd h(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd h(%s);
                 %s
                 return h;
                 }",private$dims$observations, paste(h,collapse=""))
  private$model$rcpp.function.strings$h = code

  code = sprintf("Eigen::VectorXd h_const(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec){
                 Eigen::VectorXd h(%s);
                 %s
                 return h;
                 }",private$dims$observations, paste(h,collapse=""))
  private$model$rcpp.function.strings$h_const = code

  h = sapply(seq_along(private$names$obs),
             function(i){
               term = hat2pow(private$model$obs.eqs.trans[[i]]$rhs)
               new.term = do.call(substitute, list(term, subsList))
               sprintf("out(%s) = %s;",i-1, deparse1(new.term))
             })
  code = sprintf("void h_prealoc(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec, Eigen::VectorXd& out) {
                 %s
                 }", paste(h,collapse="\n "))
  private$model$rcpp.function.strings$h_prealoc = code

  ##################################################
  # observation jacobian
  ##################################################

  # calculate all the terms and substitute variables
  dhdx = c()
  for(i in seq_along(private$names$obs)){
    for(j in seq_along(private$names$states)){
      term = hat2pow(private$model$diff.terms.obs[[i]][[j]])
      new.term = do.call(substitute, list(term, subsList))
      dhdx = c(dhdx, sprintf("dhdx(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  code = sprintf("Eigen::MatrixXd dhdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dhdx(%s,%s);
                 %s
                 return dhdx;
                 }", private$dims$observations, private$dims$states, paste(dhdx,collapse=""))
  private$model$rcpp.function.strings$dhdx = code

  code = sprintf("Eigen::MatrixXd dhdx_const(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec){
                 Eigen::MatrixXd dhdx(%s,%s);
                 %s
                 return dhdx;
                 }", private$dims$observations, private$dims$states, paste(dhdx,collapse=""))
  private$model$rcpp.function.strings$dhdx_const = code

  dhdx = c()
  for(i in seq_along(private$names$obs)){
    for(j in seq_along(private$names$states)){
      term = hat2pow(private$model$diff.terms.obs[[i]][[j]])
      new.term = do.call(substitute, list(term, subsList))
      dhdx = c(dhdx, sprintf("out(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  code = sprintf("void dhdx_prealoc(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec, Eigen::MatrixXd& out) {
                 %s
                 }", paste(dhdx,collapse=""))
  private$model$rcpp.function.strings$dhdx_prealoc = code

  ##################################################
  # observation variance
  ##################################################

  hvar = lapply(seq_along(private$model$obs.var.trans),
                function(i) {
                  term = hat2pow(private$model$obs.var.trans[[i]]$rhs)
                  new.term = do.call(substitute, list(term, subsList))
                  sprintf("hvar(%s,%s) = %s;", i-1, i-1, deparse1(new.term))
                })
  code = sprintf("Eigen::MatrixXd hvar(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd hvar(%s,%s);
                 hvar.setZero();
                 %s
                 return hvar;
                 }", private$dims$observations, private$dims$observations, paste(hvar,collapse=""))
  private$model$rcpp.function.strings$hvar = code

  code = sprintf("Eigen::MatrixXd hvar_const(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec){
                 Eigen::MatrixXd hvar(%s,%s);
                 hvar.setZero();
                 %s
                 return hvar;
                 }", private$dims$observations, private$dims$observations, paste(hvar,collapse=""))
  private$model$rcpp.function.strings$hvar_const = code

  hvar = lapply(seq_along(private$model$obs.var.trans),
                function(i) {
                  term = hat2pow(private$model$obs.var.trans[[i]]$rhs)
                  new.term = do.call(substitute, list(term, subsList))
                  sprintf("out(%s,%s) = %s;", i-1, i-1, deparse1(new.term))
                })
  code = sprintf("void hvar_prealoc(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec, Eigen::MatrixXd& out) {
                 out.setZero();
                 %s
                 }", paste(hvar,collapse=""))
  private$model$rcpp.function.strings$hvar_prealoc = code

  ##################################################
  # observation variance vector form
  ##################################################

  hvar = lapply(seq_along(private$model$obs.var.trans),
                function(i) {
                  term = hat2pow(private$model$obs.var.trans[[i]]$rhs)
                  new.term = do.call(substitute, list(term, subsList))
                  sprintf("hvar(%s) = %s;", i-1, deparse1(new.term))
                })
  code = sprintf("Eigen::ArrayXd hvar_array(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::ArrayXd hvar(%s);
                 hvar.setZero();
                 %s
                 return hvar;
                 }", private$dims$observations, paste(hvar,collapse=""))
  private$model$rcpp.function.strings$hvar_array = code

  code = sprintf("Eigen::ArrayXd hvar_array_const(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec){
                 Eigen::ArrayXd hvar(%s);
                 hvar.setZero();
                 %s
                 return hvar;
                 }", private$dims$observations, paste(hvar,collapse=""))
  private$model$rcpp.function.strings$hvar_array_const = code

  ##################################################
  # derivative w.r.t inputs (for linear kalman filter)
  ##################################################

  # We first find constant terms
  # this corresponds to an input that is always 1 (first element of U)
  zero.list <- as.list(numeric(private$dims$inputs + private$dims$states))
  names(zero.list) <-c(private$names$inputs, private$names$states)
  constant.terms <- sapply(private$model$sys.eqs.trans, function(x) Deriv::Simplify(do.call(substitute, list(x$diff.dt, zero.list))))
  if(inherits(constant.terms, "try-error")){
    constant.terms <- rep(0, private$dims$states)
  }

  jac.f.wrt.u = c()
  # 2. Now we find input-terms by differentiation
  for(i in seq_along(private$names$states)){
    #its inputs + 1 below because the first element is for constants
    for(j in 1:(private$dims$inputs+1)){
      if(j==1){
        term <- constant.terms[[i]]
      } else {
        term <- ctsmTMB_Deriv(f=private$model$sys.eqs.trans[[i]]$diff.dt, x=private$names$inputs[j-1])
      }
      # skip if zero
      if(term=="0") next
      term = hat2pow(term)
      new.term = do.call(substitute, list(term, subsList))
      jac.f.wrt.u = c(jac.f.wrt.u, sprintf("ans(%s, %s) = %s;", i-1, j-1, deparse1(new.term)))
    }
  }

  code = sprintf("Eigen::MatrixXd dfdu(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd ans(%s,%s);
                 ans.setZero();
                 %s
                 return ans;
                 }", private$dims$states, private$dims$inputs+1, paste(jac.f.wrt.u, collapse=""))

  private$model$rcpp.function.strings$dfdu = code

  code = sprintf("Eigen::MatrixXd dfdu_const(const Eigen::VectorXd& stateVec, const Eigen::VectorXd& parVec, const Eigen::VectorXd& inputVec){
                 Eigen::MatrixXd ans(%s,%s);
                 ans.setZero();
                 %s
                 return ans;
                 }", private$dims$states, private$dims$inputs+1, paste(jac.f.wrt.u, collapse=""))
  private$model$rcpp.function.strings$dfdu_const = code


}
