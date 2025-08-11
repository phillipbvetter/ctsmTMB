
create.state.space.functions.for.estimation <- function(forceAD, .envir=parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  # This "hack" where zero-matrices/vectors and initialized and not created inside the stace space functions
  # reduces MakeADFun compilation time by roughly 20%, but requires this unintended use of forcing AD which is
  # a bit unstable (but works!) and breaks for instance REPORT
  if(forceAD) {
    
    if(private$method == "ekf"){
      
      stateVec <- RTMB::AD(stateVec, force=TRUE)
      covMat <- RTMB::AD(covMat, force=TRUE)
      assign("stateVec", stateVec, envir = .envir)
      assign("covMat", covMat, envir = .envir)
      
      inputMat <- RTMB::AD(inputMat,force=TRUE)
      assign("inputMat", inputMat, envir = .envir)
      
    }
    
    if(private$method == "lkf") {
      
      stateVec <- RTMB::AD(stateVec, force=TRUE)
      covMat <- RTMB::AD(covMat, force=TRUE)
      assign("stateVec", stateVec, envir = .envir)
      assign("covMat", covMat, envir = .envir)
      
    }
    
    if(private$method == "ukf"){
      stateVec <- RTMB::AD(stateVec, force=TRUE)
      covMat <- RTMB::AD(covMat, force=TRUE)
      assign("stateVec", stateVec, envir = .envir)
      assign("covMat", covMat, envir = .envir)
      
      inputMat <- RTMB::AD(inputMat,force=TRUE)
      assign("inputMat", inputMat, envir = .envir)
      
      Fsigma <- array(list(),c(1,n.sigmapoints))
      Hsigma <- array(list(),c(1,n.sigmapoints))
      assign("Fsigma", Fsigma, envir = .envir)
      assign("Hsigma", Hsigma, envir = .envir)
      
    }
    
    # State space functions
    f_vec <- RTMB::AD(numeric(n.states), force=TRUE)
    dfdx_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.states, ncol=n.states), force=TRUE)
    g_mat <- RTMB::AD(RTMB::matrix(0,nrow=n.states, ncol=n.diffusions), force=TRUE)
    h_vec <- RTMB::AD(numeric(n.obs), force=TRUE)
    dhdx_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.obs, ncol=n.states), force=TRUE)
    hvar_mat <- RTMB::AD(RTMB::matrix(0, nrow=n.obs, ncol=n.obs), force=TRUE)
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
    
    create.function.from.string.body("f__", "f_vec", private$rtmb.function.strings.indexed2$f, .envir=.envir)
    create.function.from.string.body("dfdx__", "dfdx_mat", private$rtmb.function.strings.indexed2$dfdx, .envir=.envir)
    create.function.from.string.body("g__",  "g_mat", private$rtmb.function.strings.indexed2$g, .envir=.envir)
    create.function.from.string.body("h__", "h_vec", private$rtmb.function.strings.indexed2$h, .envir=.envir)
    create.function.from.string.body("dhdx__", "dhdx_mat", private$rtmb.function.strings.indexed2$dhdx, .envir=.envir)
    create.function.from.string.body("hvar__matrix", "hvar_mat", private$rtmb.function.strings.indexed2$hvar__matrix, .envir=.envir)
    create.function.from.string.body("dfdu__", "dfdu_mat", private$rtmb.function.strings.indexed2$dfdu, .envir=.envir)
    create.function.from.string.body("hvar__", "hvar_vec", private$rtmb.function.strings.indexed2$hvar, .envir=.envir)
    
  } else {
    
    # UKF related
    if(private$method=="ukf"){
      Fsigma <- array(list(), c(1,n.sigmapoints))
      Hsigma <- array(list(),c(1,n.sigmapoints))
      assign("Fsigma", Fsigma, envir = .envir)
      assign("Hsigma", Hsigma, envir = .envir)
    }
    
    create.function.from.string.body("f__", "ans", private$r.function.strings$f, .envir=.envir)
    create.function.from.string.body("dfdx__", "ans", private$r.function.strings$dfdx, .envir=.envir)
    create.function.from.string.body("g__",  "ans", private$r.function.strings$g, .envir=.envir)
    create.function.from.string.body("h__", "ans", private$r.function.strings$h, .envir=.envir)
    create.function.from.string.body("dhdx__", "ans", private$r.function.strings$dhdx, .envir=.envir)
    create.function.from.string.body("hvar__matrix", "ans", private$r.function.strings$hvar__matrix, .envir=.envir)
    create.function.from.string.body("dfdu__", "ans", private$r.function.strings$dfdu, .envir=.envir)
    create.function.from.string.body("hvar__", "ans", private$r.function.strings$hvar, .envir=.envir)
    
  }
  
  return(NULL)
}


create.state.space.functions.for.filtering <- function(.envir=parent.frame()){
  
  list2env(as.list(.envir), envir = environment())
  
  # UKF related
  if(private$method=="ukf"){
    Fsigma <- array(list(), c(1,n.sigmapoints))
    Hsigma <- array(list(),c(1,n.sigmapoints))
    assign("Fsigma", Fsigma, envir = .envir)
    assign("Hsigma", Hsigma, envir = .envir)
  }
  
  create.function.from.string.body("f__", "ans", private$r.function.strings$f, .envir=.envir)
  create.function.from.string.body("dfdx__", "ans", private$r.function.strings$dfdx, .envir=.envir)
  create.function.from.string.body("g__",  "ans", private$r.function.strings$g, .envir=.envir)
  create.function.from.string.body("h__", "ans", private$r.function.strings$h, .envir=.envir)
  create.function.from.string.body("dhdx__", "ans", private$r.function.strings$dhdx, .envir=.envir)
  create.function.from.string.body("hvar__matrix", "ans", private$r.function.strings$hvar__matrix, .envir=.envir)
  create.function.from.string.body("dfdu__", "ans", private$r.function.strings$dfdu, .envir=.envir)
  create.function.from.string.body("hvar__", "ans", private$r.function.strings$hvar, .envir=.envir)
  
  return(NULL)
}


create.function.from.string.body <- function(fun.name, return.name, body.string, .envir=parent.frame()){
  func_string <- sprintf("%s <- function(stateVec, parVec, inputVec) { %s; return(%s) }", fun.name, body.string, return.name)
  func <- eval(parse(text = func_string), envir=.envir)
  return(NULL)
}

create.stace.space.function.strings = function(self, private)
{
  
  obsList = lapply(seq_along(private$obs.names), function(id) substitute(obsVec[[i]],list(i=as.numeric(id))))
  parList = lapply(seq_along(private$parameter.names), function(id) substitute(parVec[[i]],list(i=as.numeric(id))))
  stateList = lapply(seq_along(private$state.names), function(id) substitute(stateVec[[i]],list(i=as.numeric(id))))
  inputList = lapply(seq_along(private$input.names), function(id) substitute(inputVec[[i]],list(i=as.numeric(id))))
  names(obsList) = private$obs.names
  names(parList) = private$parameter.names
  names(stateList) = private$state.names
  names(inputList) = private$input.names
  subsList = c(obsList, parList, stateList, inputList)
  
  ##################################################
  # drift
  ##################################################
  f.elements <- c()
  f.elements2 <- numeric(private$number.of.states)
  for(i in seq_along(private$state.names)){
    term <- deparse1(do.call(substitute, list(private$diff.terms[[i]]$dt, subsList)))
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
  private$rtmb.function.strings.indexed2$f = paste(f.elements,collapse="; ")
  private$r.function.strings$f <- paste("ans <- c(",paste(f.elements2,collapse=", "),")")
  
  ##################################################
  # drift jacobian
  ##################################################
  f.elements = c()
  f.elements2 <- 0*diag(private$number.of.states)
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term <- deparse1(do.call(substitute, list(private$diff.terms.drift[[i]][[j]], subsList)))
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
  private$rtmb.function.strings.indexed2$dfdx = paste(f.elements, collapse=";")
  private$r.function.strings$dfdx <- paste("ans <- RTMB::matrix(c(",paste(c(f.elements2),collapse=", "),"), nrow=n.states)")
  
  ##################################################
  # diffusion
  ##################################################
  f.elements = c()
  f.elements2 <- matrix(0, private$number.of.states, private$number.of.diffusions)
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term <- deparse1(do.call(substitute, list(private$diff.terms[[i]][[j+1]], subsList)))
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
  private$rtmb.function.strings.indexed2$g = paste(f.elements,collapse=";")
  private$r.function.strings$g <- paste("ans <- RTMB::matrix(c(",paste(c(f.elements2),collapse=", "),"), nrow=n.states)")
  
  ##################################################
  # observation
  ##################################################
  f.elements <- c()
  f.elements2 <- numeric(private$number.of.observations)
  for(i in seq_along(private$obs.names)){
    term <- deparse1(do.call(substitute, list(private$obs.eqs.trans[[i]]$rhs, subsList)))
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
  private$rtmb.function.strings.indexed2$h = paste(f.elements, collapse=";")
  private$r.function.strings$h <- paste("ans <- c(",paste(f.elements2,collapse=", "),")")
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  f.elements = c()
  f.elements2 <- matrix(0, private$number.of.observations, private$number.of.states)
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term <- deparse1(do.call(substitute, list(private$diff.terms.obs[[i]][[j]], subsList)))
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
  private$rtmb.function.strings.indexed2$dhdx = paste(f.elements,collapse=";")
  private$r.function.strings$dhdx <- paste("ans <- RTMB::matrix(c(",paste(c(f.elements2),collapse=", "),"), ncol=n.states)")
  
  ##################################################
  # observation variance
  ##################################################
  
  f.elements <- c()
  f2.elements <- numeric(private$number.of.observations)
  for(i in seq_along(private$obs.names)){
    term <- deparse1(do.call(substitute, list(private$obs.var.trans[[i]]$rhs, subsList)))
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
  private$rtmb.function.strings.indexed2$hvar = paste(f.elements,collapse=";")
  private$r.function.strings$hvar <- paste("ans <- c(",paste(c(f.elements2),collapse=", "),")")
  
  ### MATRIX FORM ### (for kalman)
  f.elements <- c()
  f.elements2 <- numeric(private$number.of.observations)
  for(i in seq_along(private$obs.names)){
    term <- deparse1(do.call(substitute, list(private$obs.var.trans[[i]]$rhs, subsList)))
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
  private$rtmb.function.strings.indexed2$hvar__matrix = paste(f.elements,collapse=";")
  private$r.function.strings$hvar__matrix <- paste("ans <- RTMB::diag(c(",paste(c(f.elements2),collapse=", "),"), nrow=n.obs)")
  
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
  zero.list <- as.list(numeric(private$number.of.inputs + private$number.of.states))
  names(zero.list) <-c(private$input.names, private$state.names)
  constant.terms <- try_withWarningRecovery(
    unname(sapply(
      sapply(private$sys.eqs.trans, function(x) Deriv::Simplify(do.call(substitute, list(x$diff.dt, zero.list)))),
      function(x) deparse1(do.call(substitute, list(x, subsList)))
    ))
  )
  if(inherits(constant.terms, "try-error")){
    constant.terms <- rep(0, private$number.of.states)
  }
  
  f.elements <- c()
  f.elements2 <- matrix(0, private$number.of.states, private$number.of.inputs+1)
  # 2. Now we find input-terms by differentiation
  for(i in seq_along(private$state.names)){
    #its inputs + 1 below because the first element is for constants
    for(j in 1:(private$number.of.inputs+1)){
      
      if(j==1){
        term <- constant.terms[i]
      } else {
        tempterm <- ctsmTMB.Deriv(f=private$sys.eqs.trans[[i]]$diff.dt, x=private$input.names[j-1])
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
  private$rtmb.function.strings.indexed2$dfdu = paste(f.elements,collapse=";")
  private$r.function.strings$dfdu <- paste("ans <- RTMB::matrix(c(",paste(c(f.elements2),collapse=", "),"), nrow=n.states)")
  
  return(invisible(NULL))
}

##########################################################
##########################################################
##########################################################
# USER-FUNCTION CONSTRUCTION FOR RCPP EKF
##########################################################
##########################################################
##########################################################

create.rcpp.stace.space.function.strings = function(self, private){
  
  # Create substitution translation list
  obsList = lapply(seq_along(private$obs.names), function(id) substitute(obsVec(i),list(i=as.numeric(id-1))))
  parList = lapply(seq_along(private$parameter.names), function(id) substitute(parVec(i),list(i=as.numeric(id-1))))
  stateList = lapply(seq_along(private$state.names), function(id) substitute(stateVec(i),list(i=as.numeric(id-1))))
  inputList = lapply(seq_along(private$input.names), function(id) substitute(inputVec(i),list(i=as.numeric(id-1))))
  names(obsList) = private$obs.names
  names(parList) = private$parameter.names
  names(stateList) = private$state.names
  names(inputList) = private$input.names
  subsList = c(obsList, parList, stateList, inputList)
  
  
  ##################################################
  # drift
  ##################################################
  
  # Perform substitution of parameters, inputs and states
  f = sapply( seq_along(private$state.names),
              function(i){
                drift.term = hat2pow(private$diff.terms[[i]]$dt)
                new.drift.term = do.call(substitute, list(drift.term, subsList))
                sprintf("f(%i) = %s;",i-1,deparse1(new.drift.term))
              })
  code = sprintf("Eigen::VectorXd f(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd f(%s); 
                 %s
                 return f;
                 }",private$number.of.states, paste(f,collapse=""))
  
  private$rcpp.function.strings$f = code
  
  ##################################################
  # drift jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dfdx = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      # term = hat2pow(Deriv::Deriv(private$diff.terms[[i]]$dt, x=private$state.names[j], cache.exp=F))
      term = hat2pow(ctsmTMB.Deriv(f=private$diff.terms[[i]]$dt, x=private$state.names[j]))
      new.term = do.call(substitute, list(term, subsList))
      dfdx = c(dfdx, sprintf("dfdx(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  code = sprintf("Eigen::MatrixXd dfdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dfdx(%s,%s);
                 %s
                 return dfdx;
                 }",private$number.of.states, private$number.of.states, paste(dfdx,collapse=""))
  
  private$rcpp.function.strings$dfdx = code
  
  ##################################################
  # diffusion
  ##################################################
  
  # calculate all the terms and substitute variables
  g = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term = hat2pow(private$diff.terms[[i]][[j+1]])
      new.term = do.call(substitute, list(term, subsList))
      g = c(g, sprintf("g(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  code = sprintf("Eigen::MatrixXd g(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd g(%s,%s);
                 %s
                 return g;
                 }",private$number.of.states, private$number.of.diffusions, paste(g,collapse=""))
  
  private$rcpp.function.strings$g = code
  
  ##################################################
  # observation
  ##################################################
  
  # calculate all the terms and substitute variables
  h = sapply(seq_along(private$obs.names), 
             function(i){
               term = hat2pow(private$obs.eqs.trans[[i]]$rhs)
               new.term = do.call(substitute, list(term, subsList))
               sprintf("h(%s) = %s;",i-1, deparse1(new.term))
             }) 
  code = sprintf("Eigen::VectorXd h(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd h(%s);
                 %s
                 return h;
                 }",private$number.of.observations, paste(h,collapse=""))
  
  private$rcpp.function.strings$h = code
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dhdx = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = hat2pow(private$diff.terms.obs[[i]][[j]])
      new.term = do.call(substitute, list(term, subsList))
      dhdx = c(dhdx, sprintf("dhdx(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  code = sprintf("Eigen::MatrixXd dhdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dhdx(%s,%s);
                 %s
                 return dhdx;
                 }", private$number.of.observations, private$number.of.states, paste(dhdx,collapse=""))
  
  private$rcpp.function.strings$dhdx = code
  
  
  ##################################################
  # observation variance
  ##################################################
  
  hvar = lapply(seq_along(private$obs.var.trans), 
                function(i) {
                  term = hat2pow(private$obs.var.trans[[i]]$rhs)
                  new.term = do.call(substitute, list(term, subsList))
                  sprintf("hvar(%s,%s) = %s;", i-1, i-1, deparse1(new.term))
                })
  code = sprintf("Eigen::MatrixXd hvar(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd hvar(%s,%s);
                 hvar.setZero();
                 %s
                 return hvar;
                 }", private$number.of.observations, private$number.of.observations, paste(hvar,collapse=""))
  
  private$rcpp.function.strings$hvar = code
  
  
}
