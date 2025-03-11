##########################################################
##########################################################
##########################################################
# USER-FUNCTION CONSTRUCTION FOR RTMB
##########################################################
##########################################################
##########################################################

create_rtmb_function_strings = function(self, private)
{
  
  # Create substitution translation list
  # obsList = lapply(seq_along(private$obs.names), function(id) substitute(obsVec[i],list(i=as.numeric(id))))
  # parList = lapply(seq_along(private$parameter.names), function(id) substitute(parVec[i],list(i=as.numeric(id))))
  # stateList = lapply(seq_along(private$state.names), function(id) substitute(stateVec[i],list(i=as.numeric(id))))
  # inputList = lapply(seq_along(private$input.names), function(id) substitute(inputVec[i],list(i=as.numeric(id))))
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
  f.elements = sapply( seq_along(private$state.names),
                       function(i){
                         drift.term = private$diff.terms[[i]]$dt
                         deparse1(do.call(substitute, list(drift.term, subsList)))
                       })
  
  f.function.text = paste('
  f__ = function(stateVec, parVec, inputVec){
    ans = c(F_ELEMENTS)
    return(ans)
  }')
  private$rtmb.function.strings$f = stringr::str_replace_all(f.function.text,
                                                             pattern="F_ELEMENTS",
                                                             replacement=paste(f.elements,collapse=","))
  
  ##################################################
  # drift jacobian
  ##################################################
  dfdx.elements = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term = ctsmTMB.Deriv(f=private$diff.terms[[i]]$dt, x=private$state.names[j])
      dfdx.elements = c(dfdx.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  dfdx.function.text = paste('
  dfdx__ = function(stateVec, parVec, inputVec){
    ans = RTMB::matrix(c(DFDX_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_STATES, byrow=T)
    return(ans)
  }')
  
  dfdx.function.text = stringr::str_replace_all(dfdx.function.text,
                                                pattern="NUMBER_OF_STATES",
                                                replacement=deparse(as.numeric(private$number.of.states)))
  
  dfdx.function.text = stringr::str_replace_all(dfdx.function.text,
                                                pattern="DFDX_ELEMENTS",
                                                replacement=paste(dfdx.elements,collapse=","))
  
  private$rtmb.function.strings$dfdx = dfdx.function.text
  
  ##################################################
  # diffusion
  ##################################################
  g.elements = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term = private$diff.terms[[i]][[j+1]]
      g.elements = c(g.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  g.function.text = paste('
  g__ = function(stateVec, parVec, inputVec){
    ans = RTMB::matrix(c(G_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_DIFFUSIONS, byrow=T)
    return(ans)
  }')
  
  g.function.text = stringr::str_replace_all(g.function.text,
                                             pattern="NUMBER_OF_STATES",
                                             replacement=deparse(as.numeric(private$number.of.states)))
  
  g.function.text = stringr::str_replace_all(g.function.text,
                                             pattern="NUMBER_OF_DIFFUSIONS",
                                             replacement=deparse(as.numeric(private$number.of.diffusions)))
  
  g.function.text = stringr::str_replace_all(g.function.text,
                                             pattern="G_ELEMENTS",
                                             replacement=paste(g.elements,collapse=","))
  
  private$rtmb.function.strings$g = g.function.text
  
  
  ##################################################
  # observation
  ##################################################
  h.elements = sapply(seq_along(private$obs.names),
                      function(i){
                        term = private$obs.eqs.trans[[i]]$rhs
                        deparse1(do.call(substitute, list(term, subsList)))
                      })
  
  h.function.text = paste('
  h__ = function(stateVec, parVec, inputVec){
    ans = c(H_ELEMENTS)
    return(ans)
  }')
  
  private$rtmb.function.strings$h = stringr::str_replace_all(h.function.text,
                                                             pattern="H_ELEMENTS",
                                                             replacement=paste(h.elements,collapse=","))
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dhdx.elements = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = private$diff.terms.obs[[i]][[j]]
      dhdx.elements = c(dhdx.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  dhdx.function.text = paste('
  dhdx__ = function(stateVec, parVec, inputVec){
    ans = RTMB::matrix(c(DHDX_ELEMENTS),nrow=NUMBER_OF_OBSERVATIONS, ncol=NUMBER_OF_STATES, byrow=T)
    return(ans)
  }')
  
  dhdx.function.text = stringr::str_replace_all(dhdx.function.text,
                                                pattern="NUMBER_OF_STATES",
                                                replacement=deparse(as.numeric(private$number.of.states)))
  
  dhdx.function.text= stringr::str_replace_all(dhdx.function.text,
                                               pattern="NUMBER_OF_OBSERVATIONS",
                                               replacement=deparse(as.numeric(private$number.of.observations)))
  
  dhdx.function.text = stringr::str_replace_all(dhdx.function.text,
                                                pattern="DHDX_ELEMENTS",
                                                replacement=paste(dhdx.elements,collapse=","))
  
  private$rtmb.function.strings$dhdx = dhdx.function.text
  
  ##################################################
  # observation variance
  ##################################################
  
  ### VECTOR FORM ### (for laplace)
  
  hvar.elements = sapply(seq_along(private$obs.var.trans),
                         function(i) {
                           term = private$obs.var.trans[[i]]$rhs
                           deparse1(do.call(substitute, list(term, subsList)))
                         })
  
  hvar.function.text = paste('
  hvar__ = function(stateVec, parVec, inputVec){
    ans = c(HVAR_ELEMENTS)
    return(ans)
  }')
  
  private$rtmb.function.strings$hvar = stringr::str_replace_all(hvar.function.text,
                                                                pattern="HVAR_ELEMENTS",
                                                                replacement=paste(hvar.elements,collapse=","))
  
  ### MATRIX FORM ### (for kalman)
  
  hvar.function.text = paste('
  hvar__matrix = function(stateVec, parVec, inputVec){
    ans = RTMB::matrix(c(HVAR_ELEMENTS),nrow=NUMBER_OF_OBSERVATIONS, ncol=NUMBER_OF_OBSERVATIONS, byrow=T)
    return(ans)
  }')
  
  hvar.function.text = stringr::str_replace_all(hvar.function.text,
                                                pattern="NUMBER_OF_OBSERVATIONS",
                                                replacement=deparse(as.numeric(private$number.of.observations)))
  
  hvar.function.text = stringr::str_replace_all(hvar.function.text,
                                                pattern="HVAR_ELEMENTS",
                                                replacement=paste(hvar.elements,collapse=","))
  
  private$rtmb.function.strings$hvar__matrix = hvar.function.text
  
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
  constant.terms <- unname(sapply(
    sapply(private$sys.eqs.trans, function(x) Deriv::Simplify(do.call(substitute, list(x$diff.dt, zero.list)))),
    function(x) deparse1(do.call(substitute, list(x, subsList)))
  ))
  # constant.terms <- as.vector(sapply(constant.terms, function(x) deparse1(do.call(substitute, list(x, subsList)))))
  
  dfdu.elements <- c()
  # 2. Now we find input-terms by differentiation
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$input.names)){
      if(j==1) dfdu.elements <- c(dfdu.elements, constant.terms[i])
      term <- ctsmTMB.Deriv(f = private$sys.eqs.trans[[i]]$diff.dt, x=private$input.names[j])
      dfdu.elements = c(dfdu.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  dfdu.function.text = paste('
  dfdu__ = function(stateVec, parVec, inputVec){
    ans = RTMB::matrix(c(DFDU_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_INPUTS_PLUS_ONE, byrow=T)
    return(ans)
  }')
  
  dfdu.function.text = stringr::str_replace_all(dfdu.function.text,
                                                pattern="NUMBER_OF_STATES",
                                                replacement=deparse(as.numeric(private$number.of.states)))
  
  dfdu.function.text = stringr::str_replace_all(dfdu.function.text,
                                                pattern="NUMBER_OF_INPUTS_PLUS_ONE",
                                                replacement=deparse(as.numeric(private$number.of.inputs)+1))
  
  dfdu.function.text = stringr::str_replace_all(dfdu.function.text,
                                                pattern="DFDU_ELEMENTS",
                                                replacement=paste(dfdu.elements,collapse=","))
  
  private$rtmb.function.strings$dfdu = dfdu.function.text
  
  return(invisible(NULL))
}

##########################################################
##########################################################
##########################################################
# USER-FUNCTION CONSTRUCTION FOR RTMB
##########################################################
##########################################################
##########################################################

create_rtmb_function_strings_new = function(self, private)
{
  
  # Create substitution translation list
  # Assignment with [[]] instead of [] reduces RTMB::MakeADFun compilation time
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
  for(i in seq_along(private$state.names)){
    term <- deparse1(do.call(substitute, list(private$diff.terms[[i]]$dt, subsList)))
    
    # skip if zero
    if(term=="0") next
    
    # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
    # so wrap in RTMB::AD
    out <- suppressWarnings(as.numeric(term))
    if(!is.na(out)){
      f.elements <- c(f.elements, sprintf("f_vec[[%s]] <- RTMB::AD(%s);",i,term))
    } else{
      f.elements <- c(f.elements, sprintf("f_vec[[%s]] <- %s;",i,term))
    }
  }
  ftext <- 'f__ = function(stateVec, parVec, inputVec){
  %s
  return(f_vec)
  }'
  ftext2 <- sprintf(ftext, paste(f.elements,collapse=" "))
  private$rtmb.function.strings.indexed$f = ftext2
  
  ##################################################
  # drift jacobian
  ##################################################
  dfdx.elements = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term <- deparse1(do.call(substitute, list(private$diff.terms.drift[[i]][[j]], subsList)))
      
      # skip if zero
      if(term=="0") next
      
      # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
      # so wrap in RTMB::AD
      out <- suppressWarnings(as.numeric(term))
      if(!is.na(out)){
        dfdx.elements <- c(dfdx.elements, sprintf("dfdx_mat[[%s,%s]] <- RTMB::AD(%s);",i,j,term))
      } else{
        dfdx.elements <- c(dfdx.elements, sprintf("dfdx_mat[[%s,%s]] <- %s;",i,j,term))
      }
    }
  }
  ftext <- 'dfdx__ = function(stateVec, parVec, inputVec){
  %s
  return(dfdx_mat)
  }'
  ftext2 <- sprintf(ftext, paste(dfdx.elements,collapse=" "))
  private$rtmb.function.strings.indexed$dfdx = ftext2
  
  ##################################################
  # diffusion
  ##################################################
  g.elements = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term <- deparse1(do.call(substitute, list(private$diff.terms[[i]][[j+1]], subsList)))
      
      # skip if zero
      if(term=="0") next
      
      # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
      # so wrap in RTMB::AD
      out <- suppressWarnings(as.numeric(term))
      if(!is.na(out)){
        g.elements <- c(g.elements, sprintf("g_mat[[%s,%s]] <- RTMB::AD(%s);",i,j, term))
      } else{
        g.elements <- c(g.elements, sprintf("g_mat[[%s,%s]] <- %s;",i,j, term))
      }
    }
  }
  ftext <- 'g__ = function(stateVec, parVec, inputVec){
  %s
  return(g_mat)
  }'
  ftext2 <- sprintf(ftext, paste(g.elements,collapse=" "))
  private$rtmb.function.strings.indexed$g = ftext2
  
  ##################################################
  # observation
  ##################################################
  h.elements <- c()
  for(i in seq_along(private$obs.names)){
    term <- deparse1(do.call(substitute, list(private$obs.eqs.trans[[i]]$rhs, subsList)))
    
    # skip if zero
    if(term=="0") next
    
    # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
    # so wrap in RTMB::AD
    out <- suppressWarnings(as.numeric(term))
    if(!is.na(out)){
      # hack
      h.elements <- c(h.elements, sprintf("h_vec[[%s]] <- RTMB::AD(%s);",i, out))
    } else{
      h.elements <- c(h.elements, sprintf("h_vec[[%s]] <- %s;",i, term))
    }
  }
  ftext <- 'h__ = function(stateVec, parVec, inputVec){
  %s
  return(h_vec)
  }'
  ftext2 <- sprintf(ftext, paste(h.elements,collapse=" "))
  private$rtmb.function.strings.indexed$h = ftext2
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dhdx.elements = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term <- deparse1(do.call(substitute, list(private$diff.terms.obs[[i]][[j]], subsList)))
      
      # skip if zero
      if(term=="0") next
      
      # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
      # so wrap in RTMB::AD
      out <- suppressWarnings(as.numeric(term))
      if(!is.na(out)){
        dhdx.elements <- c(dhdx.elements, sprintf("dhdx_mat[[%s,%s]] <- RTMB::AD(%s);",i,j,term))
      } else{
        dhdx.elements <- c(dhdx.elements, sprintf("dhdx_mat[[%s,%s]] <- %s;",i,j,term))
      }
    }
  }
  ftext <- 'dhdx__ = function(stateVec, parVec, inputVec){
  %s
  return(dhdx_mat)
  }'
  ftext2 <- sprintf(ftext, paste(dhdx.elements,collapse=" "))
  private$rtmb.function.strings.indexed$dhdx = ftext2
  
  ##################################################
  # observation variance
  ##################################################
  
  ### VECTOR FORM ### (for laplace)
  
  # hvar.elements = sapply(seq_along(private$obs.var.trans),
  #                        function(i) {
  #                          term = private$obs.var.trans[[i]]$rhs
  #                          deparse1(do.call(substitute, list(term, subsList)))
  #                        })
  # 
  # hvar.function.text = paste('
  # hvar__ = function(stateVec, parVec, inputVec){
  #   ans = c(HVAR_ELEMENTS)
  #   return(ans)
  # }')
  # 
  # private$rtmb.function.strings.indexed$hvar = stringr::str_replace_all(hvar.function.text,
  #                                                                       pattern="HVAR_ELEMENTS",
  #                                                                       replacement=paste(hvar.elements,collapse=","))
  
  hvar.elements <- c()
  for(i in seq_along(private$obs.names)){
    term <- deparse1(do.call(substitute, list(private$obs.var.trans[[i]]$rhs, subsList)))
    
    # skip if zero
    if(term=="0") next
    
    # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
    # so wrap in RTMB::AD
    out <- suppressWarnings(as.numeric(term))
    if(!is.na(out)){
      # hack
      hvar.elements <- c(h.elements, sprintf("hvar_vec[[%s]] <- RTMB::AD(%s);",i, out))
    } else{
      hvar.elements <- c(h.elements, sprintf("hvar_vec[[%s]] <- %s;",i, term))
    }
  }
  ftext <- 'hvar__ = function(stateVec, parVec, inputVec){
  %s
  return(hvar_vec)
  }'
  ftext2 <- sprintf(ftext, paste(hvar.elements,collapse=" "))
  private$rtmb.function.strings.indexed$hvar = ftext2
  
  ### MATRIX FORM ### (for kalman)
  hvar.elements <- c()
  for(i in seq_along(private$obs.names)){
    term <- deparse1(do.call(substitute, list(private$obs.var.trans[[i]]$rhs, subsList)))
    
    # skip if zero
    if(term=="0") next
    
    # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
    # so wrap in RTMB::AD
    out <- suppressWarnings(as.numeric(term))
    if(!is.na(out)){
      hvar.elements <- c(hvar.elements, sprintf("hvar_mat[[%s,%s]] <- RTMB::AD(%s);",i,i,term))
    } else{
      hvar.elements <- c(hvar.elements, sprintf("hvar_mat[[%s,%s]] <- %s;",i,i,term))
    }
  }
  ftext <- 'hvar__matrix = function(stateVec, parVec, inputVec){
  %s
  return(hvar_mat)
  }'
  ftext2 <- sprintf(ftext, paste(hvar.elements,collapse=" "))
  private$rtmb.function.strings.indexed$hvar__matrix = ftext2
  
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
  constant.terms <- unname(sapply(
    sapply(private$sys.eqs.trans, function(x) Deriv::Simplify(do.call(substitute, list(x$diff.dt, zero.list)))),
    function(x) deparse1(do.call(substitute, list(x, subsList)))
  ))
  
  dfdu.elements <- c()
  # 2. Now we find input-terms by differentiation
  for(i in seq_along(private$state.names)){
    # for(j in seq_along(private$input.names)){
    for(j in 1:(private$number.of.inputs+1)){
      
      # if(j==1) dfdu.elements <- c(dfdu.elements, constant.terms[i])
      # term <- ctsmTMB.Deriv(f = private$sys.eqs.trans[[i]]$diff.dt, x=private$input.names[j])
      
      if(j==1){
        term <- constant.terms[i]
      } else {
        tempterm <- ctsmTMB.Deriv(f=private$sys.eqs.trans[[i]]$diff.dt, x=private$input.names[j-1])
        term <- deparse1(do.call(substitute, list(tempterm, subsList)))
      }
      
      # skip if zero
      if(term=="0") next
      
      # check for direct numerical assignments e.g. stateVec[[1]] <- 5 which throws lost AD error
      # so wrap in RTMB::AD
      out <- suppressWarnings(as.numeric(term))
      if(!is.na(out)){
        dfdu.elements <- c(dfdu.elements, sprintf("dfdu_mat[[%s,%s]] <- RTMB::AD(%s);",i,j,term))
      } else{
        dfdu.elements <- c(dfdu.elements, sprintf("dfdu_mat[[%s,%s]] <- %s;",i,j,term))
      }
    }
  }
  ftext <- 'dfdu__ = function(stateVec, parVec, inputVec){
  %s
  return(dfdu_mat)
  }'
  ftext2 <- sprintf(ftext, paste(dfdu.elements,collapse=" "))
  private$rtmb.function.strings.indexed$dfdu = ftext2
  
  return(invisible(NULL))
}

##########################################################
##########################################################
##########################################################
# USER-FUNCTION CONSTRUCTION FOR R-IMPLEMENTATION OF EKF
##########################################################
##########################################################
##########################################################

create_rekf_function_strings = function(self, private)
{
  
  # Create substitution translation list
  obsList = lapply(seq_along(private$obs.names), function(id) substitute(obsVec[i],list(i=as.numeric(id))))
  parList = lapply(seq_along(private$parameter.names), function(id) substitute(parVec[i],list(i=as.numeric(id))))
  stateList = lapply(seq_along(private$state.names), function(id) substitute(stateVec[i],list(i=as.numeric(id))))
  inputList = lapply(seq_along(private$input.names), function(id) substitute(inputVec[i],list(i=as.numeric(id))))
  names(obsList) = private$obs.names
  names(parList) = private$parameter.names
  names(stateList) = private$state.names
  names(inputList) = private$input.names
  subsList = c(obsList, parList, stateList, inputList)
  
  ##################################################
  # drift
  ##################################################
  f.elements = sapply( seq_along(private$state.names),
                       function(i){
                         drift.term = private$diff.terms[[i]]$dt
                         deparse1(do.call(substitute, list(drift.term, subsList)))
                       })
  
  f.function.text = paste('
  f__ = function(stateVec, parVec, inputVec){
    ans = c(F_ELEMENTS)
    return(ans)
  }')
  private$rekf.function.strings$f = stringr::str_replace_all(f.function.text,
                                                             pattern="F_ELEMENTS",
                                                             replacement=paste(f.elements,collapse=","))
  
  ##################################################
  # drift jacobian
  ##################################################
  dfdx.elements = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term = ctsmTMB.Deriv(f=private$diff.terms[[i]]$dt, x=private$state.names[j])
      dfdx.elements = c(dfdx.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  dfdx.function.text = paste('
  dfdx__ = function(stateVec, parVec, inputVec){
    ans = matrix(c(DFDX_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_STATES, byrow=T)
    return(ans)
  }')
  
  dfdx.function.text = stringr::str_replace_all(dfdx.function.text,
                                                pattern="NUMBER_OF_STATES",
                                                replacement=deparse(as.numeric(private$number.of.states)))
  
  dfdx.function.text = stringr::str_replace_all(dfdx.function.text,
                                                pattern="DFDX_ELEMENTS",
                                                replacement=paste(dfdx.elements,collapse=","))
  
  private$rekf.function.strings$dfdx = dfdx.function.text
  
  ##################################################
  # diffusion
  ##################################################
  g.elements = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term = private$diff.terms[[i]][[j+1]]
      g.elements = c(g.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  g.function.text = paste('
  g__ = function(stateVec, parVec, inputVec){
    ans = matrix(c(G_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_DIFFUSIONS, byrow=T)
    return(ans)
  }')
  
  g.function.text = stringr::str_replace_all(g.function.text,
                                             pattern="NUMBER_OF_STATES",
                                             replacement=deparse(as.numeric(private$number.of.states)))
  
  g.function.text = stringr::str_replace_all(g.function.text,
                                             pattern="NUMBER_OF_DIFFUSIONS",
                                             replacement=deparse(as.numeric(private$number.of.diffusions)))
  
  g.function.text = stringr::str_replace_all(g.function.text,
                                             pattern="G_ELEMENTS",
                                             replacement=paste(g.elements,collapse=","))
  
  private$rekf.function.strings$g = g.function.text
  
  
  ##################################################
  # observation
  ##################################################
  h.elements = sapply(seq_along(private$obs.names),
                      function(i){
                        term = private$obs.eqs.trans[[i]]$rhs
                        deparse1(do.call(substitute, list(term, subsList)))
                      })
  
  h.function.text = paste('
  h__ = function(stateVec, parVec, inputVec){
    ans = c(H_ELEMENTS)
    return(ans)
  }')
  
  private$rekf.function.strings$h = stringr::str_replace_all(h.function.text,
                                                             pattern="H_ELEMENTS",
                                                             replacement=paste(h.elements,collapse=","))
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dhdx.elements = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = private$diff.terms.obs[[i]][[j]]
      dhdx.elements = c(dhdx.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  dhdx.function.text = paste('
  dhdx__ = function(stateVec, parVec, inputVec){
    ans = matrix(c(DHDX_ELEMENTS),nrow=NUMBER_OF_OBSERVATIONS, ncol=NUMBER_OF_STATES, byrow=T)
    return(ans)
  }')
  
  dhdx.function.text = stringr::str_replace_all(dhdx.function.text,
                                                pattern="NUMBER_OF_STATES",
                                                replacement=deparse(as.numeric(private$number.of.states)))
  
  dhdx.function.text= stringr::str_replace_all(dhdx.function.text,
                                               pattern="NUMBER_OF_OBSERVATIONS",
                                               replacement=deparse(as.numeric(private$number.of.observations)))
  
  dhdx.function.text = stringr::str_replace_all(dhdx.function.text,
                                                pattern="DHDX_ELEMENTS",
                                                replacement=paste(dhdx.elements,collapse=","))
  
  private$rekf.function.strings$dhdx = dhdx.function.text
  
  ##################################################
  # observation variance
  ##################################################
  
  ### VECTOR FORM ### (for laplace)
  
  hvar.elements = sapply(seq_along(private$obs.var.trans),
                         function(i) {
                           term = private$obs.var.trans[[i]]$rhs
                           deparse1(do.call(substitute, list(term, subsList)))
                         })
  
  hvar.function.text = paste('
  hvar__ = function(stateVec, parVec, inputVec){
    ans = c(HVAR_ELEMENTS)
    return(ans)
  }')
  
  private$rekf.function.strings$hvar = stringr::str_replace_all(hvar.function.text,
                                                                pattern="HVAR_ELEMENTS",
                                                                replacement=paste(hvar.elements,collapse=","))
  
  ### MATRIX FORM ### (for kalman)
  
  hvar.function.text = paste('
  hvar__matrix = function(stateVec, parVec, inputVec){
    ans = matrix(c(HVAR_ELEMENTS),nrow=NUMBER_OF_OBSERVATIONS, ncol=NUMBER_OF_OBSERVATIONS, byrow=T)
    return(ans)
  }')
  
  hvar.function.text = stringr::str_replace_all(hvar.function.text,
                                                pattern="NUMBER_OF_OBSERVATIONS",
                                                replacement=deparse(as.numeric(private$number.of.observations)))
  
  hvar.function.text = stringr::str_replace_all(hvar.function.text,
                                                pattern="HVAR_ELEMENTS",
                                                replacement=paste(hvar.elements,collapse=","))
  
  private$rekf.function.strings$hvar__matrix = hvar.function.text
  
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
  constant.terms <- unname(sapply(
    sapply(private$sys.eqs.trans, function(x) Deriv::Simplify(do.call(substitute, list(x$diff.dt, zero.list)))),
    function(x) deparse1(do.call(substitute, list(x, subsList)))
  ))
  
  dfdu.elements <- c()
  # 2. Now we find input-terms by differentiation
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$input.names)){
      if(j==1) dfdu.elements <- c(dfdu.elements, constant.terms[i])
      term <- ctsmTMB.Deriv(f = private$sys.eqs.trans[[i]]$diff.dt, x=private$input.names[j])
      dfdu.elements = c(dfdu.elements, deparse1(do.call(substitute, list(term, subsList))))
    }
  }
  
  dfdu.function.text = paste('
  dfdu__ = function(stateVec, parVec, inputVec){
    ans = matrix(c(DFDU_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_INPUTS_PLUS_ONE, byrow=T)
    return(ans)
  }')
  
  dfdu.function.text = stringr::str_replace_all(dfdu.function.text,
                                                pattern="NUMBER_OF_STATES",
                                                replacement=deparse(as.numeric(private$number.of.states)))
  
  dfdu.function.text = stringr::str_replace_all(dfdu.function.text,
                                                pattern="NUMBER_OF_INPUTS_PLUS_ONE",
                                                replacement=deparse(as.numeric(private$number.of.inputs)+1))
  
  dfdu.function.text = stringr::str_replace_all(dfdu.function.text,
                                                pattern="DFDU_ELEMENTS",
                                                replacement=paste(dfdu.elements,collapse=","))
  
  private$rekf.function.strings$dfdu = dfdu.function.text
  
  return(invisible(NULL))
}


##########################################################
##########################################################
##########################################################
# USER-FUNCTION CONSTRUCTION FOR RCPP EKF
##########################################################
##########################################################
##########################################################

create_rcpp_function_strings = function(self, private){
  
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
  code = sprintf("SEXP f(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd f(%s); 
                 %s
                 return Rcpp::wrap(f);
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
  
  code = sprintf("SEXP dfdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dfdx(%s,%s);
                 %s
                 return Rcpp::wrap(dfdx);
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
  
  code = sprintf("SEXP g(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd g(%s,%s);
                 %s
                 return Rcpp::wrap(g);
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
  
  code = sprintf("SEXP h(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd h(%s);
                 %s
                 return Rcpp::wrap(h);
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
  
  code = sprintf("SEXP dhdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dhdx(%s,%s);
                 %s
                 return Rcpp::wrap(dhdx);
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
  
  code = sprintf("SEXP hvar(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd hvar(%s,%s);
                 hvar.setZero();
                 %s
                 return Rcpp::wrap(hvar);
                 }", private$number.of.observations, private$number.of.observations, paste(hvar,collapse=""))
  
  private$rcpp.function.strings$hvar = code
  
  
}
