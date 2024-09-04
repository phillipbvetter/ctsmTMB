#######################################################
# CONSTRUCT DRIFT, DIFF, OBS FUNCTIONS FOR RTMB
#######################################################

create_rtmb_function_strings = function(self, private)
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
  private$rtmb.function.strings$f = stringr::str_replace_all(f.function.text, 
                                                             pattern="F_ELEMENTS", 
                                                             replacement=paste(f.elements,collapse=","))
  
  ##################################################
  # drift jacobian
  ##################################################
  dfdx.elements = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term = Deriv::Deriv(private$diff.terms[[i]]$dt,x=private$state.names[j], cache.exp=F)
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
  
  return(invisible(NULL))
}

create_rtmb_laplace_string = function(self, private){
  
  
  # define the negative loglikelihood 
  laplace.text = paste('
  laplace.nll = function(p){
    # set negative log-likelihood
    nll = 0
    
    # small identity matrix
    small_identity = diag(1e-8, nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_STATES)
    
    # extract state random effects
    stateMat = cbind(STATE_NAMES)
    
    # fixed effects parameter vector
    parVec = cbind(FIXED_PARAMETERS)
    
    ###### TIME LOOP START #######
    for(i in 1:(nrow(obsMat)-1)) {
    
    # Define inputs and use first order input interpolation
    inputVec = RTMB::advector(inputMat[i,])
    dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
    
    ###### BETWEEN TIME POINTS LOOP START #######
    for(j in 1:ode_timesteps[i]){
    
    # grab current and next state
    x_now = stateMat[ode_cumsum_timesteps[i]+j,]
    x_next = stateMat[ode_cumsum_timesteps[i]+j+1,]
    
    # compute drift (vector) and diffusion (matrix)
    f = f__(x_now, parVec, inputVec)
    g = g__(x_now, parVec, inputVec)
    inputVec = inputVec + dinputVec
    
    # assume multivariate gauss distribution according to euler-step
    # and calculate the likelihood
    z = x_next - (x_now + f * ode_timestep_size[i])
    v = (g %*% t(g) + small_identity) * ode_timestep_size[i]
    nll = nll - RTMB::dmvnorm(z, Sigma=v, log=TRUE)
    }
    ###### BETWEEN TIME POINTS LOOP END #######
    }
    ###### TIME LOOP END #######
    
    obsMat = RTMB::OBS(obsMat)
    ###### DATA UPDATE START #######
    for(i in 1:NUMBER_OF_OBSERVATIONS){
    for(j in 1:length(iobs[[i]])){
    
    # Get index where observation is available
    k = iobs[[i]][j]
    
    # Get corresponding input, state and observation
    inputVec = inputMat[k,]
    stateVec = stateMat[ode_cumsum_timesteps[k]+1,]
    obsScalar = obsMat[k,i]
    
    # Observation equation and variance
    h_x = h__(stateVec, parVec, inputVec)[i]
    hvar_x = hvar__(stateVec, parVec, inputVec)[i]
    
    # likelihood contribution
    nll = nll - RTMB::dnorm(obsScalar, mean=h_x, sd=sqrt(hvar_x), log=TRUE)
    }}
    ###### DATA UPDATE END #######
    
    # return
    return(invisible(nll))
}')
  
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="NUMBER_OF_STATES", 
                                          replacement=deparse(as.numeric(private$number.of.states)))
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="NUMBER_OF_DIFFUSIONS", 
                                          replacement=deparse(as.numeric(private$number.of.diffusions)))
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="NUMBER_OF_OBSERVATIONS", 
                                          replacement=deparse(as.numeric(private$number.of.observations)))
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="STATE_NAMES", 
                                          replacement=paste("p$",private$state.names,sep="",collapse=","))
  laplace.text = stringr::str_replace_all(laplace.text, 
                                          pattern="FIXED_PARAMETERS", 
                                          replacement=paste("p$",private$parameter.names,sep="",collapse=","))
  
  # Store likelihood function
  private$rtmb.nll.strings$laplace = laplace.text
  
  return(invisible(NULL))
}

create_rtmb_ekf_string = function(self, private){
  
  ################################################
  ekf.text = paste('
  ekf.nll = function(p){
    ####### INITIALIZATION #######
    # storage variables
    xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
    
    # constants
    ns <- length(stateVec)
    halflog2pi = log(2*pi)/2
    I0 <- diag(NUMBER_OF_STATES)
    E0 <- diag(NUMBER_OF_OBSERVATIONS)
    
    # set negative log-likelihood
    nll = 0
    
    # fixed effects parameter vector
    parVec = c(FIXED_PARAMETERS)
    
    # Initial Stationary Conditions
    # inputVec = inputMat[1,]
    # F <- RTMB::MakeTape(function(x) f__(x, inputVec, parVec)^2, numeric(ns))
    # X0 <- F$newton(1:ns)
    # X1 <- X0(numeric(0))
    # stateVec <- X1
    
    # A <- dfdx__(stateVec, parVec, inputVec)
    # G <- g__(stateVec, parVec, inputVec)
    # Q <- G %*% t(G)
    # I <- diag(rep(1, nrow(A)))
    # P <- kron_left(A) + kron_right(A)
    # X <- -RTMB::solve(P, as.numeric(Q))
    # covMat <- X
    
    xPrior[[1]] <- stateVec
    pPrior[[1]] <- covMat
    
    ####### First Date Update ########
    obsVec = obsMat[1,]
    inputVec = inputMat[1,]
    obsVec_bool = !is.na(obsVec)
    
    ######## DATA UPDATE - IF ANY DATA AVAILABLE ######## 
    if(any(obsVec_bool)){
      s = sum(obsVec_bool)
      y = obsVec[obsVec_bool]
      E = E0[obsVec_bool,obsVec_bool,drop=FALSE]
      H = h__(stateVec, parVec, inputVec)
      C = E %*% dhdx__(stateVec, parVec, inputVec)
      e = y - E %*% H
      V0 = hvar__matrix(stateVec, parVec, inputVec)
      # V0 = RTMB::diag(hvar__(stateVec, parVec, inputVec))
      V = E %*% V0 %*% t(E)
      R = C %*% covMat %*% t(C) + V
      Ri = RTMB::solve(R)
      K = covMat %*% t(C) %*% Ri
      stateVec = stateVec + K %*% e
      covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
      # Store innovation and covariance
      Innovation[[1]] = e
      InnovationCovariance[[1]] = R
    }
    xPost[[1]] <- stateVec
    pPost[[1]] <- covMat
    
    # ###### TIME LOOP START #######
    for(i in 1:(nrow(obsMat)-1)){
    
    # # Define inputs and use first order input interpolation
    inputVec = inputMat[i,]
    dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]
    
    ###### TIME UPDATE - ODE SOLVER #######
    for(j in 1:ode_timesteps[i]){
    sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i], ode_solver)
    stateVec = sol[[1]]
    covMat = sol[[2]]
    inptuVec = inputVec + dinputVec
    }
    xPrior[[i+1]] = stateVec
    pPrior[[i+1]] = covMat
    
    ######## DATA UPDATE - KALMAN ALGORITHM ########
    obsVec = obsMat[i+1,]
    inputVec = inputMat[i+1,]
    obsVec_bool = !is.na(obsVec)
    
    ######## DATA UPDATE - IF ANY DATA AVAILABLE ######## 
    if(any(obsVec_bool)){
    s = sum(obsVec_bool)
    y = obsVec[obsVec_bool]
    E = E0[obsVec_bool,obsVec_bool,drop=FALSE]
    H = h__(stateVec, parVec, inputVec)
    C = E %*% dhdx__(stateVec, parVec, inputVec)
    e = y - E %*% H
    V0 = hvar__matrix(stateVec, parVec, inputVec)
    # V0 = RTMB::diag(hvar__(stateVec, parVec, inputVec))
    V = E %*% V0 %*% t(E)
    R = C %*% covMat %*% t(C) + V
    Ri = RTMB::solve(R)
    K = covMat %*% t(C) %*% Ri
    
    # Update State
    stateVec = stateVec + K %*% e
    covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
    nll = nll + 0.5 * logdet(R) + 0.5 * t(e) %*% Ri %*% e + halflog2pi * s
    
    # Store innovation and covariance
    Innovation[[i+1]] = e
    InnovationCovariance[[i+1]] = R
    }
    xPost[[i+1]] = stateVec
    pPost[[i+1]] = covMat
    }
    
    # ###### MAXIMUM A POSTERIOR #######
    
    # ###### REPORT #######
    RTMB::REPORT(Innovation)
    RTMB::REPORT(InnovationCovariance)
    RTMB::REPORT(xPrior)
    RTMB::REPORT(xPost)
    RTMB::REPORT(pPrior)
    RTMB::REPORT(pPost)
    
    # ###### RETURN #######
    return(nll)
}')
  
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="NUMBER_OF_STATES", 
                                      replacement=deparse(as.numeric(private$number.of.states)))
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="NUMBER_OF_DIFFUSIONS", 
                                      replacement=deparse(as.numeric(private$number.of.diffusions)))
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="NUMBER_OF_OBSERVATIONS", 
                                      replacement=deparse(as.numeric(private$number.of.observations)))
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="STATE_NAMES", 
                                      replacement=paste("p$",private$state.names,sep="",collapse=","))
  ekf.text = stringr::str_replace_all(ekf.text, 
                                      pattern="FIXED_PARAMETERS", 
                                      replacement=paste("p$",private$parameter.names,sep="",collapse=","))
  
  private$rtmb.nll.strings$ekf = ekf.text
  
  
}
