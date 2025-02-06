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
      term = ctsmTMB.Deriv(f=private$diff.terms[[i]]$dt, x=private$state.names[j])
      # term = Deriv::Deriv(private$diff.terms[[i]]$dt,x=private$state.names[j], cache.exp=F)
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
