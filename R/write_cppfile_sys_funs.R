

get_substitution_list = function(self, private){
  
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
  
  return(subsList)
}

##################################################
# drift
##################################################

write_f = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  f <- c()
  for(i in seq_along(private$state.names)){
    drift.term <- Deriv::Simplify(private$diff.terms[[i]]$dt)
    if(!(drift.term==0)){
      drift.term = hat2pow(private$diff.terms[[i]]$dt)
      new.drift.term = do.call(substitute, list(drift.term, subsList))
      f <- c(f, sprintf("f__(%i) = %s;",i-1, deparse1(new.drift.term)))
    }
  }
  
  newtxt = "\n//////////// drift function //////////
  template<class Type>
  vector<Type> f__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> f__(%s);
    f__.setZero();
    %s
    return f__;
  }"
  
  newtxt = sprintf(newtxt, private$number.of.states, paste(f,collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# drift jacobian
##################################################

write_jac_f = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  jac.f = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term <- Deriv::Simplify(private$diff.terms.drift[[i]][[j]])
      # term <- Deriv::Simplify(ctsmTMB.Deriv(f=private$diff.terms[[i]]$dt, x=private$state.names[j]))
      # term <- Deriv::Simplify(Deriv::Deriv(private$diff.terms[[i]]$dt,x=private$state.names[j], cache.exp=F))
      if(!(term==0)){
        term = hat2pow(term)
        new.term = do.call(substitute, list(term, subsList))
        jac.f = c(jac.f, sprintf("dfdx__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
      }
    }
  }
  
  newtxt = "\n//////////// jacobian of drift function ///////////
  template<class Type>
  matrix<Type> dfdx__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> dfdx__(%s, %s);
    dfdx__.setZero();
    %s
    return dfdx__;
  }"
  
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.states, paste(jac.f, collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# diffusion
##################################################

write_g = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  g = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term <- Deriv::Simplify(private$diff.terms[[i]][[j+1]])
      if(!(term==0)){
        term = hat2pow(term)
        new.term = do.call(substitute, list(term, subsList))
        g = c(g, sprintf("g__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
      }
    }
  }
  
  newtxt = "\n//////////// diffusion function ///////////
  template<class Type>
  matrix<Type> g__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> g__(%s, %s);
    g__.setZero();
    %s
    return g__;
  }"
  
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.diffusions, paste(g,collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# observation
##################################################

write_h = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  h <- c()
  # calculate all the terms and substitute variables
  for(i in seq_along(private$obs.names)){
    term <- Deriv::Simplify(private$obs.eqs.trans[[i]]$rhs)
    if(term!=0){
      term = hat2pow(term)
      new.term = do.call(substitute, list(term, subsList))
      h <- c(h, sprintf("h__(%s) = %s;",i-1, deparse1(new.term)))
    }
  }
  
  newtxt = "\n//////////// observation function ///////////
  template<class Type>
  vector<Type> h__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> h__(%s);
    h__.setZero();
    %s
    return h__;
  }"
  
  newtxt = sprintf(newtxt, private$number.of.observations, paste(h,collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# observation jacobian
##################################################

write_jac_h = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  jac.h = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = Deriv::Simplify(private$diff.terms.obs[[i]][[j]])
      if(term!=0){
        term = hat2pow(term)
        new.term = do.call(substitute, list(term, subsList))
        jac.h = c(jac.h, sprintf("dhdx__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
      }
    }
  }
  
  newtxt = "\n//////////// jacobian of obs function ///////////
  template<class Type>
  matrix<Type> dhdx__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> dhdx__(%s, %s);
    dhdx__.setZero();
    %s
    return dhdx__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, private$number.of.states, paste(jac.h,collapse="\n\t\t"))
  
  return(newtxt)
}

##################################################
# observation variance
##################################################

write_h_var = function(self, private){
  
  subsList <- get_substitution_list(self, private)
  
  hvar <- c()
  for(i in seq_along(private$obs.var.trans)){
    term <- Deriv::Simplify(private$obs.var.trans[[i]]$rhs)
    if(term!=0){
      term <- hat2pow(term)
      new.term = do.call(substitute, list(term, subsList))
      hvar <- c(hvar, sprintf("hvar__(%s) = %s;", i-1, deparse1(new.term)))
    }
  }
  
  newtxt = "\n//////////// observation variance matrix function ///////////
  template<class Type>
  vector<Type> hvar__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> hvar__(%s);
    hvar__.setZero();
    %s
    return hvar__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(hvar,collapse="\n\t\t"))
  
  return(newtxt)
}
