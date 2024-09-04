########## UNSCENTED KALMAN FILTER FUNCTIONS ###############
############################################################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

write_laplace_functions = function(self, private){
  
  txt = c()
  
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
  txt = c(txt,newtxt)
  
  ##################################################
  # diffusion
  ##################################################
  
  # calculate all the terms and substitute variables
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
  txt = c(txt, newtxt)
  
  ##################################################
  # observation
  ##################################################
  
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
  txt = c(txt, newtxt)
  
  ##################################################
  # observation variance
  ##################################################
  
  # obs.var = lapply(seq_along(private$obs.var.trans), 
  #                  function(i) {
  #                    term = hat2pow(private$obs.var.trans[[i]]$rhs)
  #                    new.term = do.call(substitute, list(term, subsList))
  #                    sprintf("hvar__(%s) = %s;", i-1, deparse1(new.term))
  #                  })
  # newtxt = "\n//////////// observation variance matrix function ///////////
  # template<class Type>
  # vector<Type> hsqrt__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
  #   vector<Type> hvar__(%s);
  #   %s
  #   return sqrt(hvar__);
  # }"
  # newtxt = sprintf(newtxt, private$number.of.observations, paste(obs.var,collapse="\n\t\t"))
  # txt = c(txt, newtxt)
  
  ##################################################
  # observation variance
  ##################################################
  
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
  vector<Type> hsqrt__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> hvar__(%s);
    hvar__.setZero();
    %s
    return sqrt(hvar__);
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(hvar,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


write_laplace_estimate = function(self, private){
  
  txt = c()
  
  # Observation Vectors
  txt = c(txt, "\n//// observations ////")
  txt = c(txt, "DATA_MATRIX(obsMat);")
  
  # Input Vectors
  txt = c(txt, "\n//// inputs ////")
  txt = c(txt, "DATA_MATRIX(inputMat);")
  
  # State Random Effects
  txt = c(txt, "\n//// initial state ////")
  txt = c(txt, "PARAMETER_MATRIX(stateMat);")
  
  # Parameters
  txt = c(txt, "\n//// parameters ////")
  for(i in 1:length(private$parameters)){
    txt = c(txt, sprintf("PARAMETER(%s);",private$parameter.names[i]))
  }
  # Time-step
  txt = c(txt, "DATA_VECTOR(ode_timestep_size);")
  txt = c(txt, "DATA_IVECTOR(ode_timesteps);")
  txt = c(txt, "DATA_IVECTOR(ode_timesteps_cumsum);")
  
  # Iobs
  txt = c(txt, "\n//// iobs vectors ////")
  for (i in 1:private$number.of.observations) {
    nam = paste("iobs_",private$obs.names[i],sep="")
    txt = c(txt, sprintf("\t DATA_IVECTOR(%s);",nam))
  }
  
  # System size
  txt = c(txt, "\n//// system size ////")
  txt = c( txt, "DATA_INTEGER(number_of_state_eqs);")
  txt = c( txt, "DATA_INTEGER(number_of_obs_eqs);")
  txt = c( txt, "DATA_INTEGER(number_of_diffusions);")
  txt = c( txt, "int tsize = inputMat.col(0).size();")
  txt = c( txt, "int idx;")
  
  # state, par, input, obs
  txt = c(txt, "\n//// state, par, input, obs vectors ////")
  txt = c(txt, sprintf("vector<Type> inputVec(%s), dinputVec(%s), obsVec(%s), parVec(%s);", 
                       private$number.of.inputs, 
                       private$number.of.inputs,
                       private$number.of.observations,
                       private$number.of.pars
  ))
  txt = c(txt, sprintf("parVec << %s;", paste(private$parameter.names,collapse=", ")))
  
  txt = c(txt, "vector<vector<Type>> xOut(tsize);")
  
  # Initiaze variables
  txt = c(txt, "\n //////////// initialize variables ///////////")
  txt = c(txt, "vector<Type> F__(number_of_state_eqs);")
  txt = c(txt, "vector<Type> stateVec(number_of_state_eqs);")
  txt = c(txt, "vector<Type> stateVec1(number_of_state_eqs);")
  txt = c(txt, "vector<Type> Z__(number_of_state_eqs);")
  txt = c(txt, "matrix<Type> G__(number_of_state_eqs, number_of_diffusions);")
  txt = c(txt, "matrix<Type> V__(number_of_state_eqs,number_of_state_eqs);")
  txt = c(txt, "matrix<Type> I__(number_of_state_eqs,number_of_state_eqs);")
  txt = c(txt, "I__.setIdentity();")
  txt = c(txt, "I__ *= 1e-8;")
  
  # Forward simulation and likelihood contribution from states
  txt = c(txt, "\n //////////// MAIN LOOP OVER TIME POINTS ///////////")
  txt = c(txt, "for(int i=0 ; i < tsize-1 ; i++){")
  txt = c(txt, "inputVec = inputMat.row(i);")
  txt = c(txt, "dinputVec = (inputMat.row(i+1) - inputMat.row(i))/ode_timesteps(i);")
  txt = c(txt, "idx = ode_timesteps_cumsum(i);")
  # 
  txt = c(txt, "for(int j=0 ; j < ode_timesteps(i) ; j++){")
  txt = c(txt, "stateVec = stateMat.row(idx);")
  txt = c(txt, "stateVec1 = stateMat.row(idx+1);")
  txt = c(txt, "F__ = f__(stateVec, parVec, inputVec);")
  txt = c(txt, "G__ = g__(stateVec, parVec, inputVec);")
  txt = c(txt, "Z__ = stateVec1 - stateVec - F__ * ode_timestep_size(i);")
  txt = c(txt, "V__ = (G__ * G__.transpose() + I__) * ode_timestep_size(i);")
  txt = c(txt, "nll__ += MVNORM(V__)(Z__);")
  txt = c(txt, "inputVec += dinputVec;")
  txt = c(txt, "idx += 1;")
  txt = c(txt, "}")
  txt = c(txt, "xOut(i) = stateVec;")
  txt = c(txt, "}")
  txt = c(txt, "REPORT(xOut);")
  
  
  # Calculate h(x) and hvar(x) for all time-points
  txt = c(txt, "\n//////////// DATA-UPDATE ///////////")
  txt = c(txt, "matrix<Type> obsMatMean__(tsize, number_of_obs_eqs);")
  txt = c(txt, "matrix<Type> obsMatStd__(tsize, number_of_obs_eqs);")
  # FIXME: This should just be over unique(iobs) in principle, not all rows of data...
  # since some might be NA...
  txt = c(txt, "for(int i=0 ; i < tsize ; i++){")
  txt = c(txt, "idx = ode_timesteps_cumsum(i);")
  txt = c(txt, "stateVec = stateMat.row(idx);")
  txt = c(txt, "inputVec = inputMat.row(i);")
  txt = c(txt, "obsMatMean__.row(i) = h__(stateVec, parVec, inputVec);")
  txt = c(txt, "obsMatStd__.row(i) = hsqrt__(stateVec, parVec, inputVec);")
  txt = c(txt, "}")
  # 
  
  # Data Likelihood Contribution
  txt = c(txt, "Type obsScalar, obsScalarMean, obsScalarStd;")
  txt = c(txt, "int idy;")
  for(i in 1:private$number.of.observations){
    nam = paste("iobs_",private$obs.names[i], sep="")
    txt = c(txt, sprintf("for(int i=0 ; i < %s.size() ; i++){", nam))
    txt = c(txt, sprintf("idy = %s(i);", nam))
    txt = c(txt, sprintf("obsScalar = obsMat(idy, %s);",i-1))
    txt = c(txt, sprintf("obsScalarMean = obsMatMean__(idy, %s);",i-1))
    txt = c(txt, sprintf("obsScalarStd = obsMatStd__(idy, %s);", i-1))
    # 
    txt = c(txt, "nll__ -= dnorm(obsScalar, obsScalarMean, obsScalarStd, true);")
    txt = c(txt, "}")
  }
  
  return(txt)
}


