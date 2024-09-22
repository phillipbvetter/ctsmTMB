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
  
  ##################################################
  # DRIFT, DIFFUSION, OBSERVATION FUNCTIONS
  ##################################################
  txt = c(txt, write_f(self, private))
  txt = c(txt, write_g(self, private))
  txt = c(txt, write_h(self, private))
  txt = c(txt, write_h_var(self, private))
  
  return(txt)
  
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
  # txt = c(txt, "obsMatStd__.row(i) = hsqrt__(stateVec, parVec, inputVec);")
  txt = c(txt, "obsMatStd__.row(i) = sqrt(hvar__(stateVec, parVec, inputVec));")
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


