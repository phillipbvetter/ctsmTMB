
########## GENERAL HELPER FUNCTIONS ###############
############################################################
write_helper_cppfunctions = function(){
  
  txt = c()
  
  # Find function finds indices of NAs in a vector
  newtxt = "\n//////////// FIND NA INDICES IN VECTOR ///////////
  template <class Type>
  vector<Type> is_not_na(vector<Type> x){
    vector<Type> y(x.size());
    y.fill(Type(1.0));
      for(int i=0; i<x.size(); i++){
        if( R_IsNA(asDouble(x(i))) ){
          y(i) = Type(0.0);
        }
      }
    return y;
  }"
  txt = c(txt,newtxt)
  
  # This function removes NAs from a vector
  newtxt = "\n//////////// REMOVE NA'S FROM VECTOR ///////////
  template<class Type>
  vector<Type> remove_nas__(vector<Type> obsVec, int number_of_nonNA_observations, vector<Type> is_not_na_vector){
    int ii = 0;
    vector<Type> y_reduced(number_of_nonNA_observations);
      for(int i=0; i < obsVec.size(); i++){
        if(is_not_na_vector(i) == Type(1.0)){
          y_reduced(ii) = obsVec(i);
          ii++;
        }
      }
    return y_reduced;
  }"
  txt = c(txt,newtxt)
  
  # Construct Permutation Matrix E for Kalman Filter NA-removal
  newtxt = "\n//////////// helper fun: construct permutation matrix ///////////
  template <class Type>
  matrix<Type> construct_permutation_matrix(int number_of_nonNA_observations, int number_of_obs_eqs, vector<Type> is_not_na_vector){
	  matrix<Type> E(number_of_nonNA_observations, number_of_obs_eqs);
	  E.setZero();
	  int j=0;
	  for(int i=0; i < number_of_obs_eqs; i++){
      /*if p(i) is 1 then include by setting 1 in diagonal of matrix*/
		  if(is_not_na_vector(i) == Type(1.0)){
        E(j,i) = Type(1.0);
			  j += 1;
      }
	  }
	  return E;
  }"
  txt = c(txt,newtxt)
  
  # Implements Tukey and Huber loss functions
  newtxt = "\n//////////// loss function ///////////
  template<class Type>
  Type lossfunction__(Type x, vector<Type> tukeypars, Type huber_c, int lossFunc){
    Type loss;
    if(lossFunc==1){
      Type a = tukeypars(0);
      Type b = tukeypars(1);
      Type c = tukeypars(2);
      Type d = tukeypars(3);
      loss = d * ( (Type(1.0)/(Type(1.0)+exp(-a*(x-b)))) + c );
    } else if (lossFunc==2){
      Type c_squared = pow(huber_c,2);
      loss = c_squared * (sqrt(1 + (x / c_squared)) - 1);
    } else {
      loss = x;
    }
    return loss;
  }"
  txt = c(txt,newtxt)
  
  # Helper function for MAP estimation
  newtxt = "\n//////////// MAP estimation helper ///////////
  template<class Type>
  vector<Type> get_free_pars__(vector<int> mapints, int sum_mapints, vector<Type> parvec) {
	  vector<Type> ans(sum_mapints);
	  int j=0;
	  for(int i=0;i<mapints.size();i++){
		  if(mapints(i)==1){
			  ans(j) = parvec(i);
			  j += 1;
		  }
	  }
	  return ans;
  }"
  txt = c(txt, newtxt)
  
  return(txt)
}

write_ekf_functions = function(self, private){
  
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
  f = sapply( seq_along(private$state.names),
              function(i){
                drift.term = hat2pow(private$diff.terms[[i]]$dt)
                new.drift.term = do.call(substitute, list(drift.term, subsList))
                sprintf("f__(%i) = %s;",i-1,deparse1(new.drift.term))
              })
  
  newtxt = "\n//////////// drift function //////////
  template<class Type>
  vector<Type> f__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> f__(%s);
    %s
    return f__;
  }"
  newtxt = sprintf(newtxt, private$number.of.states, paste(f,collapse="\n\t\t"))
  txt = c(txt,newtxt)
  
  ##################################################
  # drift jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  jac.f = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term = hat2pow(Deriv::Deriv(private$diff.terms[[i]]$dt,x=private$state.names[j], cache.exp=F))
      new.term = do.call(substitute, list(term, subsList))
      jac.f = c(jac.f, sprintf("dfdx__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  newtxt = "\n//////////// jacobian of drift function ///////////
  template<class Type>
  matrix<Type> dfdx__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> dfdx__(%s, %s);
    %s
    return dfdx__;
  }"
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.states, paste(jac.f, collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # diffusion
  ##################################################
  
  # calculate all the terms and substitute variables
  g = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term = hat2pow(private$diff.terms[[i]][[j+1]])
      new.term = do.call(substitute, list(term, subsList))
      g = c(g, sprintf("g__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  newtxt = "\n//////////// diffusion function ///////////
  template<class Type>
  matrix<Type> g__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> g__(%s, %s);
    %s
    return g__;
  }"
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.diffusions, paste(g,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # observation
  ##################################################
  
  # calculate all the terms and substitute variables
  h = sapply(seq_along(private$obs.names), 
             function(i){
               term = hat2pow(private$obs.eqs.trans[[i]]$rhs)
               new.term = do.call(substitute, list(term, subsList))
               sprintf("h__(%s) = %s;",i-1, deparse1(new.term))
             }) 
  
  newtxt = "\n//////////// observation function ///////////
  template<class Type>
  vector<Type> h__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> h__(%s);
    %s
    return h__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(h,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  jac.h = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = hat2pow(private$diff.terms.obs[[i]][[j]])
      new.term = do.call(substitute, list(term, subsList))
      jac.h = c(jac.h, sprintf("dhdx__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  newtxt = "\n//////////// jacobian of obs function ///////////
  template<class Type>
  matrix<Type> dhdx__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> dhdx__(%s, %s);
    %s
    return dhdx__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, private$number.of.states, paste(jac.h,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # observation variance
  ##################################################
  
  obs.var = lapply(seq_along(private$obs.var.trans), 
                   function(i) {
                     term = hat2pow(private$obs.var.trans[[i]]$rhs)
                     new.term = do.call(substitute, list(term, subsList))
                     sprintf("hvar__(%s) = %s;", i-1, deparse1(new.term))
                   })
  newtxt = "\n//////////// observation variance matrix function ///////////
  template<class Type>
  vector<Type> hvar__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> hvar__(%s);
    %s
    return hvar__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(obs.var,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # 1-Step MOMENT DIFFERENTIAL EQUATIONS
  ##################################################
  
  newtxt = "\n//////////// 1-step f moment ODE ///////////
  template<class Type>
  matrix<Type> cov_ode_1step(matrix<Type> covMat, vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> cov_ode_1step = dfdx__(stateVec, parVec, inputVec) * covMat + covMat * dfdx__(stateVec, parVec, inputVec).transpose() + g__(stateVec, parVec, inputVec) * g__(stateVec, parVec, inputVec).transpose();
    return cov_ode_1step;
  }"
  txt = c(txt, newtxt)
  
  ##################################################
  # EKF ODE SOLVER
  ##################################################
  
  newtxt = "\n//////////// ODE SOLVER ///////////
  template<class Type>
  struct ode_integration {
  
	  vector<Type> X1;
	  matrix<Type> P1;
  
	  ode_integration(matrix<Type> covMat, vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec, vector<Type> dinputVec, Type dt, int ode_solver){
      /*Initial State and Cov Values*/
      vector<Type> X0 = stateVec;
      matrix<Type> P0 = covMat;
  
		  /*Forward Euler*/
		  if(ode_solver == 1){
  
		   X1 = X0 + f__(stateVec, parVec, inputVec) * dt;
       P1 = P0 + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt;
  
		  /*4th Order Runge-Kutta 4th*/
		  } else if (ode_solver == 2){
  
		   vector<Type> k1, k2, k3, k4;
		   matrix<Type> c1, c2, c3, c4;
  
		   /*1. Approx Slope at Initial Point*/
       k1 = f__(stateVec, parVec, inputVec);
       c1 = cov_ode_1step(covMat, stateVec, parVec, inputVec);
  
		   /*2. First Approx Slope at Midpoint*/
       inputVec += 0.5 * dinputVec;
       stateVec = X0 + 0.5 * dt * k1;
       covMat   = P0 + 0.5 * dt * c1;
       k2       = f__(stateVec, parVec, inputVec); 
       c2       = cov_ode_1step(covMat, stateVec, parVec, inputVec);        
  
		   /*3. Second Approx Slope at Midpoint*/
       stateVec = X0 + 0.5 * dt * k2;
       covMat   = P0 + 0.5 * dt * c2;
       k3       = f__(stateVec, parVec, inputVec); 
       c3       = cov_ode_1step(covMat, stateVec, parVec, inputVec);
  
		   /*4. Approx Slope at End Point*/
       inputVec += 0.5 * dinputVec;
       stateVec = X0 + dt * k3;
       covMat   = P0 + dt * c3;
       k4       = f__(stateVec, parVec, inputVec); 
       c4       = cov_ode_1step(covMat, stateVec, parVec, inputVec);
  
		   /*ODE UPDATE*/
		   X1 = X0 + dt * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
		   P1 = P0 + dt * (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0;
  
		 } else {
		 /*nothing*/
		 }
		}
  };"
  txt = c(txt, newtxt)
  
  # return
  return(txt)
}

########## UNSCENTED KALMAN FILTER FUNCTIONS ###############
############################################################
write_ukf_functions = function(self, private){
  
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
  f = sapply( seq_along(private$state.names),
              function(i){
                drift.term = hat2pow(private$diff.terms[[i]]$dt)
                new.drift.term = do.call(substitute, list(drift.term, subsList))
                sprintf("f__(%i) = %s;",i-1,deparse1(new.drift.term))
              })
  
  newtxt = "\n//////////// drift function //////////
  template<class Type>
  vector<Type> f__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> f__(%s);
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
      term = hat2pow(private$diff.terms[[i]][[j+1]])
      new.term = do.call(substitute, list(term, subsList))
      g = c(g, sprintf("g__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  newtxt = "\n//////////// diffusion function ///////////
  template<class Type>
  matrix<Type> g__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> g__(%s, %s);
    %s
    return g__;
  }"
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.diffusions, paste(g,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # observation
  ##################################################
  
  # calculate all the terms and substitute variables
  h = sapply(seq_along(private$obs.names), 
             function(i){
               term = hat2pow(private$obs.eqs.trans[[i]]$rhs)
               new.term = do.call(substitute, list(term, subsList))
               sprintf("h__(%s) = %s;",i-1, deparse1(new.term))
             }) 
  
  newtxt = "\n//////////// observation function ///////////
  template<class Type>
  vector<Type> h__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> h__(%s);
    %s
    return h__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(h,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # observation variance
  ##################################################
  
  obs.var = lapply(seq_along(private$obs.var.trans), 
                   function(i) {
                     term = hat2pow(private$obs.var.trans[[i]]$rhs)
                     new.term = do.call(substitute, list(term, subsList))
                     sprintf("hvar__(%s) = %s;", i-1, deparse1(new.term))
                   })
  newtxt = "\n//////////// observation variance matrix function ///////////
  template<class Type>
  vector<Type> hvar__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> hvar__(%s);
    %s
    return hvar__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(obs.var,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  # Construct sigma points function
  newtxt = "\n//////////// UKF sigma points ///////////
  template<class Type>
  matrix<Type> create_sigmapoints_from_stateVec(vector<Type> stateVec, matrix<Type> sqrt_covMat, int n, int nn){
    matrix<Type> Xsigmapoints(n, nn);
    vector<Type> Ai;
    Xsigmapoints.col(0) = stateVec;
    for(int i=1; i < n+1; i++){
      Ai = sqrt_covMat.col(i-1);
      Xsigmapoints.col(i) = stateVec + sqrt(1+n) * Ai;
      Xsigmapoints.col(i+n) = stateVec - sqrt(1+n) * Ai;
    }
    return Xsigmapoints;
  }"
  txt = c(txt, newtxt)
  
  # Construct Phi Matrix for Sigma Points
  newtxt = "\n//////////// UKF function ///////////
  template<class Type>
  matrix<Type> Phi__(matrix<Type> M){
    matrix<Type> K(M.col(0).size(), M.row(0).size());
    K.setZero();
    K.template triangularView<Eigen::Lower>() = M.template triangularView<Eigen::Lower>();
    K.diagonal() = K.diagonal()/Type(2.0);
    return K;
  }"
  txt = c(txt, newtxt)
  
  # Construct sigma points drift function
  newtxt = "\n//////////// UKF sigma points drift function ///////////
  template<class Type>
  matrix<Type> ukf_f__(matrix<Type> Xsigmapoints, vector<Type> parVec, vector<Type> inputVec, int n, int nn){
    matrix<Type> F(n,nn);
    vector<Type> stateVec;
    for(int i=0; i < nn; i++){
      stateVec = Xsigmapoints.col(i);
      F.col(i) = f__(stateVec, parVec, inputVec);
    }
    return F;
  }"
  txt = c(txt, newtxt)
  
  # Construct sigma points obs function
  newtxt = "\n//////////// UKF sigma points obs function ///////////
  template<class Type>
  matrix<Type> ukf_h__(matrix<Type> Xsigmapoints, vector<Type> parVec, vector<Type> inputVec, int m, int nn){
    matrix<Type> H(m,nn);
    vector<Type> stateVec;
    for(int i=0; i < nn; i++){
      stateVec = Xsigmapoints.col(i);
      H.col(i) = h__(stateVec, parVec, inputVec);
    }
    return H;
  }"
  txt = c(txt, newtxt)
  
  
  
  ##################################################
  # UKF 1-Step MOMENT DIFFERENTIAL EQUATIONS
  ##################################################
  
  newtxt = "\n//////////// 1-step f moment ODE ///////////
  template<class Type>
  matrix<Type> ukf_1step(matrix<Type> Xsigmapoints, matrix<Type> sqrt_covMat, matrix<Type> W, vector<Type> wm, vector<Type> parVec, vector<Type> inputVec, int n, int nn)
  {
  matrix<Type> F, G, Ainv, M, A_Phi_M, F_rhs, F_rhs0, F_rhs1(n, nn);
  F_rhs1.setZero();
  F = ukf_f__(Xsigmapoints, parVec, inputVec, n, nn);
  // G is currently not using sigma points because assumed state independent, but we hack it by just using the state vector
  G = g__(vector<Type>(Xsigmapoints.col(0)), parVec, inputVec);
  Ainv = sqrt_covMat.inverse();
  M = Ainv * (Xsigmapoints * W * F.transpose() + F * W * Xsigmapoints.transpose() + G*G.transpose()) * Ainv.transpose();
  A_Phi_M = sqrt_covMat * Phi__(M);
  F_rhs1.block(0, 1, n, n) = A_Phi_M;
  F_rhs1.block(0, n+1, n, n) = -A_Phi_M;
  F_rhs0 = (F * wm).replicate(1, nn);
  F_rhs = F_rhs0 + sqrt(3.0) * F_rhs1;
  
  return F_rhs;
  }"
  txt = c(txt, newtxt)
  
  ##################################################
  # UKF ODE SOLVER
  ##################################################
  
  newtxt = "\n//////////// ODE SOLVER ///////////
  template<class Type>
  struct ukf_ode_integration {
  
  matrix<Type> X1sigmapoints;
  
  ukf_ode_integration(matrix<Type> Xsigmapoints, matrix<Type> sqrt_covMat, matrix<Type> W, vector<Type> wm, vector<Type> parVec, vector<Type> inputVec, vector<Type> dinputVec, Type dt, int ode_solver, int n, int nn){
      /*Initial State and Cov Values*/
      matrix<Type> X0sigmapoints = Xsigmapoints;
  
		  /*Forward Euler*/
		  if(ode_solver == 1){
  
		   X1sigmapoints = X0sigmapoints + ukf_1step(Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, n, nn) * dt;
  
		  /*4th Order Runge-Kutta 4th*/
		  } else if (ode_solver == 2){
  
		   matrix<Type> c1, c2, c3, c4;
  
		   /*1. Approx Slope at Initial Point*/
       c1 = ukf_1step(Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, n, nn);
  
		   /*2. First Approx Slope at Midpoint*/
       inputVec += 0.5 * dinputVec;
       Xsigmapoints = X0sigmapoints + 0.5 * dt * c1;
       sqrt_covMat = ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt(Type(3.0))).block(0, 1, n, n);
       c2 = ukf_1step(Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, n, nn);        
  
		   /*3. Second Approx Slope at Midpoint*/
       Xsigmapoints = X0sigmapoints + 0.5 * dt * c2;
       sqrt_covMat = ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt(Type(3.0))).block(0, 1, n, n);
       c3 = ukf_1step(Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, n, nn);
  
		   /*4. Approx Slope at End Point*/
       inputVec += 0.5 * dinputVec;
       Xsigmapoints = X0sigmapoints + dt * c3;
       sqrt_covMat = ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt(Type(3.0))).block(0, 1, n, n);
       c4 = ukf_1step(Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, n, nn);
  
		   /*ODE UPDATE*/
       X1sigmapoints = X0sigmapoints + dt * (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0; 
  
		 } else {
		 /*nothing*/
		 }
		}
  };"
  txt = c(txt, newtxt)
  
  # return
  return(txt) 
}


write_ukf_estimate = function(self, private)
{
  
  txt = c()
  
  # Observation Vectors
  txt = c(txt, "\n//// observations ////")
  txt = c(txt, "DATA_MATRIX(obsMat)")
  
  # Input Vectors
  txt = c(txt, "\n//// inputs ////")
  txt = c(txt, "DATA_MATRIX(inputMat)")
  
  # Initialize State
  txt = c(txt, "\n//// initial state ////")
  txt = c(txt, "DATA_VECTOR(stateVec);")
  txt = c(txt, "DATA_MATRIX(covMat);")
  
  # Time-step
  txt = c(txt, "DATA_VECTOR(ode_timestep_size);")
  txt = c(txt, "DATA_IVECTOR(ode_timesteps);")
  
  # Loss parameters
  txt = c(txt, "\n//// loss parameters ////")
  txt = c(txt, "DATA_VECTOR(tukey_loss_parameters);")
  txt = c(txt, "DATA_INTEGER(loss_function);")
  txt = c(txt, "DATA_SCALAR(loss_threshold_value);")
  
  # UKF Transform Hyperparameters
  txt = c(txt, "\n//// unscented transform hyper-parameters ////")  
  txt = c(txt, "DATA_SCALAR(ukf_alpha);")
  txt = c(txt, "DATA_SCALAR(ukf_beta);")
  txt = c(txt, "DATA_SCALAR(ukf_kappa);")
  
  # Maximum a Posterior
  txt = c(txt, "\n//// map estimation ////")
  txt = c(txt, "DATA_INTEGER(MAP_bool);")
  
  # Parameters
  txt = c(txt, "\n//// parameters ////")
  for(i in 1:length(private$parameters)){
    txt = c(txt, sprintf("PARAMETER(%s);",private$parameter.names[i]))
  }
  
  # System size
  txt = c(txt, "\n//// system size ////")
  txt = c( txt, "DATA_INTEGER(number_of_state_eqs);")
  txt = c( txt, "DATA_INTEGER(number_of_obs_eqs);")
  txt = c( txt, "int tsize = inputMat.col(0).size();")
  txt = c( txt, "int nn = 2*number_of_state_eqs + 1;")
  
  # state, par, input, obs
  txt = c(txt, "\n//// state, par, input, obs vectors ////")
  txt = c(txt, sprintf("vector<Type> inputVec(%s), dinputVec(%s), obsVec(%s), parVec(%s);", 
                       private$number.of.inputs, 
                       private$number.of.inputs,
                       private$number.of.observations,
                       private$number.of.pars
  ))
  txt = c(txt, sprintf("parVec << %s;", paste(private$parameter.names,collapse=", ")))
  
  ##################################################
  # Storage variables
  txt = c(txt, "\n//////////// storage variables ///////////")
  txt = c(txt, "vector<vector<Type>> xPrior(tsize), xPost(tsize), Innovation(tsize);")
  txt = c(txt, "vector<matrix<Type>> pPrior(tsize), pPost(tsize), InnovationCovariance(tsize);")
  txt = c(txt, "vector<matrix<Type>> blatest1(tsize),blatest2(tsize);")
  
  txt = c(txt, "\n//////////// set initial value ///////////")
  txt = c(txt, "xPrior(0) = stateVec, xPost(0) = stateVec;")
  txt = c(txt, "pPrior(0) = covMat, pPost(0) = covMat;")
  
  txt = c(txt, "\n //////////// initialize variables ///////////")
  txt = c(txt, "int number_of_nonNA_observations;")
  txt = c(txt, "Type half_log2PI = Type(0.5) * log(2*M_PI);")
  txt = c(txt, "vector<Type> data_vector__(number_of_obs_eqs);")
  txt = c(txt, "vector<Type> y__,e__, is_not_na_obsVec;")
  txt = c(txt, "matrix<Type> E__, H__, Syy__, SyyInv__, Sxy__, K__, Xsigmapoints, sqrt_covMat;")
  
  ##################################################
  # Observation variance matrix
  txt = c(txt,"matrix<Type> V0__(number_of_obs_eqs,number_of_obs_eqs);")
  txt = c(txt,"V0__.setZero();")
  
  ##################################################
  # Weights
  txt = c(txt, "\n//////////// create weights ///////////")
  txt = c(txt, "Type ukf_lambda = pow(ukf_alpha,2)*(number_of_state_eqs + ukf_kappa) - number_of_state_eqs;")
  txt = c(txt, "Type ukf_weights = Type(1.0)/(Type(2.0)*(number_of_state_eqs + ukf_lambda));")
  txt = c(txt, "vector<Type> wm__(nn), wmC__(nn);")
  txt = c(txt, "matrix<Type> Wm__, WcDiag__(nn,nn), I__(nn,nn), W__;")
  txt = c(txt, "I__.setIdentity();")
  txt = c(txt, "WcDiag__.setZero();")
  # wm
  txt = c(txt, "wm__.fill(ukf_weights);")
  txt = c(txt, "wmC__.fill(ukf_weights);")
  txt = c(txt, "wm__(0) = ukf_lambda/(ukf_lambda + number_of_state_eqs);")
  txt = c(txt, "wmC__(0) = ukf_lambda/((number_of_state_eqs + ukf_lambda) + (1-pow(ukf_alpha,2)+ukf_beta));")
  txt = c(txt, "Wm__ = wm__.replicate(1, nn);")
  txt = c(txt, "WcDiag__.diagonal() = wmC__;")
  txt = c(txt, "W__ = (I__ - Wm__) * WcDiag__ * (I__ - Wm__).transpose();")
  
  
  
  ##################################################
  
  ########## <MAIN LOOP>  ##########
  txt = c(txt, "\n //////////// START MAIN LOOP ///////////")
  txt = c(txt, "for(int i=0 ; i < tsize - 1 ; i++){")
  
  # Get sqrt of covariance with cholesky (llt) and compute sigma points
  txt = c(txt, "\n//////////// cholesky cov and sigma points ///////////")
  txt = c(txt, "sqrt_covMat = covMat.llt().matrixL();")
  txt = c(txt, "Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, number_of_state_eqs, nn);")
  # txt = c(txt, "blatest1(i) = Xsigmapoints;")
  
  # Set inputs
  txt = c(txt, "inputVec = inputMat.row(i);")
  txt = c(txt, "dinputVec = (inputMat.row(i+1) - inputMat.row(i))/ode_timesteps(i);")
  # Solve Moment ODEs
  txt = c(txt, "\n //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////")
  txt = c(txt, "for(int j=0 ; j < ode_timesteps(i) ; j++){")
  txt = c(txt, "ukf_ode_integration<Type> odelist = {Xsigmapoints, sqrt_covMat, W__, wm__, parVec, inputVec, dinputVec, ode_timestep_size(i), ode_solver, number_of_state_eqs, nn};")
  txt = c(txt, "Xsigmapoints = odelist.X1sigmapoints;")
  txt = c(txt, "sqrt_covMat = ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt(Type(3.0))).block(0, 1, number_of_state_eqs, number_of_state_eqs );")
  txt = c(txt, "inputVec += dinputVec;")
  txt = c(txt, "}")
  txt = c(txt, "stateVec = Xsigmapoints.col(0);")
  txt = c(txt, "covMat = sqrt_covMat * sqrt_covMat.transpose();")
  txt = c(txt, "xPrior(i+1) = stateVec;")
  txt = c(txt, "pPrior(i+1) = covMat;")
  
  # txt = c(txt, "blatest2(i) = Xsigmapoints;")
  
  # Data Update
  txt = c(txt, "\n //////////// DATA-UPDATE ///////////")
  # Set observation indices (i+1)
  txt = c(txt, "obsVec = obsMat.row(i+1);")
  # Check the number of NAs in the observation vector
  txt = c(txt, "is_not_na_obsVec = is_not_na(obsVec);")
  txt = c(txt, "number_of_nonNA_observations = CppAD::Integer(sum(is_not_na_obsVec));")
  ########## <OBS IF STATEMENT>  ##########
  txt = c(txt, "if( number_of_nonNA_observations > 0 ){")
  txt = c(txt, "inputVec = inputMat.row(i+1);")
  txt = c(txt, "y__ = remove_nas__(obsVec, number_of_nonNA_observations, is_not_na_obsVec);")
  txt = c(txt, "E__ = construct_permutation_matrix(number_of_nonNA_observations, number_of_obs_eqs, is_not_na_obsVec);")
  txt = c(txt, "H__ = ukf_h__(Xsigmapoints, parVec, inputVec, number_of_obs_eqs, nn);")
  txt = c(txt, "e__  = y__ - E__ * (H__ * wm__);")
  txt = c(txt, "V0__.diagonal() << hvar__(stateVec, parVec, inputVec);")
  txt = c(txt, "Syy__  = E__ * (H__ * W__ * H__.transpose() + V0__) * E__.transpose();")
  txt = c(txt, "SyyInv__  = Syy__.inverse();")
  txt = c(txt, "Sxy__  = Xsigmapoints * W__ * H__.transpose() * E__.transpose();")
  txt = c(txt, "K__ = Sxy__ * SyyInv__;")
  txt = c(txt, "stateVec = stateVec + K__ * e__;")
  txt = c(txt, "covMat = covMat - K__ * Syy__ * K__.transpose();")
  txt = c(txt, "nll__ += Type(0.5) * atomic::logdet(Syy__) + Type(0.5) * lossfunction__((e__*(SyyInv__*e__)).sum(), tukey_loss_parameters, loss_threshold_value, loss_function) + half_log2PI * asDouble(number_of_nonNA_observations);")
  txt = c(txt, "Innovation(i+1) = e__;")
  txt = c(txt, "InnovationCovariance(i+1) = Syy__;")
  txt = c(txt, "};")
  txt = c(txt, "xPost(i+1) = stateVec;")
  txt = c(txt, "pPost(i+1) = covMat;")
  txt = c(txt, "};")
  txt = c(txt, "//////////// END MAIN LOOP ///////////")
  
  ##################################################
  # Maximum-A-Posterior
  txt = c(txt, "\n//////////// MAP CONTRIBUTION ///////////")
  txt = c(txt, "if(MAP_bool == 1){")
  txt = c(txt, "DATA_VECTOR(map_mean__);")
  txt = c(txt, "DATA_MATRIX(map_cov__);")
  txt = c(txt, "DATA_IVECTOR(map_ints__);")
  txt = c(txt, "DATA_INTEGER(sum_map_ints__);")
  txt = c(txt, "vector<Type> map_pars__;")
  txt = c(txt, "map_pars__ = get_free_pars__(map_ints__, sum_map_ints__, parVec);")
  txt = c(txt, "vector<Type> pars_eps__ = map_pars__ - map_mean__;")
  txt = c(txt, "matrix<Type> map_invcov__ = map_cov__.inverse();")
  txt = c(txt, "Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();")
  txt = c(txt, "nll__ += map_nll__;")
  txt = c(txt, "REPORT(map_nll__);")
  txt = c(txt, "REPORT(map_pars__);")
  txt = c(txt, "REPORT(pars_eps__);")
  txt = c(txt, "}")
  
  ##################################################
  # Report variables
  txt = c(txt, "\n//////////// Report //////////////")
  txt = c(txt, "REPORT(Innovation);")
  txt = c(txt, "REPORT(InnovationCovariance);")
  txt = c(txt, "REPORT(xPrior);")
  txt = c(txt, "REPORT(xPost);")
  txt = c(txt, "REPORT(pPrior);")
  txt = c(txt, "REPORT(pPost);")
  
  return(txt)
}

write_ekf_estimate = function(self, private){
  
  txt = c()
  
  # Observation Vectors
  txt = c(txt, "\n//// observations ////")
  txt = c(txt, "DATA_MATRIX(obsMat)")
  
  # Input Vectors
  txt = c(txt, "\n//// inputs ////")
  txt = c(txt, "DATA_MATRIX(inputMat)")
  
  # Initialize State
  txt = c(txt, "\n//// initial state ////")
  txt = c(txt, "DATA_VECTOR(stateVec);")
  txt = c(txt, "DATA_MATRIX(covMat);")
  
  # Time-step
  txt = c(txt, "DATA_VECTOR(ode_timestep_size);")
  txt = c(txt, "DATA_IVECTOR(ode_timesteps);")
  
  # Loss parameters
  txt = c(txt, "\n//// loss parameters ////")
  txt = c(txt, "DATA_VECTOR(tukey_loss_parameters);")
  txt = c(txt, "DATA_INTEGER(loss_function);")
  txt = c(txt, "DATA_SCALAR(loss_threshold_value);")
  
  # Maximum a Posterior
  txt = c(txt, "\n//// map estimation ////")
  txt = c(txt, "DATA_INTEGER(MAP_bool);")
  
  # Parameters
  txt = c(txt, "\n//// parameters ////")
  for(i in 1:length(private$parameters)){
    txt = c(txt, sprintf("PARAMETER(%s);",private$parameter.names[i]))
  }
  
  # system size
  txt = c(txt, "\n//// system size ////")
  txt = c( txt , "DATA_INTEGER(number_of_state_eqs);")
  txt = c( txt , "DATA_INTEGER(number_of_obs_eqs);")
  txt = c( txt , "int tsize = inputMat.col(0).size();")
  
  # state, par, input, obs
  txt = c(txt, "\n//// state, par, input, obs vectors ////")
  txt = c(txt, sprintf("vector<Type> inputVec(%s), dinputVec(%s), obsVec(%s), parVec(%s);", 
                       private$number.of.inputs, 
                       private$number.of.inputs,
                       private$number.of.observations,
                       private$number.of.pars
  ))
  txt = c(txt, sprintf("parVec << %s;", paste(private$parameter.names,collapse=", ")))
  
  # Storage variables
  txt = c(txt, "\n//////////// storage variables ///////////")
  txt = c(txt, "vector<vector<Type>> xPrior(tsize), xPost(tsize), Innovation(tsize);")
  txt = c(txt, "vector<matrix<Type>> pPrior(tsize), pPost(tsize), InnovationCovariance(tsize);")
  
  txt = c(txt, "\n//////////// set initial value ///////////")
  txt = c(txt, "xPrior(0) = stateVec, xPost(0) = stateVec;")
  txt = c(txt, "pPrior(0) = covMat, pPost(0) = covMat;")
  
  # Initiaze variables
  txt = c(txt, "\n //////////// initialize variables ///////////")
  txt = c(txt, "int number_of_nonNA_observations;")
  txt = c(txt, "Type half_log2PI = Type(0.5) * log(2*M_PI);")
  txt = c(txt, "vector<Type> is_not_na_obsVec, e__, y__, F__, H__;")
  txt = c(txt, "matrix<Type> C__, R__, K__, E__, V__, Ri__, A__, G__;")
  
  # Identity Matrix
  txt = c(txt, "\n//////////// identity matrix ///////////")
  txt = c(txt, "matrix<Type> I__(number_of_state_eqs, number_of_state_eqs), V0__(number_of_obs_eqs, number_of_obs_eqs);")
  txt = c(txt, "I__.setIdentity();")
  
  ########## <MAIN LOOP>  ##########
  txt = c(txt, "\n //////////// START MAIN LOOP ///////////")
  
  # Time for-loop
  txt = c(txt, "for(int i=0 ; i < tsize - 1 ; i++){")
  # Set inputs
  txt = c(txt, "inputVec = inputMat.row(i);")
  txt = c(txt, "dinputVec = (inputMat.row(i+1) - inputMat.row(i))/ode_timesteps(i);")
  # Solve Moment ODEs
  txt = c(txt, "\n //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////")
  txt = c(txt, "for(int j=0 ; j < ode_timesteps(i) ; j++){")
  txt = c(txt, "ode_integration<Type> odelist = {covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size(i), ode_solver};")
  txt = c(txt, "stateVec = odelist.X1;")
  txt = c(txt, "covMat = odelist.P1;")
  # txt = c(txt, "inputVec += dinputVec;")
  txt = c(txt, "}")
  txt = c(txt, "xPrior(i+1) = stateVec;")
  txt = c(txt, "pPrior(i+1) = covMat;")
  
  # Data Update
  txt = c(txt, "\n //////////// DATA-UPDATE ///////////")
  
  # Set observation indices (i+1)
  txt = c(txt, "obsVec = obsMat.row(i+1);")
  
  # Check the number of NAs in the observation vector
  txt = c(txt, "is_not_na_obsVec = is_not_na(obsVec);")
  txt = c(txt, "number_of_nonNA_observations = CppAD::Integer(sum(is_not_na_obsVec));")
  
  ########## <OBS IF STATEMENT>  ##########
  txt = c(txt, "if( number_of_nonNA_observations > 0 ){")
  txt = c(txt, "inputVec = inputMat.row(i+1);")
  txt = c(txt, "y__ = remove_nas__(obsVec, number_of_nonNA_observations, is_not_na_obsVec);")
  txt = c(txt, "E__ = construct_permutation_matrix(number_of_nonNA_observations, number_of_obs_eqs, is_not_na_obsVec);")
  txt = c(txt, "H__ = h__(stateVec, parVec, inputVec);")
  txt = c(txt, "C__ = E__ * dhdx__(stateVec, parVec, inputVec);")
  txt = c(txt, "e__ = y__ - E__ * H__;")
  txt = c(txt, "V0__.diagonal() << hvar__(stateVec, parVec, inputVec);")
  txt = c(txt, "V__ = E__ * V0__ * E__.transpose();")
  txt = c(txt, "R__ = C__ * covMat * C__.transpose() + V__;")
  txt = c(txt, "Ri__ = R__.inverse();")
  txt = c(txt, "K__ = covMat * C__.transpose() * Ri__;")
  txt = c(txt, "stateVec = stateVec + K__*e__;")
  txt = c(txt, "covMat = (I__ - K__ * C__) * covMat * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();")
  txt = c(txt, "nll__ += Type(0.5) * atomic::logdet(R__) + Type(0.5) * lossfunction__((e__*(Ri__*e__)).sum(), tukey_loss_parameters, loss_threshold_value, loss_function) + half_log2PI * asDouble(number_of_nonNA_observations);")
  txt = c(txt, "Innovation(i+1) = e__;")
  txt = c(txt, "InnovationCovariance(i+1) = R__;")
  txt = c(txt, "}")
  ########## </OBS IF STATEMENT>  ##########
  
  # Store posterior values
  txt = c(txt, "xPost(i+1) = stateVec;")
  txt = c(txt, "pPost(i+1) = covMat;")
  
  txt = c(txt, "}")
  ########## </MAIN LOOP>  ##########
  txt = c(txt, "//////////// END MAIN LOOP ///////////")
  
  # Maximum-A-Posterior
  txt = c(txt, "\n//////////// MAP CONTRIBUTION ///////////")
  txt = c(txt, "if(MAP_bool == 1){")
  txt = c(txt, "DATA_VECTOR(map_mean__);")
  txt = c(txt, "DATA_MATRIX(map_cov__);")
  txt = c(txt, "DATA_IVECTOR(map_ints__);")
  txt = c(txt, "DATA_INTEGER(sum_map_ints__);")
  txt = c(txt, "vector<Type> map_pars__;")
  txt = c(txt, "map_pars__ = get_free_pars__(map_ints__, sum_map_ints__, parVec);")
  txt = c(txt, "vector<Type> pars_eps__ = map_pars__ - map_mean__;")
  txt = c(txt, "matrix<Type> map_invcov__ = map_cov__.inverse();")
  txt = c(txt, "Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();")
  txt = c(txt, "nll__ += map_nll__;")
  txt = c(txt, "REPORT(map_nll__);")
  txt = c(txt, "REPORT(map_pars__);")
  txt = c(txt, "REPORT(pars_eps__);")
  txt = c(txt, "}")
  
  # Report variables
  txt = c(txt, "\n//////////// Report //////////////")
  txt = c(txt, "REPORT(Innovation);")
  txt = c(txt, "REPORT(InnovationCovariance);")
  txt = c(txt, "REPORT(xPrior);")
  txt = c(txt, "REPORT(xPost);")
  txt = c(txt, "REPORT(pPrior);")
  txt = c(txt, "REPORT(pPost);")
  
  # return
  return(txt)
}

write_cppfile = function(self, private) {
  
  #Initialize C++ file
  fileconn = file(paste0(private$cppfile.path,".cpp"))
  
  # Include TMB header
  txt = "#include <TMB.hpp>"
  
  # density namespace for special functions
  newtxt = "using namespace density;"
  txt = c(txt, newtxt)
  
  ##################################################
  # GENERAL HELPER FUNCTIONS
  ##################################################
  txt = c(txt, "\n//////////// HELPER FUNCTIONS ///////////")
  txt = c(txt, "//////////// HELPER FUNCTIONS ///////////")
  txt = c(txt, "//////////// HELPER FUNCTIONS ///////////")
  newtxt = write_helper_cppfunctions()
  txt = c(txt, newtxt)
  
  ####################################################################
  # FUNCTIONS FOR DRIFT, DIFFUSION, OBSERVATION, OBS VARIANCE, ETC...
  ####################################################################
  txt = c(txt, "\n//////////// EKF FUNCTIONS ///////////")
  txt = c(txt, "//////////// EKF FUNCTIONS ///////////")
  txt = c(txt, "//////////// EKF FUNCTIONS ///////////")
  newtxt = write_ekf_functions(self, private)
  txt = c(txt, newtxt)
  
  ##################################################
  # UKF FUNCTIONS
  ##################################################
  
  txt = c(txt, "\n//////////// UKF FUNCTIONS ///////////")
  txt = c(txt, "//////////// UKF FUNCTIONS ///////////")
  txt = c(txt, "//////////// UKF FUNCTIONS ///////////")
  # newtxt = write_ukf_functions(self,private)
  # txt = c(txt,newtxt)
  
  ##################################################
  # BEGIN OBJECTIVE FUNCTION
  ##################################################
  
  # Initialize objective function
  txt = c(txt, "\n//////////// OBJECTIVE FUNCTION START ///////////")
  txt = c(txt, "//////////// OBJECTIVE FUNCTION START ///////////")
  txt = c(txt, "//////////// OBJECTIVE FUNCTION START ///////////")
  
  txt = c(txt,"template<class Type>\nType objective_function<Type>::operator() ()\n{")
  
  txt = c(txt, "DATA_INTEGER(estimation_method);")
  txt = c(txt, "DATA_INTEGER(ode_solver);")
  txt = c(txt, "Type nll__ = 0;")
  
  ##################################################
  # EKF OBJECTIVE FUNCTION
  ##################################################
  
  txt = c(txt, "\n//////////// EXTENDED KALMAN FILTER ///////////")
  txt = c(txt, "//////////// EXTENDED KALMAN FILTER ///////////")
  txt = c(txt, "//////////// EXTENDED KALMAN FILTER ///////////")
  txt = c(txt, "if (estimation_method == 1){")
  
  newtxt = write_ekf_estimate(self, private)
  txt = c(txt, newtxt)
  
  txt = c(txt, "\n//////////// UNSCENTED KALMAN FILTER ///////////")
  txt = c(txt, "//////////// UNSCENTED KALMAN FILTER ///////////")
  txt = c(txt, "//////////// UNSCENTED KALMAN FILTER ///////////")
  txt = c(txt, "} else if (estimation_method == 2) {")
  
  newtxt = write_ukf_estimate(self, private)
  txt = c(txt, newtxt)
  
  txt = c(txt, "} else {")
  
  txt = c(txt, "}")
  
  # Return nll
  txt = c(txt, "\n //////////// Return //////////////")
  txt = c(txt, "return nll__;")
  txt = c(txt, "}")
  
  txt = c(txt, "\n//////////// OBJECTIVE FUNCTION END ///////////")
  txt = c(txt, "//////////// OBJECTIVE FUNCTION END ///////////")
  txt = c(txt, "//////////// OBJECTIVE FUNCTION END ///////////")
  
  # Write cpp function and close file connection
  writeLines(txt,fileconn)
  close(fileconn)
  
  # Return
  return(invisible(self))
}

write_method_cppfile = function(self, private) {
  
  #Initialize C++ file
  fileconn = file(paste0(private$cppfile.path.with.method,".cpp"))
  
  # Header etc...
  txt = "#include <TMB.hpp>"
  newtxt = "using namespace density;"
  txt = c(txt, newtxt)
  
  # Various helper functions
  newtxt = write_helper_cppfunctions()
  txt = c(txt, newtxt)
  
  # Specific method functions
  if(private$method=="ekf"){
    newtxt = write_ekf_functions(self, private)
    txt = c(txt, newtxt)
  }
  if(private$method=="ukf"){
    newtxt = write_ukf_functions(self,private)
    txt = c(txt,newtxt)
  }
  
  # Initialize TMB Objective Function
  
  txt = c(txt,"template<class Type>\nType objective_function<Type>::operator() ()\n{")
  # txt = c(txt, "DATA_INTEGER(estimation_method);")
  txt = c(txt, "DATA_INTEGER(ode_solver);")
  txt = c(txt, "Type nll__ = 0;")
  
  # Specific estimation method
  if(private$method=="ekf"){
    newtxt = write_ekf_estimate(self, private)
    txt = c(txt, newtxt)
  }
  if(private$method=="ukf"){
    newtxt = write_ukf_estimate(self, private)
    txt = c(txt, newtxt)
  }
  
  # Return statement
  txt = c(txt, "return nll__;")
  txt = c(txt, "}")
  
  # Write cpp function and close file connection
  writeLines(txt,fileconn)
  close(fileconn)
  
  # Return
  return(invisible(self))
}
