################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

write_sys_functions = function(self, private){
  
  # Perform substitution of parameters, inputs and states
  f = sapply(seq_along(private$state.names),
             function(i){
               drift.term = hat2pow(private$diff.terms[[i]]$dt)
               # new.drift.term = do.call(substitute, list(drift.term, subsList))
               # sprintf("f__(%i) = %s;",i-1,deparse1(new.drift.term))
             })
  
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
  # drift jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  jac.f = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term <- Deriv::Simplify(Deriv::Deriv(private$diff.terms[[i]]$dt,x=private$state.names[j], cache.exp=F))
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
  txt = c(txt, newtxt)
  
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
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
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
  txt = c(txt, newtxt)
  
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
  vector<Type> hvar__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> hvar__(%s);
    hvar__.setZero();
    %s
    return hvar__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(hvar,collapse="\n\t\t"))
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
  
  newtxt = "\n ////// Newton Solver for Stationary Mean //////
  template<class Type>
  struct Functor {
    vector<Type> p; // parVec
    vector<Type> u; // Inputs
    Functor(const vector<Type> &p, const vector<Type> &u) : p(p), u(u) {}
    Type operator()(const vector<Type> &s) {
  	  vector<Type> f = f__(s, p, u);
      return((f * f).sum());
    }
  };"
  txt = c(txt, newtxt)
  
  newtxt = "\n ////// Lyapunov Solve for Stationary Variance (Linearized System) ///////
  template<class Type>
  matrix<Type> LyapSolver(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    int n = stateVec.size();
    matrix<Type> I(n,n);
    I.setIdentity();
    matrix<Type> A = dfdx__(stateVec, parVec, inputVec);
    matrix<Type> GGT = g__(stateVec, parVec, inputVec) * g__(stateVec, parVec, inputVec).transpose();
    matrix<Type> GGT1d = GGT.vec().matrix();
    /* Kronecker */ 
    matrix<Type> X = tmbutils::kronecker(I, A) + tmbutils::kronecker(A, I);
    /*Solve Linear Systems */
    matrix<Type> P0 = Type(-1.0) * X.householderQr().solve(GGT1d).reshaped(n,n);
    /* Return */  
    return P0;
  }"
  txt = c(txt, newtxt)
  
  # return
  return(txt)
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


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
  txt = c(txt, "DATA_STRUCT(cfg, newton::newton_config_t);")
  
  # Parameter bounds
  txt = c(txt, "DATA_VECTOR(par_lb);")
  txt = c(txt, "DATA_VECTOR(par_ub);")
  
  # Time-step
  txt = c(txt, "DATA_VECTOR(ode_timestep_size);")
  txt = c(txt, "DATA_IVECTOR(ode_timesteps);")
  txt = c(txt, "DATA_INTEGER(ode_solver);")
  
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
  txt = c(txt, sprintf("vector<Type> inputVec(%s), dinputVec(%s), obsVec(%s), parVec0(%s), parVec(%s);", 
                       private$number.of.inputs, 
                       private$number.of.inputs,
                       private$number.of.observations,
                       private$number.of.pars,
                       private$number.of.pars
  ))
  txt = c(txt, sprintf("parVec0 << %s;", paste(private$parameter.names,collapse=", ")))
  txt = c(txt, "parVec = parVec0;")
  # txt = c(txt, "for(int i=0; i < parVec.size(); i++){")
  # txt = c(txt, "if(R_finite(asDouble(par_lb(i)))){")
  # txt = c(txt, "parVec(i) = par_ub(i)*invlogit(parVec(i))+par_lb(i);")
  # txt = c(txt ,"}")
  # txt = c(txt ,"}")
  
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
  txt = c(txt, "V0__.setZero();")
  
  # Storage variables
  txt = c(txt, "\n//////////// storage variables ///////////")
  txt = c(txt, "vector<vector<Type>> xPrior(tsize), xPost(tsize), Innovation(tsize);")
  txt = c(txt, "vector<matrix<Type>> pPrior(tsize), pPost(tsize), InnovationCovariance(tsize);")
  txt = c(txt, "vector<Type> nll_report(tsize);")
  
  # Get starting guesses as solution to stationary problem
  txt = c(txt, "\n//////////// stationary mean ///////////")
  txt = c(txt, "DATA_INTEGER(estimate_stationary_initials);")
  txt = c(txt, "DATA_SCALAR(initial_variance_scaling);")
  txt = c(txt, "if(estimate_stationary_initials == 1){")
  txt = c(txt, "inputVec = inputMat.row(0);")
  txt = c(txt, "Functor<TMBad::ad_aug> F(parVec, inputVec);")
  txt = c(txt, "stateVec = newton::Newton(F, stateVec, cfg);")
  txt = c(txt, "covMat = LyapSolver(stateVec, parVec, inputVec) * initial_variance_scaling;")
  txt = c(txt, "};")
  
  txt = c(txt, "\n//////////// set initial value ///////////")
  txt = c(txt, "xPrior(0) = stateVec;")
  txt = c(txt, "pPrior(0) = covMat;")
  
  # ######## START WITH DATA UPDATE ########
  txt = c(txt, "\n //////////// THE FIRST POINT DATA-UPDATE ///////////")
  # 
  txt = c(txt, "obsVec = obsMat.row(0);")
  txt = c(txt, "is_not_na_obsVec = is_not_na(obsVec);")
  txt = c(txt, "number_of_nonNA_observations = CppAD::Integer(sum(is_not_na_obsVec));")
  # 
  # ########## <OBS IF STATEMENT>  ##########
  txt = c(txt, "if( number_of_nonNA_observations > 0 ){")
  txt = c(txt, "inputVec = inputMat.row(0);")
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
  # 
  txt = c(txt, "// State update // ")
  txt = c(txt, "stateVec = stateVec + K__*e__;")
  txt = c(txt, "covMat = (I__ - K__ * C__) * covMat * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();")
  # 
  txt = c(txt, "// Likelihood contribution // ")
  txt = c(txt, "nll_report(0) = Type(0.5) * atomic::logdet(R__) + Type(0.5) * lossfunction__((e__*(Ri__*e__)).sum(), tukey_loss_parameters, loss_threshold_value, loss_function) + half_log2PI * asDouble(number_of_nonNA_observations);")
  txt = c(txt, "nll__ += nll_report(0);")
  #
  txt = c(txt, "// Save residual mean/covariance // ")
  txt = c(txt, "Innovation(0) = e__;")
  txt = c(txt, "InnovationCovariance(0) = R__;")
  txt = c(txt, "}")
  # ########## </OBS IF STATEMENT>  ##########
  
  txt = c(txt, "xPost(0) = stateVec;")
  txt = c(txt, "pPost(0) = covMat;")
  
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
  txt = c(txt, "inputVec += dinputVec;")
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
  # 
  txt = c(txt, "// State update // ")
  txt = c(txt, "stateVec = stateVec + K__*e__;")
  txt = c(txt, "covMat = (I__ - K__ * C__) * covMat * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();")
  # 
  txt = c(txt, "// Likelihood contribution // ")
  txt = c(txt, "nll_report(i+1) = Type(0.5) * atomic::logdet(R__) + Type(0.5) * lossfunction__((e__*(Ri__*e__)).sum(), tukey_loss_parameters, loss_threshold_value, loss_function) + half_log2PI * asDouble(number_of_nonNA_observations);")
  txt = c(txt, "nll__ += nll_report(i+1);")
  # 
  txt = c(txt, "// Save residual mean/covariance // ")
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
  txt = c(txt, "REPORT(nll_report);")
  
  # return
  return(txt)
}
