########## UNSCENTED KALMAN FILTER FUNCTIONS ###############
############################################################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

write_ukf_functions = function(self, private){
  
  txt = c()
  
  ##################################################
  # DRIFT, DIFFUSION, OBSERVATION FUNCTIONS
  ##################################################
  
  txt = c(txt, write_f(self, private))
  txt = c(txt, write_jac_f(self, private))
  txt = c(txt, write_g(self, private))
  txt = c(txt, write_h(self, private))
  txt = c(txt, write_jac_h(self, private))
  txt = c(txt, write_h_var(self, private))
  
  ##################################################
  # Construct sigma points function
  ##################################################
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
  
  ##################################################
  # Construct Phi Matrix for Sigma Points
  ##################################################
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
  
  ##################################################
  # Construct sigma points drift function
  ##################################################
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
  
  ##################################################
  # Construct sigma points obs function
  ##################################################
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
  
  
  ##################################################
  # STATIONARY SOLVERS
  ##################################################
  
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


write_ukf_estimate = function(self, private)
{
  
  txt = c()
  
  # Observation Vectors
  txt = c(txt, "\n//// observations ////")
  txt = c(txt, "DATA_MATRIX(obsMat);")
  
  # Input Vectors
  txt = c(txt, "\n//// inputs ////")
  txt = c(txt, "DATA_MATRIX(inputMat);")
  
  # Initialize State
  txt = c(txt, "\n//// initial state ////")
  txt = c(txt, "DATA_VECTOR(stateVec);")
  txt = c(txt, "DATA_MATRIX(covMat);")
  txt = c(txt, "DATA_STRUCT(cfg, newton::newton_config_t);")
  
  # Time-step
  txt = c(txt, "DATA_VECTOR(ode_timestep_size);")
  txt = c(txt, "DATA_IVECTOR(ode_timesteps);")
  txt = c(txt, "DATA_INTEGER(ode_solver);")
  
  # Loss parameters
  txt = c(txt, "\n//// loss parameters ////")
  txt = c(txt, "DATA_VECTOR(tukey_loss_parameters);")
  txt = c(txt, "DATA_INTEGER(loss_function);")
  txt = c(txt, "DATA_SCALAR(loss_threshold_value);")
  
  # UKF Transform Hyperparameters
  txt = c(txt, "\n//// unscented transform hyper-parameters ////")  
  txt = c(txt, "DATA_VECTOR(ukf_pars);")
  # txt = c(txt, "DATA_SCALAR(ukf_alpha);")
  # txt = c(txt, "DATA_SCALAR(ukf_beta);")
  # txt = c(txt, "DATA_SCALAR(ukf_kappa);")
  
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
  txt = c(txt, "Type ukf_lambda = pow(ukf_pars(0),2)*(number_of_state_eqs + ukf_pars(2)) - number_of_state_eqs;")
  txt = c(txt, "Type ukf_weights = Type(1.0)/(Type(2.0)*(number_of_state_eqs + ukf_lambda));")
  txt = c(txt, "vector<Type> wm__(nn), wmC__(nn);")
  txt = c(txt, "matrix<Type> Wm__, WcDiag__(nn,nn), I__(nn,nn), W__;")
  txt = c(txt, "I__.setIdentity();")
  txt = c(txt, "WcDiag__.setZero();")
  # wm
  txt = c(txt, "wm__.fill(ukf_weights);")
  txt = c(txt, "wmC__.fill(ukf_weights);")
  txt = c(txt, "wm__(0) = ukf_lambda/(ukf_lambda + number_of_state_eqs);")
  txt = c(txt, "wmC__(0) = ukf_lambda/((number_of_state_eqs + ukf_lambda) + (1-pow(ukf_pars(0),2)+ukf_pars(1)));")
  txt = c(txt, "Wm__ = wm__.replicate(1, nn);")
  txt = c(txt, "WcDiag__.diagonal() = wmC__;")
  txt = c(txt, "W__ = (I__ - Wm__) * WcDiag__ * (I__ - Wm__).transpose();")
  
  # Get starting guesses as solution to stationary problem
  txt = c(txt, "\n//////////// stationary mean ///////////")
  txt = c(txt, "DATA_INTEGER(estimate_stationary_initials);")
  txt = c(txt, "if(estimate_stationary_initials == 1){")
  # txt = c(txt, "inputVec = inputMat.row(0);")
  # txt = c(txt, "Functor<TMBad::ad_aug> F(parVec, inputVec);")
  # txt = c(txt, "stateVec = newton::Newton(F, stateVec, cfg);")
  # txt = c(txt, "covMat = LyapSolver(stateVec, parVec, inputVec);")
  txt = c(txt, "};")
  
  txt = c(txt, "\n//////////// set initial value ///////////")
  txt = c(txt, "xPrior(0) = stateVec;")
  txt = c(txt, "pPrior(0) = covMat;")
  
  # FIRST DATA UPDATE
  txt = c(txt, "\n //////////// FIRST DATA-UPDATE ///////////")
  
  txt = c(txt, "sqrt_covMat = covMat.llt().matrixL();")
  txt = c(txt, "Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, number_of_state_eqs, nn);")
  
  txt = c(txt, "obsVec = obsMat.row(0);")
  txt = c(txt, "is_not_na_obsVec = is_not_na(obsVec);")
  txt = c(txt, "number_of_nonNA_observations = CppAD::Integer(sum(is_not_na_obsVec));")
  
  ########## <OBS IF STATEMENT>  ##########
  txt = c(txt, "if( number_of_nonNA_observations > 0 ){")
  txt = c(txt, "inputVec = inputMat.row(0);")
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
  
  txt = c(txt, "Innovation(0) = e__;")
  txt = c(txt, "InnovationCovariance(0) = Syy__;")
  txt = c(txt, "};")
  ########## /<OBS IF STATEMENT>  ##########
  txt = c(txt, "xPost(0) = stateVec;")
  txt = c(txt, "pPost(0) = covMat;")
  
  ########## <MAIN LOOP>  ##########
  txt = c(txt, "\n //////////// START MAIN LOOP ///////////")
  txt = c(txt, "for(int i=0 ; i < tsize - 1 ; i++){")
  
  # Get sqrt of covariance with cholesky (llt) and compute sigma points
  txt = c(txt, "\n//////////// cholesky cov and sigma points ///////////")
  txt = c(txt, "sqrt_covMat = covMat.llt().matrixL();")
  txt = c(txt, "Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, number_of_state_eqs, nn);")
  
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
