// SYSINFO: NUMBER_OF_STATES
// SYSINFO: NUMBER_OF_OBS
// SYSINFO: NUMBER_OF_INPUTS
// SYSINFO: NUMBER_OF_PARS

#include <TMB.hpp>

// NAMESPACES
using namespace density;

// CONSTANTS  
const double pi = M_PI;
const double k_smooth = 5.0; // for loss function

// INSERT F

// INSERT DFDX

// INSERT G

// INSERT H

// INSERT DHDX

// INSERT HVAR

// HELPER FUNCTIONS
template<class Type>
Type erf(Type x){
  Type y = sqrt(Type(2.0)) * x;
  Type z = Type(2.0) * pnorm(y) - Type(1.0);
  return z;
}

template<class Type>
Type sigmoid_fun(Type r_squared, Type loss_c){
  Type x = 1/(1+exp(-k_smooth * (sqrt(r_squared) - loss_c)));
  return(x);
}

template<class Type>
Type lossfunction__(Type r_squared, Type loss_c, std::string loss_type){
  // Quadratic loss
  Type loss = r_squared;
  if (loss_type == "huber"){
    Type s = sigmoid_fun(r_squared, loss_c);
    loss = r_squared * (1.0 - s) + loss_c * (2.0 * sqrt(r_squared) - loss_c) * s;
  }
  if(loss_type == "tukey"){
    Type s = sigmoid_fun(r_squared, loss_c);
    loss = r_squared * (1.0 - s) + loss_c * loss_c * s;
  }
  return loss;
}

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
}

template<class Type>
vector<Type> remove_nas__(vector<Type> obsVec, int n_available_obs, vector<Type> is_not_na_vector){
  int ii = 0;
  vector<Type> y_reduced(n_available_obs);
  for(int i=0; i < obsVec.size(); i++){
    if(is_not_na_vector(i) == Type(1.0)){
      y_reduced(ii) = obsVec(i);
      ii++;
    }
  }
  return y_reduced;
}

template <class Type>
matrix<Type> construct_permutation_matrix(int n_available_obs, int n_obs, vector<Type> is_not_na_vector){
  matrix<Type> E(n_available_obs, n_obs);
  E.setZero();
  int j=0;
  for(int i=0; i < n_obs; i++){
  /*if p(i) is 1 then include by setting 1 in diagonal of matrix*/
    if(is_not_na_vector(i) == Type(1.0)){
      E(j,i) = Type(1.0);
      j += 1;
    }
  }
  return E;
}

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
}

///////////// UKF sigma points ///////////
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
}

//////////// UKF function ///////////
template<class Type>
matrix<Type> Phi__(matrix<Type> M){
  matrix<Type> K(M.col(0).size(), M.row(0).size());
  K.setZero();
  K.template triangularView<Eigen::Lower>() = M.template triangularView<Eigen::Lower>();
  K.diagonal() = K.diagonal()/Type(2.0);
  return K;
}

//////////// UKF sigma points drift function ///////////
template<class Type>
matrix<Type> ukf_f__(matrix<Type> Xsigmapoints, vector<Type> parVec, vector<Type> inputVec, int n, int nn){
  matrix<Type> F(n,nn);
  vector<Type> stateVec;
  for(int i=0; i < nn; i++){
    stateVec = Xsigmapoints.col(i);
    F.col(i) = f__(stateVec, parVec, inputVec);
  }
  return F;
}

//////////// UKF sigma points obs function ///////////
template<class Type>
matrix<Type> ukf_h__(matrix<Type> Xsigmapoints, vector<Type> parVec, vector<Type> inputVec, int m, int nn){
  matrix<Type> H(m,nn);
  vector<Type> stateVec;
  for(int i=0; i < nn; i++){
    stateVec = Xsigmapoints.col(i);
    H.col(i) = h__(stateVec, parVec, inputVec);
  }
  return H;
}

//////////// 1-step f moment ODE ///////////
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
}

//////////// ODE SOLVER ///////////
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
};

////// Newton Solver for Stationary Mean //////
template<class Type>
struct Functor {
  vector<Type> p; // parVec
  vector<Type> u; // Inputs
  Functor(const vector<Type> &p, const vector<Type> &u) : p(p), u(u) {}
  Type operator()(const vector<Type> &s) {
    vector<Type> f = f__(s, p, u);
    return((f * f).sum());
  }
};

////// Lyapunov Solve for Stationary Variance (Linearized System) ///////
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
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  Type nll__ = 0;

  //// observations ////
  DATA_MATRIX(obsMat);

  //// inputs ////
  DATA_MATRIX(inputMat);

  //// initial state ////
  DATA_VECTOR(stateVec);
  DATA_MATRIX(covMat);

  //// newton settings ////
  DATA_STRUCT(cfg, newton::newton_config_t);

  //// ode settings ////
  DATA_VECTOR(ode_timestep_size);
  DATA_IVECTOR(ode_timesteps);
  DATA_INTEGER(ode_solver);

  //// loss parameters ////
  DATA_STRING(loss_type)
  DATA_SCALAR(loss_c);

  //// unscented transform hyper-parameters ////
  DATA_VECTOR(ukf_pars);

  //// map estimation ////
  DATA_INTEGER(MAP_bool);

  //// parameters ////
  PARAMETER_VECTOR(parVec);

  //// system size ////
  DATA_INTEGER(n_states);
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_inputs);
  int tsize = inputMat.col(0).size();
  int nn = 2*n_states + 1;

  //// state, par, input, obs vectors ////
  vector<Type> inputVec(n_inputs), dinputVec(n_inputs), obsVec(n_obs);

  //////////// storage variables ///////////
  vector<vector<Type>> xPrior(tsize), xPost(tsize), Innovation(tsize);
  vector<matrix<Type>> pPrior(tsize), pPost(tsize), InnovationCovariance(tsize);
  matrix<Type> XSigmaPrior(tsize, n_states * nn);
  matrix<Type> XSigmaPost(tsize, n_states * nn);

   //////////// initialize variables ///////////
  int n_available_obs;
  Type half_log2PI = Type(0.5) * log(2*M_PI);
  vector<Type> data_vector__(n_obs);
  vector<Type> y__,e__, is_not_na_obsVec;
  matrix<Type> E__, H__, Syy__, SyyInv__, Sxy__, K__, Xsigmapoints, sqrt_covMat;
  matrix<Type> V0__(n_obs,n_obs);
  V0__.setZero();

  //////////// create weights ///////////
  Type ukf_lambda = pow(ukf_pars(0),2)*(n_states + ukf_pars(2)) - n_states;
  Type ukf_weights = Type(1.0)/(Type(2.0)*(n_states + ukf_lambda));
  vector<Type> wm__(nn), wmC__(nn);
  matrix<Type> Wm__, WcDiag__(nn,nn), I__(nn,nn), W__;
  I__.setIdentity();
  WcDiag__.setZero();
  wm__.fill(ukf_weights);
  wmC__.fill(ukf_weights);
  wm__(0) = ukf_lambda/(ukf_lambda + n_states);
  wmC__(0) = ukf_lambda/((n_states + ukf_lambda) + (1-pow(ukf_pars(0),2)+ukf_pars(1)));
  Wm__ = wm__.replicate(1, nn);
  WcDiag__.diagonal() = wmC__;
  W__ = (I__ - Wm__) * WcDiag__ * (I__ - Wm__).transpose();

  //////////// stationary solution ///////////
  DATA_INTEGER(estimate_stationary_initials);
  if(estimate_stationary_initials == 1){
    inputVec = inputMat.row(0);
    Functor<TMBad::ad_aug> F(parVec, inputVec);
    stateVec = newton::Newton(F, stateVec, cfg);
    covMat = LyapSolver(stateVec, parVec, inputVec);
  };

  //////////// set initial value ///////////
  xPrior(0) = stateVec;
  pPrior(0) = covMat;

   //////////// FIRST DATA-UPDATE ///////////
  sqrt_covMat = covMat.llt().matrixL();
  Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, n_states, nn);
  XSigmaPrior.row(0) = Xsigmapoints.reshaped();
  obsVec = obsMat.row(0);
  is_not_na_obsVec = is_not_na(obsVec);
  n_available_obs = CppAD::Integer(sum(is_not_na_obsVec));
  if( n_available_obs > 0 ){
    inputVec = inputMat.row(0);
    y__ = remove_nas__(obsVec, n_available_obs, is_not_na_obsVec);
    E__ = construct_permutation_matrix(n_available_obs, n_obs, is_not_na_obsVec);
    H__ = ukf_h__(Xsigmapoints, parVec, inputVec, n_obs, nn);
    e__  = y__ - E__ * (H__ * wm__);
    V0__.diagonal() << hvar__(stateVec, parVec, inputVec);
    Syy__  = E__ * (H__ * W__ * H__.transpose() + V0__) * E__.transpose();
    SyyInv__  = Syy__.inverse();
    Sxy__  = Xsigmapoints * W__ * H__.transpose() * E__.transpose();
    K__ = Sxy__ * SyyInv__;
    // State update // 
    stateVec = stateVec + K__ * e__;
    covMat = covMat - K__ * Syy__ * K__.transpose();
    // Likelihood contribution // 
    nll__ += Type(0.5) * atomic::logdet(Syy__) + Type(0.5) * lossfunction__((e__*(SyyInv__*e__)).sum(), loss_c, loss_type) + half_log2PI * asDouble(n_available_obs);
    // Save residual mean/covariance // 
    Innovation(0) = e__;
    InnovationCovariance(0) = Syy__;
  };
  xPost(0) = stateVec;
  pPost(0) = covMat;

   //////////// START MAIN LOOP ///////////
  for(int i=0 ; i < tsize - 1 ; i++){

    //////////// cholesky cov and sigma points ///////////
    sqrt_covMat = covMat.llt().matrixL();
    Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, n_states, nn);
    XSigmaPost.row(i+1) = Xsigmapoints.reshaped();
    inputVec = inputMat.row(i);
    dinputVec = (inputMat.row(i+1) - inputMat.row(i))/ode_timesteps(i);

   //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
    for(int j=0 ; j < ode_timesteps(i) ; j++){
      ukf_ode_integration<Type> odelist = {Xsigmapoints, sqrt_covMat, W__, wm__, parVec, inputVec, dinputVec, ode_timestep_size(i), ode_solver, n_states, nn};
      Xsigmapoints = odelist.X1sigmapoints;
      sqrt_covMat = ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt(Type(3.0))).block(0, 1, n_states, n_states );
      inputVec += dinputVec;
    }
    stateVec = Xsigmapoints.col(0);
    covMat = sqrt_covMat * sqrt_covMat.transpose();
    xPrior(i+1) = stateVec;
    pPrior(i+1) = covMat;
    XSigmaPrior.row(i+1) = Xsigmapoints.reshaped();

   //////////// DATA-UPDATE ///////////
    obsVec = obsMat.row(i+1);
    is_not_na_obsVec = is_not_na(obsVec);
    n_available_obs = CppAD::Integer(sum(is_not_na_obsVec));
    if( n_available_obs > 0 ){
      inputVec = inputMat.row(i+1);
      y__ = remove_nas__(obsVec, n_available_obs, is_not_na_obsVec);
      E__ = construct_permutation_matrix(n_available_obs, n_obs, is_not_na_obsVec);
      H__ = ukf_h__(Xsigmapoints, parVec, inputVec, n_obs, nn);
      e__  = y__ - E__ * (H__ * wm__);
      V0__.diagonal() << hvar__(stateVec, parVec, inputVec);
      Syy__  = E__ * (H__ * W__ * H__.transpose() + V0__) * E__.transpose();
      SyyInv__  = Syy__.inverse();
      Sxy__  = Xsigmapoints * W__ * H__.transpose() * E__.transpose();
      K__ = Sxy__ * SyyInv__;
      stateVec = stateVec + K__ * e__;
      covMat = covMat - K__ * Syy__ * K__.transpose();
      nll__ += Type(0.5) * atomic::logdet(Syy__) + Type(0.5) * lossfunction__((e__*(SyyInv__*e__)).sum(), loss_c, loss_type) + half_log2PI * asDouble(n_available_obs);
      Innovation(i+1) = e__;
      InnovationCovariance(i+1) = Syy__;
    };
    xPost(i+1) = stateVec;
    pPost(i+1) = covMat;
  };
  //////////// END MAIN LOOP ///////////

  //////////// MAP CONTRIBUTION ///////////
  if(MAP_bool == 1){
    DATA_VECTOR(map_mean__);
    DATA_MATRIX(map_cov__);
    DATA_IVECTOR(map_ints__);
    DATA_INTEGER(sum_map_ints__);
    vector<Type> map_pars__;
    map_pars__ = get_free_pars__(map_ints__, sum_map_ints__, parVec);
    vector<Type> pars_eps__ = map_pars__ - map_mean__;
    matrix<Type> map_invcov__ = map_cov__.inverse();
    Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();
    nll__ += map_nll__;
    REPORT(map_nll__);
    REPORT(map_pars__);
    REPORT(pars_eps__);
  }

  //////////// Report //////////////
  REPORT(Innovation);
  REPORT(InnovationCovariance);
  REPORT(xPrior);
  REPORT(xPost);
  REPORT(pPrior);
  REPORT(pPost);
  REPORT(XSigmaPrior);
  REPORT(XSigmaPost);
  return nll__;
}
