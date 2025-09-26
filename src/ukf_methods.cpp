#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_funs2.h"
#include "helpers_ukf.h"
using namespace Rcpp;
using namespace Eigen;

template <typename T1, typename T2>
List ukf_filter(
  T1 f__, 
  T2 g__,
  T2 dfdx__,
  T1 h__,
  T2 dhdx__,
  T2 hvar__,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd ode_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  Eigen::VectorXd ukf_pars,
  const int ode_solver)
{

  // constants  
  const int tsize = inputMat.col(0).size();
  const int n_inputs = inputMat.row(0).size();
  const int n_states = stateVec.size();
  const int nn = 2*n_states + 1;
  const int n_obs = obsMat.row(0).size();

  // pre-allocate and define
  VectorXd inputVec(n_inputs), dinputVec(n_inputs), obsVec(n_obs), e, y, obsVec_all;
  VectorXi bool_is_not_na_obsVec;
  MatrixXd Xsigmapoints, sqrt_covMat;;
  MatrixXd E, R, Ri, Sxy, K, V, H, Hvar;
  int n_available_obs;

  //////////// create weights ///////////
  double ukf_alpha = ukf_pars(0);
  double ukf_beta = ukf_pars(1);
  double ukf_kappa = ukf_pars(2);
  VectorXd wm(nn), wmC(nn);
  MatrixXd Wm, WcDiag(nn,nn), I(nn,nn), W;
  I.setIdentity();
  WcDiag.setZero();
  double sqrt_c = sqrt(pow(ukf_alpha,2)*(n_states + ukf_kappa));
  double ukf_lambda = pow(sqrt_c, 2) - n_states;
  double ukf_weights = 1.0/(2.0*(n_states + ukf_lambda));
  wm.fill(ukf_weights);
  wmC.fill(ukf_weights);
  wm(0) = ukf_lambda/(ukf_lambda + n_states);
  wmC(0) = ukf_lambda/((n_states + ukf_lambda) + (1-pow(ukf_alpha,2)+ukf_beta));
  Wm = wm.replicate(1, nn);
  WcDiag.diagonal() = wmC;
  W = (I - Wm) * WcDiag * (I - Wm).transpose();

  // storage
  Rcpp::List Innovation(tsize), InnovationCovariance(tsize), xPrior(tsize), xPost(tsize), pPrior(tsize), pPost(tsize);

  // store prior
  xPrior(0) = stateVec;
  pPrior(0) = covMat;

   //////////// FIRST DATA-UPDATE ///////////
  sqrt_covMat = covMat.llt().matrixL();
  Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, sqrt_c, n_states, nn);
  n_available_obs = number_of_available_obs(0);
  if( n_available_obs > 0 ){
    inputVec = inputMat.row(0);
    obsVec_all = obsMat.row(0);
    bool_is_not_na_obsVec = bool_is_not_na_obsMat.row(0);
    // remove potential NA entries in obsVec and construct permutation matrix
    obsVec = remove_NAs2(obsVec_all, n_available_obs, bool_is_not_na_obsVec);
    E = construct_permutation_matrix2(n_available_obs, n_obs, bool_is_not_na_obsVec);
    // Call funs
    H = ukf_h(h__, Xsigmapoints, parVec, inputVec, n_obs, nn);
    Hvar = hvar__(stateVec, parVec, inputVec);
    // Kalman Filter
    e  = obsVec - E * H * wm;
    V = E * Hvar * E.transpose();
    R  = E * H * W * H.transpose() * E.transpose() + V;
    Ri  = R.inverse();
    K = Xsigmapoints * W * H.transpose() * E.transpose() * Ri;
    stateVec = stateVec + K * e;
    covMat = covMat - K * R * K.transpose();
    // Store innovations
    Innovation(0) = e;
    InnovationCovariance(0) = R;
  };

  // store posterior
  xPost(0) = stateVec;
  pPost(0) = covMat;

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < (tsize-1) ; i++){

    sqrt_covMat = covMat.llt().matrixL();
    Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, sqrt_c, n_states, nn);
    inputVec = inputMat.row(i);
    dinputVec = (inputMat.row(i+1) - inputMat.row(i))/ode_timesteps(i);

    //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
    for(int j=0 ; j < ode_timesteps(i) ; j++){
      Xsigmapoints = ukf_ode_integration(f__, g__, Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, dinputVec, ode_timestep_size(i), sqrt_c, ode_solver, n_states, nn);
      sqrt_covMat = sigma2chol(Xsigmapoints, sqrt_c, nn, n_states);  
      // sqrt_covMat = ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt(3.0)).block(0, 1, n_states, n_states );
      inputVec += dinputVec;
    }
    stateVec = Xsigmapoints.col(0);
    covMat = sqrt_covMat * sqrt_covMat.transpose();
    xPrior(i+1) = stateVec;
    pPrior(i+1) = covMat;

    sqrt_covMat = covMat.llt().matrixL();
    Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, sqrt_c, n_states, nn);
    n_available_obs = number_of_available_obs(i+1);
    if( n_available_obs > 0 ){
      inputVec = inputMat.row(i+1);
      obsVec_all = obsMat.row(i+1);
      bool_is_not_na_obsVec = bool_is_not_na_obsMat.row(i+1);
      // remove potential NA entries in obsVec and construct permutation matrix
      obsVec = remove_NAs2(obsVec_all, n_available_obs, bool_is_not_na_obsVec);
      E = construct_permutation_matrix2(n_available_obs, n_obs, bool_is_not_na_obsVec);
      // Call funs
      H = ukf_h(h__, Xsigmapoints, parVec, inputVec, n_obs, nn);
      Hvar = hvar__(stateVec, parVec, inputVec);
      // Kalman Filter
      e  = obsVec - E * H * wm;
      V = E * Hvar * E.transpose();
      R  = E * H * W * H.transpose() * E.transpose() + V;
      Ri  = R.inverse();
      K = Xsigmapoints * W * H.transpose() * E.transpose() * Ri;
      stateVec = stateVec + K * e;
      covMat = covMat - K * R * K.transpose();
      // Store innovations
      Innovation(i+1) = e;
      InnovationCovariance(i+1) = R;
    };

    xPost(i+1) = stateVec;
    pPost(i+1) = covMat;
  }
  
  return List::create(
    Named("xPrior") = xPrior,
    Named("xPost") = xPost,
    Named("pPrior") = pPrior,
    Named("pPost") = pPost,
    Named("Innovation") = Innovation,
    Named("InnovationCovariance") = InnovationCovariance
    );
}

// Function typedefs
typedef Eigen::VectorXd (*funPtr_vec)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
typedef Eigen::MatrixXd (*funPtr_mat)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

// Function exported to R that performs (Extended Kalman) filtering
// [[Rcpp::export]]
List ukf_filter_rcpp(
  SEXP f__R, 
  SEXP g__R,
  SEXP dfdx__R,
  SEXP h__R,
  SEXP dhdx__R,
  SEXP hvar__R,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd ode_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  Eigen::VectorXd ukf_pars,
  int ode_solver)
{

  funPtr_vec     f__ = *XPtr<funPtr_vec>(f__R);
  funPtr_mat     g__ = *XPtr<funPtr_mat>(g__R);
  funPtr_mat  dfdx__ = *XPtr<funPtr_mat>(dfdx__R);
  funPtr_vec     h__ = *XPtr<funPtr_vec>(h__R);
  funPtr_mat  dhdx__ = *XPtr<funPtr_mat>(dhdx__R);
  funPtr_mat  hvar__ = *XPtr<funPtr_mat>(hvar__R);

  return ukf_filter<funPtr_vec, funPtr_mat>(
    f__, 
    g__, 
    dfdx__, 
    h__, 
    dhdx__,
    hvar__,
    obsMat,
    inputMat,
    parVec,
    covMat, 
    stateVec,
    ode_timestep_size,
    ode_timesteps,
    bool_is_not_na_obsMat,
    number_of_available_obs,
    ukf_pars,
    ode_solver);

}

template <typename T1, typename T2>
List ukf_predict(
  T1 f__, 
  T2 g__,
  T2 dfdx__,
  T1 h__,
  T2 dhdx__,
  T2 hvar__,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd ode_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  Eigen::VectorXd ukf_pars,
  const int last_pred_id,
  const int k_step_ahead,
  const int ode_solver)
{

  Rcpp::List filt = ukf_filter(f__, 
                              g__, 
                              dfdx__, 
                              h__, 
                              dhdx__,
                              hvar__,
                              obsMat,
                              inputMat,
                              parVec,
                              covMat, 
                              stateVec,
                              ode_timestep_size,
                              ode_timesteps,
                              bool_is_not_na_obsMat,
                              number_of_available_obs,
                              ukf_pars,
                              ode_solver);

  Rcpp::List xPost = filt["xPost"];
  Rcpp::List pPost = filt["pPost"];

  // constants  
  const int n_states = stateVec.size();
  const int n_inputs = inputMat.row(0).size();
  const int nn = 2*n_states + 1;
  const int n_squared = n_states*n_states;

  // pre-allocate and define
  VectorXd inputVec(n_inputs), dinputVec(n_inputs);
  MatrixXd Xsigmapoints, sqrt_covMat;;
  Eigen::MatrixXd predMat(k_step_ahead+1, n_states + n_squared);
  predMat.setZero();
  Rcpp::List xk(last_pred_id);

  //////////// create weights ///////////
  double ukf_alpha = ukf_pars(0);
  double ukf_beta = ukf_pars(1);
  double ukf_kappa = ukf_pars(2);
  VectorXd wm(nn), wmC(nn);
  MatrixXd Wm, WcDiag(nn,nn), I(nn,nn), W;
  I.setIdentity();
  WcDiag.setZero();
  double sqrt_c = sqrt(pow(ukf_alpha,2)*(n_states + ukf_kappa));
  double ukf_lambda = pow(sqrt_c, 2) - n_states;
  double ukf_weights = 1.0/(2.0*(n_states + ukf_lambda));
  wm.fill(ukf_weights);
  wmC.fill(ukf_weights);
  wm(0) = ukf_lambda/(ukf_lambda + n_states);
  wmC(0) = ukf_lambda/((n_states + ukf_lambda) + (1-pow(ukf_alpha,2)+ukf_beta));
  Wm = wm.replicate(1, nn);
  WcDiag.diagonal() = wmC;
  W = (I - Wm) * WcDiag * (I - Wm).transpose();

   //////////// FIRST DATA-UPDATE ///////////
  sqrt_covMat = covMat.llt().matrixL();
  Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, sqrt_c, n_states, nn);
 
  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){
    stateVec = xPost[i];
    covMat = pPost[i];
    predMat.row(0).head(n_states) = stateVec;
    predMat.row(0).tail(n_squared) = covMat.reshaped();

    sqrt_covMat = covMat.llt().matrixL();
    Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, sqrt_c, n_states, nn);
    for(int k=0 ; k < k_step_ahead ; k++){
      inputVec = inputMat.row(i+k);
      dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k))/ode_timesteps(i+k);

      //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
      for(int j=0 ; j < ode_timesteps(i+k) ; j++){
        Xsigmapoints = ukf_ode_integration(f__, g__, Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, dinputVec, ode_timestep_size(i+k), sqrt_c, ode_solver, n_states, nn);
        sqrt_covMat = sigma2chol(Xsigmapoints, sqrt_c, nn, n_states);  
        inputVec += dinputVec;
      }

      // Save the prediction
      stateVec = Xsigmapoints.col(0);
      covMat = sqrt_covMat * sqrt_covMat.transpose();
      predMat.row(k+1).head(n_states) = stateVec;
      predMat.row(k+1).tail(n_squared) = covMat.reshaped();
    }

    xk[i] = predMat;

  }
  
  return List::create(
    Named("predMats") = xk
    );
}

// Function typedefs
typedef Eigen::VectorXd (*funPtr_vec)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
typedef Eigen::MatrixXd (*funPtr_mat)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

// Function exported to R that performs (Extended Kalman) filtering
// [[Rcpp::export]]
List ukf_predict_rcpp(
  SEXP f__R, 
  SEXP g__R,
  SEXP dfdx__R,
  SEXP h__R,
  SEXP dhdx__R,
  SEXP hvar__R,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd ode_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  Eigen::VectorXd ukf_pars,
  const int last_pred_id,
  const int k_step_ahead,
  const int ode_solver)
{

  funPtr_vec     f__ = *XPtr<funPtr_vec>(f__R);
  funPtr_mat     g__ = *XPtr<funPtr_mat>(g__R);
  funPtr_mat  dfdx__ = *XPtr<funPtr_mat>(dfdx__R);
  funPtr_vec     h__ = *XPtr<funPtr_vec>(h__R);
  funPtr_mat  dhdx__ = *XPtr<funPtr_mat>(dhdx__R);
  funPtr_mat  hvar__ = *XPtr<funPtr_mat>(hvar__R);

  return ukf_predict<funPtr_vec, funPtr_mat>(
    f__, 
    g__, 
    dfdx__, 
    h__, 
    dhdx__,
    hvar__,
    obsMat,
    inputMat,
    parVec,
    covMat, 
    stateVec,
    ode_timestep_size,
    ode_timesteps,
    bool_is_not_na_obsMat,
    number_of_available_obs,
    ukf_pars,
    last_pred_id,
    k_step_ahead,
    ode_solver);

}

template <typename T1, typename T2>
List ukf_simulate(
  T1 f__, 
  T2 g__,
  T2 dfdx__,
  T1 h__,
  T2 dhdx__,
  T2 hvar__,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd ode_timesteps,
  Eigen::VectorXd simulation_timestep_size,
  Eigen::VectorXd simulation_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  Eigen::VectorXd ukf_pars,
  const int ng,
  const int last_pred_id,
  const int k_step_ahead,
  const int ode_solver,
  const int nsims)
{

  Rcpp::List filt = ukf_filter(f__, 
                              g__, 
                              dfdx__, 
                              h__, 
                              dhdx__,
                              hvar__,
                              obsMat,
                              inputMat,
                              parVec,
                              covMat, 
                              stateVec,
                              ode_timestep_size,
                              ode_timesteps,
                              bool_is_not_na_obsMat,
                              number_of_available_obs,
                              ukf_pars,
                              ode_solver);

  Rcpp::List xPost = filt["xPost"];
  Rcpp::List pPost = filt["pPost"];

// misc  
  const int n = stateVec.size();
  const int ni = inputMat.row(0).size();
  VectorXd inputVec(ni), dinputVec(ni);
  MatrixXd stateMat(nsims, n), randN(n, nsims);

  // storage for predictions
  List xk_simulate(last_pred_id), xk_simulate_temp(k_step_ahead+1);
  
  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){

    stateVec = xPost[i];
    covMat = pPost[i];

    /* 
    We draw from a multivariate normal z = u + A*dB
    where u is the mean (stateVec) and A is cholesky factor of covariance matrix (sqrt(covMat))
    and dB is i.d.d normal vector
    We do simultaneously for all #nsims simulations, so dB is here a matrix of #nsims i.d.d vectors
    and similarly u is repeated with replicate
    */
    for(int j=0; j < n; j++){
      for(int k=0; k < nsims; k++){
        // randN(j,k) = zigg.norm();
        randN(j,k) = ziggurat.rnorm();
      }
    }
    stateMat = (stateVec.replicate(1, nsims) + covMat.llt().matrixL() * randN).transpose(); 
    xk_simulate_temp[0] = stateMat;

    for(int k=0 ; k < k_step_ahead ; k++){
      inputVec = inputMat.row(i+k);
      dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k))/simulation_timesteps(i+k);

      for(int j=0 ; j < simulation_timesteps(i+k) ; j++){
        stateMat = euler_maruyama_simulation2(f__, g__, stateMat, parVec, inputVec, simulation_timestep_size(i+k), nsims, n, ng);
        inputVec += dinputVec;
      }

      xk_simulate_temp[k+1] = stateMat;
    }

    // Save a clone of the temporary list (must use clone)
    xk_simulate[i] = clone(xk_simulate_temp);
  }

  // Return
  return xk_simulate;
}

// Function typedefs
typedef Eigen::VectorXd (*funPtr_vec)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
typedef Eigen::MatrixXd (*funPtr_mat)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

// Function exported to R that performs (Extended Kalman) filtering
// [[Rcpp::export]]
List ukf_simulate_rcpp(
  SEXP f__R, 
  SEXP g__R,
  SEXP dfdx__R,
  SEXP h__R,
  SEXP dhdx__R,
  SEXP hvar__R,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd ode_timesteps,
  Eigen::VectorXd simulation_timestep_size,
  Eigen::VectorXd simulation_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  Eigen::VectorXd ukf_pars,
  const int ng,
  const int last_pred_id,
  const int k_step_ahead,
  const int ode_solver,
  const int nsims)
{

  funPtr_vec     f__ = *XPtr<funPtr_vec>(f__R);
  funPtr_mat     g__ = *XPtr<funPtr_mat>(g__R);
  funPtr_mat  dfdx__ = *XPtr<funPtr_mat>(dfdx__R);
  funPtr_vec     h__ = *XPtr<funPtr_vec>(h__R);
  funPtr_mat  dhdx__ = *XPtr<funPtr_mat>(dhdx__R);
  funPtr_mat  hvar__ = *XPtr<funPtr_mat>(hvar__R);

  return ukf_simulate<funPtr_vec, funPtr_mat>(
    f__, 
    g__, 
    dfdx__, 
    h__, 
    dhdx__,
    hvar__,
    obsMat,
    inputMat,
    parVec,
    covMat, 
    stateVec,
    ode_timestep_size,
    ode_timesteps,
    simulation_timestep_size,
    simulation_timesteps,
    bool_is_not_na_obsMat,
    number_of_available_obs,
    ukf_pars,
    ng,
    last_pred_id,
    k_step_ahead,
    ode_solver,
    nsims);

}