#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_funs2.h"
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

//  This is the predict kalman filter function
template <typename T1, typename T2>
List lkf_filter(
  T1 f__, 
  T2 g__,
  T2 dfdx__,
  T1 h__,
  T2 dhdx__,
  T2 hvar__,
  T2 dfdu__,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs)
{

  // constants  
  const int tsize = inputMat.col(0).size();
  const int ni = inputMat.row(0).size();
  const int n = stateVec.size();
  const int m = obsMat.row(0).size();

  // pre-allocate and define
  Eigen::VectorXd inputVec(ni), dinputVec(ni), obsVec(m), e, y, obsVec_all;
  Eigen::VectorXi bool_is_not_na_obsVec;
  Eigen::MatrixXd C, R, K, E, V, Ri, I(n,n);
  I.setIdentity();
  Eigen::VectorXd H;
  Eigen::MatrixXd Hvar, dHdX;

  // create matrix exponential
  inputVec = inputMat.row(0);
  double dt_timestep = ode_timestep_size(0);
  Eigen::MatrixXd A = dfdx__(stateVec, parVec, inputVec);
  Eigen::MatrixXd B = dfdu__(stateVec, parVec, inputVec);
  Eigen::MatrixXd G = g__(stateVec, parVec, inputVec);
  // # [A B \\ 0 0]
  Eigen::MatrixXd Phi1(n+ni+1, n+ni+1);
  Phi1.setZero();
  Phi1.block(0, 0, n, n) = A;
  Phi1.block(0, n, n, ni+1) = B;
  // # [-A GG^T \\ 0 A^t]
  Eigen::MatrixXd Phi2(2*n,2*n);
  Phi2.setZero();
  Phi2.block(0, 0, n, n)  = -A;
  Phi2.block(0, n, n, n)  = G*G.transpose();
  Phi2.block(n, n, n, n)  = A.transpose();
  // Calculate matrix exponentials
  Eigen::MatrixXd ePhi1 = (Phi1 * dt_timestep).exp();
  Eigen::MatrixXd ePhi2 = (Phi2 * dt_timestep).exp();
  // 
  Eigen::MatrixXd Ahat = ePhi1.block(0, 0, n, n);
  Eigen::MatrixXd Ahat_T = Ahat.transpose();
  Eigen::MatrixXd Bhat = ePhi1.block(0, n, n, ni+1);
  Eigen::MatrixXd Q12 = ePhi2.block(0, n, n, n);
  Eigen::MatrixXd Q22 = ePhi2.block(n, n, n ,n);
  Eigen::MatrixXd Vhat = Q22.transpose() * Q12;

  Eigen::VectorXd constant_and_inputVec(ni+1);
  constant_and_inputVec(0) = 1.0;

  // storage
  Rcpp::List ode_1step_integration(2);
  Rcpp::List Innovation(tsize), InnovationCovariance(tsize), xPrior(tsize), xPost(tsize), pPrior(tsize), pPost(tsize);

  // store prior
  xPrior(0) = stateVec;
  pPrior(0) = covMat;

  //////////// INITIAL DATA-UPDATE ///////////
  if( number_of_available_obs(0) > 0 ){
    inputVec = inputMat.row(0);
    obsVec_all = obsMat.row(0);
    bool_is_not_na_obsVec = bool_is_not_na_obsMat.row(0);
    // remove potential NA entries in obsVec and construct permutation matrix
    obsVec = remove_NAs2(obsVec_all, number_of_available_obs(0), bool_is_not_na_obsVec);
    E = construct_permutation_matrix2(number_of_available_obs(0), m, bool_is_not_na_obsVec);
    // Call funs
    H = h__(stateVec, parVec, inputVec);
    dHdX = dhdx__(stateVec, parVec, inputVec);
    Hvar = hvar__(stateVec, parVec, inputVec);
    // Kalman Filter
    C = E * dHdX;
    e = obsVec - E * H;
    V = E * Hvar * E.transpose();
    R = C * covMat * C.transpose() + V;
    Ri = R.inverse();
    K = covMat * C.transpose() * Ri;
    stateVec = stateVec + K * e;
    covMat = (I - K * C) * covMat * (I - K * C).transpose() + K*V*K.transpose();
    // Store innovations
    Innovation(0) = e;
    InnovationCovariance(0) = R;
  }

  // store posterior
  xPost(0) = stateVec;
  pPost(0) = covMat;

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < (tsize-1) ; i++){

    //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
    constant_and_inputVec.tail(ni) = inputMat.row(i);
    stateVec = Ahat * stateVec + Bhat * constant_and_inputVec;
    covMat = Ahat * covMat * Ahat_T + Vhat;
    xPrior(i+1) = stateVec;
    pPrior(i+1) = covMat;

    //////////// DATA-UPDATE ///////////
    if( number_of_available_obs(i+1) > 0 ){
      inputVec = inputMat.row(i+1);
      obsVec_all = obsMat.row(i+1);
      bool_is_not_na_obsVec = bool_is_not_na_obsMat.row(i+1);
      // remove potential NA entries in obsVec and construct permutation matrix
      obsVec = remove_NAs2(obsVec_all, number_of_available_obs(i+1), bool_is_not_na_obsVec);
      E = construct_permutation_matrix2(number_of_available_obs(i+1), m, bool_is_not_na_obsVec);
      // Call funs
      H = h__(stateVec, parVec, inputVec);
      dHdX = dhdx__(stateVec, parVec, inputVec);
      Hvar = hvar__(stateVec, parVec, inputVec);
      // Kalman Filter
      C = E * dHdX;
      e = obsVec - E * H;
      V = E * Hvar * E.transpose();
      R = C * covMat * C.transpose() + V;
      Ri = R.inverse();
      K = covMat * C.transpose() * Ri;
      stateVec = stateVec + K * e;
      covMat = (I - K * C) * covMat * (I - K * C).transpose() + K*V*K.transpose();
      // Store innovations
      Innovation(i+1) = e;
      InnovationCovariance(i+1) = R;
    }
    xPost(i+1) = stateVec;
    pPost(i+1) = covMat;
  }
  
  return List::create(
    Named("xPrior") = xPrior,
    Named("xPost") = xPost,
    Named("pPrior") = pPrior,
    Named("pPost") = pPost,
    Named("Innovation") = Innovation,
    Named("InnovationCovariance") = InnovationCovariance,
    Named("Ahat") = Ahat,
    Named("Bhat") = Bhat,
    Named("Vhat") = Vhat
    );
}

// Function typedefs
typedef Eigen::VectorXd (*funPtr_vec)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
typedef Eigen::MatrixXd (*funPtr_mat)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

// Function exported to R that performs linear kalman filtering
// [[Rcpp::export]]
List lkf_filter_rcpp(
  SEXP f__R, 
  SEXP g__R,
  SEXP dfdx__R,
  SEXP h__R,
  SEXP dhdx__R,
  SEXP hvar__R,
  SEXP dfdu__R,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs)
{

  funPtr_vec     f__ = *XPtr<funPtr_vec>(f__R);
  funPtr_mat     g__ = *XPtr<funPtr_mat>(g__R);
  funPtr_mat  dfdx__ = *XPtr<funPtr_mat>(dfdx__R);
  funPtr_vec     h__ = *XPtr<funPtr_vec>(h__R);
  funPtr_mat  dhdx__ = *XPtr<funPtr_mat>(dhdx__R);
  funPtr_mat  hvar__ = *XPtr<funPtr_mat>(hvar__R);
  funPtr_mat  dfdu__ = *XPtr<funPtr_mat>(dfdu__R);

  return lkf_filter<funPtr_vec, funPtr_mat>(
    f__, 
    g__, 
    dfdx__, 
    h__, 
    dhdx__,
    hvar__, 
    dfdu__,
    obsMat,
    inputMat,
    parVec,
    covMat, 
    stateVec,
    ode_timestep_size,
    bool_is_not_na_obsMat,
    number_of_available_obs);
}

//  This is the predict kalman filter function
template <typename T1, typename T2>
List lkf_predict(
  T1 f__, 
  T2 g__,
  T2 dfdx__,
  T1 h__,
  T2 dhdx__,
  T2 hvar__,
  T2 dfdu__,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  const int last_pred_id,
  const int k_step_ahead)
{

  Rcpp::List filt = lkf_filter(f__, 
                               g__,
                               dfdx__,
                               h__,
                               dhdx__,
                               hvar__,
                               dfdu__,
                               obsMat,
                               inputMat,
                               parVec,
                               covMat, 
                               stateVec,
                               ode_timestep_size,
                               bool_is_not_na_obsMat,
                               number_of_available_obs);
  Rcpp::List xPost = filt["xPost"];
  Rcpp::List pPost = filt["pPost"];
  Eigen::MatrixXd Ahat = filt["Ahat"];
  Eigen::MatrixXd Ahat_T = Ahat.transpose();
  Eigen::MatrixXd Bhat = filt["Bhat"];
  Eigen::MatrixXd Vhat = filt["Vhat"];

  // constants 
  const int ni = inputMat.row(0).size();
  const int n = stateVec.size();
  const int n_squared = n*n;

  // pre-allocate and define
  Eigen::VectorXd inputVec(ni), dinputVec(ni);
  Eigen::MatrixXd predMat(k_step_ahead+1, n + n_squared);
  predMat.setZero();
  Rcpp::List xk(last_pred_id);
  Eigen::VectorXd constant_and_inputVec(ni+1);
  constant_and_inputVec(0) = 1.0;

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){
    stateVec = xPost[i];
    covMat = pPost[i];
    predMat.row(0).head(n) = stateVec;
    predMat.row(0).tail(n_squared) = covMat.reshaped();

    for(int k=0 ; k < k_step_ahead ; k++){
      constant_and_inputVec.tail(ni) = inputMat.row(i+k);
      stateVec = Ahat * stateVec + Bhat * constant_and_inputVec;
      covMat = Ahat * covMat * Ahat_T + Vhat;
      predMat.row(k+1).head(n) = stateVec;
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

// Function exported to R that performs linear kalman filtering
// [[Rcpp::export]]
List lkf_predict_rcpp(
  SEXP f__R, 
  SEXP g__R,
  SEXP dfdx__R,
  SEXP h__R,
  SEXP dhdx__R,
  SEXP hvar__R,
  SEXP dfdu__R,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  const int last_pred_id,
  const int k_step_ahead)
{

  funPtr_vec     f__ = *XPtr<funPtr_vec>(f__R);
  funPtr_mat     g__ = *XPtr<funPtr_mat>(g__R);
  funPtr_mat  dfdx__ = *XPtr<funPtr_mat>(dfdx__R);
  funPtr_vec     h__ = *XPtr<funPtr_vec>(h__R);
  funPtr_mat  dhdx__ = *XPtr<funPtr_mat>(dhdx__R);
  funPtr_mat  hvar__ = *XPtr<funPtr_mat>(hvar__R);
  funPtr_mat  dfdu__ = *XPtr<funPtr_mat>(dfdu__R);

  return lkf_predict<funPtr_vec, funPtr_mat>(
    f__, 
    g__, 
    dfdx__, 
    h__, 
    dhdx__,
    hvar__, 
    dfdu__,
    obsMat,
    inputMat,
    parVec,
    covMat, 
    stateVec,
    ode_timestep_size,
    bool_is_not_na_obsMat,
    number_of_available_obs,
    last_pred_id,
    k_step_ahead);
}

template <typename T1, typename T2>
List lkf_simulate(
  T1 f__, 
  T2 g__,
  T2 dfdx__,
  T1 h__,
  T2 dhdx__,
  T2 hvar__,
  T2 dfdu__,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd simulation_timestep_size,
  Eigen::VectorXd simulation_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  int ng,
  int last_pred_id,
  int k_step_ahead,
  int nsims)
{

  Rcpp::List filt = lkf_filter(f__, 
                               g__,
                               dfdx__,
                               h__,
                               dhdx__,
                               hvar__,
                               dfdu__,
                               obsMat,
                               inputMat,
                               parVec,
                               covMat, 
                               stateVec,
                               ode_timestep_size,
                               bool_is_not_na_obsMat,
                               number_of_available_obs);

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

// Function exported to R that performs linear kalman filtering
// [[Rcpp::export]]
List lkf_simulate_rcpp(
  SEXP f__R, 
  SEXP g__R,
  SEXP dfdx__R,
  SEXP h__R,
  SEXP dhdx__R,
  SEXP hvar__R,
  SEXP dfdu__R,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd simulation_timestep_size,
  Eigen::VectorXd simulation_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  int ng,
  int last_pred_id,
  int k_step_ahead,
  int nsims)
{

  funPtr_vec     f__ = *XPtr<funPtr_vec>(f__R);
  funPtr_mat     g__ = *XPtr<funPtr_mat>(g__R);
  funPtr_mat  dfdx__ = *XPtr<funPtr_mat>(dfdx__R);
  funPtr_vec     h__ = *XPtr<funPtr_vec>(h__R);
  funPtr_mat  dhdx__ = *XPtr<funPtr_mat>(dhdx__R);
  funPtr_mat  hvar__ = *XPtr<funPtr_mat>(hvar__R);
  funPtr_mat  dfdu__ = *XPtr<funPtr_mat>(dfdu__R);

  return lkf_simulate<funPtr_vec, funPtr_mat>(
    f__, 
    g__, 
    dfdx__, 
    h__, 
    dhdx__,
    hvar__, 
    dfdu__,
    obsMat,
    inputMat,
    parVec,
    covMat, 
    stateVec,
    ode_timestep_size,
    simulation_timestep_size,
    simulation_timesteps,
    bool_is_not_na_obsMat,
    number_of_available_obs,
    ng,
    last_pred_id,
    k_step_ahead,
    nsims);
}