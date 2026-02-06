#include <Rcpp.h>
#include <RcppEigen.h>

#include "function_typedefs.h"
#include "misc_helpers.h"
#include "sde_solvers.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
List lkf_filter_rcpp(
  List funPtrs,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  LogicalVector any_available_obs,
  List non_na_ids)
{

  auto h__ = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");
  auto dfdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dfdx_const");
  auto dhdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar__ = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");
  auto dfdu__ = get_funptr<funPtr_mat_const>(funPtrs, "dfdu_const");

  // constants  
  const int tsize = inputMat.col(0).size();
  const int ni = inputMat.row(0).size();
  const int n = stateVec.size();
  const int m = obsMat.row(0).size();
  int idx, n_available_obs;

  // pre-allocate and define
  Eigen::VectorXd H, inputVec(ni), dinputVec(ni), obsVec(m), e(m);
  Eigen::MatrixXd C, R, K, E, V, KC, IKC, Hvar, dHdX;
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
  Eigen::VectorXi obs_ids;

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
  if(any_available_obs(0)){
    inputVec = inputMat.row(0);
    obsVec = obsMat.row(0);
    obs_ids = Rcpp::as<Eigen::VectorXi>(non_na_ids(0));
    n_available_obs = obs_ids.size();
    // Calculcate H
    H = h__(stateVec, parVec, inputVec);
    dHdX = dhdx__(stateVec, parVec, inputVec);
    Hvar = hvar__(stateVec, parVec, inputVec);
    // Extract to reduce dimensions to fit number of observations
    C = dHdX.topRows(n_available_obs);
    V = Hvar.topLeftCorner(n_available_obs, n_available_obs);
    E = e.head(n_available_obs);
    // Compute innovation and remove rows/cols from dHdX and V
    for(int j=0; j < n_available_obs; j++){
      // Grab indices where obsVec has actual (non-NA) entries
      idx = obs_ids[j];
      // Innovations
      E(j) = obsVec(idx) - H(idx);
      // Obs Jacobian
      C.row(j) = dHdX.row(idx);
      // Variance
      V(j,j) = Hvar(idx, idx);
    }
    // Kalman Gain
    R = C * covMat * C.transpose() + V;
    Eigen::LLT<Eigen::MatrixXd> llt(R);
    K.transpose() = llt.solve(C * covMat);
    /*Eigen::LDLT<Eigen::MatrixXd> ldlt(R);
    K.transpose() = ldlt.solve(C * covMat);*/
    // State Update
    stateVec += K*E;
    // Covariance Update - Joseph Form
    IKC = I - K * C;
    covMat = IKC * covMat * IKC.transpose() + K * V * K.transpose();
    // Store innovations
    Innovation(0) = E;
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
    if(any_available_obs(i+1)){
      inputVec = inputMat.row(i+1);
      obsVec = obsMat.row(i+1);
      obs_ids = Rcpp::as<Eigen::VectorXi>(non_na_ids(i+1));
      n_available_obs = obs_ids.size();
      // Calculcate H
      H = h__(stateVec, parVec, inputVec);
      dHdX = dhdx__(stateVec, parVec, inputVec);
      Hvar = hvar__(stateVec, parVec, inputVec);
      C = dHdX.topRows(n_available_obs);
      V = Hvar.topLeftCorner(n_available_obs, n_available_obs);
      E = e.head(n_available_obs);
      // Compute innovation and remove rows/cols from dHdX and V
      for(int j=0; j < n_available_obs; j++){
        idx = obs_ids[j];
        // Innovations
        E(j) = obsVec(idx) - H(idx);
        // Obs Jacobian
        C.row(j) = dHdX.row(idx);
        // Variance
        V(j,j) = Hvar(idx,idx);
      }
      // Kalman Gain
      R = C * covMat * C.transpose() + V;
      Eigen::LLT<Eigen::MatrixXd> llt(R);
      K.transpose() = llt.solve(C * covMat);
      /*Eigen::LDLT<Eigen::MatrixXd> ldlt(R);
      K.transpose() = ldlt.solve(C * covMat);*/
      // State Update
      stateVec += K*E;
      // Covariance Update - Joseph Form
      IKC = I - K * C;
      covMat = IKC * covMat * IKC.transpose() + K * V * K.transpose();
      // Store innovations
      Innovation(i+1) = E;
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

// [[Rcpp::export]]
List lkf_predict_rcpp(
  List funPtrs,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  LogicalVector any_available_obs,
  List non_na_ids,
  const int last_pred_id,
  const int k_step_ahead)
{
  List filt = lkf_filter_rcpp(
    funPtrs,
    obsMat,
    inputMat,
    parVec,
    covMat, 
    stateVec,
    ode_timestep_size,
    any_available_obs,
    non_na_ids
    );

  List xPost = filt["xPost"];
  List pPost = filt["pPost"];
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
  List xk(last_pred_id);
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

//  This is the predict kalman filter function
// [[Rcpp::export]]
List lkf_simulate_rcpp(
  List funPtrs,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd simulation_timestep_size,
  Eigen::VectorXd simulation_timesteps,
  LogicalVector any_available_obs,
  List non_na_ids,
  int ng,
  int last_pred_id,
  int k_step_ahead,
  int nsims,
  Nullable<int> seed)
{

  // Set simulating seed if seed is not NULL
  set_simulation_seed(seed, ziggurat_states);

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");

  List filt = lkf_filter_rcpp(
    funPtrs,
    obsMat,
    inputMat,
    parVec,
    covMat, 
    stateVec,
    ode_timestep_size,
    any_available_obs,
    non_na_ids
    );
  List xPost = filt["xPost"];
  List pPost = filt["pPost"];

  // misc  
  const int n = stateVec.size();
  const int ni = inputMat.row(0).size();
  VectorXd inputVec(ni), dinputVec(ni);
  MatrixXd stateMat(nsims, n), randN(n, nsims);
  VectorXd simulation_timesteps_inv = simulation_timesteps.cwiseInverse();

  // storage for predictions
  List outer_simulate_list(last_pred_id);
  
  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){

    List inner_simulate_list(k_step_ahead + 1);

    // Extract posteriors as starting points
    stateVec = xPost(i);
    covMat = pPost(i);

    /* 
    1. We draw from a multivariate normal by z = u + A * dB where A = chol(covMat), dB is i.d.d normal vector
    2. We do simultaneously for all simulations i.e. dB is matrix of #nsims i.d.d vectors and u is repeated means
    */
    for(int j=0; j < n; j++){
      for(int k=0; k < nsims; k++){
        randN(j,k) = ziggurat_states.rnorm();
      }
    }
    stateMat = (stateVec.replicate(1, nsims) + covMat.llt().matrixL() * randN).transpose();
    inner_simulate_list(0) = stateMat;

    /* For each prediction horizon k: */
    for(int k=0 ; k < k_step_ahead ; k++){
      inputVec = inputMat.row(i+k);
      dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k)) * simulation_timesteps_inv(i+k);

      for(int j=0 ; j < simulation_timesteps(i+k) ; j++){
        euler_maruyama_simulation_inplace(
          f__, g__, 
          stateMat,
          parVec, inputVec, 
          simulation_timestep_size(i+k), 
          nsims, n, ng
          );
        inputVec += dinputVec;
      }

      inner_simulate_list(k+1) = stateMat;
    }

    outer_simulate_list(i) = inner_simulate_list;
  }

  // Return
  return outer_simulate_list;
}

// [[Rcpp::export]]
Rcpp::List test_return_list(Eigen::MatrixXd a){
  Rcpp::List A(3);
  A(0) = a;
  A(1) = a;
  return A;
}