#include <Rcpp.h>
#include <RcppEigen.h>
#include <zigg/header>

#include "function_typedefs.h"
#include "misc_helpers.h"
#include "ode_solvers.h"
#include "sde_solvers.h"

using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

zigg::Ziggurat ziggurat_states;

// [[Rcpp::export]]
List ekf_filter_rcpp(
  List funPtrs, 
  const Eigen::MatrixXd & obsMat,
  const Eigen::MatrixXd & inputMat,
  const Eigen::VectorXd & parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  const Eigen::VectorXd & ode_timestep_size,
  const Eigen::VectorXd & ode_timesteps,
  LogicalVector any_available_obs,
  List non_na_ids,
  const int ode_solver)
{

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto h__ = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");
  auto dfdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dfdx_const");
  auto dhdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar__ = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");

  // constants  
  const int tsize = inputMat.col(0).size();
  const int ni = inputMat.row(0).size();
  const int n = stateVec.size();
  const int m = obsMat.row(0).size();
  int idx, n_available_obs;

  // pre-allocate and define
  Eigen::VectorXd H, inputVec(ni), dinputVec(ni), obsVec(m), e(m);
  Eigen::VectorXd E;
  Eigen::MatrixXd C, R, K, V, KC, IKC, Hvar, dHdX;
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
  Eigen::VectorXd inv_ode_timesteps = ode_timesteps.cwiseInverse();
  Eigen::VectorXi obs_ids;
  // Pre-allocate for ODE solver
  Eigen::VectorXd k1(n), k2(n), k3(n), k4(n);
  Eigen::MatrixXd c1(n,n), c2(n,n), c3(n,n), c4(n,n);
  // Pre-allocate storage for output
  Rcpp::List Innovation(tsize), InnovationCovariance(tsize), xPrior(tsize), xPost(tsize), pPrior(tsize), pPost(tsize);

  // store prior
  xPrior(0) = stateVec;
  pPrior(0) = covMat;

   //////////// INITIAL DATA-UPDATE ///////////
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

    inputVec = inputMat.row(i);
    dinputVec = (inputMat.row(i+1) - inputMat.row(i)) * inv_ode_timesteps(i);

    //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
    for(int j=0 ; j < ode_timesteps(i) ; j++){
      ode_integrator_inplace(
        f__, g__, dfdx__, 
        covMat, stateVec, 
        parVec, 
        inputVec, dinputVec, 
        ode_timestep_size(i), ode_solver, 
        k1, k2, k3, k4, c1, c2, c3, c4
        );
      inputVec += dinputVec;
    }
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
    Named("InnovationCovariance") = InnovationCovariance
    );
}

//  This is the predict kalman filter function
// [[Rcpp::export]]
List ekf_predict_rcpp(
  List funPtrs,
  const Eigen::MatrixXd & obsMat,
  const Eigen::MatrixXd & inputMat,
  const Eigen::VectorXd & parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  const Eigen::VectorXd & ode_timestep_size,
  const Eigen::VectorXd & ode_timesteps,
  Rcpp::LogicalVector any_available_obs,
  List non_na_ids,
  const int ode_solver,
  const int last_pred_id,
  const int k_step_ahead)
{

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");
  auto dfdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dfdx_const");

  List filt = ekf_filter_rcpp(
    funPtrs,
    obsMat,
    inputMat,
    parVec,
    covMat,
    stateVec,
    ode_timestep_size,
    ode_timesteps,
    any_available_obs,
    non_na_ids,
    ode_solver
    );
  List xPost = filt["xPost"];
  List pPost = filt["pPost"];

  // pre-allocate and misc  
  const int n = stateVec.size();
  const int ni = inputMat.row(0).size();
  const int n_squared = n*n;
  Eigen::VectorXd inputVec(ni), dinputVec(ni);
  Eigen::MatrixXd predMat(k_step_ahead+1, n + n_squared);
  predMat.setZero();
  Rcpp::List xk(last_pred_id), ode_1step_integration(2);
  Eigen::VectorXd inv_ode_timesteps = ode_timesteps.cwiseInverse();

  // Pre-allocate for ODE solver
  Eigen::VectorXd k1(n), k2(n), k3(n), k4(n);
  Eigen::MatrixXd c1(n,n), c2(n,n), c3(n,n), c4(n,n);

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){
    stateVec = xPost(i);
    covMat = pPost(i);

    predMat.row(0).head(n) = stateVec;
    predMat.row(0).tail(n_squared) = covMat.reshaped();

    //////////// K-STEP-AHEAD LOOP ///////////
    for(int k=0 ; k < k_step_ahead ; k++){
      inputVec = inputMat.row(i+k);
      dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k)) * inv_ode_timesteps(i+k);

      //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
      for(int j=0 ; j < ode_timesteps(i+k) ; j++){
      ode_integrator_inplace(
        f__, g__, dfdx__, 
        covMat, stateVec, 
        parVec, 
        inputVec, dinputVec, 
        ode_timestep_size(i+k), ode_solver, 
        k1, k2, k3, k4, c1, c2, c3, c4
        );
      inputVec += dinputVec;
    }
      // Save the prediction
      predMat.row(k+1).head(n) = stateVec;
      predMat.row(k+1).tail(n_squared) = covMat.reshaped();
    }

    xk(i) = predMat;
    }

  return List::create(
    Named("predMats") = xk
    );
}

//  This is the predict kalman filter function
// [[Rcpp::export]]
List ekf_simulate_rcpp(
  List funPtrs,
  Eigen::MatrixXd obsMat,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd ode_timestep_size,
  Eigen::VectorXd ode_timesteps, 
  Eigen::VectorXd simulation_timestep_size,
  Eigen::VectorXd simulation_timesteps,
  Rcpp::LogicalVector any_available_obs,
  List non_na_ids,
  const int ode_solver,
  const int last_pred_id,
  const int k_step_ahead,
  const int ng,
  const int nsims,
  Nullable<int> seed)
{

  // Set simulating seed if seed is not NULL
  set_simulation_seed(seed, ziggurat_states);

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");

  List filt = ekf_filter_rcpp(
    funPtrs,
    obsMat,
    inputMat,
    parVec,
    covMat,
    stateVec,
    ode_timestep_size,
    ode_timesteps,
    any_available_obs,
    non_na_ids,
    ode_solver
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
