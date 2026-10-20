#include <Rcpp.h>
#include <RcppEigen.h>
#include <zigg/header>

#include "function_typedefs.h"
#include "misc_helpers.h"
#include "ode_solvers.h"
#include "sde_solvers.h"
#include "extra_utils.h"

using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

zigg::Ziggurat ziggurat_states;

// This structure holds the results of the kalman filter
struct EKF_structure {
  Rcpp::List xPrior, xPost, pPrior, pPost, cPrior, Innovation, InnovationCovariance; // remove these eventually
  Eigen::MatrixXd priorMat, postMat, innoMat;
};

EKF_structure ekf_forward_pass(
    Rcpp::List funPtrs,
    Eigen::MatrixXd obsMat,
    Eigen::MatrixXd inputMat,
    Eigen::VectorXd parVec,
    Eigen::MatrixXd covMat,
    Eigen::VectorXd stateVec,
    Eigen::VectorXd ode_timestep_size,
    Eigen::VectorXi ode_timesteps,
    Rcpp::LogicalVector any_available_obs,
    Rcpp::List non_na_ids,
    Rcpp::CharacterVector ode_solver,
    bool first_order_input_hold,
    bool store_cross_covariance)
{

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto h__ = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");
  auto dfdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dfdx_const");
  auto dhdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar__ = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");

  // constants
  const int tsize = inputMat.rows();
  const int ni = inputMat.cols();
  const int n = stateVec.size();
  const int m = obsMat.cols();
  const int n_squared = n*n;
  const int m_squared = m*m;

  // ODE solver
  const int ODE_solver = choose_solver(ode_solver);

  // pre-allocate and define
  VectorXd inputVec, obsVec;
  VectorXd H, E;
  MatrixXd C, R, K, V, IKC, Hvar, dHdX;
  MatrixXd I = MatrixXd::Identity(n, n);

  VectorXd inv_ode_timesteps = ode_timesteps.cast<double>().cwiseInverse();
  VectorXi obs_ids;
  VectorXd dinputVec = VectorXd::Zero(ni);

  // Pre-allocate for ODE solver
  VectorXd k1,k2,k3,k4;
  MatrixXd c1,c2,c3,c4;

  // Pre-allocate storage for output
  List Innovation(tsize), InnovationCovariance(tsize), xPrior(tsize), xPost(tsize), pPrior(tsize), pPost(tsize);
  MatrixXd priorMat(tsize, n + n_squared), postMat(tsize, n + n_squared), innoMat(tsize, m + m_squared);
  innoMat.setConstant(NA_REAL);
  MatrixXd Rtemp = MatrixXd::Constant(m, m, NA_REAL);

  // Smoothing allocations
  Eigen::MatrixXd cross_covMat, d1, d2, d3, d4;
  Rcpp::List cPrior(store_cross_covariance ? tsize : 0);

  // storage functions for posterior and priors
  auto store_prior = [&](int k) {
    xPrior(k) = stateVec;
    pPrior(k) = covMat;
    priorMat.row(k).head(n) = stateVec;
    priorMat.row(k).tail(n_squared) = covMat.reshaped();
  };

  auto store_posterior = [&](int k) {
    xPost(k) = stateVec;
    pPost(k) = covMat;
    postMat.row(k).head(n) = stateVec;
    postMat.row(k).tail(n_squared) = covMat.reshaped();
  };

  auto store_innovations = [&](int k) {
    Innovation(k) = E;
    InnovationCovariance(k) = R;
    Rtemp.setConstant(NA_REAL);
    Rtemp(obs_ids, obs_ids) = R;
    innoMat(k, obs_ids) = E;
    innoMat.row(k).tail(m_squared) = Rtemp.reshaped();
  };

  auto data_update = [&](int k) {
    inputVec = inputMat.row(k);
    obsVec   = obsMat.row(k);
    obs_ids  = Rcpp::as<Eigen::VectorXi>(non_na_ids(k));

    H    = h__(stateVec, parVec, inputVec);
    dHdX = dhdx__(stateVec, parVec, inputVec);
    Hvar = hvar__(stateVec, parVec, inputVec);

    // Reduce to available observations
    E = obsVec(obs_ids) - H(obs_ids);
    C = dHdX(obs_ids, Eigen::all);
    V = Hvar(obs_ids, obs_ids);

    // Kalman gain
    R = C * covMat * C.transpose() + V;
    LLT<MatrixXd> llt(R);
    if (llt.info() != Eigen::Success)
      Rcpp::stop("Innovation covariance not positive definite at time index %d", k + 1);
    K = llt.solve(C * covMat).transpose();

    // State and covariance update (Joseph form)
    stateVec += K * E;
    IKC = I - K * C;
    covMat = IKC * covMat * IKC.transpose() + K * V * K.transpose();
    covMat = 0.5 * (covMat + covMat.transpose());

    // Store innovations
    store_innovations(k);
  };


  // Main Script Start
  // ----------------------------------------

  store_prior(0);

  // Initial data update
  if (any_available_obs(0)) {
    data_update(0);
  }
  store_posterior(0);

  // Main loop
  for (int i = 0; i < tsize - 1; i++) {
    if (store_cross_covariance) cross_covMat = covMat; // Integrating starts from covMat
    // Time update
    inputVec = inputMat.row(i);
    if (first_order_input_hold)
      dinputVec = (inputMat.row(i+1) - inputMat.row(i)) * inv_ode_timesteps(i);

    for (int j = 0; j < ode_timesteps(i); j++) {
      if (store_cross_covariance)
        ode_integrator_with_cross_cov_inplace(
          f__, g__, dfdx__, covMat, cross_covMat, stateVec,
          parVec, inputVec, dinputVec, ode_timestep_size(i), ODE_solver,
          k1, k2, k3, k4, c1, c2, c3, c4, d1, d2, d3, d4);
      else
        ode_integrator_inplace(
          f__, g__, dfdx__, covMat, stateVec, parVec,
          inputVec, dinputVec, ode_timestep_size(i), ODE_solver,
          k1, k2, k3, k4, c1, c2, c3, c4);
      inputVec += dinputVec;
    }
    store_prior(i+1);
    if (store_cross_covariance) cPrior(i+1) = cross_covMat;

    // Data update
    if (any_available_obs(i+1)) {
      data_update(i+1);
    }
    store_posterior(i+1);
  }
  // ----------------------------------------
  // Main Script End

  // Return
  return {
    xPrior, xPost, pPrior, pPost, cPrior, Innovation, InnovationCovariance,
    priorMat, postMat, innoMat
  };

}

// [[Rcpp::export]]
List ekf_filter_rcpp(
    Rcpp::List funPtrs,
    Eigen::MatrixXd obsMat,
    Eigen::MatrixXd inputMat,
    Eigen::VectorXd parVec,
    Eigen::MatrixXd covMat,
    Eigen::VectorXd stateVec,
    Eigen::VectorXd ode_timestep_size,
    Eigen::VectorXi ode_timesteps,
    Rcpp::LogicalVector any_available_obs,
    Rcpp::List non_na_ids,
    Rcpp::CharacterVector ode_solver,
    bool first_order_input_hold)
{
  EKF_structure fw = ekf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                   ode_timestep_size, ode_timesteps, any_available_obs, non_na_ids,
                                   ode_solver, first_order_input_hold, false);
  return List::create(
    Named("priorMat") = fw.priorMat,
    Named("postMat") = fw.postMat,
    Named("innoMat") = fw.innoMat
  );
}

// [[Rcpp::export]]
List ekf_smooth_rcpp(
    Rcpp::List funPtrs,
    Eigen::MatrixXd obsMat,
    Eigen::MatrixXd inputMat,
    Eigen::VectorXd parVec,
    Eigen::MatrixXd covMat,
    Eigen::VectorXd stateVec,
    Eigen::VectorXd ode_timestep_size,
    Eigen::VectorXi ode_timesteps,
    Rcpp::LogicalVector any_available_obs,
    Rcpp::List non_na_ids,
    Rcpp::CharacterVector ode_solver,
    bool first_order_input_hold)
{

  EKF_structure fw = ekf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                   ode_timestep_size, ode_timesteps, any_available_obs, non_na_ids,
                                   ode_solver, first_order_input_hold, true);

  // constants
  const int tsize = inputMat.rows();
  const int n = stateVec.size();
  const int n_squared = n*n;

  // pre-allocations
  MatrixXd smoothMat(tsize, n + n_squared);
  VectorXd xPrior, xPost, xSmooth;
  MatrixXd G, pPrior, pPost, pSmooth, cross_covMat;

  // The posterior and smooth coincide at the last time point
  smoothMat.row(tsize-1) = fw.postMat.row(tsize-1);

  xSmooth = smoothMat.row(tsize-1).head(n);
  pSmooth = smoothMat.row(tsize-1).tail(n_squared).reshaped(n,n);
  // Main loop back in time
  for(int i = tsize - 1; i > 0; i--) {
    xPrior = fw.priorMat.row(i).head(n);
    pPrior = fw.priorMat.row(i).tail(n_squared).reshaped(n,n);
    xPost = fw.postMat.row(i-1).head(n);
    pPost = fw.postMat.row(i-1).tail(n_squared).reshaped(n,n);
    // Solve for G:
    // G = C P^{-1}
    // G^{T} = P^{-1} C^{T}
    // P G^{T} = C^{T}
    // G = solve(P, C^{T}).transpose()
    Eigen::LLT<Eigen::MatrixXd> llt(pPrior);
    if (llt.info() != Eigen::Success) {
      Rcpp::print(wrap(pPrior));
      Rcpp::stop("Prior-covariance not positive definite at time index %d", i+1);
    }
    cross_covMat = fw.cPrior(i);
    G = llt.solve(cross_covMat.transpose()).transpose();

    // Compute "next" smoothed state
    xSmooth = xPost + G * (xSmooth - xPrior);
    pSmooth = pPost + G * (pSmooth - pPrior) * G.transpose();

    // store
    smoothMat.row(i-1).head(n) = xSmooth;
    smoothMat.row(i-1).tail(n_squared) = pSmooth.reshaped();
  }

  return List::create(
    Named("priorMat") = fw.priorMat,
    Named("postMat") = fw.postMat,
    Named("smoothMat") = smoothMat
  );

}

//  This is the predict kalman filter function
// [[Rcpp::export]]
List ekf_predict_rcpp(
    Rcpp::List funPtrs,
    Eigen::MatrixXd obsMat,
    Eigen::MatrixXd inputMat,
    Eigen::VectorXd parVec,
    Eigen::MatrixXd covMat,
    Eigen::VectorXd stateVec,
    Eigen::VectorXd ode_timestep_size,
    Eigen::VectorXi ode_timesteps,
    Rcpp::LogicalVector any_available_obs,
    Rcpp::List non_na_ids,
    Rcpp::CharacterVector ode_solver,
    int last_pred_id,
    int k_step_ahead,
    bool first_order_input_hold)
{

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");
  auto dfdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dfdx_const");

  EKF_structure fw = ekf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                   ode_timestep_size, ode_timesteps, any_available_obs, non_na_ids,
                                   ode_solver, first_order_input_hold, false);
  const MatrixXd& postMat = fw.postMat;

  // pre-allocate and misc
  const int n = stateVec.size();
  const int n_squared = n*n;
  const int ni = inputMat.cols();

  VectorXd inputVec, dinputVec = VectorXd::Zero(ni);

  List xk(last_pred_id);
  MatrixXd predMat = MatrixXd::Zero(k_step_ahead+1, n + n_squared);
  VectorXd inv_ode_timesteps = ode_timesteps.cast<double>().cwiseInverse();

  // ODE solver
  const int ODE_solver = choose_solver(ode_solver);

  // Pre-allocate for ODE solver
  VectorXd k1(n), k2(n), k3(n), k4(n);
  MatrixXd c1(n,n), c2(n,n), c3(n,n), c4(n,n);

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){
    stateVec = postMat.row(i).head(n);
    covMat = postMat.row(i).tail(n_squared).reshaped(n,n);

    predMat.row(0).head(n) = stateVec;
    predMat.row(0).tail(n_squared) = covMat.reshaped();

    //////////// K-STEP-AHEAD LOOP ///////////
    for(int k=0 ; k < k_step_ahead ; k++){

      inputVec = inputMat.row(i+k);
      if(first_order_input_hold){
        dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k)) * inv_ode_timesteps(i+k);
      }

      //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
      for(int j=0 ; j < ode_timesteps(i+k) ; j++){
        ode_integrator_inplace(
          f__, g__, dfdx__,
          covMat, stateVec,
          parVec,
          inputVec, dinputVec,
          ode_timestep_size(i+k), ODE_solver,
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

// [[Rcpp::export]]
List ekf_simulate_rcpp(
    List funPtrs,
    Eigen::MatrixXd obsMat,
    Eigen::MatrixXd inputMat,
    Eigen::VectorXd parVec,
    Eigen::MatrixXd covMat,
    Eigen::VectorXd stateVec,
    Eigen::VectorXd ode_timestep_size,
    Eigen::VectorXi ode_timesteps,
    Eigen::VectorXd simulation_timestep_size,
    Eigen::VectorXi simulation_timesteps,
    Rcpp::LogicalVector any_available_obs,
    List non_na_ids,
    CharacterVector ode_solver,
    int last_pred_id,
    int k_step_ahead,
    int ng,
    int nsims,
    Nullable<int> seed,
    bool first_order_input_hold)
{

  // Set simulating seed if seed is not NULL
  set_simulation_seed(seed, ziggurat_states);

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");

  EKF_structure fw = ekf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                   ode_timestep_size, ode_timesteps, any_available_obs, non_na_ids,
                                   ode_solver, first_order_input_hold, false);
  const MatrixXd& postMat = fw.postMat;

  // misc
  const int n = stateVec.size();
  const int n_squared = n*n;
  const int ni = inputMat.cols();
  VectorXd inputVec(ni);
  VectorXd dinputVec = VectorXd::Zero(ni);
  MatrixXd stateMat(nsims, n), randN(n, nsims);
  VectorXd simulation_timesteps_inv = simulation_timesteps.cast<double>().cwiseInverse();

  // storage for predictions
  List outer_simulate_list(last_pred_id);

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){

    List inner_simulate_list(k_step_ahead + 1);

    // Extract posteriors as starting points
    stateVec = postMat.row(i).head(n);
    covMat = postMat.row(i).tail(n_squared).reshaped(n,n);

    /*
     1. We draw from a multivariate normal by z = u + A * dB where A = chol(covMat), dB is i.d.d normal vector
     2. We do simultaneously for all simulations i.e. dB is matrix of #nsims i.d.d vectors and u is repeated means
     */
    for(int j=0; j < n; j++){
      for(int k=0; k < nsims; k++){
        randN(j,k) = ziggurat_states.rnorm();
      }
    }

    LLT<MatrixXd> llt(covMat);
    if (llt.info() != Eigen::Success) {
      Rcpp::stop("Posterior covariance not positive definite at time index %d", i+1);
    }
    stateMat = (stateVec.replicate(1, nsims) + llt.matrixL() * randN).transpose();
    inner_simulate_list(0) = stateMat;

    /* For each prediction horizon k: */
    for(int k=0 ; k < k_step_ahead ; k++){

      inputVec = inputMat.row(i+k);
      if(first_order_input_hold){
        dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k)) * simulation_timesteps_inv(i+k);
      }

      for(int j=0; j < simulation_timesteps(i+k); j++){
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
