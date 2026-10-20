#include <Rcpp.h>
#include <RcppEigen.h>

#include "function_typedefs.h"
#include "misc_helpers.h"
#include "extra_utils.h"
#include "sde_solvers.h"
#include "helpers_ukf.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// This structure holds the UKF weights
struct UKF_weights {
  Eigen::VectorXd wm;
  Eigen::MatrixXd W;
  double sqrt_c;
};

UKF_weights ukf_create_weights(const Eigen::VectorXd& ukf_pars, const int n)
{
  const int nn = 2*n + 1;
  const double alpha = ukf_pars(0);
  const double beta  = ukf_pars(1);
  const double kappa = ukf_pars(2);

  const double sqrt_c = std::sqrt(alpha*alpha*(n + kappa));
  const double lambda = sqrt_c*sqrt_c - n;

  VectorXd wm = VectorXd::Constant(nn, 1.0/(2.0*(n + lambda)));
  VectorXd wc = wm;
  wm(0) = lambda/(lambda + n);
  wc(0) = lambda/(n + lambda) + (1 - alpha*alpha + beta);

  MatrixXd IWm = MatrixXd::Identity(nn, nn) - wm.replicate(1, nn);
  MatrixXd W = IWm * wc.asDiagonal() * IWm.transpose();

  return {wm, W, sqrt_c};
}

// Lower Cholesky factor, stopping with an informative error if it fails
MatrixXd chol_or_stop(const MatrixXd& P, const char* what, const int k)
{
  LLT<MatrixXd> llt(P);
  if (llt.info() != Eigen::Success) {
    Rcpp::stop("%s covariance not positive definite at time index %d", what, k + 1);
  }
  return llt.matrixL();
}

// This structure holds the results of the unscented kalman filter
struct UKF_structure {
  Rcpp::List xPrior, xPost, pPrior, pPost, Innovation, InnovationCovariance; // remove these eventually
  Eigen::MatrixXd priorMat, postMat, innoMat;
};

UKF_structure ukf_forward_pass(
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
    Eigen::VectorXd ukf_pars,
    Rcpp::CharacterVector ode_solver,
    bool first_order_input_hold)
{

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto h__ = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");
  auto dhdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar__ = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");

  // constants
  const int tsize = inputMat.rows();
  const int ni = inputMat.cols();
  const int n = stateVec.size();
  const int nn = 2*n + 1;
  const int m = obsMat.cols();
  const int n_squared = n*n;
  const int m_squared = m*m;

  // ODE solver
  const int ODE_solver = choose_solver(ode_solver);

  // UKF weights
  const UKF_weights uw = ukf_create_weights(ukf_pars, n);
  const VectorXd& wm = uw.wm;
  const MatrixXd& W = uw.W;
  const double sqrt_c = uw.sqrt_c;

  // pre-allocate and define
  VectorXd inputVec, obsVec, E;
  VectorXi obs_ids;
  MatrixXd Xsigmapoints, sqrt_covMat, H, Hvar, V, R, Sxy, K, C, dHdX, IKC;
  MatrixXd I = MatrixXd::Identity(n, n);
  VectorXd inv_ode_timesteps = ode_timesteps.cast<double>().cwiseInverse();
  VectorXd dinputVec = VectorXd::Zero(ni);

  // Pre-allocate storage for output
  List Innovation(tsize), InnovationCovariance(tsize), xPrior(tsize), xPost(tsize), pPrior(tsize), pPost(tsize);
  MatrixXd priorMat(tsize, n + n_squared), postMat(tsize, n + n_squared), innoMat(tsize, m + m_squared);
  innoMat.setConstant(NA_REAL);
  MatrixXd Rtemp = MatrixXd::Constant(m, m, NA_REAL);

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

    // Observation function of sigma points, reduced to available observations
    H = ukf_h(h__, Xsigmapoints, parVec, inputVec, m, nn)(obs_ids, Eigen::all);
    Hvar = hvar__(stateVec, parVec, inputVec);
    V = Hvar(obs_ids, obs_ids);

    // Innovation, innovation covariance and state-observation cross covariance
    E = obsVec(obs_ids) - H * wm;
    R = H * W * H.transpose() + V;
    Sxy = Xsigmapoints * W * H.transpose();

    // Kalman gain
    LLT<MatrixXd> llt(R);
    if (llt.info() != Eigen::Success)
      Rcpp::stop("Innovation covariance not positive definite at time index %d", k + 1);
    K = llt.solve(Sxy.transpose()).transpose();

    // State and covariance update
    stateVec += K * E;
    covMat -= K * R * K.transpose();
    // IKC = I - K * C;
    // covMat = IKC * covMat * IKC.transpose() + K * V * K.transpose();

    // Store innovations
    store_innovations(k);
  };


  // Main Script Start
  // ----------------------------------------

  store_prior(0);

  // Initial data update
  if (any_available_obs(0)) {
    sqrt_covMat = chol_or_stop(covMat, "Prior", 0);
    Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, sqrt_c, n, nn);
    data_update(0);
  }
  store_posterior(0);

  // Main loop
  for (int i = 0; i < tsize - 1; i++) {
    // Time update
    inputVec = inputMat.row(i);
    if (first_order_input_hold)
      dinputVec = (inputMat.row(i+1) - inputMat.row(i)) * inv_ode_timesteps(i);

    sqrt_covMat = chol_or_stop(covMat, "Posterior", i);
    Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, sqrt_c, n, nn);
    for (int j = 0; j < ode_timesteps(i); j++) {
      Xsigmapoints = ukf_ode_integration(
        f__, g__, Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, dinputVec,
        ode_timestep_size(i), sqrt_c, ODE_solver, n, nn);
      sqrt_covMat = sigma2chol(Xsigmapoints, sqrt_c, nn, n);
      inputVec += dinputVec;
    }
    stateVec = Xsigmapoints.col(0);
    covMat = sqrt_covMat * sqrt_covMat.transpose();
    store_prior(i+1);

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
    xPrior, xPost, pPrior, pPost, Innovation, InnovationCovariance,
    priorMat, postMat, innoMat
  };

}

// [[Rcpp::export]]
List ukf_filter_rcpp(
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
    Eigen::VectorXd ukf_pars,
    Rcpp::CharacterVector ode_solver,
    bool first_order_input_hold)
{
  UKF_structure fw = ukf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                      ode_timestep_size, ode_timesteps, any_available_obs, non_na_ids,
                                      ukf_pars, ode_solver, first_order_input_hold);
  return List::create(
    Named("priorMat") = fw.priorMat,
    Named("postMat") = fw.postMat,
    Named("innoMat") = fw.innoMat
  );
}

// [[Rcpp::export]]
List ukf_predict_rcpp(
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
    Eigen::VectorXd ukf_pars,
    int last_pred_id,
    int k_step_ahead,
    Rcpp::CharacterVector ode_solver,
    bool first_order_input_hold)
{

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");

  UKF_structure fw = ukf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                      ode_timestep_size, ode_timesteps, any_available_obs, non_na_ids,
                                      ukf_pars, ode_solver, first_order_input_hold);
  const MatrixXd& postMat = fw.postMat;

  // constants
  const int n = stateVec.size();
  const int nn = 2*n + 1;
  const int n_squared = n*n;
  const int ni = inputMat.cols();

  // ODE solver
  const int ODE_solver = choose_solver(ode_solver);

  // UKF weights
  const UKF_weights uw = ukf_create_weights(ukf_pars, n);
  const VectorXd& wm = uw.wm;
  const MatrixXd& W = uw.W;
  const double sqrt_c = uw.sqrt_c;

  // pre-allocate and define
  VectorXd inputVec, dinputVec = VectorXd::Zero(ni);
  VectorXd inv_ode_timesteps = ode_timesteps.cast<double>().cwiseInverse();
  MatrixXd predMat = MatrixXd::Zero(k_step_ahead+1, n + n_squared);
  MatrixXd Xsigmapoints, sqrt_covMat;
  List xk(last_pred_id);

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){
    stateVec = postMat.row(i).head(n);
    covMat = postMat.row(i).tail(n_squared).reshaped(n,n);

    predMat.row(0).head(n) = stateVec;
    predMat.row(0).tail(n_squared) = covMat.reshaped();

    sqrt_covMat = chol_or_stop(covMat, "Posterior", i);
    Xsigmapoints = create_sigmapoints_from_stateVec(stateVec, sqrt_covMat, sqrt_c, n, nn);

    //////////// K-STEP-AHEAD LOOP ///////////
    for(int k=0 ; k < k_step_ahead ; k++){

      inputVec = inputMat.row(i+k);
      if(first_order_input_hold){
        dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k)) * inv_ode_timesteps(i+k);
      }

      //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
      for(int j=0 ; j < ode_timesteps(i+k) ; j++){
        Xsigmapoints = ukf_ode_integration(
          f__, g__, Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, dinputVec,
          ode_timestep_size(i+k), sqrt_c, ODE_solver, n, nn);
        sqrt_covMat = sigma2chol(Xsigmapoints, sqrt_c, nn, n);
        inputVec += dinputVec;
      }

      // Save the prediction
      stateVec = Xsigmapoints.col(0);
      covMat = sqrt_covMat * sqrt_covMat.transpose();
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
List ukf_simulate_rcpp(
    Rcpp::List funPtrs,
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
    Rcpp::List non_na_ids,
    Eigen::VectorXd ukf_pars,
    int ng,
    int last_pred_id,
    int k_step_ahead,
    Rcpp::CharacterVector ode_solver,
    int nsims,
    Nullable<int> seed,
    bool first_order_input_hold)
{

  // Set simulating seed if seed is not NULL
  set_simulation_seed(seed, ziggurat_states);

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");

  UKF_structure fw = ukf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                      ode_timestep_size, ode_timesteps, any_available_obs, non_na_ids,
                                      ukf_pars, ode_solver, first_order_input_hold);
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

    stateMat = (stateVec.replicate(1, nsims) + chol_or_stop(covMat, "Posterior", i) * randN).transpose();
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
