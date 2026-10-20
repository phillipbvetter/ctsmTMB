#include <Rcpp.h>
#include <RcppEigen.h>

#include "function_typedefs.h"
#include "misc_helpers.h"
#include "sde_solvers.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// This structure holds the results of the linear kalman filter
struct LKF_structure {
  Rcpp::List xPrior, xPost, pPrior, pPost, Innovation, InnovationCovariance; // remove these eventually
  Eigen::MatrixXd priorMat, postMat, innoMat;
  Eigen::MatrixXd Ahat, Bhat, BhatFOH, Vhat;
  double dt_timestep;
};

LKF_structure lkf_forward_pass(
    Rcpp::List funPtrs,
    Eigen::MatrixXd obsMat,
    Eigen::MatrixXd inputMat,
    Eigen::VectorXd parVec,
    Eigen::MatrixXd covMat,
    Eigen::VectorXd stateVec,
    Eigen::VectorXd ode_timestep_size,
    Rcpp::LogicalVector any_available_obs,
    Rcpp::List non_na_ids,
    bool first_order_input_hold)
{

  auto h__ = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");
  auto dfdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dfdx_const");
  auto dhdx__ = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar__ = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");
  auto dfdu__ = get_funptr<funPtr_mat_const>(funPtrs, "dfdu_const");

  // constants
  const int tsize = inputMat.rows();
  const int ni = inputMat.cols();
  const int ni1 = ni + 1;
  const int n = stateVec.size();
  const int m = obsMat.cols();
  const int n_squared = n*n;
  const int m_squared = m*m;

  // pre-allocate and define
  VectorXd inputVec, obsVec;
  VectorXd H, E;
  VectorXi obs_ids;
  MatrixXd I = MatrixXd::Identity(n, n);
  MatrixXd C, R, K, V, IKC, Hvar, dHdX;
  MatrixXd A, B, G, Phi1, ePhi1, Ahat, Ahat_T, Bhat, Phi2, ePhi2, Q12, Q22, Vhat;

  inputVec = inputMat.row(0);
  // We assume all time-steps identical
  const double dt_timestep = inputMat(1, 0) - inputMat(0, 0);
  Eigen::VectorXd data_dt = inputMat.col(0).tail(tsize - 1) - inputMat.col(0).head(tsize - 1);
  if ((data_dt.array() - dt_timestep).abs().maxCoeff() > 1e-10 * dt_timestep)
    Rcpp::stop("The LKF requires equidistant observation times.");

  // Construct one-step system matrices
  // ---------------------------------------------

  A = dfdx__(stateVec, parVec, inputVec);
  B = dfdu__(stateVec, parVec, inputVec);
  G = g__(stateVec, parVec, inputVec);

  // Phi1 - Deterministic Component
  // ZOH: [A B \\ 0 0]
  // FOH: [A B 0, 0 0 I, 0 0 0]
  const int size1 = first_order_input_hold ? n + 2*ni1 : n + ni1;
  Phi1 = MatrixXd::Zero(size1, size1);
  Phi1.block(0, 0, n, n) = A;
  Phi1.block(0, n, n, ni1) = B;
  if (first_order_input_hold) {
    Phi1.block(n, n + ni1, ni1, ni1) = MatrixXd::Identity(ni1, ni1);
  }
  // Compute matrix exponential
  ePhi1 = (Phi1 * dt_timestep).exp();
  // Extract
  Ahat = ePhi1.block(0, 0, n, n);
  Ahat_T = Ahat.transpose();
  Bhat = ePhi1.block(0, n, n, ni1);
  MatrixXd BhatFOH = MatrixXd::Zero(n, ni1);
  if (first_order_input_hold) {
    BhatFOH = ePhi1.block(0, n + ni1, n, ni1);
  }

  // Phi2 - Noise Component
  // [-A GG^T \\ 0 A^t]
  Phi2 = MatrixXd::Zero(2*n, 2*n);
  Phi2.block(0, 0, n, n) = -A;
  Phi2.block(0, n, n, n) = G*G.transpose();
  Phi2.block(n, n, n, n) = A.transpose();
  // Matrix exponential
  ePhi2 = (Phi2 * dt_timestep).exp();
  // Extract
  Vhat = ePhi2.block(n, n, n, n).transpose() * ePhi2.block(0, n, n, n);
  // ---------------------------------------------

  VectorXd constant_and_inputVec = VectorXd::Ones(ni1);   // First element must be 1
  VectorXd constant_and_dinputVec = VectorXd::Zero(ni1);  // First element must be 0

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
  for (int i = 0; i < tsize - 1; i++) {
    // Time update
    constant_and_inputVec.tail(ni) = inputMat.row(i);
    stateVec = Ahat * stateVec + Bhat * constant_and_inputVec;
    if (first_order_input_hold) {
      constant_and_dinputVec.tail(ni) = (inputMat.row(i+1) - inputMat.row(i)) / dt_timestep;
      stateVec += BhatFOH * constant_and_dinputVec;
    }
    covMat = Ahat * covMat * Ahat_T + Vhat;
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
    priorMat, postMat, innoMat,
    Ahat, Bhat, BhatFOH, Vhat,
    dt_timestep
  };

}

// [[Rcpp::export]]
List lkf_filter_rcpp(
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
    bool first_order_input_hold)
{
  LKF_structure fw = lkf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                      ode_timestep_size, any_available_obs, non_na_ids,
                                      first_order_input_hold);
  return List::create(
    Named("priorMat") = fw.priorMat,
    Named("postMat") = fw.postMat,
    Named("innoMat") = fw.innoMat,
    Named("xPrior") = fw.xPrior,
    Named("xPost") = fw.xPost,
    Named("pPrior") = fw.pPrior,
    Named("pPost") = fw.pPost,
    Named("Innovation") = fw.Innovation,
    Named("InnovationCovariance") = fw.InnovationCovariance,
    Named("Ahat") = fw.Ahat,
    Named("Bhat") = fw.Bhat,
    Named("BhatFOH") = fw.BhatFOH,
    Named("Vhat") = fw.Vhat
  );
}

// [[Rcpp::export]]
List lkf_smooth_rcpp(
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
    bool first_order_input_hold)
{

  LKF_structure fw = lkf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                      ode_timestep_size, any_available_obs, non_na_ids,
                                      first_order_input_hold);

  MatrixXd Ahat = fw.Ahat;

  // constants
  const int tsize = inputMat.rows();
  const int n = stateVec.size();
  const int n_squared = n*n;

  // pre-allocations
  MatrixXd smoothMat(tsize, n + n_squared);
  VectorXd xPrior, xPost, xSmooth;
  MatrixXd G, pPrior, pPost, pSmooth, cross_covMat_transpose;

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
    // G = P_post * Ahat^{T} * P_prior^{-1}
    // G P_prior = P_post * Ahat^{T}
    // P_prior * G^T = Ahat * P_post
    // G = solve(P_prior, Ahat * P_post).transpose()
    Eigen::LLT<Eigen::MatrixXd> llt(pPrior);
    if (llt.info() != Eigen::Success)
      Rcpp::stop("Prior-covariance not positive definite at time index %d", i+1);
    cross_covMat_transpose = Ahat * pPost;
    G = llt.solve(cross_covMat_transpose).transpose();

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

// [[Rcpp::export]]
List lkf_predict_rcpp(
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
    int last_pred_id,
    int k_step_ahead,
    bool first_order_input_hold)
{

  LKF_structure fw = lkf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                      ode_timestep_size, any_available_obs, non_na_ids,
                                      first_order_input_hold);
  const MatrixXd& postMat = fw.postMat;
  const MatrixXd& Ahat = fw.Ahat;
  const MatrixXd& Bhat = fw.Bhat;
  const MatrixXd& BhatFOH = fw.BhatFOH;
  const MatrixXd& Vhat = fw.Vhat;
  const MatrixXd Ahat_T = Ahat.transpose();
  const double dt_timestep = fw.dt_timestep;

  // constants
  const int ni = inputMat.cols();
  const int ni1 = ni + 1;
  const int n = stateVec.size();
  const int n_squared = n*n;

  // pre-allocate
  List xk(last_pred_id);
  MatrixXd predMat = MatrixXd::Zero(k_step_ahead+1, n + n_squared);
  VectorXd constant_and_inputVec = VectorXd::Ones(ni1);   // First element must be 1
  VectorXd constant_and_dinputVec = VectorXd::Zero(ni1);  // First element must be 0

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){

    // Get posteriors as starting point for k-ahead predictions
    stateVec = postMat.row(i).head(n);
    covMat = postMat.row(i).tail(n_squared).reshaped(n,n);

    // Store first prediction
    predMat.row(0).head(n) = stateVec;
    predMat.row(0).tail(n_squared) = covMat.reshaped();

    // k-step-ahead loop
    for(int k=0 ; k < k_step_ahead ; k++){
      constant_and_inputVec.tail(ni) = inputMat.row(i+k);

      // Compute next state value with zero order hold and add first order term if needed
      stateVec = Ahat * stateVec + Bhat * constant_and_inputVec;
      if(first_order_input_hold){
        constant_and_dinputVec.tail(ni) = (inputMat.row(i+k+1) - inputMat.row(i+k)) / dt_timestep;
        stateVec += BhatFOH * constant_and_dinputVec;
      }
      // Compute next covariance
      covMat = Ahat * covMat * Ahat_T + Vhat;

      // Save the k'th prediction
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
List lkf_simulate_rcpp(
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
    int ng,
    int last_pred_id,
    int k_step_ahead,
    int nsims,
    Nullable<int> seed,
    bool first_order_input_hold)
{

  // Set simulating seed if seed is not NULL
  set_simulation_seed(seed, ziggurat_states);

  auto f__ = get_funptr<funPtr_vec_const>(funPtrs, "f_const");
  auto g__ = get_funptr<funPtr_mat_const>(funPtrs, "g_const");

  LKF_structure fw = lkf_forward_pass(funPtrs, obsMat, inputMat, parVec, covMat, stateVec,
                                      ode_timestep_size, any_available_obs, non_na_ids,
                                      first_order_input_hold);
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
