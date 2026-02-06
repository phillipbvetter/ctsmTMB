#include <Rcpp.h>
#include <RcppEigen.h>

#include "function_typedefs.h"
#include "misc_helpers.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

/*MatrixXd calculate_prediction_observations(
  const MatrixXd predMat,
  const List funPtrs,
  Eigen::MatrixXd & inputMat,
  Eigen::VectorXd & parVec,
  const int n_states,
  const int n_obs,
  const int last_pred_id,
  const int k_ahead)
{

  auto h = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto dhdx = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");

  const int n_states_squared = n_states * n_states;
  const int n_obs_squared = n_obs * n_obs;
  MatrixXd outMat(predMat.rows(), n_obs+n_obs_squared);

  for(int i=0; i < last_pred_id; i++){
    for(int k=0; k < k_ahead, k++){
      auto inputVec = inputMat.row(i+k);

      auto stateVec = predMat.row(i+k).head(n_states);
      auto covMat = predMat.row(i+k).tail(n_states_squared);

      auto H = h(stateVec, parVec, inputVec);
      auto C = dhdx(stateVec, parVec, inputVec);
      auto V = hvar(stateVec, parVec, inputVec);
      auto R = C * covMat * C.transpose() + V;

      outMat.row(i+k).head(n_obs) = H;
      outMat.row(i+k).tail(n_obs_squared) = R.reshaped();
    }
  }

  return outMat;
}*/

// [[Rcpp::export]]
Eigen::MatrixXd calculate_prediction_observations(
  const Eigen::MatrixXd& predMat,
  const List& funPtrs,
  const Eigen::MatrixXd& inputMat,
  const Eigen::VectorXd& parVec,
  const int n_states,
  const int n_obs,
  const int last_pred_id,
  const int k_ahead,
  bool compute_dispersion)
{
  auto h    = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto dhdx = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");

  const int n_states_sq = n_states * n_states;
  const int n_obs_sq    = n_obs * n_obs;

  VectorXd inputVec, stateVec, covVec, H;
  MatrixXd C, V, R;
  Eigen::MatrixXd outMat(predMat.rows(), n_obs + (compute_dispersion ? n_obs_sq : 0));

  for(int i = 0; i < last_pred_id; ++i){
    for(int k = 0; k < k_ahead; ++k){

      int idx = i + k;

      inputVec = inputMat.row(idx);
      stateVec = predMat.row(idx).head(n_states);
      covVec = predMat.row(idx).tail(n_states_sq);
      Eigen::Map<const Eigen::MatrixXd> covMat(covVec.data(), n_states, n_states);

      // Mean
      auto H = h(stateVec, parVec, inputVec);
      outMat.row(idx).head(n_obs) = H;
      if(compute_dispersion){
        // Covariance
        C = dhdx(stateVec, parVec, inputVec);
        V = hvar(stateVec, parVec, inputVec);
        R = C * covMat * C.transpose() + V;
        outMat.row(idx).tail(n_obs_sq) = R.reshaped();
      }
    }
  }

  return outMat;
}