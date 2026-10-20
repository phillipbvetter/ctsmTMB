#include <Rcpp.h>
#include <RcppEigen.h>

#include "function_typedefs.h"
#include "misc_helpers.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

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
  const int k_aheads = k_ahead + 1;
  int idx;

  VectorXd inputVec, stateVec, covVec;
  MatrixXd C, V, R;
  Eigen::MatrixXd outMat(predMat.rows(), n_obs + (compute_dispersion ? n_obs_sq : 0));
  outMat.setZero();

  for(int i = 0; i < last_pred_id; i++){
    for(int k = 0; k < k_aheads; k++){
      // Index and input
      idx = k + i*k_aheads;
      inputVec = inputMat.row(i+k);

      // Mean
      stateVec = predMat.row(idx).head(n_states);
      outMat.row(idx).head(n_obs) = h(stateVec, parVec, inputVec);

      // Covariance
      if(compute_dispersion){
        covVec = predMat.row(idx).tail(n_states_sq);
        Eigen::Map<const Eigen::MatrixXd> covMat(covVec.data(), n_states, n_states);
        C = dhdx(stateVec, parVec, inputVec);
        V = hvar(stateVec, parVec, inputVec);
        R = C * covMat * C.transpose() + V;
        outMat.row(idx).tail(n_obs_sq) = R.reshaped();
      }
    }
  }

  return outMat;
}
