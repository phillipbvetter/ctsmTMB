// #include <Rcpp.h>
#include <RcppEigen.h>

#include "function_typedefs.h"
#include "misc_helpers.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// // [[Rcpp::export]]
// List calculate_filtering_observations(
//   List funPtrs,
//   Eigen::MatrixXd priorMat,
//   Eigen::MatrixXd postMat,
//   Eigen::MatrixXd inputMat,
//   Eigen::VectorXd parVec,
//   int m,
//   int n,
//   bool observation_noise)
// {
//
//   auto h = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
//   auto dhdx = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
//   auto hvar = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");
//
//   // dimensions
//   const int n_squared = n * n;
//   const int m_squared = m * m;
//   const int tsize = inputMat.rows();
//
//   // Pre-allocate
//   VectorXd stateVec, inputVec, H;
//   Eigen::MatrixXd covMat, C, V, R;
//   MatrixXd priorObsMat(tsize, m + m_squared), postObsMat(tsize, m + m_squared);
//
//   // compute observation mean and covariance for every row of a state matrix
//   auto observation_moments = [&](const MatrixXd& stateMat, MatrixXd& obsMat) {
//     for (int i = 0; i < tsize; i++) {
//       // Arguments
//       stateVec = stateMat.row(i).head(n);
//       covMat = stateMat.row(i).tail(n_squared).reshaped(n, n);
//       inputVec = inputMat.row(i);
//       // Mean
//       H = h(stateVec, parVec, inputVec);
//       // Covariance
//       C = dhdx(stateVec, parVec, inputVec);
//       R = C * covMat * C.transpose();
//       if(observation_noise) {
//         V = hvar(stateVec, parVec, inputVec);
//         R += V;
//       }
//       // Fill
//       obsMat.row(i).head(m) = H;
//       obsMat.row(i).tail(m_squared) = R.reshaped();
//     }
//   };
//
//   // Calculate prior and posterior observations
//   observation_moments(priorMat, priorObsMat);
//   observation_moments(postMat, postObsMat);
//
//   // Return
//   return List::create(
//     Named("priorObsMat") = priorObsMat,
//     Named("postObsMat") = postObsMat
//   );
//
// }

// [[Rcpp::export]]
List from_states_to_observations(
    List funPtrs,
    Eigen::MatrixXd stateMat,
    Eigen::MatrixXd inputMat,
    Eigen::VectorXd parVec,
    int m,
    int n,
    bool observation_noise)
{

  auto h = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto dhdx = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");

  // dimensions
  const int n_squared = n * n;
  const int m_squared = m * m;
  const int tsize = inputMat.rows();

  // Pre-allocate
  VectorXd stateVec, inputVec, H;
  MatrixXd covMat, C, V, R;
  MatrixXd obsMat(tsize, m + m_squared);

  // compute observation mean and covariance for every row of a state matrix
  auto observation_moments = [&](const MatrixXd& stateMat, MatrixXd& obsMat) {
    for (int i = 0; i < tsize; i++) {
      // Arguments
      stateVec = stateMat.row(i).head(n);
      covMat = stateMat.row(i).tail(n_squared).reshaped(n, n);
      inputVec = inputMat.row(i);
      // Mean
      H = h(stateVec, parVec, inputVec);
      // Covariance
      C = dhdx(stateVec, parVec, inputVec);
      R = C * covMat * C.transpose();
      if(observation_noise) {
        V = hvar(stateVec, parVec, inputVec);
        R += V;
      }
      // Fill
      obsMat.row(i).head(m) = H;
      obsMat.row(i).tail(m_squared) = R.reshaped();
    }
  };

  // Calculate prior and posterior observations
  observation_moments(stateMat, obsMat);

  // Return
  return List::create(
    Named("obsMat") = obsMat
  );

}

// [[Rcpp::export]]
List from_states_to_observations_no_state_cross_covariances(
    List funPtrs,
    Eigen::MatrixXd stateMat,
    Eigen::MatrixXd inputMat,
    Eigen::VectorXd parVec,
    int m,
    int n,
    bool observation_noise)
{

  auto h = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto dhdx = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");

  // dimensions
  const int m_squared = m*m;
  const int tsize = inputMat.rows();

  // Pre-allocate
  VectorXd stateVec, inputVec, H;
  MatrixXd C, V, R;
  MatrixXd diag_covMat;
  MatrixXd obsMat(tsize, m + m_squared);

  // compute observation mean and covariance for every row of a state matrix
  auto observation_moments = [&](const MatrixXd& stateMat, MatrixXd& obsMat) {
    for (int i = 0; i < tsize; i++) {
      // Arguments
      stateVec = stateMat.row(i).head(n);
      diag_covMat = stateMat.row(i).tail(n).asDiagonal();
      inputVec = inputMat.row(i);
      // Mean
      H = h(stateVec, parVec, inputVec);
      // Covariance
      C = dhdx(stateVec, parVec, inputVec);
      R = C * diag_covMat * C.transpose();
      if(observation_noise) {
        V = hvar(stateVec, parVec, inputVec);
        R += V;
      }
      // Fill
      obsMat.row(i).head(m) = H;
      obsMat.row(i).tail(m_squared) = R.reshaped();
    }
  };

  // Calculate prior and posterior observations
  observation_moments(stateMat, obsMat);

  // Return
  return List::create(
    Named("obsMat") = obsMat
  );

}
