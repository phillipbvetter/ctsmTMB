#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_funs2.h"
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

//  This is the predict kalman filter function
template <typename T1, typename T2>
List ekf_filter(
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
  const int n,
  const int m,
  const int ode_solver)
{

  // misc  
  const int tsize = inputMat.col(0).size();
  const int ni = inputMat.row(0).size();
  // const int n_squared = n*n;
  Eigen::VectorXd inputVec(ni), dinputVec(ni), obsVec(m), e, y, obsVec_all;
  Eigen::VectorXi bool_is_not_na_obsVec;
  Eigen::MatrixXd C, R, K, E, V, Ri, I(n,n);
  I.setIdentity();
  Eigen::VectorXd H;
  Eigen::MatrixXd Hvar, dHdX;

  // storage
  Rcpp::List ode_1step_integration(2);
  Rcpp::List Innovation(tsize), InnovationCovariance(tsize), xPrior(tsize), xPost(tsize), pPrior(tsize), pPost(tsize);
  //Eigen::MatrixXd xPrior(tsize, n), xPost(tsize,n);
  //Eigen::MatrixXd pPrior(tsize, n_squared), pPost(tsize, n_squared);

  /*xPrior.row(0) = stateVec;*/
  /*pPrior.row(0) = covMat.reshaped();*/
  xPrior(0) = stateVec;
  pPrior(0) = covMat;

  //////////// INITIAL DATA-UPDATE ///////////
  // Only update if there is any available data
  if( number_of_available_obs(0) > 0 ){
    // Get i+1 entries
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
  /*xPost.row(0) = stateVec;*/
  /*pPost.row(0) = covMat.reshaped();*/
  xPost(0) = stateVec;
  pPost(0) = covMat;

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < (tsize-1) ; i++){

    //////////// K-STEP-AHEAD LOOP ///////////
    inputVec = inputMat.row(i);
    dinputVec = (inputMat.row(i+1) - inputMat.row(i))/ode_timesteps(i);

    //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
    // Integrate moments between two time points
    for(int j=0 ; j < ode_timesteps(i) ; j++){
      ode_1step_integration = ode_integrator2(f__, g__, dfdx__, covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size(i), ode_solver);
      stateVec = ode_1step_integration["X1"];
      covMat = ode_1step_integration["P1"];
      inputVec += dinputVec;
    }
    /*xPrior.row(i+1) = stateVec;*/
    /*pPrior.row(i+1) = covMat.reshaped();*/
    xPrior(i+1) = stateVec;
    pPrior(i+1) = covMat;

    //////////// DATA-UPDATE ///////////
    // Only update if there is any available data
    if( number_of_available_obs(i+1) > 0 ){
      // Get i+1 entries
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
    /*xPost.row(i+1) = stateVec;*/
    /*pPost.row(i+1) = covMat.reshaped();*/
    xPost(i+1) = stateVec;
    pPost(i+1) = covMat;
  }
  
  /*return List::create(
    Named("Xpred")=xk , 
    Named("Ppred")=pk
    );*/
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
List ekf_filter_cpp(
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
  int n,
  int m,
  int ode_solver)
{

  funPtr_vec     f__ = *XPtr<funPtr_vec>(f__R);
  funPtr_mat     g__ = *XPtr<funPtr_mat>(g__R);
  funPtr_mat  dfdx__ = *XPtr<funPtr_mat>(dfdx__R);
  funPtr_vec     h__ = *XPtr<funPtr_vec>(h__R);
  funPtr_mat  dhdx__ = *XPtr<funPtr_mat>(dhdx__R);
  funPtr_mat  hvar__ = *XPtr<funPtr_mat>(hvar__R);

  return ekf_filter<funPtr_vec, funPtr_mat>(
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
    n,
    m,
    ode_solver);
}
