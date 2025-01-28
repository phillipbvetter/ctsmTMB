#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_funs.h"
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

//  This is the predict kalman filter function
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
List ekf_likelihood(
  T1 f__, 
  T2 g__,
  T3 dfdx__,
  T4 h__,
  T5 dhdx__,
  T6 hvar__,
  const Eigen::MatrixXd obsMat,
  const Eigen::MatrixXd inputMat,
  const Eigen::VectorXd parVec,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  const Eigen::VectorXd ode_timestep_size,
  const Eigen::VectorXd ode_timesteps,
  const Eigen::MatrixXi bool_is_not_na_obsMat,
  const Eigen::VectorXi number_of_available_obs,
  int n,
  int m,
  int ode_solver)
{

  // misc  
  int ni = inputMat.row(0).size();
  int tsize = inputMat.col(0).size();
  double nll = 0;
  double half_log2PI = 0.5 * log(2*M_PI);
  double numberOfObs;
  Eigen::VectorXd inputVec(ni), dinputVec(ni), obsVec(m), e, y, obsVec_all;
  Eigen::VectorXi bool_is_not_na_obsVec;
  Eigen::MatrixXd C, R, K, E, V, Ri, I(n,n);
  double Rdet;
  I.setIdentity();
  Rcpp::NumericVector H;
  Rcpp::NumericMatrix Hvar, dHdX;

  // storage for predictions
  //Rcpp::List xPrior(tsize),xPost(tsize),pPrior(tsize),pPost(tsize);
  Rcpp::List ode_1step_integration;

  //xPrior[0] = stateVec;
  //pPrior[0] = covMat;

  //////////// DATA-UPDATE ///////////
  // Only update if there is any available data
  if( number_of_available_obs(0) > 0 ){
    // Get i+1 entries
    inputVec = inputMat.row(0);
    obsVec_all = obsMat.row(0);
    bool_is_not_na_obsVec = bool_is_not_na_obsMat.row(0);

    // remove potential NA entries in obsVec and construct permutation matrix
    obsVec = remove_NAs(obsVec_all, number_of_available_obs(0), bool_is_not_na_obsVec);
    E = construct_permutation_matrix(number_of_available_obs(0), m, bool_is_not_na_obsVec);

    // Call funs and transform from Rcpp to Eigen
    H = h__(stateVec, parVec, inputVec);
    dHdX = dhdx__(stateVec, parVec, inputVec);
    Hvar = hvar__(stateVec, parVec, inputVec);
    Eigen::Map<Eigen::VectorXd> H_eigen = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(H);
    Eigen::Map<Eigen::MatrixXd> dHdX_eigen = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(dHdX);
    Eigen::Map<Eigen::MatrixXd> Hvar_eigen = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Hvar);
    
    // Kalman Filter
    C = E * dHdX_eigen;
    e = obsVec - E * H_eigen;
    V = E * Hvar_eigen * E.transpose();
    R = C * covMat * C.transpose() + V;
    Ri = R.inverse();
    K = covMat * C.transpose() * Ri;
    stateVec = stateVec + K * e;
    covMat = (I - K * C) * covMat * (I - K * C).transpose() + K*V*K.transpose();

    // Likelihood Contribution
    Rdet = R.llt().matrixL().determinant();
    numberOfObs = number_of_available_obs(0);
    nll += 0.5 * log(Rdet) + 0.5 * e.transpose() * Ri * e + half_log2PI * numberOfObs;
  }

  //xPost[0] = stateVec;
  //pPost[0] = covMat;

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < tsize-1 ; i++){
    inputVec = inputMat.row(i);
    dinputVec = (inputMat.row(i+1) - inputMat.row(i))/ode_timesteps(i);

    //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
    // Integrate moments between two time points
    for(int j=0 ; j < ode_timesteps(i) ; j++){
      ode_1step_integration = ode_integrator(f__, g__, dfdx__, covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size(i), ode_solver);
      stateVec = ode_1step_integration["X1"];
      covMat = ode_1step_integration["P1"];
      inputVec += dinputVec;
    }

    //xPrior[i+1] = stateVec; 
    //pPrior[i+1] = covMat;

    //////////// DATA-UPDATE ///////////
    // Only update if there is any available data
    if( number_of_available_obs(i+1) > 0 ){
      // Get i+1 entries
      inputVec = inputMat.row(i+1);
      obsVec_all = obsMat.row(i+1);
      bool_is_not_na_obsVec = bool_is_not_na_obsMat.row(i+1);

      // remove potential NA entries in obsVec and construct permutation matrix
      obsVec = remove_NAs(obsVec_all, number_of_available_obs(i+1), bool_is_not_na_obsVec);
      E = construct_permutation_matrix(number_of_available_obs(i+1), m, bool_is_not_na_obsVec);

      // Call funs and transform from Rcpp to Eigen
      H = h__(stateVec, parVec, inputVec);
      dHdX = dhdx__(stateVec, parVec, inputVec);
      Hvar = hvar__(stateVec, parVec, inputVec);
      Eigen::Map<Eigen::VectorXd> H_eigen = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(H);
      Eigen::Map<Eigen::MatrixXd> dHdX_eigen = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(dHdX);
      Eigen::Map<Eigen::MatrixXd> Hvar_eigen = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Hvar);
      
      // Kalman Filter
      C = E * dHdX_eigen;
      e = obsVec - E * H_eigen;
      V = E * Hvar_eigen * E.transpose();
      R = C * covMat * C.transpose() + V;
      Ri = R.inverse();
      K = covMat * C.transpose() * Ri;
      stateVec = stateVec + K * e;
      covMat = (I - K * C) * covMat * (I - K * C).transpose() + K*V*K.transpose();

      // Likelihood Contribution
    Rdet = R.llt().matrixL().determinant();
    numberOfObs = number_of_available_obs(i+1);
    nll += 0.5 * log(Rdet) + 0.5 * e.transpose() * Ri * e + half_log2PI * numberOfObs;
    }

    //xPost[i+1] = stateVec;
    //pPost[i+1] = covMat;
  }
  
  return List::create(
    //Named("xPrior") = xPrior, 
    //Named("xPost") = xPost, 
    //Named("pPrior") = pPrior,
    //Named("pPost") = pPost,
    Named("nll") = nll
    );
}

// Function typedefs
typedef SEXP (*funPtr)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

// Function exported to R that performs (Extended Kalman) filtering
// [[Rcpp::export]]
List ekf_rcpp_likelihood(
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

  funPtr     f__ = *XPtr<funPtr>(f__R);
  funPtr     g__ = *XPtr<funPtr>(g__R);
  funPtr  dfdx__ = *XPtr<funPtr>(dfdx__R);
  funPtr     h__ = *XPtr<funPtr>(h__R);
  funPtr  dhdx__ = *XPtr<funPtr>(dhdx__R);
  funPtr  hvar__ = *XPtr<funPtr>(hvar__R);

  return ekf_likelihood<funPtr, funPtr, funPtr, funPtr, funPtr, funPtr>(
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