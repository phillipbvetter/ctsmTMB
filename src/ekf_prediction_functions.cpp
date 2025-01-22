#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_funs.h"
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

Eigen::MatrixXd construct_permutation_matrix(int number_of_available_obs, int number_of_obs_eqs, Eigen::VectorXi bool_is_not_na_obsVec){
  Eigen::MatrixXd E(number_of_available_obs, number_of_obs_eqs);
  E.setZero();
  int j=0;
  for(int i=0; i < number_of_obs_eqs; i++){
    /*if p(i) is 1 then include by setting 1 in diagonal of matrix*/
    if(bool_is_not_na_obsVec(i) == 1){
      E(j,i) = 1.0;
      j += 1;
    }
  }
  return E;
}

Eigen::VectorXd remove_NAs(Eigen::VectorXd obsVec, int number_of_available_obs, Eigen::VectorXi bool_is_not_na_obsVec){
  // Initialize
  int ii = 0;
  Eigen::VectorXd obsVec_without_NAs(number_of_available_obs);
  // For loop to remove NA observations
  for(int i=0; i < obsVec.size(); i++){
    if( bool_is_not_na_obsVec(i) == 1){
      obsVec_without_NAs(ii) = obsVec(i);
      ii++;
    }
  }
  // Return
  return obsVec_without_NAs;
}


//  This is the predict kalman filter function
template <typename T>
List ekf_prediction(
  T f__, 
  T g__,
  T dfdx__,
  T h__,
  T dhdx__,
  T hvar__,
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
  int last_pred_id,
  int k_step_ahead,
  int ode_solver)
{

  // misc  
  int ni = inputMat.row(0).size();
  Eigen::VectorXd inputVec(ni), dinputVec(ni), obsVec(m), e, y, obsVec_all;
  Eigen::VectorXi bool_is_not_na_obsVec;
  Eigen::MatrixXd C, R, K, E, V, Ri, I(n,n);
  I.setIdentity();
  Rcpp::NumericVector H;
  Rcpp::NumericMatrix Hvar, dHdX;

  // storage for predictions
  Rcpp::List xk_temp(k_step_ahead+1), pk_temp(k_step_ahead+1), xk(last_pred_id), pk(last_pred_id);
  Rcpp::List ode_1step_integration;

  //////////// INITIAL DATA-UPDATE ///////////
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
  }

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){
    xk_temp[0] = stateVec;
    pk_temp[0] = covMat;

    //////////// K-STEP-AHEAD LOOP ///////////
    for(int k=0 ; k < k_step_ahead ; k++){
      inputVec = inputMat.row(i+k);
      dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k))/ode_timesteps(i+k);

      //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
      // Integrate moments between two time points
      for(int j=0 ; j < ode_timesteps(i+k) ; j++){
        ode_1step_integration = ode_integrator(f__, g__, dfdx__, covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size(i+k), ode_solver);
        stateVec = ode_1step_integration["X1"];
        covMat = ode_1step_integration["P1"];
        inputVec += dinputVec;
      }

      // Save the prediction
      xk_temp[k+1] = stateVec;
      pk_temp[k+1] = covMat;
    }

    // Save a clone of the lists (if we dont use clone then all entries in xk and pk are identical,
    // because of the way pointers work here)
    xk[i] = clone(xk_temp);
    pk[i] = clone(pk_temp);

    //////////// DATA-UPDATE ///////////
    // Grab the one-step-ahead prediction (prior prediction)
    stateVec = xk_temp[1];
    covMat = pk_temp[1];

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
    }
  }
  
  return List::create(
    Named("Xpred")=xk , 
    Named("Ppred")=pk
    );
}

// Function typedefs
typedef SEXP (*funPtr)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

//'@ export
// [[Rcpp::export]]
List execute_ekf_prediction(
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
  int last_pred_id,
  int k_step_ahead,
  int ode_solver)
{

  funPtr     f__ = *XPtr<funPtr>(f__R);
  funPtr     g__ = *XPtr<funPtr>(g__R);
  funPtr  dfdx__ = *XPtr<funPtr>(dfdx__R);
  funPtr     h__ = *XPtr<funPtr>(h__R);
  funPtr  dhdx__ = *XPtr<funPtr>(dhdx__R);
  funPtr  hvar__ = *XPtr<funPtr>(hvar__R);

  return ekf_prediction<funPtr>(
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
    last_pred_id,
    k_step_ahead,
    ode_solver);
}