#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_funs2.h"
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

Eigen::MatrixXd construct_permutation_matrix2(int number_of_available_obs, int number_of_obs_eqs, Eigen::VectorXi bool_is_not_na_obsVec){
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

Eigen::VectorXd remove_NAs2(Eigen::VectorXd obsVec, int number_of_available_obs, Eigen::VectorXi bool_is_not_na_obsVec){
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
template <typename T1, typename T2>
List ekf_prediction2(
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
  const int last_pred_id,
  const int k_step_ahead,
  const int ode_solver)
{

  // misc  
  const int ni = inputMat.row(0).size();
  const int n_squared = n*n;
  Eigen::VectorXd inputVec(ni), dinputVec(ni), obsVec(m), e, y, obsVec_all;
  Eigen::VectorXi bool_is_not_na_obsVec;
  Eigen::MatrixXd C, R, K, E, V, Ri, I(n,n);
  I.setIdentity();
  Eigen::VectorXd H;
  Eigen::MatrixXd Hvar, dHdX;

  // storage for predictions
  //Rcpp::List xk_temp(k_step_ahead+1), pk_temp(k_step_ahead+1), xk(last_pred_id), pk(last_pred_id);
  Rcpp::List xk(last_pred_id), ode_1step_integration(2);

  Eigen::MatrixXd predMat(k_step_ahead+1, n + n_squared);
  predMat.setZero();

  /*We can save the predictions in a matrix of columns:
  n.states + n.states^2, and k.step.ahead+1 rows
  */
  //Eigen::MatrixXd xk_temp_2(k_step_ahead+1, n + n*n);

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
  }

  //////////// MAIN LOOP OVER TIME POINTS ///////////
  for(int i=0 ; i < last_pred_id ; i++){
    //xk_temp[0] = stateVec;
    //pk_temp[0] = covMat;
    predMat.row(0).head(n) = stateVec;
    predMat.row(0).tail(n_squared) = covMat.reshaped();

    //////////// K-STEP-AHEAD LOOP ///////////
    for(int k=0 ; k < k_step_ahead ; k++){
      inputVec = inputMat.row(i+k);
      dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k))/ode_timesteps(i+k);

      //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
      // Integrate moments between two time points
      for(int j=0 ; j < ode_timesteps(i+k) ; j++){
        ode_1step_integration = ode_integrator2(f__, g__, dfdx__, covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size(i+k), ode_solver);
        stateVec = ode_1step_integration["X1"];
        covMat = ode_1step_integration["P1"];
        inputVec += dinputVec;
      }

      // Save the prediction
      //xk_temp[k+1] = stateVec;
      //pk_temp[k+1] = covMat;
      predMat.row(k+1).head(n) = stateVec;
      predMat.row(k+1).tail(n_squared) = covMat.reshaped();
    }

    // Save a clone of the lists (if we dont use clone then all entries in xk and pk are identical,
    // because of the way pointers work here)
    //xk[i] = clone(xk_temp);
    //pk[i] = clone(pk_temp);
    xk[i] = predMat;

    //////////// DATA-UPDATE ///////////
    // Grab the one-step-ahead prediction (prior prediction)
    //stateVec = xk_temp[1];
    //covMat = pk_temp[1];
    stateVec = predMat.row(1).head(n);
    covMat = predMat.row(1).tail(n_squared).reshaped(n,n);

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
    }
  }
  
  /*return List::create(
    Named("Xpred")=xk , 
    Named("Ppred")=pk
    );*/
  return List::create(
    Named("predMats") = xk
    );
}

// Function typedefs
typedef Eigen::VectorXd (*funPtr_vec)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
typedef Eigen::MatrixXd (*funPtr_mat)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

// Function exported to R that performs (Extended Kalman) filtering
// [[Rcpp::export]]
List execute_ekf_prediction2(
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

  funPtr_vec     f__ = *XPtr<funPtr_vec>(f__R);
  funPtr_mat     g__ = *XPtr<funPtr_mat>(g__R);
  funPtr_mat  dfdx__ = *XPtr<funPtr_mat>(dfdx__R);
  funPtr_vec     h__ = *XPtr<funPtr_vec>(h__R);
  funPtr_mat  dhdx__ = *XPtr<funPtr_mat>(dhdx__R);
  funPtr_mat  hvar__ = *XPtr<funPtr_mat>(hvar__R);

  return ekf_prediction2<funPtr_vec, funPtr_mat>(
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