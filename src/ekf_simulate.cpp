#include <Rcpp.h>
#include <RcppEigen.h>
#include <Ziggurat.h>
#include "helper_funs2.h"
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

//  This is the predict kalman filter function
template<typename T1, typename T2>
List ekf_simulation2(
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
  Eigen::VectorXd simulation_timestep_size,
  Eigen::VectorXd simulation_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  const int n,
  const int m,
  const int ng,
  const int last_pred_id,
  const int k_step_ahead,
  const int ode_solver,
  const int nsims)
{

  // misc  
  const int ni = inputMat.row(0).size();
  VectorXd inputVec(ni), dinputVec(ni), obsVec(m), e, y, obsVec_all;
  VectorXi bool_is_not_na_obsVec;
  MatrixXd C, R, K, E, V, Ri, I(n,n), stateMat(nsims,n), randN(n, nsims);
  I.setIdentity();
  VectorXd H;
  MatrixXd Hvar, dHdX;

  // storage for predictions
  List xk_simulate(last_pred_id), xk_simulate_temp(k_step_ahead+1), ode_1step_integration(2);

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

    // Call funs and transform from Rcpp to Eigen
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

    //////////// SIMULATION START ///////////
    //////////// SIMULATION START ///////////
    //////////// SIMULATION START ///////////

    /* 
    We draw from a multivariate normal i.e
    z = u + A*dB
    where u is the mean (stateVec) and A is cholesky factor of covariance matrix (sqrt(covMat))
    and dB is i.d.d normal vector
    We do simultaneously for all #nsims simulations, so dB is here a matrix of #nsims i.d.d vectors
    and similarly u is repeated with replicate
    */
    for(int j=0; j < n; j++){
      for(int k=0; k < nsims; k++){
        randN(j,k) = zigg.norm();
      }
    }
    stateMat = (stateVec.replicate(1, nsims) + covMat.llt().matrixL() * randN).transpose(); 
    xk_simulate_temp[0] = stateMat;

    // K-Step-Simulation - Euler Maruyama Scheme - First order input interpolation
    for(int k=0 ; k < k_step_ahead ; k++){
      inputVec = inputMat.row(i+k);
      dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k))/simulation_timesteps(i+k);

      for(int j=0 ; j < simulation_timesteps(i+k) ; j++){
        stateMat = euler_maruyama_simulation2(f__, g__, stateMat, parVec, inputVec, simulation_timestep_size(i+k), nsims, n, ng);
        inputVec += dinputVec;
      }

      // store the simulation
      xk_simulate_temp[k+1] = stateMat;
    }

    // Save a clone of the temporary list (must use clone)
    xk_simulate[i] = clone(xk_simulate_temp);

    //////////// SIMULATION END ///////////
    //////////// SIMULATION END ///////////
    //////////// SIMULATION END ///////////

    //////////// 1 STEP PREDICTION START ///////////
    //////////// 1 STEP PREDICTION START ///////////
    //////////// 1 STEP PREDICTION START ///////////

    //////////// TIME-UPDATE ///////////
    // 1-Step-Prediction - Solve Moment ODEs - First-Order Input Interpolation
    inputVec = inputMat.row(i);
    dinputVec = (inputMat.row(i+1) - inputMat.row(i))/ode_timesteps(i);
    for(int j=0 ; j < ode_timesteps(i) ; j++){
      ode_1step_integration = ode_integrator2(f__, g__, dfdx__, covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size(i), ode_solver);
      stateVec = ode_1step_integration["X1"];
      covMat = ode_1step_integration["P1"];
      inputVec += dinputVec;
    }

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

      // Call funs and transform from Rcpp to Eigen
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

      //////////// 1 STEP PREDICTION END ///////////
      //////////// 1 STEP PREDICTION END ///////////
      //////////// 1 STEP PREDICTION END ///////////
  }

  // Return
  return xk_simulate;
}

// Function typedefs
typedef Eigen::VectorXd (*funPtr_vec)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
typedef Eigen::MatrixXd (*funPtr_mat)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

// Function exported to R that performs stochastic simulations with (Extended Kalman) filtering updates
// [[Rcpp::export]]
List execute_ekf_simulation2(
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
  Eigen::VectorXd simulation_timestep_size,
  Eigen::VectorXd simulation_timesteps,
  Eigen::MatrixXi bool_is_not_na_obsMat,
  Eigen::VectorXi number_of_available_obs,
  int n,
  int m,
  int ng,
  int last_pred_id,
  int k_step_ahead,
  int ode_solver,
  int nsims)
{

  funPtr_vec     f__ = *XPtr<funPtr_vec>(f__R);
  funPtr_mat     g__ = *XPtr<funPtr_mat>(g__R);
  funPtr_mat  dfdx__ = *XPtr<funPtr_mat>(dfdx__R);
  funPtr_vec     h__ = *XPtr<funPtr_vec>(h__R);
  funPtr_mat  dhdx__ = *XPtr<funPtr_mat>(dhdx__R);
  funPtr_mat  hvar__ = *XPtr<funPtr_mat>(hvar__R);

  return ekf_simulation2<funPtr_vec, funPtr_mat>(
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
    simulation_timestep_size,
    simulation_timesteps,
    bool_is_not_na_obsMat,
    number_of_available_obs,
    n,
    m,
    ng,
    last_pred_id,
    k_step_ahead,
    ode_solver,
    nsims);
}
