#include <Rcpp.h>
#include <RcppEigen.h>
#include <Ziggurat.h>
#include "helper_funs.h"
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

//  This is the predict kalman filter function
template<typename T>
List ekf_simulation(
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

  // misc  
  int ni = inputMat.row(0).size();
  VectorXd inputVec(ni), dinputVec(ni), obsVec(m), e, y, obsVec_all;
  VectorXi bool_is_not_na_obsVec;
  MatrixXd C, R, K, E, V, Ri, I(n,n), stateMat(nsims,n), randN(n, nsims);
  I.setIdentity();
  NumericVector H;
  NumericMatrix Hvar, dHdX;

  // storage for predictions
  List xk_simulate(last_pred_id), xk_simulate_temp(k_step_ahead+1), ode_1step_integration;

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
        stateMat = euler_maruyama_simulation(f__, g__, stateMat, parVec, inputVec, simulation_timestep_size(i+k), nsims, n, ng);
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
      ode_1step_integration = ode_integrator(f__, g__, dfdx__, covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size(i), ode_solver);
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
      obsVec = remove_NAs(obsVec_all, number_of_available_obs(i+1), bool_is_not_na_obsVec);
      E = construct_permutation_matrix(number_of_available_obs(i+1), m, bool_is_not_na_obsVec);

      // Call funs and transform from Rcpp to Eigen
      H = h__(stateVec, parVec, inputVec);
      dHdX = dhdx__(stateVec, parVec, inputVec);
      Hvar = hvar__(stateVec, parVec, inputVec);
      Map<VectorXd> H_eigen = as<Map<VectorXd> >(H);
      Map<MatrixXd> dHdX_eigen = as<Map<MatrixXd> >(dHdX);
      Map<MatrixXd> Hvar_eigen = as<Map<MatrixXd> >(Hvar);

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

      //////////// 1 STEP PREDICTION END ///////////
      //////////// 1 STEP PREDICTION END ///////////
      //////////// 1 STEP PREDICTION END ///////////
  }

  // Return
  return xk_simulate;
}

// Function typedefs
typedef SEXP (*funPtr)(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

//'@ export
// [[Rcpp::export]]
List execute_ekf_simulation(
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

  funPtr     f__ = *XPtr<funPtr>(f__R);
  funPtr     g__ = *XPtr<funPtr>(g__R);
  funPtr  dfdx__ = *XPtr<funPtr>(dfdx__R);
  funPtr     h__ = *XPtr<funPtr>(h__R);
  funPtr  dhdx__ = *XPtr<funPtr>(dhdx__R);
  funPtr  hvar__ = *XPtr<funPtr>(hvar__R);

  return ekf_simulation<funPtr>(
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
