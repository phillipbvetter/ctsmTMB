#ifndef ODE_SOLVERS_H
#define ODE_SOLVERS_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <zigg/header>

#include "function_typedefs.h"

using namespace Rcpp;
using namespace Eigen;

// Solving the variance moment differential equation 1-step forward
template <typename MatrixTypeDef>
void cov_ode_1step_inplace(
  MatrixTypeDef g__,
  MatrixTypeDef dfdx__,
  const Eigen::MatrixXd& covMat, 
  const Eigen::VectorXd& stateVec,
  const Eigen::VectorXd& parVec,
  const Eigen::VectorXd& inputVec,
  Eigen::MatrixXd& c)
{

  // Compute drift jacobian and diffusion 
  Eigen::MatrixXd A_covMat = dfdx__(stateVec, parVec, inputVec) * covMat;
  Eigen::MatrixXd G = g__(stateVec, parVec, inputVec);

  // Calculate the right-hand side of covariance moment differential equation
  c = A_covMat + A_covMat.transpose() + G * G.transpose();

}

template <typename VectorTypeDef, typename MatrixTypeDef>
void ode_integrator_inplace(
  VectorTypeDef f__,
  MatrixTypeDef g__,
  MatrixTypeDef dfdx__,
  Eigen::MatrixXd& covMat,
  Eigen::VectorXd& stateVec,
  Eigen::VectorXd parVec,
  Eigen::VectorXd inputVec,
  Eigen::VectorXd dinputVec,
  const double dt,
  const int ode_solver,
  Eigen::VectorXd& k1,
  Eigen::VectorXd& k2,
  Eigen::VectorXd& k3,
  Eigen::VectorXd& k4,
  Eigen::MatrixXd& c1,
  Eigen::MatrixXd& c2,
  Eigen::MatrixXd& c3,
  Eigen::MatrixXd& c4)
{

  if (ode_solver==1) {
  // Forward Euler

    k1 = f__(stateVec, parVec, inputVec);
    cov_ode_1step_inplace(g__, dfdx__, covMat, stateVec, parVec, inputVec, c1);

    stateVec += dt * k1;
    covMat   += dt * c1;

    return;

  } else if (ode_solver==2){
  // Runge Kutta 4th Order
    
    const Eigen::VectorXd X0 = stateVec;
    const Eigen::MatrixXd P0 = covMat;

    // RK4
    // step 1
    k1 = f__(stateVec, parVec, inputVec);
    cov_ode_1step_inplace(g__, dfdx__, covMat, stateVec, parVec, inputVec, c1);

    // step 2
    inputVec += 0.5 * dinputVec;
    stateVec = X0 + 0.5 * dt * k1;
    covMat   = P0 + 0.5 * dt * c1;
    k2 = f__(stateVec, parVec, inputVec);
    cov_ode_1step_inplace(g__, dfdx__, covMat, stateVec, parVec, inputVec, c2);

    // step 3
    stateVec = X0 + 0.5 * dt * k2;
    covMat   = P0 + 0.5 * dt * c2;
    k3 = f__(stateVec, parVec, inputVec);
    cov_ode_1step_inplace(g__, dfdx__, covMat, stateVec, parVec, inputVec, c3);

    // step 4
    inputVec += 0.5 * dinputVec;
    stateVec = X0 + dt * k3;
    covMat   = P0 + dt * c3;
    k4 = f__(stateVec, parVec, inputVec);
    cov_ode_1step_inplace(g__, dfdx__, covMat, stateVec, parVec, inputVec, c4);

    // solution
    stateVec = X0 + dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    covMat   = P0 + dt * (c1 + 2*c2 + 2*c3 + c4) / 6.0;
  } else {
    // Do nothing
  }
}

#endif