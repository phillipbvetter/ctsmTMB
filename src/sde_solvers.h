#ifndef SDE_SOLVERS_H
#define SDE_SOLVERS_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <zigg/header>

#include "function_typedefs.h"
#include "ziggurat_seeder.h"

using namespace Rcpp;
using namespace Eigen;

// This ziggurat instance is defined in ekf_methods.cpp and used in lkf_methods.cpp and ukf_methods.cpp
extern zigg::Ziggurat ziggurat_states;

// Euler-maruyama simulation scheme for 1-step for all simulations
template<typename VectorTypeDef, typename MatrixTypeDef>
void euler_maruyama_simulation_inplace(
  VectorTypeDef f__,
  MatrixTypeDef g__,
  Eigen::MatrixXd& stateMat,
  Eigen::VectorXd& parVec, 
  Eigen::VectorXd& inputVec, 
  const double timestep,
  const int nsims,
  const int n,
  const int ng)
{
  // We perform one step euler-maruyama for each simulation vector (i.e each row in stateMat)
  // This produces the next step of nsims vectors (1 step ahead) called stateMat_next

  // Pre-allocate
  Eigen::MatrixXd G(ng,n);
  Eigen::VectorXd F(n), dW(ng);
  double sqrt_timestep = sqrt(timestep);
  Eigen::VectorXd stateVec;

  for(int i=0; i < nsims; i++){

    // Grab brownian increments
    for(int j=0; j < ng; j++){
      dW(j) = ziggurat_states.rnorm();
    }

    // Grab initial state vector
    stateVec = stateMat.row(i);

    // Calculate drift and diffusion
    F = f__(stateVec, parVec, inputVec);
    G = g__(stateVec, parVec, inputVec);

    // Perform euler-maruyama step'
    stateVec += F * timestep;
    stateVec += G * dW * sqrt_timestep;
    stateMat.row(i) = stateVec;
    //stateMat.row(i) = stateVec + F * timestep + G * dW * sqrt_timestep;
  }

}

// Euler-maruyama simulation scheme for 1-step for all simulations
/*inline Eigen::MatrixXd euler_maruyama_simulation(
  funPtr_vec f__,
  funPtr_mat g__,
  Eigen::MatrixXd stateMat, 
  Eigen::VectorXd parVec, 
  Eigen::VectorXd inputVec, 
  double timestep, 
  int nsims,
  int n,
  int ng)
{
  // Returns a matrix with nsims rows and columns equal to the number of system states - each row corresponds to taking a single
  // Euler-Maruyama step
  Eigen::MatrixXd stateMat_next(nsims, n), G;
  Eigen::VectorXd stateVec, dW(ng), F;
  double sqrt_timestep = sqrt(timestep);
  
  // Perform one-step simulation for each row in stateMat
  for(int i=0; i < nsims; i++){

    // extract state values
    stateVec = stateMat.row(i);

    // Calculate drift and diffusion
    F = f__(stateVec, parVec, inputVec);
    G = g__(stateVec, parVec, inputVec);

    // Generate dW vector by sampling from standard normal
    for(int i=0; i < ng; i++){
      dW(i) = ziggurat_states.rnorm();
    }

    // Perform euler-maruyama step
    stateMat_next.row(i) =  stateVec + F * timestep + G * dW * sqrt_timestep;
  }

  // Return
  return stateMat_next;
}*/

#endif