#ifndef HELPER_FUNS2_H
#define HELPER_FUNS2_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <zigg/header>  
static zigg::Ziggurat ziggurat;
const double pi = M_PI;

// Inverse logit function
// double invlogit2(double x);

// Note to self:
// The use of 'inline' and header guards allow (recommended only for small utility) functions to be placed in this header file
// as opposed to placing the full function in a .cpp file and placing a function declaration here (as commented out below)
// function for constructing permutation matrix that removes NA data
// Eigen::MatrixXd construct_permutation_matrix2(int number_of_available_obs, int number_of_obs_eqs, Eigen::VectorXi bool_is_not_na_obsVec);
// Function for removing NA entries from data vector
// Eigen::VectorXd remove_NAs2(Eigen::VectorXd obsVec, int number_of_available_obs, Eigen::VectorXi bool_is_not_na_obsVec);
inline Eigen::MatrixXd construct_permutation_matrix2(int number_of_available_obs, int number_of_obs_eqs, Eigen::VectorXi bool_is_not_na_obsVec){
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

inline Eigen::VectorXd remove_NAs2(Eigen::VectorXd obsVec, int number_of_available_obs, Eigen::VectorXi bool_is_not_na_obsVec){
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

// Solving the variance moment differential equation 1-step forward
template<typename T2>
Eigen::MatrixXd cov_ode_1step2(
  T2 g__,
  T2 dfdx__,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd parVec,
  Eigen::VectorXd inputVec)
{

  // Compute drift jacobian and diffusion 
  Eigen::MatrixXd A = dfdx__(stateVec, parVec, inputVec);
  Eigen::MatrixXd G = g__(stateVec, parVec, inputVec);

  // Calculate the right-hand side of covariance moment differential equation
  Eigen::MatrixXd cov_ode_1step = A * covMat + covMat * A.transpose() + G * G.transpose();

  // Return
  return cov_ode_1step;
}

// forward-euler and rk4 ode integrator function
template <typename T1, typename T2>
Rcpp::List ode_integrator2(
  T1 f__, 
  T2 g__,
  T2 dfdx__, 
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec, 
  Eigen::VectorXd parVec, 
  Eigen::VectorXd inputVec,
  Eigen::VectorXd dinputVec, 
  double dt, 
  int ode_solver)
{

  // Initialize return values
  Eigen::VectorXd X1;
  Eigen::MatrixXd P1;

  // Forward Euler
  if(ode_solver == 1)
  {
    Eigen::VectorXd F = f__(stateVec, parVec, inputVec);
    Eigen::MatrixXd P_rhs = cov_ode_1step2(g__, dfdx__, covMat, stateVec, parVec, inputVec);

    X1 = stateVec + F * dt;
    P1 = covMat + P_rhs * dt;
  }

  // 4th Order Runge-Kutta
  if(ode_solver==2)
  {

    // Initialize
    Eigen::VectorXd X0 = stateVec;
    Eigen::MatrixXd P0 = covMat;
    Eigen::VectorXd k1, k2, k3, k4;
    Eigen::MatrixXd c1,c2,c3,c4;

    // /*1. Approx Slope at Initial Point*/
    k1 = f__(stateVec, parVec, inputVec);
    c1 = cov_ode_1step2(g__, dfdx__, covMat, stateVec, parVec, inputVec);

    /*2. First Approx Slope at Midpoint*/
    inputVec += 0.5 * dinputVec;
    stateVec = X0 + 0.5 * dt * k1;
    covMat   = P0 + 0.5 * dt * c1;
    k2  = f__(stateVec, parVec, inputVec); 
    c2  = cov_ode_1step2(g__, dfdx__, covMat, stateVec, parVec, inputVec);

    /*3. Second Approx Slope at Midpoint*/
    stateVec = X0 + 0.5 * dt * k2;
    covMat   = P0 + 0.5 * dt * c2;
    k3 = f__(stateVec, parVec, inputVec);
    c3 = cov_ode_1step2(g__, dfdx__, covMat, stateVec, parVec, inputVec);

    /*4. Approx Slope at End Point*/
    inputVec += 0.5 * dinputVec;
    stateVec = X0 + dt * k3;
    covMat   = P0 + dt * c3;
    k4 = f__(stateVec, parVec, inputVec); 
    c4 = cov_ode_1step2(g__, dfdx__, covMat, stateVec, parVec, inputVec);

    /*ODE UPDATE*/
    X1 = X0 + dt * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    P1 = P0 + dt * (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0;
  }

  // Return
  return Rcpp::List::create(Rcpp::Named("X1")=X1, Rcpp::Named("P1")=P1);
}

// Euler-maruyama simulation scheme for 1-step for all simulations
template<typename T1, typename T2>
Eigen::MatrixXd euler_maruyama_simulation2(
  T1 f__, 
  T2 g__, 
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

    // Convert F and G from Rcpp::Numeric to Eigen
    F = f__(stateVec, parVec, inputVec);
    G = g__(stateVec, parVec, inputVec);

    // Generate dW vector by sampling from standard normal
    for(int i=0; i < ng; i++){
      // dW(i) = zigg.norm();
      dW(i) = ziggurat.rnorm();
    }

    // Perform euler-maruyama step
    stateMat_next.row(i) =  stateVec + F * timestep + G * dW * sqrt_timestep;
  }

  // Return
  return stateMat_next;
}

#endif