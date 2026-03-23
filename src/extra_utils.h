#ifndef EXTRA_UTILS_H
#define EXTRA_UTILS_H

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

inline int choose_solver(CharacterVector ode_solver) {
  CharacterVector x1 = "euler";
  CharacterVector x2 = "rk4";
  CharacterVector x3 = "implicit_euler";

  bool y1 = Rcpp::is_true(Rcpp::all(ode_solver == x1));
  bool y2 = Rcpp::is_true(Rcpp::all(ode_solver == x2));
  bool y3 = Rcpp::is_true(Rcpp::all(ode_solver == x3));

  if(y1){
    // euler method
    return 1;
  } else if (y2){
    // runge-kutta 4th order method
    return 2;
  } else if (y3){
    // implicit euler method requires a newton's solver, and AD for drift function f(stateVec<AD>, parVec, inputVec)
    return 3;
  } else {
    Rcpp::stop("The implicit euler solver is not implemented");
  }
}

inline Eigen::MatrixXd construct_permutation_matrix(
  int number_of_available_obs, 
  int number_of_obs_eqs, 
  Eigen::VectorXi bool_is_not_na_obsVec)
{
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

inline Eigen::VectorXd remove_NAs(
  Eigen::VectorXd obsVec, 
  int number_of_available_obs, 
  Eigen::VectorXi bool_is_not_na_obsVec)
{
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

#endif