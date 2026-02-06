#ifndef EXTRA_UTILS_H
#define EXTRA_UTILS_H

#include <Rcpp.h>
#include <RcppEigen.h>

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