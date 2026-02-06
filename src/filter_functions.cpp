// #include <Rcpp.h>
#include <RcppArmadillo.h> // Must be included before RcppEigen and must not include Rcpp.h
#include <RcppEigen.h>

#include "function_typedefs.h"
#include "misc_helpers.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

MatrixXd cbind_time_column(const VectorXd& time, const MatrixXd& A){
  MatrixXd B(A.rows(), A.cols()+1);
  B.rightCols(A.cols()) = A;
  B.col(0) = time;
  return B;
}

// [[Rcpp::export]]
List calculate_filtering_observations(
  const List filtration_raw,
  const List funPtrs,
  Eigen::MatrixXd inputMat,
  Eigen::VectorXd parVec,
  int n_obs)
{

  auto h = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto dhdx = get_funptr<funPtr_mat_const>(funPtrs, "dhdx_const");
  auto hvar = get_funptr<funPtr_mat_const>(funPtrs, "hvar_const");

  int N = inputMat.rows();
  Eigen::VectorXd stateVec, inputVec;
  Eigen::MatrixXd covMat;
  Eigen::VectorXd H;
  Eigen::MatrixXd C,V,R;

  // Pre-allocate lists and their names
  List outlist = List::create(
    Named("mean") = R_NilValue, 
    Named("sd") = R_NilValue,
    Named("cov") = R_NilValue
    );
  List mean_list = List::create(Named("prior") = R_NilValue, Named("posterior") = R_NilValue);
  List sd_list = List::create(Named("prior") = R_NilValue, Named("posterior") = R_NilValue);
  List cov_list = List::create(Named("prior") = R_NilValue, Named("posterior") = R_NilValue);

  // Pre-allocate storage for posterior / prior
  MatrixXd mean_prior(N, n_obs), sd_prior(N, n_obs), mean_post(N, n_obs), sd_post(N, n_obs);
  List cov_prior(N), cov_post(N);

  List stateList_prior = filtration_raw["xPrior"];
  List covList_prior = filtration_raw["pPrior"];
  List stateList_post = filtration_raw["xPost"];
  List covList_post = filtration_raw["pPost"];

  for(int i=0; i<N; i++){
    inputVec = inputMat.row(i);
    stateVec = stateList_prior(i);
    covMat = covList_prior(i);
    H = h(stateVec, parVec, inputVec);
    C = dhdx(stateVec, parVec, inputVec);
    V = hvar(stateVec, parVec, inputVec);
    R = C * covMat * C.transpose() + V;
    mean_prior.row(i) = H;
    sd_prior.row(i) = R.diagonal().array().sqrt();
    cov_prior(i) = R;

    // posteriors
    stateVec = stateList_post(i);
    covMat = covList_post(i);
    H = h(stateVec, parVec, inputVec);
    C = dhdx(stateVec, parVec, inputVec);
    V = hvar(stateVec, parVec, inputVec);
    R = C * covMat * C.transpose() + V;
    // 
    mean_post.row(i) = H;
    sd_post.row(i) = R.diagonal().array().sqrt();
    cov_post(i) = R;
  }
  mean_list(0) = cbind_time_column(inputMat.col(0), mean_prior);
  sd_list(0) = cbind_time_column(inputMat.col(0), sd_prior);
  cov_list(0) = cov_prior;
  mean_list(1) = cbind_time_column(inputMat.col(0), mean_post);
  sd_list(1) = cbind_time_column(inputMat.col(0), sd_post);
  cov_list(1) = cov_post;

  outlist(0) = mean_list;
  outlist(1) = sd_list;
  outlist(2) = cov_list;

  return outlist;
}
