#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

/* Constants */
const double pi = M_PI;

/* Custom Functions*/
double invlogit(double x){return 1.0/(1.0+exp(-x));}

/* --------------------------------------------
   Typedefs for function pointer signatures
---------------------------------------------*/
typedef Eigen::VectorXd (*funPtr_vec)(Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd);
typedef Eigen::MatrixXd (*funPtr_mat)(Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd);

/* --------------------------------------------
   User-generated functions
---------------------------------------------*/

// INSERT F

// INSERT DFDX

// INSERT G

// INSERT H

// INSERT DHDX

// INSERT HVAR

// INSERT DFDU

/* --------------------------------------------
   Dispatcher returning XPtrs
---------------------------------------------*/
// [[Rcpp::export]]
SEXP get_sysfun_cpp_function_ptrs()
{
  return List::create(
    Named("f") = XPtr<funPtr_vec>(new funPtr_vec(&f), false),
    Named("dfdx") = XPtr<funPtr_mat>(new funPtr_mat(&dfdx), false),
    Named("g") = XPtr<funPtr_mat>(new funPtr_mat(&g), false),
    Named("h") = XPtr<funPtr_vec>(new funPtr_vec(&h), false),
    Named("dhdx") = XPtr<funPtr_mat>(new funPtr_mat(&dhdx), false),
    Named("hvar") = XPtr<funPtr_mat>(new funPtr_mat(&hvar), false),
    Named("dfdu") = XPtr<funPtr_mat>(new funPtr_mat(&dfdu), false)
    );
}
