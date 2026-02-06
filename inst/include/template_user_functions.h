#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace RcppEigen;
// [[Rcpp::depends(RcppEigen)]]

/* Constants */
constexpr double pi = 3.14159265358979323846;

/* Custom Functions*/
inline double invlogit(double x){
   return 1.0/(1.0+exp(-x));
}

typedef Eigen::VectorXd (*funPtr_vec_const)(
   const Eigen::VectorXd&, 
   const Eigen::VectorXd&, 
   const Eigen::VectorXd&);
typedef Eigen::MatrixXd (*funPtr_mat_const)(
   const Eigen::VectorXd&, 
   const Eigen::VectorXd&, 
   const Eigen::VectorXd&);
typedef Eigen::ArrayXd (*funPtr_array_const)(
   const Eigen::VectorXd&, 
   const Eigen::VectorXd&, 
   const Eigen::VectorXd&);

/* --------------------------------------------
   User-generated functions
---------------------------------------------*/

// INSERT F_CONST

// INSERT DFDX_CONST

// INSERT H_CONST

// INSERT G_CONST

// INSERT HVAR_CONST

// INSERT DHDX_CONST

// INSERT HVAR_ARRAY_CONST

// INSERT DFDU_CONST

/* --------------------------------------------
   Dispatcher returning XPtrs
---------------------------------------------*/
// [[Rcpp::export]]
SEXP get_sysfun_cpp_function_ptrs()
{
  return List::create(
    Named("f_const") = XPtr<funPtr_vec_const>(new funPtr_vec_const(&f_const), false),
    Named("h_const") = XPtr<funPtr_vec_const>(new funPtr_vec_const(&h_const), false),
    Named("dfdx_const") = XPtr<funPtr_mat_const>(new funPtr_mat_const(&dfdx_const), false),
    Named("dhdx_const") = XPtr<funPtr_mat_const>(new funPtr_mat_const(&dhdx_const), false),
    Named("g_const") = XPtr<funPtr_mat_const>(new funPtr_mat_const(&g_const), false),
    Named("hvar_const") = XPtr<funPtr_mat_const>(new funPtr_mat_const(&hvar_const), false),
    Named("hvar_array_const") = XPtr<funPtr_array_const>(new funPtr_array_const(&hvar_array_const), false),
    Named("dfdu_const") = XPtr<funPtr_mat_const>(new funPtr_mat_const(&dfdu_const), false)
    );
}
