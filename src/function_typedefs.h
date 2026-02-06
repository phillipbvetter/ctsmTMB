#ifndef FUNCTION_TYPEDEFS_H
#define FUNCTION_TYPEDEFS_H

/* --------------------------------------------
   Typedefs for function pointer signatures
---------------------------------------------*/

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

#endif