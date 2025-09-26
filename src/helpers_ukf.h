#ifndef HELPERS_UKF_H
#define HELPERS_UKF_H

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

MatrixXd sigma2chol(MatrixXd Xsigmapoints, double sqrt_c, int nn, int n){
  return ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt_c).block(0, 1, n, n);
}

MatrixXd create_sigmapoints_from_stateVec(VectorXd stateVec, MatrixXd sqrt_covMat, double sqrt_c, int n, int nn){
  MatrixXd Xsigmapoints(n, nn);
  VectorXd Ai;
  Xsigmapoints.col(0) = stateVec;
  for(int i=1; i < n+1; i++){
    Ai = sqrt_covMat.col(i-1);
    //Xsigmapoints.col(i) = stateVec + sqrt(1+n) * Ai;
    //Xsigmapoints.col(i+n) = stateVec - sqrt(1+n) * Ai;
    Xsigmapoints.col(i) = stateVec + sqrt_c * Ai;
    Xsigmapoints.col(i+n) = stateVec - sqrt_c * Ai;
  }
  return Xsigmapoints;
}

//////////// UKF function ///////////
MatrixXd Phi__(MatrixXd M){
  MatrixXd K(M.col(0).size(), M.row(0).size());
  K.setZero();
  K.template triangularView<Eigen::Lower>() = M.template triangularView<Eigen::Lower>();
  K.diagonal() = K.diagonal()/2.0;
  return K;
}

//////////// UKF sigma points drift function ///////////
template<typename T1>
MatrixXd ukf_f(T1 f__, MatrixXd Xsigmapoints, VectorXd parVec, VectorXd inputVec, int n, int nn){
  MatrixXd F(n,nn);
  VectorXd stateVec;
  for(int i=0; i < nn; i++){
    stateVec = Xsigmapoints.col(i);
    F.col(i) = f__(stateVec, parVec, inputVec);
  }
  return F;
}

//////////// UKF sigma points obs function ///////////
template<typename T1>
MatrixXd ukf_h(T1 h__, MatrixXd Xsigmapoints, VectorXd parVec, VectorXd inputVec, int m, int nn){
  MatrixXd H(m, nn);
  VectorXd stateVec;
  for(int i=0; i < nn; i++){
    stateVec = Xsigmapoints.col(i);
    H.col(i) = h__(stateVec, parVec, inputVec);
  }
  return H;
}

//////////// 1-step f moment ODE ///////////
template<typename T1, typename T2>
MatrixXd ukf_1step(T1 f__, T2 g__, MatrixXd Xsigmapoints, MatrixXd sqrt_covMat, MatrixXd W, VectorXd wm, VectorXd parVec, VectorXd inputVec, double sqrt_c, int n, int nn)
{
  MatrixXd F, G, Ainv, M, A_Phi_M, F_rhs, F_rhs0, F_rhs1(n, nn);
  F_rhs1.setZero();
  F = ukf_f(f__, Xsigmapoints, parVec, inputVec, n, nn);
  // G is currently not using sigma points because assumed state independent, but we hack it by just using the state vector
  G = g__(VectorXd(Xsigmapoints.col(0)), parVec, inputVec);
  Ainv = sqrt_covMat.inverse();
  M = Ainv * (Xsigmapoints * W * F.transpose() + F * W * Xsigmapoints.transpose() + G*G.transpose()) * Ainv.transpose();
  A_Phi_M = sqrt_covMat * Phi__(M);
  F_rhs1.block(0, 1, n, n) = A_Phi_M;
  F_rhs1.block(0, n+1, n, n) = -A_Phi_M;
  F_rhs0 = (F * wm).replicate(1, nn);
  // F_rhs = F_rhs0 + sqrt(3.0) * F_rhs1;
  F_rhs = F_rhs0 + sqrt_c * F_rhs1;
  return F_rhs;
}

//////////// ODE SOLVER ///////////
template<typename T1, typename T2>
MatrixXd ukf_ode_integration(
  T1 f__, 
  T2 g__, 
  MatrixXd Xsigmapoints, 
  MatrixXd sqrt_covMat, 
  MatrixXd W, 
  VectorXd wm, 
  VectorXd parVec, 
  VectorXd inputVec, 
  VectorXd dinputVec, 
  double dt,
  double sqrt_c, 
  int ode_solver, 
  int n, 
  int nn){

    /*Initial State and Cov Values*/
    MatrixXd X1sigmapoints;
    MatrixXd X0sigmapoints = Xsigmapoints;

    /*Forward Euler*/
    if(ode_solver == 1){

     X1sigmapoints = X0sigmapoints + ukf_1step(f__, g__, Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, sqrt_c, n, nn) * dt;

    }

    /*4th Order Runge-Kutta 4th*/
    if (ode_solver == 2){

     MatrixXd c1, c2, c3, c4;

     /*1. Approx Slope at Initial Point*/
     c1 = ukf_1step(f__, g__,Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, sqrt_c, n, nn);

     /*2. First Approx Slope at Midpoint*/
     inputVec += 0.5 * dinputVec;
     Xsigmapoints = X0sigmapoints + 0.5 * dt * c1;
     // sqrt_covMat = ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt(3.0)).block(0, 1, n, n);
     sqrt_covMat = sigma2chol(Xsigmapoints, sqrt_c, nn, n);
     c2 = ukf_1step(f__, g__,Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, sqrt_c, n, nn);        

     /*3. Second Approx Slope at Midpoint*/
     Xsigmapoints = X0sigmapoints + 0.5 * dt * c2;
     //sqrt_covMat = ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt(3.0)).block(0, 1, n, n);
     sqrt_covMat = sigma2chol(Xsigmapoints, sqrt_c, nn, n);
     c3 = ukf_1step(f__, g__,Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, sqrt_c, n, nn);

     /*4. Approx Slope at End Point*/
     inputVec += 0.5 * dinputVec;
     Xsigmapoints = X0sigmapoints + dt * c3;
     //sqrt_covMat = ((Xsigmapoints - Xsigmapoints.col(0).replicate(1, nn))/sqrt(3.0)).block(0, 1, n, n);
     sqrt_covMat = sigma2chol(Xsigmapoints, sqrt_c, nn, n);
     c4 = ukf_1step(f__, g__,Xsigmapoints, sqrt_covMat, W, wm, parVec, inputVec, sqrt_c, n, nn);

     /*ODE UPDATE*/
     X1sigmapoints = X0sigmapoints + dt * (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0; 

   }

   return X1sigmapoints;
  }

#endif
