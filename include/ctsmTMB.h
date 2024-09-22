#include <Rcpp.h>

const double pi = M_PI;

//////////////////// inverse logit function ////////////////////
double invlogit(double x){
  return 1/(1 + exp(-x));
};
