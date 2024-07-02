#include <Rcpp.h>

//////////////////// inverse logit function ////////////////////
double invlogit(double x){
  return 1/(1 + exp(-x));
};
