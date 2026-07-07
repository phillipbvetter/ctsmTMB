#include <Rcpp.h>

// [[Rcpp::export]]
bool check_if_funptr_is_null(const Rcpp::List& x, const char* name){
  if (R_ExternalPtrAddr(x[name]) == nullptr){
    return true;
  }
  return false;
}