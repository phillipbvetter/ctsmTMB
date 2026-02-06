#ifndef MISC_HELPERS_H
#define MISC_HELPERS_H

#include <Rcpp.h>

// This function just makes for less code when function pointers are extracted
template <typename T>
T get_funptr(const Rcpp::List& x, const char* name) {
  return *Rcpp::as<Rcpp::XPtr<T>>(x[name]);
}

#endif