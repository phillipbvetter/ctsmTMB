#ifndef ZIGGURAT_SEEDER_H
#define ZIGGURAT_SEEDER_H

#include <zigg/header>

// Set seed for RNG draws using ziggurat
inline void set_simulation_seed(Rcpp::Nullable<int> seed, zigg::Ziggurat& ziggurat) {
  if (!seed.isNull()) {
    uint32_t su = static_cast<uint32_t>(Rcpp::as<int>(seed));
    ziggurat.setSeed(su);
  }
}

#endif