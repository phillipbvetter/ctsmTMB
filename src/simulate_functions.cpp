// #include <Rcpp.h>
#include <RcppArmadillo.h> // Must be included before RcppEigen and must not include Rcpp.h
#include <RcppEigen.h>
#include <zigg/header>

#include "function_typedefs.h"
#include "misc_helpers.h"
#include "ziggurat_seeder.h"

// This instance of ziggurat is only defined here - it is different from the instance used for state simulations
// This allows for two different seeds, such that two random sources can be differentiated.
static zigg::Ziggurat ziggurat_obs;

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
List calculate_simulation_observations(
  const List simulation_raw_states,
  const List funPtrs,
  const Eigen::Map<Eigen::MatrixXd>  & inputMat_T,
  const Eigen::Map<Eigen::VectorXd>  & parVec,
  int n_states, 
  int n_obs,
  int k_ahead, 
  int n_sims,
  Nullable<int> seed)
{

  // Set simulating seed if seed is not NULL
  set_simulation_seed(seed, ziggurat_obs);
  
  auto h = get_funptr<funPtr_vec_const>(funPtrs, "h_const");
  auto hvar_array = get_funptr<funPtr_array_const>(funPtrs, "hvar_array_const");

  // transpose to access columns instead of rows
  Eigen::ArrayXd  epsilon(n_obs); // noise vector

  int n_forecasts = simulation_raw_states.length();
  int k_aheads = k_ahead + 1;

  // Create matrix to store "simulated" observations
  Eigen::MatrixXd simulated_obs_matrix(n_obs, n_sims);

  List simulation_raw_obs(n_forecasts);

  for(int i=0; i < n_forecasts; i++){

    // Allocate
    List inner_output_list(k_aheads);

    // Extract the inner list (contains k_aheads matrices, one per prediction time step)
    const List inner_list = simulation_raw_states[i];

    for(int j=0; j < k_aheads; j++ ){

      // Extract the simulation matrix (each row is a simulation, each column is a state)
      Eigen::MatrixXd inner_state_matrix = inner_list[j];
      Eigen::MatrixXd inner_state_matrix_T = inner_state_matrix.transpose();

      // Read the relevant inputs from data
      const Eigen::Ref<const Eigen::VectorXd> inputVec = inputMat_T.col(i+j);

      for(int k=0; k < n_sims; k++){

        // Extract simulated state vector
        const Eigen::Ref<const Eigen::VectorXd> stateVec = inner_state_matrix_T.col(k);

        // Generate standard normal noise vector
        for(int ii=0; ii < n_obs; ii++){
          epsilon(ii) = ziggurat_obs.rnorm();
        }

        // Calculcate observation vector and assign to matrix
        simulated_obs_matrix.col(k) = h(stateVec, parVec, inputVec) + (hvar_array(stateVec, parVec, inputVec).sqrt() * epsilon).matrix(); 

      }

      // Store - tranpose back to original format
      inner_output_list(j) = simulated_obs_matrix.transpose();
    }
    // Store
    simulation_raw_obs(i) = inner_output_list;
  }

  return simulation_raw_obs;
}

// [[Rcpp::export]]
List build_simulation_returnlist(List simulation_raw, 
                                 NumericVector time_vector,
                                 int n_states, 
                                 int k_ahead, 
                                 int n_sims)
{
  int n_forecasts = simulation_raw.length();
  int k_ahead_plus_1 = k_ahead + 1;

  List output_list(n_states);

  for(int i = 0; i < n_states; i++){

    List horizon_list(n_forecasts);

    for(int k = 0; k < n_forecasts; k++){

      List inner_list = simulation_raw[k];

      // Store sims x horizon (column-major, efficient)
      arma::mat outmat(n_sims, k_ahead_plus_1, arma::fill::none);

      for(int j = 0; j < k_ahead_plus_1; j++){

        // Convert once per matrix
        arma::mat inner = as<arma::mat>(inner_list[j]);

        // Fast column copy (no loops, vectorized, cache-friendly)
        outmat.col(j) = inner.col(i);
      }

      // Transpose once at the end to match your expected shape
      horizon_list[k] = wrap(outmat.t());
    }

    output_list[i] = horizon_list;
  }

  return output_list;
}

// [[Rcpp::export]]
List build_simulation_timelists(NumericVector time_vector,
                                int n_forecasts,
                                int k_ahead)
{
  int k_ahead_plus_1 = k_ahead + 1;
  int ii;

  List times_list(n_forecasts);

  CharacterVector mat_names = CharacterVector::create("i","j","t.i","t.j","k.ahead");

  // Pre-build the k-ahead sequence once
  arma::vec ka = arma::regspace(0, k_ahead);

  for(int k = 0; k < n_forecasts; k++){

    arma::mat premat(k_ahead_plus_1, 5, arma::fill::none);

    premat.col(0).fill(k);                  // i
    premat.col(1) = k + ka;                 // j
    premat.col(2).fill(time_vector(k));     // t.i
    for(ii=0;ii<k_ahead_plus_1;ii++){
     premat(ii, 3) = time_vector(k+ii);     // t.j
    }
    premat.col(4) = ka;                     // k.ahead

    NumericMatrix out = wrap(premat);  // Must wrap arma::mat can't be inside List
    colnames(out) = mat_names;
    times_list[k] = out;
  }

  return times_list;
}
