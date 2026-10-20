// #include <Rcpp.h>
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
      MatrixXd outmat(k_ahead_plus_1, n_sims);
      for(int j = 0; j < k_ahead_plus_1; j++){
        Map<MatrixXd> inner = as<Map<MatrixXd>>(inner_list[j]);
        outmat.row(j) = inner.col(i).transpose();
      }
      horizon_list[k] = wrap(outmat);
    }
    output_list[i] = horizon_list;
  }
  return output_list;
}
