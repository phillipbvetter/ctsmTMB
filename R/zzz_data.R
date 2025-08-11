## roxygen documentation for data objects

#' Sample from a simulated Ornstein-Uhlenbeck process with time-dependent mean
#' 
#' @description
#' The data was simulated using a standard Euler-Maruyama method. 
#' 
#' The simulated process is governed by the SDE
#' #' dx ~ theta * (mu + u - x) * dt + sigma_x * dw
#' 
#' The parameters used for simulation were
#' theta = 2, mu = 0.5, sigma_x = 1.358, sigma_y = 1e-8
#' 
#' The simulation time-step was 1e-3, and observation time-step 1e-1. 
#' The simulation was taken from t = 0..20
#' 
#' @format
#' A data frame of 201 rows and 3 columns. The columns represent the 
#' variables: `t` (time), `y` (observation) and `u` (input).
#' 
#' @usage Ornstein
#' 
#' @name Ornstein
#' @docType data
#' @keywords data
"Ornstein"
