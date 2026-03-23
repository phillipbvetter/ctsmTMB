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
#' theta = 5, mu = 3, sigma_x = 1, sigma_y = 0.1
#' 
#' The simulation time-step was 1e-3, and observation time-step 1e-1. 
#' The simulation was taken from t = 0..20
#' 
#' The simulated input was given by \code{u.sim = cumsum(rnorm(length(t.sim),sd=0.05))}
#' where \code{t.sim} is the simulated time vector.
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

#' Sample from a simulated Ornstein-Uhlenbeck process with time-dependent mean
#' 
#' @description
#' The data was simulated using a standard Euler-Maruyama method. 
#' 
#' The simulated process is governed by the SDE
#' #' dx ~ theta * (mu + u - x) * dt + sigma_x * dw
#' 
#' The parameters used for simulation were
#' theta = 10, mu = 1, sigma_x = 1, sigma_y = 0.1
#' 
#' The simulation time-step was 1e-3, and observation time-step 1e-2. 
#' The simulation was taken from t = 0..5 seconds
#' 
#' @format
#' A data frame of 501 rows and 3 columns. The columns represent the 
#' variables: `t` (time), `y` (observation) and `u` (input).
#' 
#' @usage Ornstein500
#' 
#' @name Ornstein500
#' @docType data
#' @keywords data
"Ornstein500"

#' Sample from a simulated two-state Ornstein-Uhlenbeck process
#'
#' @description
#' The data was simulated using a standard Euler-Maruyama method.
#'
#' The simulated process is governed by the coupled SDEs
#'
#' dx1 ~ theta * (mu + u - x1) * dt + sigma_x1 * dw1
#' dx2 ~ alpha * (x1 - x2)    * dt + sigma_x2 * dw2
#'
#' with observations
#'
#' y1 ~ x1
#' y2 ~ x2
#'
#' x2 mean-reverts toward x1, acting as a lagged/smoothed version of x1.
#' The dataset is intended for testing multi-dimensional (> 1 state) inference.
#'
#' The parameters used for simulation were
#' theta = 5, mu = 3, alpha = 2, sigma_x1 = 1, sigma_x2 = 0.5,
#' sigma_y1 = 0.1, sigma_y2 = 0.1
#'
#' The simulation time-step was 1e-3, and observation time-step 1e-1.
#' The simulation was taken from t = 0..20
#'
#' The simulated input was given by \code{u.sim = cumsum(rnorm(length(t.sim), sd=0.05))}
#' where \code{t.sim} is the simulated time vector.
#'
#' @format
#' A data frame of 201 rows and 4 columns. The columns represent the
#' variables: `t` (time), `y1` (observation of x1), `y2` (observation of x2),
#' and `u` (input).
#'
#' @usage Ornstein_augmented
#'
#' @name Ornstein_augmented
#' @docType data
#' @keywords data
"Ornstein_augmented"
