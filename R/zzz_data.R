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
