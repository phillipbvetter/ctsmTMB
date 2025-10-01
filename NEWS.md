# ctsmTMB 1.0.2 (????-??-??)

* Previously only 'ekf' had support for prediction and simulation. Nowe 'lkf' and 'ukf' are also supported in both R and with Rcpp.
* New method 'setAdvancedSettings' for various extra features (force.ad and tape configurations TMB::config and RTMB::TapeConfig)
* New method 'setTrainingMethod' which enables training the likelihood as an ODE minimizing distance to full prediction - still need support for k-step.

# ctsmTMB 1.0.1 (2025-08-27)

* Internal large code-clean-up
* Rcpp RNG simulations now use 'zigg' instead of 'RcppZiggurator' thereby avoiding need for 'SystemRequirements: GNU GSL'.
* Removed RTMBode possibilities.
* Renamed smoothing method 'laplace2' to 'laplace.thygesen'.
* Fixed a simple typo bug in the results reporting functions.
* Added an the implicit euler ODE solver to 'ode.solver' options.


# ctsmTMB 1.0.0 (2025-04-08)

* Initial release
