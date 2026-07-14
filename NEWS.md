# ctsmTMB 1.1.1 (2026-12-07)

* The methods supports zero and first order input hold for all RTMB methods, via the function argument 'first.order.input.hold'.
* The plot.ctsmTMB.fit method now shows standardized residuals in the residual plot by default, which may be changed via the 'standardized.residuals' argument.
* Added some author-only tests, including data sets and functions for quickly defining 1D and 2D Ornstein Uhlenbeck models.
* Added better tests for all methods, in particular to make sure we also catch bugs in models where n.states != n.obs.
* Internal code clean-up and some refactoring, specifically bunching together several private fields i.e. private$model now contains all model equations etc.
* Readme edits.


# ctsmTMB 1.1.0 (2026-01-02)

* System C++ functions are now compiled on the first call to a method, and are only recompiled if the model equation changes.
* The C++ functionality for 'filter', 'predict' and 'simulate' is now greatly improved:
  * The compilation time is reduced from 5-10 seconds to 2-3 seconds.
  * The computations are faster.
  * The report functionality is now carried out in C++ (rather than pure R code), and is faster.
  * The three methods now also return forecasted observations, and associated dispersion.
* The reporting for the 'estimate' method is now faster, due to reporting improvements (see above)
* Previously the package returned data.frames as outputs in most cases. This has now been changed to matrices since these 
are usually faster to work with internally. Users must therefore use e.g. [,"t"] instead of $t when accessing columns.


# ctsmTMB 1.0.2 (2025-11-01)

* Previously only 'ekf' had support for prediction and simulation. Now 'lkf' and 'ukf' are also supported in both R and with Rcpp.
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
