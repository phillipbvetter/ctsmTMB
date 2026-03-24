data(Ornstein)
df <- Ornstein[1:10,]
df$z = df$y

model = ctsmTMB$new()
model$addSystem(
  dx ~ theta * (mu*u - x) * dt + sigma_x*dw,
  dx2 ~ sigma_x * dw1
)
model$addObs(
  y ~ x,
  z ~ x^2
)
model$setVariance(
  y ~ sigma_y^2,
  z ~ sigma_y^2 * x^2
)
model$setAlgebraics(
  theta   ~ exp(logtheta),
  sigma_x ~ exp(logsigma_x),
  sigma_y ~ exp(logsigma_y)
)
model$addInput(u)
model$setParameter(
  logtheta   = log(c(initial = 1, lower = 1e-5, upper = 50)),
  mu         = c(initial = 1.5, lower = 0, upper = 5),
  logsigma_x = log(c(initial = 1, lower = 1e-10, upper = 2)),
  logsigma_y = log(c(initial = 1, lower = 1e-10, upper = 2)),
  logsigma_y = log(0.1)
)
model$setInitialState(list(rep(1, 2), 0.656 * diag(2)))
model$estimate(df, method = "ekf", silent = TRUE, trace = 0)

########################################################################
# set_parameters: private$pars is always named with correct length
########################################################################
testthat::test_that("set_parameters produces named private$pars with correct length", {

  prv = model$.__enclos_env__$private
  n   = length(prv$parameters)

  check_pars = function(pars_arg) {
    model$filter(df, method = "ekf", pars = pars_arg, use.cpp = FALSE, silent = TRUE)
    p = prv$pars
    testthat::expect_equal(length(p), n)
    testthat::expect_false(is.null(names(p)))
    testthat::expect_false(any(is.na(names(p))))
    testthat::expect_equal(names(p), prv$parameter.names)
  }

  # NULL path: best-available values (estimated > initial)
  check_pars(NULL)

  # Named vector path: partial override by name
  check_pars(c(mu = 1.0))

  # Unnamed full-length path
  check_pars(unname(model$getParameters(value = "initial")))

  # Unnamed free-only path
  check_pars(unname(model$getParameters(type = "free", value = "initial")))
})

########################################################################
# setParameter round-trip via getParameters() data.frame
########################################################################
testthat::test_that("getParameters() data.frame can be passed back to setParameter()", {

  before = model$getParameters()

  # Modify one initial value in the returned data.frame and feed it back
  modified = before
  modified["mu", "initial"] = 2.5
  testthat::expect_no_error(model$setParameter(modified))

  after = model$getParameters()

  # The changed parameter reflects the new initial value
  testthat::expect_equal(after["mu", "initial"], 2.5)

  # All other parameters are unchanged
  other = setdiff(rownames(before), "mu")
  testthat::expect_equal(after[other, "initial"], before[other, "initial"])
  testthat::expect_equal(after[other, "lower"],   before[other, "lower"])
  testthat::expect_equal(after[other, "upper"],   before[other, "upper"])
  testthat::expect_equal(after[other, "type"],    before[other, "type"])

  # Names and length are intact
  testthat::expect_equal(rownames(after), rownames(before))
  testthat::expect_equal(nrow(after), nrow(before))
})

########################################################################
# Partial named pars override: only the named entries change
########################################################################
testthat::test_that("partial named pars only modifies the specified entries of private$pars", {

  prv = model$.__enclos_env__$private

  run_and_capture = function(method, call_fn) {
    call_fn(pars = NULL)
    baseline = prv$pars

    override = c(mu = 3.0, logtheta = -1.0)
    call_fn(pars = override)
    result = prv$pars

    # Overridden entries match the supplied values
    testthat::expect_equal(result[names(override)], override,
                           info = paste("method:", method))

    # All other entries are unchanged from the NULL baseline
    others = setdiff(names(baseline), names(override))
    testthat::expect_equal(result[others], baseline[others],
                           info = paste("method:", method))
  }

  run_and_capture("filter",
    function(pars) model$filter(df, method = "ekf", pars = pars,
                                use.cpp = FALSE, silent = TRUE))

  run_and_capture("predict",
    function(pars) model$predict(df, method = "ekf", pars = pars,
                                 use.cpp = FALSE, silent = TRUE))

  run_and_capture("simulate",
    function(pars) model$simulate(df, method = "ekf", pars = pars,
                                  use.cpp = FALSE, silent = TRUE))
})
