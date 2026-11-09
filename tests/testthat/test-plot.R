
# Tests for the plot.ctsmTMB.fit S3 method. In particular that the smoothed
# state estimates produced by the Laplace methods in the fit object can be
# plotted via state.type="smoothed".

# ----------------------------------------
# Data and model
# ----------------------------------------

data(Ornstein)

# How many rows we want (maximum 201)
how.many.rows <- 10

df <- Ornstein[1:how.many.rows,]
model <- create.Ornstein1D.model()

# We see warnings from the Laplace methods, but this does not reflect a real issue.
my_expect_no_error <- function(object, ..., message = NULL, class = NULL){
  suppressWarnings(testthat::expect_no_error(object=object,...,message=message,class=class))
}

# ----------------------------------------
# Kalman filter fits: prior/posterior states
# ----------------------------------------

testthat::test_that("plot prior/posterior states for Kalman filter fit", {
  fit <- model$estimate(df, method="ekf", silent=TRUE, trace=0)
  p.prior <- plot(fit, type="states", state.type="prior", against.obs="y", print.plot=0)
  p.post  <- plot(fit, type="states", state.type="posterior", print.plot=0)
  testthat::expect_s3_class(p.prior[[1]], "ggplot")
  testthat::expect_s3_class(p.post[[1]], "ggplot")
  # smoothed states are not available for the Kalman methods
  testthat::expect_error(
    plot(fit, type="states", state.type="smoothed", print.plot=0),
    "only available for the Laplace"
  )
})

# ----------------------------------------
# Laplace fits: smoothed states
# ----------------------------------------

testthat::test_that("plot smoothed states for Laplace fits", {
  for (m in c("laplace", "laplace.thygesen")) {
    fit <- suppressWarnings(model$estimate(df, method=m, silent=TRUE, trace=0))

    # the fit object contains the smoothed estimates
    testthat::expect_true(all(c("mean","sd") %in% names(fit$states)))
    testthat::expect_true(!is.null(fit$states$mean$smoothed))

    # explicit state.type="smoothed"
    p <- plot(fit, type="states", state.type="smoothed", against.obs="y", print.plot=0)
    testthat::expect_s3_class(p[[1]], "ggplot")
    testthat::expect_match(p[[1]]$labels$title, "Smoothed")

    # the plotted line should be the smoothed state mean from the fit object
    line.data <- ggplot2::ggplot_build(p[[1]])$data[[2]]
    testthat::expect_equal(line.data$y, fit$states$mean$smoothed$x)

    # the default state.type should automatically switch to smoothed
    p2 <- suppressMessages(plot(fit, type="states", print.plot=0))
    testthat::expect_s3_class(p2[[1]], "ggplot")
    testthat::expect_match(p2[[1]]$labels$title, "Smoothed")
  }
})
