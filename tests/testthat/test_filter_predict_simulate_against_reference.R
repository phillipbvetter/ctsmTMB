# Tests if the outputs from filter, predict and simulate against the Ornstein2D data, matches
# the results previously obtained and stored in the OutputReferenceData.

library(ctsmTMB)
data(Ornstein2D)
how.many.rows <- 50
df <- Ornstein2D[1:how.many.rows,]
data(OutputReferenceData)
od <- OutputReferenceData

model <- create.Ornstein2D.model()
kalman.methods <- c("ekf","lkf","ukf")

testthat::test_that("Testing filter, predict and simulate",{
  # Skip if test is not the author
  testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - skipping due to minor numerical inconsistencies across operating systems and versions")

  # ZERO ORDER HOLD
  for(m in kalman.methods) testthat::expect_equal(model$filter(df, method = m, silent=TRUE), od$filter[[m]])
  for(m in kalman.methods) testthat::expect_equal(model$predict(df, method = m, silent=TRUE), od$predict[[m]])
  for(m in kalman.methods) testthat::expect_equal(model$simulate(df, method = m, silent=TRUE, cpp.seeds = c(123, 456)), od$simulate[[m]])

  # FIRST ORDER HOLD
  for(m in kalman.methods) testthat::expect_equal(model$filter(df, method = m, silent=TRUE, first.order.input.hold = TRUE), od$filter[[paste0(m,"_foh")]])
  for(m in kalman.methods) testthat::expect_equal(model$predict(df, method = m, silent=TRUE, first.order.input.hold = TRUE), od$predict[[paste0(m,"_foh")]])
  for(m in kalman.methods) testthat::expect_equal(model$simulate(df, method = m, silent=TRUE, first.order.input.hold = TRUE, cpp.seeds = c(123, 456)), od$simulate[[paste0(m,"_foh")]])
})
