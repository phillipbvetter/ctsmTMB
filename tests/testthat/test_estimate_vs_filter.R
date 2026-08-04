#### Setup ####
data(Ornstein2D)
how.many.rows <- 10
df <- Ornstein2D[1:how.many.rows,]
data(EstimateReferenceData)
ed <- EstimateReferenceData


model <- create.Ornstein2D.model()
methods = c("ekf", "lkf", "ukf")

testthat::test_that("Check if estimation filter matches filter (ZOH)",{

  for(m in methods){
    fit <- model$estimate(df, method=m, silent=T, trace=0, first.order.input.hold = FALSE)
    rfit <- fit[c("states","residuals","observations")]
    filt <- model$filter(df, method=m, silent=T, pars=fit$par.fixed, first.order.input.hold = FALSE)
    testthat::expect_equal(rfit, filt)
  }
})

testthat::test_that("Check if estimation filter matches filter (FOH",{

  for(m in methods){
    fit <- model$estimate(df, method=m, silent=T, trace=0, first.order.input.hold = TRUE)
    rfit <- fit[c("states","residuals","observations")]
    filt <- model$filter(df, method=m, silent=T, pars=fit$par.fixed, first.order.input.hold = TRUE)
    testthat::expect_equal(rfit, filt)
  }

})

