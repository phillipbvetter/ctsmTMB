# Load data
# ------------------
data(Ornstein2D)
how.many.rows <- 50
df <- Ornstein2D[1:how.many.rows,]
# add missing data to check those kind of errors
df[1,c("y1","y2")] <- NA
df[2,"y1"] <- NA
df[3,"y2"] <- NA

model <- create.Ornstein2D.model()
methods = c("ekf", "lkf", "ukf")
fields <- c("states","observations","residuals")

testthat::test_that("Check if estimation filter matches filter (ZOH)",{

  for(m in methods){
    fit <- model$estimate(df, method=m, silent=T, trace=0, first.order.input.hold = FALSE)
    rfit <- fit[fields]
    filt <- model$filter(df, method=m, silent=T, pars=fit$par.fixed, first.order.input.hold = FALSE)
    testthat::expect_equal(rfit, filt)
  }
})

testthat::test_that("Check if estimation filter matches filter (FOH",{

  for(m in methods){
    fit <- model$estimate(df, method=m, silent=T, trace=0, first.order.input.hold = TRUE)
    rfit <- fit[fields]
    filt <- model$filter(df, method=m, silent=T, pars=fit$par.fixed, first.order.input.hold = TRUE)
    testthat::expect_equal(rfit, filt)
  }

})

