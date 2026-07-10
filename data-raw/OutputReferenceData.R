## Generate reference outputs for C++ filter / predict / simulate tests using the Ornstein2D data.

library(ctsmTMB)
data(Ornstein2D)
how.many.rows <- 50
df <- Ornstein2D[1:how.many.rows,]

# -----------------------------------------------------------------------
# 2D Ornstein model
# -----------------------------------------------------------------------

model <- create.Ornstein2D.model()

# -----------------------------------------------------------------------
# Filter Reference outputs
# -----------------------------------------------------------------------

filter.ekf   = model$filter(df,
                            method  = "ekf",
                            use.cpp = TRUE,
                            silent  = TRUE)
filter.lkf   = model$filter(df,
                            method  = "lkf",
                            use.cpp = TRUE,
                            silent  = TRUE)
filter.ukf   = model$filter(df,
                            method  = "ukf",
                            use.cpp = TRUE,
                            silent  = TRUE)

# -----------------------------------------------------------------------
# Predict Reference outputs
# -----------------------------------------------------------------------

predict.ekf   = model$predict(df,
                              method  = "ekf",
                              use.cpp = TRUE,
                              silent  = TRUE)
predict.lkf  = model$predict(df,
                             method  = "lkf",
                             use.cpp = TRUE,
                             silent  = TRUE)
predict.ukf   = model$predict(df,
                              method  = "ukf",
                              use.cpp = TRUE,
                              silent  = TRUE)

# -----------------------------------------------------------------------
# Simulate Reference outputs
# -----------------------------------------------------------------------

simulate.ekf   = model$simulate(df,
                                method  = "ekf",
                                use.cpp = TRUE,
                                cpp.seeds = c(123, 456),
                                silent  = TRUE)
simulate.lkf   = model$simulate(df,
                                method  = "lkf",
                                use.cpp = TRUE,
                                cpp.seeds = c(123, 456),
                                silent  = TRUE)
simulate.ukf   = model$simulate(df,
                                method  = "ukf",
                                use.cpp = TRUE,
                                cpp.seeds = c(123, 456),
                                silent  = TRUE)

# -----------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------

OutputReferenceData = list(
  filters = list(
    ekf = filter.ekf,
    lkf = filter.lkf,
    ukf = filter.ukf
  ),
  predicts = list(
    ekf = predict.ekf,
    lkf = predict.lkf,
    ukf = predict.ukf
  ),
  simulate = list(
    ekf = simulate.ekf,
    lkf = simulate.lkf,
    ukf = simulate.ukf
  )
)

usethis::proj_set("~/github/ctsmTMB")
usethis::use_data(OutputReferenceData, overwrite = TRUE)
