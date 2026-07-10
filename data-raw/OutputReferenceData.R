## Generate reference outputs for C++ filter / predict / simulate tests using the Ornstein2D data.

library(ctsmTMB)
data(Ornstein2D)
how.many.rows <- 50
df <- Ornstein2D[1:how.many.rows,]

# -----------------------------------------------------------------------
# 2D Ornstein model
# -----------------------------------------------------------------------

model <- create.Ornstein2D.model()
kalman.methods <- c("ekf","lkf","ukf")
OutputReferenceData <- list()

# ZERO ORDER HOLD
for(m in kalman.methods) OutputReferenceData$filter[[m]] <- model$filter(df, method  = m, silent  = TRUE)
for(m in kalman.methods) OutputReferenceData$predict[[m]] <- model$predict(df, method  = m, silent  = TRUE)
for(m in kalman.methods) OutputReferenceData$simulate[[m]] <- model$simulate(df, method  = m, silent  = TRUE, cpp.seeds = c(123, 456))

# FIRST ORDER HOLD
for(m in kalman.methods) OutputReferenceData$filter[[paste0(m,"_foh")]] <- model$filter(df, method  = m, silent  = TRUE, first.order.input.hold = TRUE )
for(m in kalman.methods) OutputReferenceData$predict[[paste0(m,"_foh")]] <- model$predict(df, method  = m, silent  = TRUE, first.order.input.hold = TRUE)
for(m in kalman.methods) OutputReferenceData$simulate[[paste0(m,"_foh")]] <- model$simulate(df, method  = m, silent  = TRUE, first.order.input.hold = TRUE, cpp.seeds = c(123, 456))

usethis::proj_set("~/github/ctsmTMB")
usethis::use_data(OutputReferenceData, overwrite = TRUE)
