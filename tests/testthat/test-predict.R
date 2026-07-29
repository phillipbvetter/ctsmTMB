# In this script we run the estimate function on four different models with
# number of states and observations as follows:
# States / Observations:
# 1 / 1
# 2 / 2
# 2 / 3
# 3 / 2
# It is of particular important to run models with different number of states
# to make sure matrix dimensions are correctly coded etc.

# ----------------------------------------
# Data
# ----------------------------------------

data(Ornstein)
data(Ornstein2D)

# How many rows we want (maximum 201)
how.many.rows <- 10

# We need a data-set for each of the model
df11 <- Ornstein[1:how.many.rows,] # 1 / 1
df22 <- Ornstein2D[1:how.many.rows,] # 2 / 2
df23 <- df22; df23$y3 <- df23$y1 + df23$y2 # 2 / 3
df32 <- df22 # 3 / 2


# ----------------------------------------
# Models
# ----------------------------------------

model11 <- create.Ornstein1D.model() # 1 / 1
model22 <- create.Ornstein2D.model() # 2 / 2
{ # 2 / 3
  model23 <- create.Ornstein2D.model()
  model23$addObs(y3 ~ x1 + x2)
  model23$setVariance(y3 ~ sigma_y1^2 + sigma_y2^2)
}
{ # 3 /2
  model32 <- create.Ornstein2D.model()
  model32$addSystem(dx3 ~ (theta+alpha)*(x2-x3) * dt + sigma_x1 * dw1 + sigma_x2 * dw2 + dw3)
  model32$setInitialState(list(rep(1,3), 1e-1*diag(3)))
}
# ----------------------------------------
# Estimation
# ----------------------------------------

# Modify testthat::expect_no_error with warning supression. We see warnings
# from the Laplace methods, but this does not reflect a real issue.
my_expect_no_error <- function(object, ..., message = NULL, class = NULL){
  suppressWarnings(testthat::expect_no_error(object=object,...,message=message,class=class))
}

caller <- "predict"
methods <- c("ekf", "lkf", "ukf")

testthat::test_that("Checking for 'estimate' errors in 1/1 model",{
  # 1 / 1
  for (m in methods) {
    my_expect_no_error(model11[[caller]](df11, method=m, silent = TRUE, trace=0))
  }
})

testthat::test_that("Checking for 'estimate' errors in 2/2 model",{
  # 2 / 2
  for (m in methods) {
    my_expect_no_error(model22[[caller]](df22, method=m, silent = TRUE, trace=0))
  }
})

testthat::test_that("Checking for 'estimate' errors in 2/3 model",{
  # 2 / 3
  for (m in methods) {
    my_expect_no_error(model23[[caller]](df23, method=m, silent = TRUE, trace=0))
  }
})

testthat::test_that("Checking for 'estimate' errors in 3/2 model",{
  # 3 / 2
  for (m in methods) {
    my_expect_no_error(model32[[caller]](df32, method=m, silent = TRUE, trace=0))
  }
})

