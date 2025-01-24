
# Hacks to solve issues with R CMD check



# 1. Prevent 'no visile binding for global variables' for functions defined
# by eval(parse(text=function_as_string) type constructions used in RTMB

f__ <- function(stateVec, parVec, inputVec) NULL
dfdx__ <- function(stateVec, parVec, inputVec) NULL
g__ <- function(stateVec, parVec, inputVec) NULL
ekf.nll <- function(p) NULL
laplace.nll <- function(p) NULL
x <- NULL



# The vignettes generate cpp and dynamic libraries, but this gives
# errors on github windows. We remove the files with the following 
# function 
