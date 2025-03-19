
# Hacks to solve issues with R CMD check

# 1. Prevent 'no visile binding for global variables' for functions defined
# by eval(parse(text=function_as_string) type constructions used in RTMB

f__ <- function(stateVec, parVec, inputVec) NULL
dfdx__ <- function(stateVec, parVec, inputVec) NULL
g__ <- function(stateVec, parVec, inputVec) NULL
dhdx__ <- function(stateVec, parVec, inputVec) NULL
dfdu__ <- function(stateVec, parVec, inputVec) NULL
h__ <- function(stateVec, parVec, inputVec) NULL
hvar__ <- function(stateVec, parVec, inputVec) NULL
hvar__matrix <- function(stateVec, parVec, inputVec) NULL
x <- NULL
