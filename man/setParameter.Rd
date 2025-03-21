% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method_documentation.R
\name{setParameter}
\alias{setParameter}
\title{Specify parameters in the model object}
\usage{
setParameter(...)
}
\arguments{
\item{...}{a comma-seperated series of vectors/matrix entries}
}
\description{
Declare which variables that are (fixed effects) parameters in
the specified model, and specify the initial optimizer values, as well as
lower / upper bounds. Parameters can be declared either as vectors or as
matrices. The first entry is the initial value, the second entry is the lower
bound and the third entry is the upper bound. Providing only a first entry
fixes the particular parameter to that value.
}
\section{Example 2}{
parmat = matrix(rep(c(0,-10,10),times=3),ncol=3,byrow=T)
rownames(parmat) = c("a","b","c")
colnames(parmat) = c("initial","lower","upper")
setParameter(parmat)
}

