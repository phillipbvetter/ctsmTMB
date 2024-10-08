% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method_documentation.R
\name{addInput}
\alias{addInput}
\title{Specify input variables in the model object.}
\usage{
addInput(...)
}
\arguments{
\item{...}{a series of variable names (unquouted) that match variable names
in the defined system which should be treated as input variables.}
}
\description{
Declare whether a variable contained in system, observation or observation
variance equations is an input variable.
}
\examples{
# Example 1
addInput(u)

# Example 2
addInput(u1, u2, u3)

}
