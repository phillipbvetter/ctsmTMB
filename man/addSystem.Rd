% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method_documentation.R
\name{addSystem}
\alias{addSystem}
\title{Add state equations to model object}
\usage{
addSystem(form,
            ...)
}
\arguments{
\item{form}{formula specifying the stochastic differential equation to be
added to the system.}

\item{...}{formulas similar to \code{form}, used to allow specifying
multiple formulas by comma-seperation rather than providing a list.}
}
\description{
Add stochastic differential equation to the \code{ctsmTMB} model-object that
governs the differential evolution of states in the specified model.
}
\details{
\subsection{Usable functions}{

The formulas can contain most elementary functions such as \code{log},
\code{exp}, \code{logit} and \code{invlogit}. In general the supported
functions are only those that are both 1) defined in the derivative table
of the \code{Deriv} package, and 2) undestood by \emph{TMB} on the \emph{C++} side.
}
}
\examples{
# Example 1 - Linear System
addSystem(dx ~ theta * (mu - x + u) * dt + sigma * dw)

# Example 2 - Non-Linear System
addSystem(dx ~ theta * (mu - exp(x)^2 + u) * dt + sigma * x * (1-x) * dw)

}
