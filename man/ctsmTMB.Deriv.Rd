% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deriv_funs.R
\name{ctsmTMB.Deriv}
\alias{ctsmTMB.Deriv}
\title{Modify Deriv::Deriv}
\usage{
ctsmTMB.Deriv(
  f,
  x = if (is.function(f)) NULL else all.vars(if (is.character(f)) parse(text = f) else f),
  env = get.Deriv.env(),
  use.D = FALSE,
  cache.exp = FALSE,
  nderiv = NULL,
  combine = "c",
  drule = get.Deriv.drule()
)
}
\description{
We create our own Deriv function based on Deriv::Deriv with custom environment and
drules to support new functions e.g erf
}
