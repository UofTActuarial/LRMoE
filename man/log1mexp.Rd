% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{log1mexp}
\alias{log1mexp}
\title{Compute f(a) = log(1 +/- exp(-a)) Numerically Optimally}
\usage{
log1mexp(a, cutoff = log(2))
}
\arguments{
\item{a}{\cr
numeric vector of positive values}

\item{cutoff}{\cr
log(2) is “optimal”, but the exact value is unimportant, and anything in (0.5, 1) is fine.}
}
\description{
Compute f(a) = log(1 +/- exp(-a)) Numerically Optimally
This helper function was a directly copy from \code{copula} package. It was initialized here intend to maintain the minimum package complexity
}
