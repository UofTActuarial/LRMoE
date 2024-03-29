% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation.R
\name{sim_logit_gating}
\alias{sim_logit_gating}
\title{Simulation based on Logit Gating Function}
\usage{
sim_logit_gating(probs)
}
\arguments{
\item{probs}{The probability matrix with N x P dimentions.}
}
\value{
result The result matrix with N X P results
}
\description{
Given alpha and x, generate a probability matrix. Give the simulation matrix with each row
of the probability matrix to be multinomial probability.
}
\keyword{internal}
