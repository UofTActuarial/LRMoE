% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Predict.R
\name{predict_class_prior}
\alias{predict_class_prior}
\title{predict_class_prior}
\usage{
predict_class_prior(X, alpha)
}
\arguments{
\item{X}{A matrix of covariates.}

\item{alpha}{A matrix of logit regression coefficients.}
}
\value{
prob
A matrix of latent class probabilities.

max_prob_idx
A matrix of the most likely latent class for each observation.
}
\description{
: Predicts the latent class probabilities, given covariates \code{X} and logit regression coefficients \code{alpha}.
}
