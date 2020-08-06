#' Predict the latent class probabilities, given a fixed covariate matrix X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param prob If \code{TRUE}, the function returns a matrix of latent class probabilities. If \code{FALSE},
#'        it returns the index of the most likely latent class.
#'
#'
#' @return A matrix of latent class probabilities by observation and by component.
#'
#' @importFrom matrixStats rowLogSumExps
#'
#' @export PredictClassPrior
PredictClassPrior = function(X, alpha, prob = TRUE)
{
  weighting = GateLogit(X, alpha)
  if(prob==TRUE){
    return(exp(weighting))
  }else{
    return(ColMaxIdx(weighting))
  }
}
