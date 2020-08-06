#' Predict the posterior latent class probabilities, given a fixed covariate matrix X, Y and a model.
#'
#' @param Y A matrix of observed responses for \code{X}.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param prob If \code{TRUE}, the function returns a matrix of latent class probabilities. If \code{FALSE},
#'        it returns the index of the most likely latent class.

#'
#' @return A matrix of latent class probabilities by observation and by component.
#'
#' @importFrom matrixStats rowLogSumExps
#'
#'
#' @export PredictClassPosterior
PredictClassPosterior = function(Y, X, alpha, comp.dist, zero.prob, params.list, prob = TRUE)
{
  gate.ll = GateLogit(X, alpha)
  expert.list = DimCompExpertLL(Y, comp.dist, zero.prob, params.list)
  ll.list = GateExpertLL(alpha, gate.ll, expert.list, penalty=FALSE) # , hyper.alpha, hyper.params)
  weighting = EMEzkz(gate.ll, expert.list, ll.list)$z.e.obs

  if(prob==TRUE){
    return(weighting)
  }else{
    return(ColMaxIdx(weighting))
  }
}
