#' A function to calculate the mean of a matrix of given distributions (positive part only), by dimension and by component.
#'
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param paramas.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#' @return A matrix of mean values by dimension and by component.
#'
#' @keywords internal
#'
# #' @export dim.comp.mean.y.pos
DimCompPosYMean = function(comp.dist, params.list)
{
  n.comp = ncol(comp.dist)
  dim.m = nrow(comp.dist)

  result = matrix(0, nrow = dim.m, ncol = n.comp)

  for(k in 1:dim.m)
  {
    for(j in 1:n.comp)
    {
      result[k, j] = PosMean(comp.dist[k,j], params.list[[k]][[j]])
    }
  }

  return(result)
}


#' A function to calculate the mean of a matrix of given distributions (with zero inflation), by dimension and by component.
#'
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param zero.prob A d*g matrix of numbers between 0 and 1, describing zero probability masses by dimension and by component.
#' @param paramas.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#'
#' @return A matrix of mean values by dimension and by component.
#'
#' @keywords internal
#'
# #' @export DimCompYMean
DimCompYMean = function(comp.dist, zero.prob, params.list)
{
  temp = DimCompPosYMean(comp.dist, params.list)
  result = (1-zero.prob)*temp

  return(result)
}


# #' Predict the mean of y, given a fixed covariate vector X and a model.
# #'
# #' @param X A matrix of covariates.
# #' @param alpha A matrix of logit regression coefficients.
# #' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
# #'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
# #' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
# #' @param params.list A list of length d, where each element is a sublist of length g.
# #'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
# #'
# #'
# #' @return A vector of mean values by dimension.
# #'
# #' @keywords internal
# #'
# #' @export ind.predict.mean
# ind.predict.mean = function(X, alpha, comp.dist, zero.prob, params.list)
# {
#   n.comp = nrow(alpha)
#   dim.m = nrow(comp.dist)
#
#   weighting = exp(gate.logit(X, alpha))
#   temp = dim.comp.mean.y(comp.dist, zero.prob, params.list)
#   result = temp%*%t(weighting)
#
#   return(t(result))
# }


#' Predict the mean of y, given a fixed covariate matrix X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#'
#' @return A matrix of mean values by observation and by dimension.
#'
#' @export PredictMeanPrior
PredictMeanPrior = function(X, alpha, comp.dist, zero.prob, params.list)
{
  X.size = nrow(X)
  dim.m = nrow(comp.dist)
  result = array(0, dim = c(X.size, dim.m))

  meantmp = DimCompYMean(comp.dist, zero.prob, params.list)
  weighting = PredictClassPrior(X, alpha, TRUE)

  for(d in 1:dim.m){
    result[,d] = tcrossprod(weighting, t(meantmp[d,]))
  }

  return(result)

}


#' Predict the mean of y, given a fixed covariate matrix X and a model.
#'
#' @param Y A matrix of observed responses for \code{X}.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#'
#' @return A matrix of mean values by observation and by dimension.
#'
#' @export PredictMeanPosterior
PredictMeanPosterior = function(Y, X, alpha, comp.dist, zero.prob, params.list)
{
  X.size = nrow(X)
  dim.m = nrow(comp.dist)
  result = array(0, dim = c(X.size, dim.m))

  meantmp = DimCompYMean(comp.dist, zero.prob, params.list)
  weighting = PredictClassPosterior(Y, X, alpha, comp.dist, zero.prob, params.list, TRUE)

  for(d in 1:dim.m){
    result[,d] = tcrossprod(weighting, t(meantmp[d,]))
  }

  return(result)

}
