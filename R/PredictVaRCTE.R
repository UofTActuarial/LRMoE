#' A function to calculate the cdf of a matrix of given distributions (positive part only), by dimension and by component.
#'
#' @param y A numeric vector by dimension.
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param paramas.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#'
#' @return A matrix of cdf values by dimension and by component.
#'
#' @keywords internal
#'
# #' @export DimCompPosYCDF
DimCompPosYCDF = function(y, comp.dist, params.list)
{
  n.comp = ncol(comp.dist)
  dim.m = nrow(comp.dist)

  result = matrix(0, nrow = dim.m, ncol = n.comp)

  for(k in 1:dim.m)
  {
    for(j in 1:n.comp)
    {
      result[k, j] = PosCDF(y[k], comp.dist[k,j], params.list[[k]][[j]])
    }
  }

  return(result)
}

#' A function to calculate the cdf of a matrix of given distributions (with zero inflation), by dimension and by component.
#'
#' @param y A numeric vector by dimension.
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param zero.prob A d*g matrix of numbers between 0 and 1, describing zero probability masses by dimension and by component.
#' @param paramas.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#'
#'
#' @return A matrix of cdf values by dimension and by component.
#'
#' @keywords internal
#'
# #' @export DimCompYCDF
DimCompYCDF = function(y, comp.dist, zero.prob, params.list)
{
  temp = DimCompPosYCDF(y, comp.dist, params.list)
  result = zero.prob*as.numeric(y>=0) + (1-zero.prob)*temp

  return(result)
}


#' Calculate the cdf of y, given a fixed covariate vector X and a model.
#'
#' @param y A numeric vector by dimension.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#'
#' @return A vector of cdf values by dimension.
#'
#' @keywords internal
#'
# #' @export indCDF
indCDFPrior = function(y, X, alpha, comp.dist, zero.prob, params.list)
{
  n.comp = nrow(alpha)
  dim.m = nrow(comp.dist)

  # X.alpha = X %*% t(alpha)
  weighting = weighting = PredictClassPrior(X, alpha, TRUE)
  temp = DimCompYCDF(y, comp.dist, zero.prob, params.list)
  result = tcrossprod(weighting, temp) #  temp%*%t(weighting)

  # return(t(result))
  return(result)
}

#' Predict the VaR of y, given a fixed covariate vector X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param prob A vector of probabilities. Default is a vector of length d of 0.95, if no value is provided.
#'
#' @return A vector of VaR by dimension.
#'
#' @keywords internal
#'
#' @importFrom stats uniroot
#'
# #' @export ind.predict.quantile
#'
indSolveQuantilePrior = function(X, alpha, comp.dist, zero.prob, params.list, prob = NULL)
{
  dim.m = nrow(comp.dist)
  result = rep(NA, dim.m)

  if(is.null(prob))
  {
    prob = rep(0.95, dim.m)
  }

  for(k in 1:dim.m)
  {

    temp.func = function(y, X, alpha, comp.dist, zero.prob, params.list, quant)
    {
      return( indCDFPrior(rep(y, dim.m), X, alpha, comp.dist, zero.prob, params.list)[1,k] - quant )
    }

    result[k] = tryCatch(
      {uniroot(f = temp.func, interval = c(0, 99999999999),
               tol = .Machine$double.eps^0.25, maxiter = 1000,
               X = X, alpha = alpha,
               comp.dist = comp.dist, zero.prob = zero.prob, params.list = params.list,
               quant = prob[k])$root},
      error = function(e){return(NA)})
  }
  return(result)
}

#' Predict the VaR of y, given a fixed covariate matrix X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param prob A vector of probabilities. Default is a vector of length d of 0.95, if no value is provided.
#'
#'
#' @return A matrix of VaR.
#'
#'
#' @export PredictVaRPrior
PredictVaRPrior = function(X, alpha, comp.dist, zero.prob, params.list, prob = NULL)
{
  X.size = nrow(X)
  dim.m = nrow(comp.dist)

  if(is.null(prob))
  {
    prob = rep(0.95, dim.m)
  }

  if(is.null(X.size))
  {
    return( indSolveQuantilePrior(prob = prob, X, alpha, comp.dist, zero.prob, params.list)  )
  }

  VaR = array(0, dim = c(X.size, dim.m))

  VaR = apply(X, MARGIN = 1, FUN = indSolveQuantilePrior,
                 prob = prob, alpha = alpha, comp.dist = comp.dist, zero.prob = zero.prob, params.list = params.list)

  return(VaR)

}


#' Predict the CTE of y, given a fixed covariate matrix X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param prob A vector of probabilities. Default is a vector of length d of 0.95, if no value is provided.
#'
#'
#' @return A vector of CTE by dimension.
#'
#' @keywords internal
#'
# #' @export indSolveCTEPrior
#'
indSolveCTEPrior = function(X, alpha, comp.dist, zero.prob, params.list, prob = NULL)
{
  dim.m = nrow(comp.dist)
  if(is.null(prob))
  {
    prob = rep(0.95, dim.m)
  }

  if(is.null(dim(X))){
    X = t(data.matrix(X))
  }

  value.at.risk = indSolveQuantilePrior(X, alpha, comp.dist, zero.prob, params.list, prob) #ind.predict.quantile(prob, X, alpha, comp.dist, zero.prob, params.list)
  expect.excess = PredictLimExPrior(X, alpha, comp.dist, zero.prob, params.list, value.at.risk)[[2]]
  result = value.at.risk + expect.excess / (1-prob)

  return(result)
}

#' Predict the CTE of y, given a fixed covariate matrix X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param prob A vector of probabilities. Default is a vector of length d of 0.95, if no value is provided.
#'
#'
#' @return A matrix of CTE.
#'
#'
#' @export PredictCTEPrior
PredictCTEPrior = function(X, alpha, comp.dist, zero.prob, params.list, prob = NULL)
{
  X.size = nrow(X)
  dim.m = nrow(comp.dist)

  if(is.null(prob))
  {
    prob = rep(0.95, dim.m)
  }

  if(is.null(X.size))
  {
    return( indSolveCTEPrior(X, alpha, comp.dist, zero.prob, params.list, prob)  )
  }

  CTE = array(0, dim = c(X.size, dim.m))

  CTE = apply(X, MARGIN = 1, FUN = indSolveCTEPrior,
                 prob = prob, alpha = alpha, comp.dist = comp.dist, zero.prob = zero.prob, params.list = params.list)
  return(t(CTE))

}
