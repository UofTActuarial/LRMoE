## Simulate y, given a fixed covariate matrix X and a model
#' Simulate y, given a fixed covariate matrix X and a model
#'
#' @param X A N*P matrix of of covariates. The first column must be 1. Each row may be different.
#' @param alpha A g*P matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param zero.prob A d*g matrix of numbers between 0 and 1, describing zero probability masses by dimension and by component.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the initial parameter guess for the corresponding \code{comp.dist}.
#'
#' @return A matrix of simulated values, where each row represents a policyholder and each column a dimension of the response variable.
#'
#' @export SimYSet
SimYSet = function(X, alpha, comp.dist, zero.prob, params.list)
{
  X.size = nrow(X)
  dim.m = nrow(comp.dist)

  # y.sim = array(0, dim = c(X.size, dim.m))
  # for(i in 1:(X.size) )
  # {
  #   y.sim[i,] = ind.data.simulator(X[i,], alpha, comp.dist, zero.prob, params.list, sample.size = 1)
  # }

  y.sim = apply(X, MARGIN = 1, FUN = SimYInd,
                alpha = alpha, comp.dist = comp.dist,
                zero.prob = zero.prob, params.list = params.list,
                sample.size = 1)

  # return(t(y.sim))
  # return(t(t(y.sim)))
  return(y.sim)
}
