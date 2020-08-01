## Simulate y, given a fixed covariate vector X and a model
#' Simulate y, given a fixed covariate vector X and a model
#'
#' @param X A P-element vector of covariates. The first element must be 1.
#' @param alpha A g*P matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param zero.prob A d*g matrix of numbers between 0 and 1, describing zero probability masses by dimension and by component.
#' @param paramas.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the initial parameter guess for the corresponding \code{comp.dist}.
#' @param sample.size Sample size to simulate. Default is 1.
#'
#' @return A matrix of simulated values.
#'
#' @keywords internal
#'
#' @importFrom stats rmultinom rbinom
#'
#' @export SimYInd
SimYInd = function(X, alpha, comp.dist, zero.prob, params.list,
                              sample.size = 1)
{
  n.comp = nrow(alpha)
  dim.m = nrow(comp.dist)

  weighting = exp(GateLogit(X, alpha))

  # sim.pi = array(0, dim = c(sample.size, n.comp))
  # for(i in c(1:sample.size))
  # {
  #   sim.pi[i,] = t(rmultinom(1, size = 1, prob = c( exp(X.alpha) / sum(exp(X.alpha))) ))
  # }

  sim.pi = matrix(apply(weighting, 1, rmultinom, n = sample.size, size = 1), nrow = sample.size, byrow = TRUE)

  y.comp.sim = array(0, dim = c(sample.size, n.comp))
  y.sim = array(0, dim = c(sample.size, dim.m))

  for(k in 1:dim.m)
  {
    for(j in 1:n.comp)
    {
      # Zero mass
      zero.ind = rbinom(sample.size, 1, zero.prob[k,j])
      # Positive part
      y.pos = SimPosY(sample.size, comp.dist[k,j], params.list[[k]][[j]])
      # Combine two parts
      y.comp.sim[,j] = zero.ind * 0 + (1-zero.ind) * y.pos
    }
    y.sim[,k] = apply(sim.pi*y.comp.sim, 1, FUN = sum)
  }

  return(y.sim)
}
