## Computes the entire loglikelihood
#' Computes the entire loglikelihood, including gating function and experts.
#'
#' @param alpha A \code{g*P} matrix, where \code{g} is the number of components and P is the number of covariates.
#' @param gate.ll An object returned by function \code{\link{GateLogit}}.
#' @param expert.list An object returned by function \code{\link{DimCompExpertLL}}.
#' @param penalty TRUE/FALSE, which indicates whether parameter penalty should be applied. Default (and recommended) is TRUE.
#' @param hyper.alpha A numeric, which penalizes the magnitude of \code{alpha}.
#' @param hyper.params A list of length \code{d}. Each element is a sublist of length \code{g}.
#'                     Each element of a sublist is a vector of numerics, which penalizes expert parameters. See also \code{\link{DimCompExpertPenalty}}.
#'
#' @importFrom matrixStats rowLogSumExps
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
#' @export GateExpertLL
GateExpertLL = function(alpha, gate.ll, expert.list, penalty = TRUE, hyper.alpha = NULL, hyper.params = NULL)
{
  sample.size.n = nrow(expert.list[[1]]$expert.ll) # sample size
  n.comp = length(expert.list[[1]]$dim.comp.dist) # no. of experts
  dim.m = length(expert.list) # no. of dimensions of y
  # Note: dimension of gate.func is: sample.size.n * n.comp

  # First, aggregate the loglikelihood of experts BY DIMENSION of y
  dim.aggre.pos.expert.ll     = array(0, dim = c(sample.size.n, n.comp))
  dim.aggre.pos.expert.tn     = array(0, dim = c(sample.size.n, n.comp))
  dim.aggre.pos.expert.tn.bar = array(0, dim = c(sample.size.n, n.comp))
  dim.aggre.expert.ll         = array(0, dim = c(sample.size.n, n.comp))
  dim.aggre.expert.tn         = array(0, dim = c(sample.size.n, n.comp))
  dim.aggre.expert.tn.bar     = array(0, dim = c(sample.size.n, n.comp))

  for(k in 1:dim.m)
  {
    dim.aggre.pos.expert.ll     = dim.aggre.pos.expert.ll     + expert.list[[k]]$pos.expert.ll
    dim.aggre.pos.expert.tn     = dim.aggre.pos.expert.tn     + expert.list[[k]]$pos.expert.tn
    dim.aggre.pos.expert.tn.bar = dim.aggre.pos.expert.tn.bar + expert.list[[k]]$pos.expert.tn.bar
    dim.aggre.expert.ll         = dim.aggre.expert.ll         + expert.list[[k]]$expert.ll
    dim.aggre.expert.tn         = dim.aggre.expert.tn         + expert.list[[k]]$expert.tn
    dim.aggre.expert.tn.bar     = dim.aggre.expert.tn.bar     + expert.list[[k]]$expert.tn.bar
  }

  # Calculate loglik (without penalty). Dimension: sample.size.n * n.comp
  ll.pos.ind        = gate.ll + dim.aggre.pos.expert.ll
  ll.pos.ind.tn     = gate.ll + dim.aggre.pos.expert.tn
  # ll.pos.ind.tn.bar = gate.ll + dim.aggre.pos.expert.tn.bar
  ll.pos.ind.tn.bar = gate.ll + log1mexp(-dim.aggre.pos.expert.tn)

  ll.ind            = gate.ll + dim.aggre.expert.ll
  ll.ind.tn         = gate.ll + dim.aggre.expert.tn
  # ll.ind.tn.bar     = gate.ll + dim.aggre.expert.tn.bar
  ll.ind.tn.bar     = gate.ll + log1mexp(-dim.aggre.expert.tn)

  # Aggregate loglik by observation. Dimension: sample.size.n * 1
  comp.aggre.ll.pos.ind         = rowLogSumExps(ll.pos.ind)
  comp.aggre.ll.pos.ind.tn      = rowLogSumExps(ll.pos.ind.tn)
  comp.aggre.ll.pos.ind.tn.bar  = rowLogSumExps(ll.pos.ind.tn.bar)
  comp.aggre.ll.ind             = rowLogSumExps(ll.ind)
  comp.aggre.ll.ind.tn          = rowLogSumExps(ll.ind.tn)
  comp.aggre.ll.ind.tn.bar      = rowLogSumExps(ll.ind.tn.bar)

  # Normalize by the truncation factor. Dimension: sample.size.n * 1
  aggre.ll.pos = comp.aggre.ll.pos.ind  - comp.aggre.ll.pos.ind.tn
  aggre.ll     = comp.aggre.ll.ind      - comp.aggre.ll.ind.tn

  # Calculate penalty terms
  ll.penalty = 0
  if(penalty == TRUE)
  {
    # Penalty on alpha
    ll.penalty.alpha = - sweep(as.matrix(alpha^2), 2, 1/(2*hyper.alpha^2), FUN="*", check.margin = FALSE)

    # Penalty on component parameters
    penalty.mat = array(0, dim = c(dim.m, n.comp))
    for(k in 1:dim.m)
    {
      for(j in 1:n.comp)
      {
        penalty.mat[k,j] = DimCompExpertPenalty(expert.list[[k]]$dim.comp.dist[j], expert.list[[k]]$pos.params[[j]], hyper.params[[k]][[j]])
      }
    }
    ll.penalty = sum(ll.penalty.alpha) + sum(penalty.mat)
  }

  # Calculate entire loglik
  ll.np = sum(aggre.ll)
  ll = sum(aggre.ll) + ll.penalty

  return(list(# Constants
    sample.size.n = sample.size.n, # sample size
    n.comp = n.comp, # no. of experts
    dim.m = dim.m, # no. of dimensions of y
    # Numbers
    ll = ll, ll.np = ll.np, ll.penalty = ll.penalty,
    # Dimension: sample.size.n * 1
    aggre.ll = aggre.ll, aggre.ll.pos = aggre.ll.pos,
    # Dimension: sample.size.n * 1
    comp.aggre.ll.ind = comp.aggre.ll.ind, comp.aggre.ll.ind.tn = comp.aggre.ll.ind.tn, comp.aggre.ll.ind.tn.bar = comp.aggre.ll.ind.tn.bar,
    comp.aggre.ll.pos.ind = comp.aggre.ll.pos.ind, comp.aggre.ll.pos.ind.tn = comp.aggre.ll.pos.ind.tn, comp.aggre.ll.pos.ind.tn.bar = comp.aggre.ll.pos.ind.tn.bar,
    # Dimension: sample.size.n * n.comp. WITH gate function
    ll.ind = ll.ind, ll.ind.tn = ll.ind.tn, ll.ind.tn.bar = ll.ind.tn.bar,
    ll.pos.ind = ll.pos.ind, ll.pos.ind.tn = ll.pos.ind.tn, ll.pos.ind.tn.bar = ll.pos.ind.tn.bar,
    # Dimension: sample.size.n * n.comp. NO gate function. Expert ONLY
    dim.aggre.expert.ll = dim.aggre.expert.ll, dim.aggre.expert.tn = dim.aggre.expert.tn, dim.aggre.expert.tn.bar = dim.aggre.expert.tn.bar,
    dim.aggre.pos.expert.ll = dim.aggre.pos.expert.ll, dim.aggre.pos.expert.tn = dim.aggre.pos.expert.tn, dim.aggre.pos.expert.tn.bar = dim.aggre.pos.expert.tn.bar
  )
  )

}
