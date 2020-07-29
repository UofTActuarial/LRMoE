## Computation of z
#' ECM: E-Step for \code{z.zero.e.obs}.
#'
#' @param yl.k A vector of length N: lower bounds of censoring.
#' @param comp.kj.zero.inflation TRUE/FALSE. See returned value \code{zero.inflation} in \code{\link{DimCompExpertLL}}.
#' @param comp.kj.zero.prob.old A numeric between 0 and 1. See returned value \code{zero.prob} in \code{\link{DimCompExpertLL}}.
#' @param comp.kj.pos.expert.ll See returned value \code{pos.expert.ll} in \code{\link{DimCompExpertLL}}
#'
#' @return \code{z.zero.e.obs} Numerical vectors of length N.
#'
#' @keywords internal
#'
#' @importFrom copula log1pexp
#' @export EMEzzeroobs
EMEzzeroobs = function(yl.k, comp.kj.zero.inflation, comp.kj.zero.prob.old, comp.kj.pos.expert.ll)
{
  sample.size.n = length(yl.k)
  temp = array(0, dim = c(sample.size.n, 1))

  if(comp.kj.zero.inflation == TRUE) # Otherwise, there is no possibility of zero component.
  {
    obs.y.zero.idx = (yl.k==0) # yl=0: possible for observation to be from zero component
    # Update z.zero only for those indices. Otherwise, they stay zero.
    temp[obs.y.zero.idx] = # exp( - log1pexp( log(1/comp.kj.zero.prob.old-1) + comp.kj.pos.expert.ll[obs.y.zero.idx]) )
      comp.kj.zero.prob.old / (comp.kj.zero.prob.old + (1-comp.kj.zero.prob.old)*exp(comp.kj.pos.expert.ll[obs.y.zero.idx]))
  }

  return(temp)
}
