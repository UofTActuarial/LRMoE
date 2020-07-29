## Computation of z
#' ECM: E-Step for \code{z.zero.e.lat}.
#'
#' @param tl.k A vector of length N: lower bounds of truncation.
#' @param comp.kj.zero.inflation TRUE/FALSE. See returned value \code{zero.inflation} in \code{\link{DimCompExpertLL}}.
#' @param comp.kj.zero.prob.old A numeric between 0 and 1. See returned value \code{zero.prob} in \code{\link{DimCompExpertLL}}.
#' @param comp.kj.pos.expert.tn.bar See returned value \code{pos.expert.tn.bar} in \code{\link{DimCompExpertLL}}
#'
#' @return \code{z.zero.e.lat} Numerical vectors of length N.
#'
#' @keywords internal
#'
#' @importFrom copula log1pexp
#' @export EMEzzerolat
EMEzzerolat = function(tl.k, comp.kj.zero.inflation, comp.kj.zero.prob.old, comp.kj.pos.expert.tn.bar)
{
  sample.size.n = length(tl.k)
  temp = array(0, dim = c(sample.size.n, 1))

  if(comp.kj.zero.inflation == TRUE) # Otherwise, there is no possibility of zero component.
  {
    lat.y.zero.idx = (tl.k>0)
    # Note that the expert.pos part already contains (1-zero.mass)
    temp[lat.y.zero.idx] = # exp( - log1pexp( log(1/comp.kj.zero.prob.old-1) + comp.kj.pos.expert.tn.bar[lat.y.zero.idx]) )
      comp.kj.zero.prob.old / (comp.kj.zero.prob.old + (1-comp.kj.zero.prob.old)*exp(comp.kj.pos.expert.tn.bar[lat.y.zero.idx]))
  }

  return(temp)
}
