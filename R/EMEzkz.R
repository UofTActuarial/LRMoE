## Computation of z, k, z
#' ECM: E-Step for \code{z.e.obs}, \code{k.e} and \code{z.e.lat}.
#'
#' @param gate.ll An object returned by \code{\link{GateLogit}}.
#' @param expert.list An object returned by \code{\link{DimCompExpertLL}}.
#' @param ll.list An object returned by \code{\link{GateExpertLL}}.
#'
#' @return \code{z.e.obs},\code{k.e},\code{z.e.lat} Numerical vectors of length N.
#'
#' @keywords internal
#'
#' @importFrom matrixStats rowLogSumExps
#'
#' @export EMEzkz
EMEzkz = function(gate.ll, expert.list, ll.list)
{

  sample.size.n = ll.list$sample.size.n # sample size
  n.comp = ll.list$n.comp # no. of experts
  dim.m = ll.list$dim.m # no. of dimensions of y

  # Note: For the purpose of optimizing wrt gating weights alpha's,
  #       there is no need to consider z.zero and z.pos yet.
  z.e.obs       = array(0, dim = c(sample.size.n, n.comp))
  z.e.lat       =  array(0, dim = c(sample.size.n, n.comp))
  k.e           = array(0, dim = c(sample.size.n, 1))

  z.e.obs = exp(XColMinusY(ll.list$ll.ind, ll.list$comp.aggre.ll.ind))
    # exp( sweep(ll.list$ll.ind, 1, ll.list$comp.aggre.ll.ind, FUN = "-", check.margin = FALSE) )

  z.e.lat = exp(XColMinusY(ll.list$ll.ind.tn.bar, ll.list$comp.aggre.ll.ind.tn.bar))
    # exp( sweep(ll.list$ll.ind.tn.bar, 1, ll.list$comp.aggre.ll.ind.tn.bar, FUN = "-", check.margin = FALSE) )
  z.e.lat[is.na(z.e.lat)] = 1/n.comp

  k.e     = expm1(- ll.list$comp.aggre.ll.ind.tn)
    # exp( - ll.list$comp.aggre.ll.ind.tn ) - 1

  return(list(z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e))

}
