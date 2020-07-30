## Computation of z
#' ECM: M-Step for \code{zero.prob}.
#'
#' @param z.zero.e.obs An object returned by \code{\link{EMEzzeroobs}}.
#' @param z.pos.e.obs An object returned by \code{\link{EMEzzeroobs}}.
#' @param z.zero.e.lat An object returned by \code{\link{EMEzzerolat}}.
#' @param z.pos.e.lat An object returned by \code{\link{EMEzzerolat}}.
#' @param k.e An object returned by \code{\link{EMEzkz}}.
#'
#' @return \code{zero.prob} Updated zero.prob.
#'
#' @keywords internal
#'
#' @export EMMzero
EMMzero = function(z.zero.e.obs, z.pos.e.obs, z.zero.e.lat, z.pos.e.lat, k.e)
{
  term.zero = XPlusYColTimesZ(z.zero.e.obs, z.zero.e.lat, k.e)
    # sweep(matrix(z.zero.e.obs), 1, sweep(matrix(z.zero.e.lat), 1, matrix(k.e), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE)
    # z.zero.e.obs  + sweep(z.zero.e.lat, 1, k.e, FUN = "*", check.margin = FALSE)
  term.pos  = XPlusYColTimesZ(z.pos.e.obs, z.pos.e.lat, k.e)
    #  sweep(matrix(z.pos.e.obs), 1, sweep(matrix(z.pos.e.lat), 1, matrix(k.e), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE)
    # z.pos.e.obs   + sweep(z.pos.e.lat, 1, k.e, FUN = "*", check.margin = FALSE)

  numerator   = sum(term.zero) # apply(term.zero, 2, sum)
  denominator = numerator + sum(term.pos) # apply(term.pos, 2, sum) # apply(term.zero, 2, sum) + apply(term.pos, 2, sum)
  return( numerator / denominator )
}
