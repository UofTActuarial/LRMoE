#' A switch to calculate the limited mean of a given distribution.
#'
#' @param comp.dist A string which indicates the expert function.
#' \itemize{
#'     \item \code{gamma}: Gamma
#'     \item \code{lnorm}: Log Normal
#'     \item \code{invgauss}: Inverse Gaussian
#'     \item \code{weibull}: Weibull
#'     \item \code{burr}: Burr
#'     \item \code{poisson}: Poisson
#'     \item \code{nbinom}: Negative Binomial
#'     \item \code{binom}: Binomial
#'     \item \code{gammacount}: Gamma Count
#'     \item \code{ZI-root}: Zero-inflated versions of the distributions above, e.g. \code{ZI-gamma}.
#' }
#' @param params A vector of parameter values.
#' \itemize{
#'     \item \code{gamma}: \code{(m, theta)}
#'     \item \code{lnorm}: \code{(meanlog, sdlog)}
#'     \item \code{invgauss}: \code{(mean.mu, scale.lambda)}
#'     \item \code{weibull}: \code{(shape.k, scale.lambda)}
#'     \item \code{burr}: \code{(shape1.k, shape2.c, scale.lambda)}
#'     \item \code{poisson}: \code{(mean.theta)}
#'     \item \code{nbinom}: \code{(size.n, prob.p)}
#'     \item \code{binom}: \code{(size.n, prob.p)}
#'     \item \code{gammacount}: \code{(m, s)}
#' }
#'
#' @param limit A numeric of limit value.
#'
#' @importFrom actuar levgamma levinvgauss levlnorm levweibull levburr
#'
#' @keywords internal
#'
#' @export PosLimit
PosLimit = function(comp.dist, params, limit)
{
  temp = NULL
  switch (comp.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = actuar::levgamma(limit = limit, shape = params[1], scale = params[2], order = 1) },
          "ZI-gamma"    = {temp = actuar::levgamma(limit = limit, shape = params[1], scale = params[2], order = 1) },
          "invgauss"    = {temp = actuar::levinvgauss(limit = limit, mean = params[1], shape = params[2], order = 1) },
          "ZI-invgauss" = {temp = actuar::levinvgauss(limit = limit, mean = params[1], shape = params[2], order = 1) },
          "lnorm"       = {temp = actuar::levlnorm(limit = limit, meanlog = params[1], sdlog = params[2], order = 1) },
          "ZI-lnorm"    = {temp = actuar::levlnorm(limit = limit, meanlog = params[1], sdlog = params[2], order = 1) },
          "weibull"     = {temp = actuar::levweibull(limit = limit, shape = params[1], scale = params[2], order = 1) },
          "ZI-weibull"  = {temp = actuar::levweibull(limit = limit, shape = params[1], scale = params[2], order = 1) },
          "burr"        = {temp = actuar::levburr(limit = limit, shape1 = params[1], shape2 = params[2], scale = params[3], order = 1) },
          "ZI-burr"     = {temp = actuar::levburr(limit = limit, shape1 = params[1], shape2 = params[2], scale = params[3], order = 1) },
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = NA },
          "ZI-poisson"  = {temp = NA },
          "nbinom"      = {temp = NA },
          "ZI-nbinom"   = {temp = NA },
          "binom"       = {temp = NA },
          "ZI-binom"    = {temp = NA },
          "gammacount"  = {temp = NA },
          "ZI-gammacount"  = {temp = NA }
  )
  return(temp)
}
