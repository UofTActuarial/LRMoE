#' A switch to calculate the mean of a given distribution.
#'
#' @param comp.dist A string which indicates the expert function.
#' \itemize{
#'     \item \code{gamma}: Gamma
#'     \item \code{lnorm}: Log Normal
#'     \item \code{invgauss}: Inverse Gaussian
#'     \item \code{weibull}: Weibull
#'     \item \code{burr}: Burr
#'     \item \code{poisson}: Poisson
#'     \item \code{ztpoisson}: Zero-Truncated Poisson
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
#'      \item \code{ztpoisson}: \code{(mean.theta)}
#'     \item \code{nbinom}: \code{(size.n, prob.p)}
#'     \item \code{binom}: \code{(size.n, prob.p)}
#'     \item \code{gammacount}: \code{(m, s)}
#' }
#'
#'
#' @importFrom actuar mburr
#'
#' @keywords internal
#'
#' @export PosMean
PosMean = function(comp.dist, params)
{
  temp = NULL
  switch (comp.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = params[1]*params[2] },
          "ZI-gamma"    = {temp = params[1]*params[2] },
          "invgauss"    = {temp = params[1] },
          "ZI-invgauss" = {temp = params[1] },
          "lnorm"       = {temp = exp(params[1] + 0.5*params[2]*params[2]) },
          "ZI-lnorm"    = {temp = exp(params[1] + 0.5*params[2]*params[2]) },
          "weibull"     = {temp = params[2] * gamma(1+1/params[1]) },
          "ZI-weibull"  = {temp = params[2] * gamma(1+1/params[1]) },
          "burr"        = {temp = actuar::mburr(order = 1, shape1 = params[1], shape2 = params[2], scale = params[3]) },
          "ZI-burr"     = {temp = actuar::mburr(order = 1, shape1 = params[1], shape2 = params[2], scale = params[3]) },
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = params[1] },
          "ZI-poisson"  = {temp = params[1] },
          "ztpoisson"     = {temp = params[1]/(1-exp(-params[1])) },
          "ZI-ztpoisson"  = {temp = params[1]/(1-exp(-params[1])) },
          "nbinom"      = {temp = params[1] * (1-params[2])/params[2] },
          "ZI-nbinom"   = {temp = params[1] * (1-params[2])/params[2] },
          "binom"       = {temp = params[1] * params[2] },
          "ZI-binom"    = {temp = params[1] * params[2] },
          "gammacount"  = {temp = mgammacount(1, m = params[1], s = params[2]) },
          "ZI-gammacount"  = {temp = mgammacount(1, m = params[1], s = params[2]) }
  )
  return(temp)
}
