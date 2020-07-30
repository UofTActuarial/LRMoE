#' Generate a series of density values, given a model
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
#' @param y.series A vector of y values.
#'
#' @return A vector of density values evaluatd at y given.
#'
#' @importFrom stats dgamma dlnorm dweibull dpois dnbinom dbinom
#' @importFrom statmod dinvgauss
#' @importFrom actuar dburr
#'
#' @keywords internal
#'
#' @export PosDensity
PosDensity = function(comp.dist, params, y.series)
{
  temp = NULL
  switch (comp.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = dgamma(y.series, shape = params[1], scale = params[2], log = FALSE) },
          "ZI-gamma"    = {temp = dgamma(y.series, shape = params[1], scale = params[2], log = FALSE) },
          "invgauss"    = {temp = statmod::dinvgauss(y.series, mean = params[1], shape = params[2], log = FALSE) },
          "ZI-invgauss" = {temp = statmod::dinvgauss(y.series, mean = params[1], shape = params[2], log = FALSE) },
          "lnorm"       = {temp = dlnorm(y.series, meanlog = params[1], sdlog = params[2], log = FALSE) },
          "ZI-lnorm"    = {temp = dlnorm(y.series, meanlog = params[1], sdlog = params[2], log = FALSE) },
          "weibull"     = {temp = dweibull(y.series, shape = params[1], scale = params[2], log = FALSE) },
          "ZI-weibull"  = {temp = dweibull(y.series, shape = params[1], scale = params[2], log = FALSE) },
          "burr"        = {temp = actuar::dburr(y.series, shape1 = params[1], shape2 = params[2], scale = params[3], log = FALSE) },
          "ZI-burr"     = {temp = actuar::dburr(y.series, shape1 = params[1], shape2 = params[2], scale = params[3], log = FALSE) },
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = dpois(y.series, lambda = params[1], log = FALSE) },
          "ZI-poisson"  = {temp = dpois(y.series, lambda = params[1], log = FALSE) },
          "nbinom"      = {temp = dnbinom(y.series, size = params[1], prob = params[2], log = FALSE) },
          "ZI-nbinom"   = {temp = dnbinom(y.series, size = params[1], prob = params[2], log = FALSE) },
          "binom"       = {temp = dbinom(y.series, size = params[1], prob = params[2], log = FALSE) },
          "ZI-binom"    = {temp = dbinom(y.series, size = params[1], prob = params[2], log = FALSE) },
          "gammacount"  = {temp = dgammacount(y.series, m = params[1], s = params[2], log = FALSE) },
          "ZI-gammacount"  = {temp = dgammacount(y.series, m = params[1], s = params[2], log = FALSE) }
  )
  return(temp)
}
