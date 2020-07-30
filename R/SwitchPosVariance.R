#' A switch to calculate the variance of a given distribution.
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
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @importFrom actuar mburr
#'
#' @keywords internal
#'
#' @export PosVariance
PosVariance = function(comp.dist, params)
{
  temp = NULL
  switch (comp.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = params[1]*params[2]*params[2] },
          "ZI-gamma"    = {temp = params[1]*params[2]*params[2] },
          "invgauss"    = {temp = params[1]^3/params[2] },
          "ZI-invgauss" = {temp = params[1]^3/params[2] },
          "lnorm"       = {temp = exp(2*params[1] + params[2]*params[2]) * (exp(params[2]*params[2])-1) },
          "ZI-lnorm"    = {temp = exp(2*params[1] + params[2]*params[2]) * (exp(params[2]*params[2])-1) },
          "weibull"     = {temp = params[2]^2 * ( gamma(1+2/params[1]) - (gamma(1+1/params[1]))^2 ) },
          "ZI-weibull"  = {temp = params[2]^2 * ( gamma(1+2/params[1]) - (gamma(1+1/params[1]))^2 ) },
          "burr"        = {temp = actuar::mburr(order = 2, shape1 = params[1], shape2 = params[2], scale = params[3]) -
            (actuar::mburr(order = 1, shape1 = params[1], shape2 = params[2], scale = params[3]))^2 },
          "ZI-burr"     = {temp = actuar::mburr(order = 2, shape1 = params[1], shape2 = params[2], scale = params[3]) -
            (actuar::mburr(order = 1, shape1 = params[1], shape2 = params[2], scale = params[3]))^2 },
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = params[1] },
          "ZI-poisson"  = {temp = params[1] },
          "nbinom"      = {temp = params[1] * (1-params[2])/(params[2])^2 },
          "ZI-nbinom"   = {temp = params[1] * (1-params[2])/(params[2])^2 },
          "binom"       = {temp = params[1] * params[2] * (1-params[2]) },
          "ZI-binom"    = {temp = params[1] * params[2] * (1-params[2]) },
          "gammacount"  = {temp = mgammacount(2, m = params[1], s = params[2])[1] - (mgammacount(1, m = params[1], s = params[2])[1])^2 },
          "ZI-gammacount"  = {temp = mgammacount(2, m = params[1], s = params[2])[1] - (mgammacount(1, m = params[1], s = params[2])[1])^2 }
  )
  return(temp)
}
