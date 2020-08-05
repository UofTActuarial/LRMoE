#' A switch to calculate the cdf of a given distribution.
#'
#' @param y A numeric vector at which the CDF is to be calculated.
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
#'     \item \code{ztpoisson}: \code{(mean.theta)}
#'     \item \code{nbinom}: \code{(size.n, prob.p)}
#'     \item \code{binom}: \code{(size.n, prob.p)}
#'     \item \code{gammacount}: \code{(m, s)}
#' }
#'
#' @importFrom stats pgamma plnorm pweibull ppois pnbinom pbinom
#' @importFrom statmod pinvgauss
#' @importFrom actuar pburr
#' @importFrom countreg pztpois
#'
#' @keywords internal
#'
#' @export PosCDF
PosCDF = function(y, comp.dist, params)
{
  temp = NULL
  switch (comp.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = ifelse(y<=0, 0, pgamma(y, shape = params[1], scale = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "ZI-gamma"    = {temp = ifelse(y<=0, 0, pgamma(y, shape = params[1], scale = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "invgauss"    = {temp = ifelse(y<=0, 0, statmod::pinvgauss(y, mean = params[1], shape = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "ZI-invgauss" = {temp = ifelse(y<=0, 0, statmod::pinvgauss(y, mean = params[1], shape = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "lnorm"       = {temp = ifelse(y<=0, 0, plnorm(y, meanlog = params[1], sdlog = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "ZI-lnorm"    = {temp = ifelse(y<=0, 0, plnorm(y, meanlog = params[1], sdlog = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "weibull"     = {temp = ifelse(y<=0, 0, pweibull(y, shape = params[1], scale = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "ZI-weibull"  = {temp = ifelse(y<=0, 0, pweibull(y, shape = params[1], scale = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "burr"        = {temp = ifelse(y<=0, 0, actuar::pburr(y, shape1 = params[1], shape2 = params[2], scale = params[3], lower.tail = TRUE, log.p = FALSE)) },
          "ZI-burr"     = {temp = ifelse(y<=0, 0, actuar::pburr(y, shape1 = params[1], shape2 = params[2], scale = params[3], lower.tail = TRUE, log.p = FALSE)) },
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = ifelse(y<0, 0, ppois(y, lambda = params[1], lower.tail = TRUE, log.p = FALSE)) },
          "ZI-poisson"  = {temp = ifelse(y<0, 0, ppois(y, lambda = params[1], lower.tail = TRUE, log.p = FALSE)) },
          "ztpoisson"     = {temp = ifelse(y<0, 0, pztpois(y, lambda = params[1], lower.tail = TRUE, log.p = FALSE)) },
          "ZI-ztpoisson"  = {temp = ifelse(y<0, 0, pztpois(y, lambda = params[1], lower.tail = TRUE, log.p = FALSE)) },
          "nbinom"      = {temp = ifelse(y<0, 0, pnbinom(y, size = params[1], prob = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "ZI-nbinom"   = {temp = ifelse(y<0, 0, pnbinom(y, size = params[1], prob = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "binom"       = {temp = ifelse(y<0, 0, pbinom(y, size = params[1], prob = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "ZI-binom"    = {temp = ifelse(y<0, 0, pbinom(y, size = params[1], prob = params[2], lower.tail = TRUE, log.p = FALSE)) },
          "gammacount"  = {temp = ifelse(y<0, 0, pgammacount(y, m = params[1], s = params[2], log.p = FALSE)) },
          "ZI-gammacount"  = {temp = ifelse(y<0, 0, pgammacount(y, m = params[1], s = params[2], log.p = FALSE)) }
  )
  return(temp)
}
