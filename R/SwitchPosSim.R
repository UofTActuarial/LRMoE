## A function to simulate positive part of observations y
#' A function to simulate positive part of observations y
#'
#' @param sample.size.n The size of y to simulate.
#' @param comp.dist A string which indicates the expert function.
#' @param params A vector of parameters for \code{ind.dist}.
#'
#' @importFrom stats rgamma rlnorm rweibull rpois rnbinom rbinom
#' @importFrom statmod rinvgauss
#' @importFrom actuar rburr
#' @importFrom rmutil rgammacount
#'
#'
#' @keywords internal
#'
#' @export SimPosY
SimPosY = function(sample.size.n, comp.dist, params)
{
  temp = NULL
  switch (comp.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = rgamma(sample.size.n, shape = params[1], scale = params[2])},
          "ZI-gamma"    = {temp = rgamma(sample.size.n, shape = params[1], scale = params[2])},
          "invgauss"    = {temp = statmod::rinvgauss(sample.size.n, mean = params[1], shape = params[2])},
          "ZI-invgauss" = {temp = statmod::rinvgauss(sample.size.n, mean = params[1], shape = params[2])},
          "lnorm"       = {temp = rlnorm(sample.size.n, meanlog = params[1], sdlog = params[2])},
          "ZI-lnorm"    = {temp = rlnorm(sample.size.n, meanlog = params[1], sdlog = params[2])},
          "weibull"     = {temp = stats::rweibull(sample.size.n, shape = params[1], scale = params[2])},
          "ZI-weibull"  = {temp = stats::rweibull(sample.size.n, shape = params[1], scale = params[2])},
          "burr"        = {temp = actuar::rburr(sample.size.n, shape1 = params[1], shape2 = params[2], scale = params[3])},
          "ZI-burr"     = {temp = actuar::rburr(sample.size.n, shape1 = params[1], shape2 = params[2], scale = params[3])},
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = rpois(sample.size.n, lambda = params[1])},
          "ZI-poisson"  = {temp = rpois(sample.size.n, lambda = params[1])},
          "nbinom"      = {temp = rnbinom(sample.size.n, size = params[1], prob = params[2])},
          "ZI-nbinom"   = {temp = rnbinom(sample.size.n, size = params[1], prob = params[2])},
          "binom"       = {temp = rbinom(sample.size.n, size = params[1], prob = params[2])},
          "ZI-binom"    = {temp = rbinom(sample.size.n, size = params[1], prob = params[2])},
          "gammacount"  = {temp = rmutil::rgammacount(sample.size.n, m = params[1], s = params[2])},
          "ZI-gammacount"  = {temp = rmutil::rgammacount(sample.size.n, m = params[1], s = params[2])},
          # Error
          stop("Invalid distribution!")
  )
  return(temp)
}
