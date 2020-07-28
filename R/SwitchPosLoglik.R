## Positive part of Expert functions
#' A switch to calculate the loglikelihood of expert functions.
#'
#' @param ind.dist A string which indicates the expert function.
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
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param params A vector of parameters for \code{ind.dist}.
#'
#' @return A list of matrices of expert loglikelihood for the chosen distribution \code{ind.dist}.
#' \itemize{
#'     \item \code{expert.ll}: An N*1 matrix of loglikelihood: density for exact observation; probability for censored observation.
#'     \item \code{expert.tn}: An N*1 matrix of loglikelihood within truncation limits.
#'     \item \code{expert.tn.bar}: An N*1 matrix of loglikelihood outside truncation limits, which is \code{log(1-exp(expert.tn))}.
#' }
#'
#' @seealso \code{\link{ExpertGamma}}, \code{\link{ExpertLognormal}}, \code{\link{ExpertInvgauss}}, \code{\link{ExpertWeibull}}, \code{\link{ExpertBurr}},
#'          \code{\link{ExpertPoisson}}, \code{\link{ExpertNbinom}}, \code{\link{ExpertBinom}}, \code{\link{ExpertGammaCount}}
#'
#' @keywords internal
#'
#' @export PosExpertLL
PosExpertLL = function(ind.dist, tl, yl, yu, tu, params)
{
  temp = NULL
  switch (ind.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = ExpertGamma(tl, yl, yu, tu, params[1], params[2])},
          "ZI-gamma"    = {temp = ExpertGamma(tl, yl, yu, tu, params[1], params[2])},
          "invgauss"    = {temp = ExpertInvgauss(tl, yl, yu, tu, params[1], params[2])},
          "ZI-invgauss" = {temp = ExpertInvgauss(tl, yl, yu, tu, params[1], params[2])},
          "lnorm"       = {temp = ExpertLognormal(tl, yl, yu, tu, params[1], params[2])},
          "ZI-lnorm"    = {temp = ExpertLognormal(tl, yl, yu, tu, params[1], params[2])},
          "weibull"     = {temp = ExpertWeibull(tl, yl, yu, tu, params[1], params[2])},
          "ZI-weibull"  = {temp = ExpertWeibull(tl, yl, yu, tu, params[1], params[2])},
          "burr"        = {temp = ExpertBurr(tl, yl, yu, tu, params[1], params[2], params[3])},
          "ZI-burr"     = {temp = ExpertBurr(tl, yl, yu, tu, params[1], params[2], params[3])},
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = ExpertPoisson(tl, yl, yu, tu, params[1])},
          "ZI-poisson"  = {temp = ExpertPoisson(tl, yl, yu, tu, params[1])},
          "nbinom"      = {temp = ExpertNbinom(tl, yl, yu, tu, params[1], params[2])},
          "ZI-nbinom"   = {temp = ExpertNbinom(tl, yl, yu, tu, params[1], params[2])},
          "binom"       = {temp = ExpertBinom(tl, yl, yu, tu, params[1], params[2])},
          "ZI-binom"    = {temp = ExpertBinom(tl, yl, yu, tu, params[1], params[2])},
          "gammacount"  = {temp = ExpertGammaCount(tl, yl, yu, tu, params[1], params[2])},
          "ZI-gammacount"  = {temp = ExpertGammaCount(tl, yl, yu, tu, params[1], params[2])},
          # Error
          stop("Invalid distribution!")
  )
  return(temp)
}
