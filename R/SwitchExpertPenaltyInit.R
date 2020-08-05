## Initialize enalty of parameters
#' Initialize penalty of parameters for expert functions by dimension and by component.
#'
#' @param ind.dist A string which indicates the expert function.
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
#'
#' @return A list of hyper parameters by dimension and by component.
#'
#' @seealso \code{\link{DimCompExpertLL}}.
#'
#'
#' @keywords internal
#'
#' @export DimCompExpertPenaltyInit
DimCompExpertPenaltyInit = function(ind.dist)
{
  temp = 0
  switch (ind.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = c(2, 1, 2, 1)  },
          "ZI-gamma"    = {temp = c(2, 1, 2, 1)  },
          "invgauss"    = {temp = c(1, Inf, 1, Inf) },
          "ZI-invgauss" = {temp = c(1, Inf, 1, Inf) },
          "lnorm"       = {temp = c(Inf, 1, Inf) },
          "ZI-lnorm"    = {temp = c(Inf, 1, Inf) },
          "weibull"     = {temp = c(2, 1) },
          "ZI-weibull"  = {temp = c(2, 1) },
          "burr"        = {temp = c(2, 1, 2, 1, 2, 1)  },
          "ZI-burr"     = {temp = c(2, 1, 2, 1, 2, 1)  },
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = c(2, 1) },
          "ZI-poisson"  = {temp = c(2, 1) },
          "ztpoisson"     = {temp = c(2, 1) },
          "ZI-ztpoisson"  = {temp = c(2, 1) },
          "nbinom"      = {temp = c(2, 1) },
          "ZI-nbinom"   = {temp = c(2, 1) },
          "binom"      = {temp = c() },
          "ZI-binom"   = {temp = c() },
          "gammacount"  = {temp = c(2, 1, 2, 1) },
          "ZI-gammacount" = {temp = c(2, 1, 2, 1) },
          # Error
          stop("Invalid distribution!")
  )
  return(temp)

}
