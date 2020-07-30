## A function to format model output: add names for parameters
#' A function to format model output: add names for parameters
#'
#' @seealso \code{\link{pos.expert.loglik.calc}}
#'
#' @keywords internal
#'
#' @export PosName
PosName = function(comp.dist)
{
  temp = NULL
  switch (comp.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = c("shape", "scale") },
          "ZI-gamma"    = {temp = c("shape", "scale") },
          "invgauss"    = {temp = c("mean", "scale") },
          "ZI-invgauss" = {temp = c("mean", "scale") },
          "lnorm"       = {temp = c("meanlog", "sdlog") },
          "ZI-lnorm"    = {temp = c("meanlog", "sdlog") },
          "weibull"     = {temp = c("shape", "scale") },
          "ZI-weibull"  = {temp = c("shape", "scale") },
          "burr"        = {temp = c("shape1", "shape2", "scale") },
          "ZI-burr"     = {temp = c("shape1", "shape2", "scale") },
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = c("mean") },
          "ZI-poisson"  = {temp = c("mean") },
          "nbinom"      = {temp = c("size", "prob") },
          "ZI-nbinom"   = {temp = c("size", "prob") },
          "binom"      = {temp = c("size", "prob") },
          "ZI-binom"   = {temp = c("size", "prob") },
          "gammacount"  = {temp = c("m", "s") },
          "ZI-gammacount"  = {temp = c("m", "s") },
          # Error
          stop("Invalid distribution!")
  )
  return(temp)
}
