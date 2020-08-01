#' Count the number of parameters in the positive part.
#'
#' @param comp.dist A d*g matrix of strings describing the component distributions.
#'
#' @return The number of parameters in the positive part.
#'
#'
#' @keywords internal
#'
#' @export PosCountParams
PosCountParams = function(comp.dist)
{
  result = 0

  for(k in 1:(nrow(comp.dist)))
  {
    for(j in 1:(ncol(comp.dist)))
    {
      switch (comp.dist[k,j],
              # Severity distributions & their zero-inflation
              "gamma"       = {result = result + 2},
              "ZI-gamma"    = {result = result + 2},
              "invgauss"    = {result = result + 2},
              "ZI-invgauss" = {result = result + 2},
              "lnorm"       = {result = result + 2},
              "ZI-lnorm"    = {result = result + 2},
              "weibull"     = {result = result + 2},
              "ZI-weibull"  = {result = result + 2},
              "burr"        = {result = result + 3},
              "ZI-burr"     = {result = result + 3},
              # Frequency distributions & their zero-inflation
              "poisson"     = {result = result + 1},
              "ZI-poisson"  = {result = result + 1},
              "nbinom"      = {result = result + 2},
              "ZI-nbinom"   = {result = result + 2},
              "binom"       = {result = result + 1},
              "ZI-binom"    = {result = result + 1},
              "gammacount"  = {result = result + 2},
              "ZI-gammacount"  = {result = result + 2}
      )
    }
  }
  return(result)
}
