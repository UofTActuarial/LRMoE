## Penalty of parameters
#' Computes penalty of parameters for expert functions by dimension and by component.
#'
#' @param ind.dist A string which indicates the expert function.
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
#' @param hyper.params A vector of parameters penalization for \code{params}.
#'     For simplicity, positive parameters are given Gamma priors, real parameters are
#'     given normal priors, while parameters with bounded ranges are not penalized at all.
#' \itemize{
#'     \item \code{gamma}: \code{(hyper.m.1, hyper.m.2, hyper.theta.1, hyper.theta.2)}.
#'           We assume \code{m} has a Gamma prior with given shape and scale parameters. Similar for \code{theta}.
#'     \item \code{lnorm}: \code{(hyper.meanlog.1, hyper.sdlog.1, hyper.sdlog.2)}
#'     \item \code{invgauss}: \code{(hyper.mean.mu.1, hyper.scale.lambda.1, hyper.scale.lambda.2)}
#'     \item \code{weibull}: \code{(hyper.shape.k.1, hyper.shape.k.2, hyper.scale.lambda.1, hyper.scale.lambda.2)}
#'     \item \code{burr}: \code{(hyper.shape1.k.1, hyper.shape1.k.2, hyper.shape2.c.1, hyper.shape2.c.2, hyper.scale.lambda.1, hyper.scale.lambda.2)}
#'     \item \code{poisson}: \code{(hyper.mean.theta.1)}
#'     \item \code{nbinom}: \code{(hyper.size.n.1, hyper.size.n.2)}
#'     \item \code{binom}: \code{()}
#'     \item \code{gammacount}: \code{(hyper.m,1, hyper.m.2, hyper.s.1, hyper.s.2)}
#' }
#'
#' @return A negative numerical value of penalty for a particular dimension and component.
#'
#' @seealso \code{\link{DimCompExpertLL}}.
#'
#'
#' @keywords internal
#'
#' @export DimCompExpertPenalty
DimCompExpertPenalty = function(ind.dist, params, hyper.params)
{
  temp = 0
  switch (ind.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] +
            (hyper.params[3]-1)*log(params[2]) - params[2]/hyper.params[4]  },
          "ZI-gamma"    = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] +
            (hyper.params[3]-1)*log(params[2]) - params[2]/hyper.params[4]  },
          "invgauss"    = {temp = 0 },
          "ZI-invgauss" = {temp = 0 },
          "lnorm"       = {temp = 0 },
          "ZI-lnorm"    = {temp = 0 },
          "weibull"     = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] },
          "ZI-weibull"  = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] },
          "burr"        = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] +
            (hyper.params[3]-1)*log(params[2]) - params[2]/hyper.params[4] +
            (hyper.params[5]-1)*log(params[3]) - params[3]/hyper.params[6]  },
          "ZI-burr"     = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] +
            (hyper.params[3]-1)*log(params[2]) - params[2]/hyper.params[4] +
            (hyper.params[5]-1)*log(params[3]) - params[3]/hyper.params[6]  },
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] },
          "ZI-poisson"  = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] },
          "nbinom"      = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] },
          "ZI-nbinom"   = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] },
          "binom"       = {temp = 0 },
          "ZI-binom"    = {temp = 0 },
          "gammacount"  = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] +
            (hyper.params[3]-1)*log(params[2]) - params[2]/hyper.params[4] },
          "ZI-gammacount" = {temp = (hyper.params[1]-1)*log(params[1]) - params[1]/hyper.params[2] +
            (hyper.params[3]-1)*log(params[2]) - params[2]/hyper.params[4] },
          # Error
          stop("Invalid distribution!")
  )

  return(temp)

}
