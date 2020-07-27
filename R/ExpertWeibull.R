## Expert Function: Weibull
#' Expert Function: Weibull.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param shape.k A vector of length \code{1}: Weibull shape parameters.
#' @param scale.lambda A vector of length \code{1}: Weibull scale parameters.
#' @return A list of matrices of expert loglikelihood for Weibull.
#'
#' @seealso \code{\link[stats]{Weibull}}.
#'
#' @importFrom stats pweibull dweibull
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
#' @export ExpertWeibull
ExpertWeibull = function(tl, yl, yu, tu, shape.k, scale.lambda)
{
  # Initialization: return value are N * 1 matrices
  expert.ll = expert.tn = expert.tn.bar = array(-Inf, dim=c(length(yu),1))

  # Find indexes of unequal yl & yu: Not exact observations, but censored
  censor.idx = (yl!=yu)
  prob.log.yu = stats::pweibull(yu[censor.idx], shape = shape.k, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)
  prob.log.yl = stats::pweibull(yl[censor.idx], shape = shape.k, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)

  # Compute loglikelihood for expert j, first for y
  expert.ll[censor.idx,1] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl) # likelihood of censored interval: some easy algebra
  expert.ll[!censor.idx,1] = stats::dweibull(yu[!censor.idx], shape = shape.k, scale = scale.lambda, log = TRUE) # exact likelihood

  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu = stats::pweibull(tu, shape = shape.k, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)
  prob.log.tl = stats::pweibull(tl, shape = shape.k, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)

  # Normalizing factor for truncation limits, in log
  expert.tn[,1]= prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)

  ###################################################################
  # Deal with no truncation case
  no.trunc.idx = (tl==tu)
  expert.tn[no.trunc.idx,1] = stats::dweibull(tu[no.trunc.idx], shape = shape.k, scale = scale.lambda, log = TRUE)
  ###################################################################

  ###################################################################
  # Deal with exact zero case
  zero.idx = (tu==0)
  expert.ll[zero.idx,1]=(-Inf)
  expert.tn[zero.idx,1]=(-Inf)
  ###################################################################

  # Log of Pr(outside of truncation interval)
  expert.tn.bar[!no.trunc.idx,1] = log1mexp(-expert.tn[!no.trunc.idx,1])
  expert.tn.bar[no.trunc.idx,1] = 0

  # Return values
  list(expert.ll = expert.ll, expert.tn = expert.tn, expert.tn.bar = expert.tn.bar)
}
