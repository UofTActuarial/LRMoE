## Expert Function: Log Normal
#' Expert Function: Log Normal.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param g A numeric which is set to 1. See Note.
#' @param meanlog A vector of length \code{g}: Lognormal mean parameters.
#' @param sdlog A vector of length \code{g}: Lognormal sd parameters.
#' @return A list of matrices of expert loglikelihood for Log Normal.
#'
#' @seealso \code{\link[stats]{Lognormal}}.
#'
#' @importFrom stats plnorm dlnorm
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
#' @export ExpertLognormal
ExpertLognormal = function(tl, yl, yu, tu, meanlog, sdlog)
{
  # Initialization: return value are N * g matrices
  expert.ll = expert.tn = expert.tn.bar = array(-Inf, dim=c(length(yu),1))

  # Find indexes of unequal yl & yu: Not exact observations, but censored
  censor.idx = (yl!=yu)
  prob.log.yu = plnorm(yu[censor.idx], meanlog = meanlog, sdlog = sdlog, log.p=TRUE)
  prob.log.yl = plnorm(yl[censor.idx], meanlog = meanlog, sdlog = sdlog, log.p=TRUE)

  # Compute loglikelihood for expert j, first for y
  expert.ll[censor.idx,1] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl) # likelihood of censored interval: some easy algebra
  expert.ll[!censor.idx,1] = dlnorm(yu[!censor.idx], meanlog = meanlog, sdlog = sdlog, log=TRUE) # exact likelihood

  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu=plnorm(tu, meanlog = meanlog, sdlog = sdlog, log.p=TRUE)
  prob.log.tl=plnorm(tl, meanlog = meanlog, sdlog = sdlog, log.p=TRUE)

  # Normalizing factor for truncation limits, in log
  expert.tn[,1]= prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)

  ###################################################################
  # Deal with no truncation case
  no.trunc.idx = (tl==tu)
  expert.tn[no.trunc.idx,1] = dlnorm(tu[no.trunc.idx], meanlog = meanlog, sdlog = sdlog, log=TRUE)
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
