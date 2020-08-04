## Expert Function: Gamma Count
#' Expert Function: Gamma Count.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param m A vector of length \code{1}: Gamma Count shape parameters.
#' @param s A vector of length \code{1}: Gamma Count dispersion parameters.
#' @return A list of matrices of expert loglikelihood for Gamma Count.
#'
#' @seealso \code{\link[rmutil]{GammaCount}}.
#'
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
#' @export ExpertGammaCount
ExpertGammaCount = function(tl, yl, yu, tu, m, s)
{
  # Initialization: return value are N * 1 matrices
  expert.ll = expert.tn = expert.tn.bar = array(-Inf, dim=c(length(yu),1))

  # Find indexes of unequal yl & yu: Not exact observations, but censored
  censor.idx = (yl!=yu)
  prob.log.yu = ifelse(yu[censor.idx]==Inf, 0, pgammacount(yu[censor.idx], m = m, s = s, log.p = TRUE))
  prob.log.yl = pgammacount(ceiling(yl[censor.idx])-1, m = m, s = s, log.p = TRUE)

  # Compute loglikelihood for expert j, first for y
  expert.ll[censor.idx,1] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl) # likelihood of censored interval: some easy algebra
  expert.ll[!censor.idx,1] = dgammacount(yl[!censor.idx], m = m, s = s, log = TRUE)  # exact likelihood

  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu = ifelse(tu==Inf, 0, pgammacount(tu, m = m, s = s, log.p = TRUE))
  prob.log.tl = pgammacount(ceiling(tl)-1, m = m, s = s, log.p = TRUE)

  # Normalizing factor for truncation limits, in log
  expert.tn[,1]= prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)

  ###################################################################
  # Deal with no truncation case
  no.trunc.idx = (tl==tu)
  expert.tn[no.trunc.idx,1] = dgammacount(tu[no.trunc.idx], m =  m, s = s, log = TRUE)
  ###################################################################

  ###################################################################
  # Deal with exact zero case: The following code treating zeros is NOT applicable for frequency distributions!!!
  # zero.idx = (tu==0)
  # expert.ll[zero.idx,1]=(-Inf)
  # expert.tn[zero.idx,1]=(-Inf)
  ###################################################################

  # Log of Pr(outside of truncation interval)
  expert.tn.bar[!no.trunc.idx,1] = log1mexp(-expert.tn[!no.trunc.idx,1])
  expert.tn.bar[no.trunc.idx,1] = log1mexp(-expert.tn[no.trunc.idx,1])

  # Return values
  list(expert.ll = expert.ll, expert.tn = expert.tn, expert.tn.bar = expert.tn.bar)
}
