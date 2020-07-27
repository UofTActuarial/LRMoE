## Expert Function: Burr
#' Expert Function: Burr.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param shape1.k A vector of length \code{1}: Burr shape1 parameters.
#' @param shape2.c A vector of length \code{1}: Burr shape2 parameters.
#' @param scale.lambda A vector of length \code{1}: Burr scale parameters.
#' @return A list of matrices of expert loglikelihood for Burr.
#'
#' @seealso \code{\link[actuar]{Burr}}.
#'
#' @importFrom actuar pburr dburr
#' @importFrom copula log1mexp log1pexp
#'
#' @keywords internal
#'
#' @export ExpertBurr
ExpertBurr = function(tl, yl, yu, tu, shape1.k, shape2.c, scale.lambda)
{
  # Initialization: return value are N * 1 matrices
  expert.ll = expert.tn = expert.tn.bar = array(-Inf, dim=c(length(yu),1))

  # Find indexes of unequal yl & yu: Not exact observations, but censored
  censor.idx = (yl!=yu)
  prob.log.yu = actuar::pburr(yu[censor.idx], shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)
  prob.log.yl = actuar::pburr(yl[censor.idx], shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)

  # Compute loglikelihood for expert j, first for y
  expert.ll[censor.idx,1] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl) # likelihood of censored interval: some easy algebra
  expert.ll[!censor.idx,1] = actuar::dburr(yu[!censor.idx], shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, log = TRUE)  # exact likelihood

  ###################################################################
  # Deal with numerical underflow: prob.log.yu and prob.log.yl can both be -Inf
  # NA.idx = which(is.na(expert.ll[,j]))
  expert.ll[which(is.na(expert.ll[,1])), 1] = -Inf
  ###################################################################

  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu = actuar::pburr(tu, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)
  prob.log.tl = actuar::pburr(tl, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)

  # Normalizing factor for truncation limits, in log
  expert.tn[,1]= prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)

  ###################################################################
  # Deal with numerical underflow: prob.log.tu and prob.log.tl can both be -Inf
  # NA.idx = which(is.na(expert.tn[,j]))
  expert.tn[which(is.na(expert.tn[,1])), 1] = -Inf
  ###################################################################

  ###################################################################
  # Deal with no truncation case
  no.trunc.idx = (tl==tu)
  expert.tn[no.trunc.idx,1] = actuar::dburr(tu[no.trunc.idx], shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, log = TRUE)
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

