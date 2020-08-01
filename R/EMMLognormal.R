## ECM algorithm of LogNormal expert
#' ECM: M-Step for LognOrmal expert.
#'
#' @importFrom stats pnorm
#'
#' @keywords internal
#'
#' @export EMMLognormal
EMMLognormal = function(params.old,
                        tl, yl, yu, tu,
                        expert.ll, expert.tn, expert.tn.bar,
                        z.e.obs, z.e.lat, k.e,
                        penalty, hyper.params)
{
  # Constants
  sample.size.n = length(tl)
  # Take old parameters
  meanlog = params.old[1]
  sdlog   = params.old[2]
  # Take hyper parameters
  hyper.meanlog = hyper.params[1]
  hyper.sdlog.1   = hyper.params[2]
  hyper.sdlog.2   = hyper.params[3]
  # Value to return
  params.new = params.old


  # E-Step: conditional expectations for log(y), (log(y))^2
  censor.idx = (yl!=yu)
  y.log.e.obs = rep(0, sample.size.n) # array(0, dim = c(sample.size.n, 1))
  y.log.e.lat = rep(0, sample.size.n) # array(0, dim = c(sample.size.n, 1))
  y.log.sq.e.obs = rep(0, sample.size.n) # array(0, dim = c(sample.size.n, 1))
  y.log.sq.e.lat = rep(0, sample.size.n) # array(0, dim = c(sample.size.n, 1))

  # Conditional expectation of log(y): untruncated and uncensored case.
  y.log.e.obs[!censor.idx] = log(yl[!censor.idx])
  # Conditional expectation of log(y): untruncated but censored case.
  diff.dens.untrunc = exp(-0.5 * (( log(yl[censor.idx])-meanlog )/sdlog)^2 ) - exp(-0.5 * (( log(yu[censor.idx])-meanlog )/sdlog)^2 )
  diff.dist.untrunc = pnorm(log(yu[censor.idx]), meanlog, sdlog, lower.tail=TRUE, log.p=FALSE) - pnorm(log(yl[censor.idx]), meanlog, sdlog, lower.tail=TRUE, log.p=FALSE)
  y.log.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * ( (sdlog) * (1/sqrt(2*pi)) * (0.5) * diff.dens.untrunc + meanlog * diff.dist.untrunc)
  # Conditional expectation of log(y): truncated case.
  diff.dens.trunc = exp(-0.5 * (( log(tl)-meanlog )/sdlog)^2 ) - exp(-0.5 * (( log(tu)-meanlog )/sdlog)^2 )
  diff.dist.trunc = pnorm(log(tu), meanlog, sdlog, lower.tail=TRUE, log.p=FALSE) - pnorm(log(tl), meanlog, sdlog, lower.tail=TRUE, log.p=FALSE)
  y.log.e.lat = exp(-expert.tn.bar) * ( meanlog - ( (sdlog) * (1/sqrt(2*pi)) * (0.5) * diff.dens.trunc + meanlog * diff.dist.trunc) )
  y.log.e.lat[is.na(y.log.e.lat)] = 0 # Hardcode: to prevent NaN

  # Conditional expectation of (log(y))^2: untruncated and uncensored case.
  y.log.sq.e.obs[!censor.idx] = (log(yl[!censor.idx]))^2
  # Conditional expectation of (log(y))^2: untruncated but censored case.
  z.dens.z.func = function(z)
  {
    temp = z * exp( -0.5 * (z^2) )
    temp[which(z==Inf)] = 0
    temp[which(z==-Inf)] = 0
    temp[which(is.na(temp))] = 0
    return(temp)
  }
  if(sum(censor.idx!=0))
  {
    diff.ydensy.untrunc = sapply( (log(yl[censor.idx])-meanlog)/sdlog, FUN = z.dens.z.func ) - sapply( (log(yu[censor.idx])-meanlog)/sdlog, FUN = z.dens.z.func )
    y.log.sq.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * ( (meanlog^2 + sdlog^2) * diff.dist.untrunc + meanlog*sdlog*(1/sqrt(2*pi)) * diff.dens.untrunc + (sdlog^2)*(1/sqrt(2*pi)) * diff.ydensy.untrunc )
  }
  # Conditional expectation of (log(y))^2: truncated case.
  diff.ydensy.trunc = sapply( (log(tl)-meanlog)/sdlog, FUN = z.dens.z.func ) - sapply( (log(tu)-meanlog)/sdlog, FUN = z.dens.z.func )
  y.log.sq.e.lat = exp(-expert.tn.bar) * ( (meanlog^2+sdlog^2) - ( (meanlog^2 + sdlog^2) * diff.dist.trunc + meanlog*sdlog*(1/sqrt(2*pi)) * diff.dens.trunc + (sdlog^2)*(1/sqrt(2*pi)) * diff.ydensy.trunc ) )
  y.log.sq.e.lat[is.na(y.log.sq.e.lat)] = 0 # Hardcode: to prevent NaN

  # M-Step: maximization of loglik
  # I should only use those y's that are POSITIVE for severity distributions
  pos.idx = (yu!=0)
  # Assuming no penalty, the paramater updates have closed-form solutions.
  # term.zkz = z.e.obs + k.e * z.e.lat
  # term.zkz.logy = z.e.obs * y.log.e.obs + k.e * z.e.lat * y.log.e.lat
  # term.zkz.logy.sq = z.e.obs * y.log.sq.e.obs + k.e * z.e.lat * y.log.sq.e.lat
  #
  # meanlog.new = sum(term.zkz.logy[pos.idx]) / sum(term.zkz[pos.idx])
  # sdlog.new = sqrt( 1/sum(term.zkz[pos.idx]) * ( sum(term.zkz.logy.sq[pos.idx]) - 2*meanlog.new*sum(term.zkz.logy[pos.idx]) + (meanlog.new^2)*sum(term.zkz[pos.idx]) ) )

  term.zkz = XPlusYZ(z.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx])
    # XPlusYColTimesZ(z.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx])
    # sweep(matrix(z.e.obs[pos.idx]), 1,
    #               sweep(matrix(k.e[pos.idx]), 1, matrix(z.e.lat[pos.idx]), FUN = "*", check.margin = FALSE),
    #               FUN = "+", check.margin = FALSE)
    # z.e.obs[pos.idx] + k.e[pos.idx] * z.e.lat[pos.idx]
  term.zkz.logy = XAPlusYZB(z.e.obs[pos.idx], y.log.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx], y.log.e.lat[pos.idx])
    # sweep(
    # sweep(matrix(z.e.obs[pos.idx]), 1, matrix(y.log.e.obs[pos.idx]), FUN = "*", check.margin = FALSE), 1,
    # sweep(matrix(k.e[pos.idx]), 1,
    #       sweep(matrix(z.e.lat[pos.idx]), 1, matrix(y.log.e.lat[pos.idx]), FUN = "*", check.margin = FALSE),
    #       FUN = "*", check.margin = FALSE),
    # FUN = "+", check.margin = FALSE)
    # z.e.obs[pos.idx] * y.log.e.obs[pos.idx] + k.e[pos.idx] * z.e.lat[pos.idx] * y.log.e.lat[pos.idx]
  term.zkz.logy.sq = XAPlusYZB(z.e.obs[pos.idx], y.log.sq.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx], y.log.sq.e.lat[pos.idx])
    # sweep(
    # sweep(matrix(z.e.obs[pos.idx]), 1, matrix(y.log.sq.e.obs[pos.idx]), FUN = "*", check.margin = FALSE), 1,
    # sweep(matrix(k.e[pos.idx]), 1,
    #       sweep(matrix(z.e.lat[pos.idx]), 1, matrix(y.log.sq.e.lat[pos.idx]), FUN = "*", check.margin = FALSE),
    #       FUN = "*", check.margin = FALSE),
    # FUN = "+", check.margin = FALSE)
    # z.e.obs[pos.idx] * y.log.sq.e.obs[pos.idx] + k.e[pos.idx] * z.e.lat[pos.idx] * y.log.sq.e.lat[pos.idx]

  meanlog.new = sum(term.zkz.logy) / sum(term.zkz)
    # apply(term.zkz.logy, 2, FUN = "sum") / apply(term.zkz, 2, FUN = "sum")
    # sum(term.zkz.logy) / sum(term.zkz)
  sdlog.new = sqrt( 1/sum(term.zkz) * ( sum(term.zkz.logy.sq) - 2*meanlog.new*sum(term.zkz.logy) + (meanlog.new^2)*sum(term.zkz) ) )
    # sqrt( 1/apply(term.zkz, 2, FUN = "sum") * ( apply(term.zkz.logy.sq, 2, FUN = "sum") - 2*meanlog.new*apply(term.zkz.logy, 2, FUN = "sum") + (meanlog.new^2)*apply(term.zkz, 2, FUN = "sum") ) )
    # sqrt( 1/sum(term.zkz) * ( sum(term.zkz.logy.sq) - 2*meanlog.new*sum(term.zkz.logy) + (meanlog.new^2)*sum(term.zkz) ) )

  # Update the parameters and return the result
  params.new[1] = meanlog.new
  params.new[2] = sdlog.new

  # print(meanlog.new)
  # print(sdlog.new)

  return(params.new)

}
