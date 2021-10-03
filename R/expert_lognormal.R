######################################################################################
# Initialize the parameters of the distribuion
######################################################################################
lognormal.params_init <- function(y) {
  log_y = log(y)
  meanlog = mean(log_y)
  sdlog = sd(log_y)
  return( list(meanlog = meanlog, sdlog = sdlog) )
}

lognormal.exposurize <- function(params, exposure) {
  return( list(meanlog = params[["meanlog"]], sdlog = params[["sdlog"]]) )
}

lognormal.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
lognormal.expert_ll_exact <- function(y, params) {
  return( lognormal.logpdf(params = params, x = y) )
}

lognormal.expert_ll_not_exact = function(tl, yl, yu, tu, params)
{
  # Check dim tl == yl == tu == yu
  # Check value

  # Find indexes of unequal yl & yu: Not exact observations, but censored
  expert.ll = expert.tn = expert.tn.bar = rep(-Inf, length(yu))
  censor.idx = (yl!=yu)
  no.trunc.idx = (tl==tu)

  # Expert LL calculation
  ######################################################################################
  prob.log.yu = lognormal.logcdf(params, yu[censor.idx])
  prob.log.yl = lognormal.logcdf(params, yl[censor.idx])

  # Compute loglikelihood for expert j, first for y
  # likelihood of censored interval: some easy algebra
  expert.ll[censor.idx] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl)
  # exact likelihood
  expert.ll[!censor.idx] = lognormal.logpdf(params,yu[!censor.idx])
  ######################################################################################

  # Expert TN Calculation
  ######################################################################################
  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu=lognormal.logcdf(params, tu)
  prob.log.tl=lognormal.logcdf(params, tl)

  # Normalizing factor for truncation limits, in log
  expert.tn= prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)

  # Deal with no truncation case
  expert.tn[no.trunc.idx] = lognormal.logpdf(params, tu[no.trunc.idx])

  # Deal with exact zero case
  zero.idx = (tu==0)
  expert.ll[zero.idx]=(-Inf)
  expert.tn[zero.idx]=(-Inf)
  ######################################################################################

  # Expert TN Bar Calculation
  ######################################################################################
  expert.tn.bar[!no.trunc.idx] = log1mexp(-expert.tn[!no.trunc.idx])
  expert.tn.bar[no.trunc.idx] = 0
  ######################################################################################

  # Return values
  return( list(expert_ll = expert.ll, expert_tn = expert.tn, expert_tn_bar = expert.tn.bar) )
}

lognormal.penalty = function(params, penalty_params){
  if(!length(penalty_params)) { penalty_params = c(Inf, 1.0, Inf) }
  # return( (params[["meanlog"]]/penalty_params[1])^2 + (penalty_params[2] - 1)*log(params[["sdlog"]]) - params[["sdlog"]]/penalty_params[3] )
  return(0)
}

lognormal.default_penalty <- function() {
  return(c(Inf, 1.0, Inf))
}
######################################################################################
# ddistribution, pdistribution, qdistribution and rdistribution implementations.
######################################################################################
lognormal.simulation <- function(params, n) {
  return( rlnorm(n = n, meanlog = params[["meanlog"]], sdlog = params[["sdlog"]]) )
}

lognormal.mean <- function(params) {
  mean = exp(params[["meanlog"]] + (params[["sdlog"]]^2)/2)
  return( mean )
}

lognormal.variance <- function(params) {
  variance = (exp(params[["sdlog"]]^2) - 1) * exp(2*params[["meanlog"]] + params[["sdlog"]]^2)
  return( variance )
}

lognormal.logpdf <- function(params, x) {
  return( dlnorm(x = x, meanlog = params[["meanlog"]], sdlog = params[["sdlog"]], log = T) )
}

lognormal.pdf <- function(params, x) {
  return( dlnorm(x = x, meanlog = params[["meanlog"]], sdlog = params[["sdlog"]], log = F) )
}

lognormal.logcdf <- function(params, q) {
  return( plnorm(q = q, meanlog = params[["meanlog"]], sdlog = params[["sdlog"]], log.p = T) )
}

lognormal.cdf <- function(params, q) {
  return( plnorm(q = q, meanlog = params[["meanlog"]], sdlog = params[["sdlog"]], log.p = F) )
}

lognormal.quantile <- function(params, p) {
  return( qlnorm(p = p, meanlog = params[["meanlog"]], sdlog = params[["sdlog"]]) )
}

lognormal.lev <- function(params, u) {
  mu = params[["meanlog"]]
  s = params[["sdlog"]]

  return(exp(mu + 0.5 * s^2) * pnorm(log(u), mean = mu + s^2, sd = s) + u * (1-pnorm(log(u), mean = mu, sd = s)))
}
######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
lognormal.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # What is the usage of expert old, penalty and pen_params
  logY_e_obs = log(ye)
  logY_sq_e_obs = logY_e_obs^2

  pos_idx = (ye != 0)
  term_zkz = z_e_obs[pos_idx]
  term_zkz_logY = (z_e_obs[pos_idx] * logY_e_obs[pos_idx])
  term_zkz_logY_sq = (z_e_obs[pos_idx] * logY_sq_e_obs[pos_idx])

  meanlog_new = sum(term_zkz_logY) / sum(term_zkz)
  sdlog_new = sqrt( 1/sum(term_zkz) *
                      (sum(term_zkz_logY_sq) - 2*meanlog_new*sum(term_zkz_logY)+ meanlog_new^2*sum(term_zkz)) )

  return(list(meanlog = meanlog_new, sdlog = sdlog_new))
}

lognormal.EM_notexact <- function(expert_old,
                              tl, yl, yu, tu,
                              exposure,
                              z_e_obs, z_e_lat, k_e,
                              penalty, pen_params) {

  # Old parameters
  meanlog_old = expert_old$get_params()$meanlog
  sdlog_old = expert_old$get_params()$sdlog

  # Old loglikelihoods
  expert_ll = rep(-Inf, length(yl))
  expert_tn = rep(-Inf, length(yl))
  expert_tn_bar = rep(-Inf, length(yl))

  for(i in 1:length(yl)){
    expert_expo = expert_old$exposurize(exposure[i])
    result_set = expert_expo$ll_not_exact(tl[i], yl[i], yu[i], tu[i])
    expert_ll[i] = result_set[["expert_ll"]]
    expert_tn[i] = result_set[["expert_tn"]]
    expert_tn_bar[i] = result_set[["expert_tn_bar"]]
  }

  # Additional E-step
  cencor_idx = (yl!=yu)
  logY_e_obs = rep(0, length(yl))
  logY_e_lat = rep(0, length(yl))
  logY_sq_e_obs = rep(0, length(yl))
  logY_sq_e_lat = rep(0, length(yl))

  # Uncensored observations
  logY_e_obs[!cencor_idx] = log(yl[!cencor_idx])
  logY_sq_e_obs[!cencor_idx] = (log(yl[!cencor_idx]))^2

  # Conditional expectation of log(y): untruncated but censored case.
  diff.dens.untrunc = exp(-0.5 * (( log(yl[cencor_idx])-meanlog_old )/sdlog_old)^2 ) -
    exp(-0.5 * (( log(yu[cencor_idx])-meanlog_old )/sdlog_old)^2 )
  diff.dist.untrunc = pnorm(log(yu[cencor_idx]), meanlog_old, sdlog_old, lower.tail=TRUE, log.p=FALSE) -
    pnorm(log(yl[cencor_idx]), meanlog_old, sdlog_old, lower.tail=TRUE, log.p=FALSE)
  logY_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) *
    ( (sdlog_old) * (1/sqrt(2*pi)) * (0.5) * diff.dens.untrunc + meanlog_old * diff.dist.untrunc)
  # Conditional expectation of log(y): truncated case.
  diff.dens.trunc = exp(-0.5 * (( log(tl)-meanlog_old )/sdlog_old)^2 ) - exp(-0.5 * (( log(tu)-meanlog_old )/sdlog_old)^2 )
  diff.dist.trunc = pnorm(log(tu), meanlog_old, sdlog_old, lower.tail=TRUE, log.p=FALSE) -
    pnorm(log(tl), meanlog_old, sdlog_old, lower.tail=TRUE, log.p=FALSE)
  logY_e_lat = exp(-expert_tn_bar) * ( meanlog_old - ( (sdlog_old) * (1/sqrt(2*pi)) * (0.5) * diff.dens.trunc + meanlog_old * diff.dist.trunc) )
  logY_e_lat[is.na(logY_e_lat)] = 0 # Hardcode: to prevent NaN

  # Conditional expectation of (log(y))^2: untruncated but censored case.
  z.dens.z.func = function(z)
  {
    temp = z * exp( -0.5 * (z^2) )
    temp[is.infinite(z)] = 0
    temp[is.na(temp)] = 0
    return(temp)
  }
  if(sum(cencor_idx!=0))
  {
    diff.ydensy.untrunc = sapply( (log(yl[cencor_idx])-meanlog_old)/sdlog_old, FUN = z.dens.z.func ) -
      sapply( (log(yu[cencor_idx])-meanlog_old)/sdlog_old, FUN = z.dens.z.func )
    logY_sq_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) *
      ( (meanlog_old^2 + sdlog_old^2) * diff.dist.untrunc + meanlog_old*sdlog_old*(1/sqrt(2*pi)) * diff.dens.untrunc +
          (sdlog_old^2)*(1/sqrt(2*pi)) * diff.ydensy.untrunc )
  }
  # Conditional expectation of (log(y))^2: truncated case.
  diff.ydensy.trunc = sapply( (log(tl)-meanlog_old)/sdlog_old, FUN = z.dens.z.func ) -
    sapply( (log(tu)-meanlog_old)/sdlog_old, FUN = z.dens.z.func )
  logY_sq_e_lat = exp(-expert_tn_bar) *
    ( (meanlog_old^2+sdlog_old^2) -
        ( (meanlog_old^2 + sdlog_old^2) * diff.dist.trunc + meanlog_old*sdlog_old*(1/sqrt(2*pi)) *
            diff.dens.trunc + (sdlog_old^2)*(1/sqrt(2*pi)) * diff.ydensy.trunc ) )
  logY_sq_e_lat[is.na(logY_sq_e_lat)] = 0 # Hardcode: to prevent NaN


  # Update parameters
  pos_idx = (yu!=0)

  term_zkz = XPlusYZ(z_e_obs[pos_idx], z_e_lat[pos_idx], k_e[pos_idx])
  term_zkz_logY = XAPlusYZB(z_e_obs[pos_idx], logY_e_obs[pos_idx], z_e_lat[pos_idx], k_e[pos_idx], logY_e_lat[pos_idx])
  term_zkz_logY_sq = XAPlusYZB(z_e_obs[pos_idx], logY_sq_e_obs[pos_idx], z_e_lat[pos_idx], k_e[pos_idx], logY_sq_e_lat[pos_idx])

  meanlog_new = sum(term_zkz_logY) / sum(term_zkz)

  sdlog_new = sqrt( 1/sum(term_zkz) *
                      ( sum(term_zkz_logY_sq) - 2*meanlog_new*sum(term_zkz_logY) + (meanlog_new^2)*sum(term_zkz) ) )

  return(list(meanlog = meanlog_new, sdlog = sdlog_new))
}
######################################################################################
# Register the distribution at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################
