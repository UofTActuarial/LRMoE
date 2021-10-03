# gammacount Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the gammacount
######################################################################################
gammacount.params_init <- function(y){
  # Calculate all the parameters that needed for further calculation
  mu = mean(y)
  variance = var(y)
  s_init = variance / mu
  m_init = mu/(s_init)

  return( list(m = m_init, s = s_init) )
}

gammacount.exposurize <- function(params, exposure){
  # Calculate the exposures
  return( list(m = params[["m"]], s = params[["s"]]/exposure) )
}

gammacount.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
gammacount.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  return( gammacount.logpdf(params = params, x = y) )
}

gammacount.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # Check dim tl == yl == tu == yu
  # Check value

  # Find indexes of unequal yl & yu: Not exact observations, but censored
  expert.ll = expert.tn = expert.tn.bar = rep(-Inf, length(yu))
  censor.idx = (yl!=yu)
  no.trunc.idx = (tl==tu)

  # Expert LL calculation
  ######################################################################################
  prob.log.yu = ifelse(yu[censor.idx]==Inf, 0, gammacount.logcdf(params, yu[censor.idx]))
  prob.log.yl = gammacount.logcdf(params, ceiling(yl[censor.idx])-1)
  # Compute loglikelihood for expert j, first for y
  # likelihood of censored interval: some easy algebra
  expert.ll[censor.idx] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl)
  # exact likelihood
  expert.ll[!censor.idx] = gammacount.logpdf(params, yl[!censor.idx])
  ######################################################################################

  # Expert TN Calculation
  ######################################################################################
  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu = ifelse(tu==Inf, 0, gammacount.logcdf(params, tu))
  prob.log.tl = gammacount.logcdf(params, ceiling(tl)-1)

  # Normalizing factor for truncation limits, in log
  expert.tn = prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)

  # Deal with no truncation case
  expert.tn[no.trunc.idx] = gammacount.logpdf(params, tu[no.trunc.idx])
  ######################################################################################

  # Expert TN Bar Calculation
  ######################################################################################
  expert.tn.bar[!no.trunc.idx] = log1mexp(-expert.tn[!no.trunc.idx])
  expert.tn.bar[no.trunc.idx] = log1mexp(-expert.tn[no.trunc.idx])
  ######################################################################################

  # Return values
  return( list(expert_ll = expert.ll, expert_tn = expert.tn, expert_tn_bar = expert.tn.bar) )
}

gammacount.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  if(!length(penalty_params)) { penalty_params = c(2.0, 10.0, 2.0, 10.0) }
  p = penalty_params
  d.m = params[["m"]]
  d.s = params[["s"]]
  return( (p[1]-1)*log(d.m) - d.m/p[2] + (p[3]-1)*log(d.s) - d.s/p[4] )
}

gammacount.default_penalty <- function() {
  return(c(2.0, 10.0, 2.0, 10.0))
}

######################################################################################
# dgammacount, pgammacount, qgammacount and rgammacount implementations.
######################################################################################
gammacount.simulation <- function(params, n) {
  return( rgammacount(n, params[["m"]], params[["s"]]) )
}

gammacount.mean <- function(params) {
  return( mgammacount(1, m = params[["s"]], s = params[["s"]]) )
}

gammacount.variance <- function(params) {
  # Calculate the variance based on the params
  return( mgammacount(2, m = params[["s"]], s = params[["s"]]) -
            (mgammacount(1, m = params[["s"]], s = params[["s"]]))^2 )
}

gammacount.logpdf <- function(params, x) {
  return( dgammacount(y = x, m = params[["m"]], s = params[["s"]], log = TRUE) )
}

gammacount.pdf <- function(params, x) {
  # Return the pdf based on the input x
  return( dgammacount(y = x, m = params[["m"]], s = params[["s"]]) )
}

gammacount.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  return( pgammacount(q, m = params[["m"]], s = params[["s"]], log.p = TRUE) )
}

gammacount.cdf <- function(params, q) {
  # return the cdf based on the input x
  return( pgammacount(q, m = params[["m"]], s = params[["s"]]) )
}

gammacount.quantile <- function(params, p) {
  return( qgammacount(p, m = params[["m"]], s = params[["s"]]) )
}

######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
# gammacount._EStep <- function() {
#   # Perform the E step
#   NULL
# }
#
# gammacount._MStep <- function() {
#   # Perform the M step
#   NULL
# }
#
# gammacount.compute_EM <- function() {
#   # Perform the EM optimization
#   NULL
# }

gammacount.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations

  m_old = expert_old$get_params()$m
  s_old = expert_old$get_params()$s

  lognew = stats::optim(par = c(log(m_old), log(s_old)),
                             fn = gammacount_optim_params_exact,
                             expert_old = expert_old,
                             ye = ye,
                             exposure = exposure,
                             z_e_obs = z_e_obs,
                             penalty = penalty,
                             pen_params = pen_params,
                             # lower =  max(log(k_old)-2.0, 0.0),
                             # upper = log(k_old)+2.0,
                             method = "Nelder-Mead", #  "L-BFGS-B",
                             # lower = c(log(m_old)-2, log(s_old)-2),
                             # upper = c(log(m_old)+2, log(s_old)+2),
                             control = list(fnscale = +1))$par
  m_new = exp(lognew[1])
  s_new = exp(lognew[2])

  if( (1-pgamma(m_new*s_new, s_new, 1)) > 0.99999 | (is.na(1-pgamma(m_new*s_new, s_new, 1)))){
    m_new = m_old
    s_new = s_old
  }

  return(list(m = m_new, s = s_new))
}

gammacount.EM_notexact <- function(expert_old,
                                         tl, yl, yu, tu,
                                         exposure,
                                         z_e_obs, z_e_lat, k_e,
                                         penalty, pen_params) {

  m_old = expert_old$get_params()$m
  s_old = expert_old$get_params()$s

  lognew = stats::optim(par = c(log(m_old), log(s_old)),
                        fn = gammacount_optim_params,
                        expert_old = expert_old,
                        tl = tl, yl = yl, yu = yu, tu = tu,
                        exposure = exposure,
                        z_e_obs = z_e_obs, z_e_lat = z_e_lat, k_e = k_e,
                        penalty = penalty,
                        pen_params = pen_params,
                        # lower =  max(log(k_old)-2.0, 0.0),
                        # upper = log(k_old)+2.0,
                        method = "Nelder-Mead", #  "L-BFGS-B",
                        # lower = c(log(m_old)-2, log(s_old)-2),
                        # upper = c(log(m_old)+2, log(s_old)+2),
                        control = list(fnscale = +1))$par
  m_new = exp(lognew[1])
  s_new = exp(lognew[2])

  if( (1-pgamma(m_new*s_new, s_new, 1)) > 0.99999 | (is.na(1-pgamma(m_new*s_new, s_new, 1)))){
    m_new = m_old
    s_new = s_old
  }

  return(list(m = m_new, s = s_new))
}

gammacount_optim_params_exact <- function(lognew, expert_old,
                                          ye, exposure, z_e_obs,
                                          penalty = TRUE, pen_params = c(1, Inf, 1, Inf)) {
  m_tmp = exp(lognew[1])
  s_tmp = exp(lognew[2])
  expert_tmp = ExpertFunction$new("gammacount",
                                  list(m = m_tmp, s = s_tmp),
                                  pen_params)

  densY_e_obs = rep(0.0, length(ye))
  for(i in 1:length(ye)){
    expert_expo = expert_tmp$exposurize(exposure[i])
    densY_e_obs[i] = expert_expo$ll_exact(ye[i])
  }

  densY_e_obs[is.na(densY_e_obs)] = 0.0

  term_zkz_Y = z_e_obs * densY_e_obs

  obj = sum(term_zkz_Y)
  p = ifelse(penalty,
             (pen_params[1]-1)*lognew[1] - m_tmp/pen_params[2] +
               (pen_params[3]-1)*lognew[2] - s_tmp/pen_params[4],
             0)
  return((obj+p)*(-1))
}

gammacount_optim_params <- function(lognew,
                                     expert_old,
                                     tl, yl, yu, tu,
                                     exposure,
                                     z_e_obs, z_e_lat, k_e,
                                     y_e_obs, y_e_lat,
                                     penalty = TRUE, pen_params = c(1, Inf, 1, Inf)) {

  m_tmp = exp(lognew[1])
  s_tmp = exp(lognew[2])
  expert_tmp = ExpertFunction$new("gammacount",
                                  list(m = m_tmp, s = s_tmp),
                                  pen_params)

  densY_e_obs = rep(0.0, length(yl))
  densY_e_lat = rep(0.0, length(yl))
  for(i in 1:length(yl)){
    # expert_expo = expert_old$exposurize(exposure[i])
    expert_expo_new = expert_tmp$exposurize(exposure[i])
    # result_set = expert_expo$ll_not_exact(tl[i], yl[i], yu[i], tu[i])
    # expert_ll = result_set[["expert_ll"]]
    # expert_tn_bar = result_set[["expert_tn_bar"]]
    result_set_new = expert_expo_new$ll_not_exact(tl[i], yl[i], yu[i], tu[i])
    expert_ll_new = result_set_new[["expert_ll"]]
    expert_tn_bar_new = result_set_new[["expert_tn_bar"]]
    densY_e_obs[i] = expert_ll_new # exp(-expert_ll) * expert_ll_new
    densY_e_lat[i] = expert_tn_bar_new # exp(-expert_tn_bar) * expert_tn_bar_new
  }

  densY_e_obs[is.na(densY_e_obs)] = 0.0
  densY_e_lat[is.na(densY_e_lat)] = 0.0
  densY_e_obs[is.infinite(densY_e_obs)] = 0.0
  densY_e_lat[is.infinite(densY_e_lat)] = 0.0

  term_zkz_Y = XAPlusYZB(z_e_obs, densY_e_obs, z_e_lat, k_e, densY_e_lat)

  obj = sum(term_zkz_Y)
  p = ifelse(penalty,
             (pen_params[1]-1)*lognew[1] - m_tmp/pen_params[2] +
               (pen_params[3]-1)*lognew[2] - s_tmp/pen_params[4],
             0)
  return((obj+p)*(-1))

}


######################################################################################
# Register the gammacount at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################
