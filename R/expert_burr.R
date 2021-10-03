# burr Expert Function
######################################################################################
# Initialize the parameters of the distribuion
######################################################################################
burr.params_init <- function(y){
  # There is no clear way we could estimate burr distribution params from observations, so we would use MLE to estimate
  # Calculate all the parameters that needed for further calculation
  y = y[which(y>0)]
  init_obj <- function(logparams) {
    n = length(y)
    mean_tmp = exp(logparams[1])
    shape_tmp = exp(logparams[2])
    k_tmp = n / (sum(log1p((y / mean_tmp)^shape_tmp)))
    result = -1 * (n*log(mean_tmp*k_tmp) - n*(mean_tmp-1)*log(shape_tmp) + mean_tmp*sum(log(y)) -
                     (k_tmp+1)*sum(log1p((y / shape_tmp)^mean_tmp)) +
                     log(mean_tmp) - mean_tmp/100 + log(shape_tmp) - shape_tmp/100)  # to avoid spurious models
    return(result)
  }
  logparams_init = optim(par = c(0,0), init_obj)

  c_init = exp(logparams_init$par[1])
  lambda_init = exp(logparams_init$par[2])
  k_init = length(y) / (sum(log1p((y / lambda_init)^c_init)))

  return( list(shape1 = k_init, shape2 = c_init, scale = lambda_init) )
}

burr.exposurize <- function(params, exposure){
  return( params )
}

burr.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}

######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
burr.expert_ll_exact <- function(y, params){
  return( burr.logpdf(params = params, x = y) )
}

burr.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # Check dim tl == yl == tu == yu
  # Check value

  # Find indexes of unequal yl & yu: Not exact observations, but censored
  expert.ll = expert.tn = expert.tn.bar = rep(-Inf, length(yu))
  censor.idx = (yl!=yu)
  no.trunc.idx = (tl==tu)

  # Expert LL calculation
  ######################################################################################
  prob.log.yu = burr.logcdf(params, yu[censor.idx])
  prob.log.yl = burr.logcdf(params, yl[censor.idx])

  # Compute loglikelihood for expert j, first for y
  # likelihood of censored interval: some easy algebra
  expert.ll[censor.idx] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl)
  # exact likelihood
  expert.ll[!censor.idx] = burr.logpdf(params, yu[!censor.idx])
  ######################################################################################

  ######################################################################################
  # Deal with numerical underflow: prob.log.yu and prob.log.yl can both be -Inf
  # NA.idx = which(is.na(expert.ll[,j]))
  expert.ll[which(is.na(expert.ll))] = -Inf
  ######################################################################################

  # Expert TN Calculation
  ######################################################################################
  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu = burr.logcdf(params, tu)
  prob.log.tl = burr.logcdf(params, tl)

  # Normalizing factor for truncation limits, in log
  expert.tn = prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)

  ###################################################################
  # Deal with numerical underflow: prob.log.tu and prob.log.tl can both be -Inf
  # NA.idx = which(is.na(expert.tn[,j]))
  expert.tn[which(is.na(expert.tn))] = -Inf
  ###################################################################

  # Deal with no truncation case
  expert.tn[no.trunc.idx] = burr.logpdf(params, tu[no.trunc.idx])

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

burr.penalty <- function(params, penalty_params) {
  if(!length(penalty_params)) { penalty_params = c(2.0, 10.0, 2.0, 10.0, 2.0, 10.0) }
  c = params[["shape1"]]
  k = params[["shape2"]]
  s = params[["scale"]]
  p = penalty_params

  result = (p[1] - 1) * log(k) - k/p[2] + (p[3] - 1)*log(c) - c/p[4] + (p[5] -1)*log(s) - s/p[6]
  return( result )
}

burr.default_penalty <- function() {
  return(c(2.0, 10.0, 2.0, 10.0, 2.0, 10.0))
}

######################################################################################
# ddistribution, pdistribution, qdistribution and rdistribution implementations.
######################################################################################
burr.simulation <- function(params, n) {
  # simulated n points based on the distribution params
  return( rburr(n, shape1 = params[["shape1"]], shape2 = params[["shape2"]], scale = params[["scale"]]) )
}

burr.mean <- function(params) {
  a = params[["shape1"]]
  b = params[["shape2"]]
  return( b*beta(b - 1/a, b + 1/a) )
}

burr.variance <- function(params) {
  a = params[["shape1"]]
  b = params[["shape2"]]
  variance = b*beta(b - 1/a, b + 1/a) * b*beta(b - 2/a, b + 2/a)
  return( variance )
}

burr.logpdf <- function(params, x) {
  return( dburr(x, shape1 = params[["shape1"]], shape2 = params[["shape2"]], scale = params[["scale"]], log = T) )
}

burr.pdf <- function(params, x) {
  return( dburr(x, shape1 = params[["shape1"]], shape2 = params[["shape2"]], scale = params[["scale"]], log = F) )
}

burr.logcdf <- function(params, q) {
  return( pburr(q, shape1 = params[["shape1"]], shape2 = params[["shape2"]], scale = params[["scale"]], log.p = T) )
}

burr.cdf <- function(params, q) {
  return( pburr(q, shape1 = params[["shape1"]], shape2 = params[["shape2"]], scale = params[["scale"]], log.p = F) )
}

burr.quantile <- function(params, p) {
  return( qburr(p, shape1 = params[["shape1"]], shape2 = params[["shape2"]], scale = params[["scale"]]) )
}

burr.lev <- function(params, u) {
  k = params[["shape1"]]
  c = params[["shape2"]]
  lambda = params[["scale"]]

  uu = 1/(1+(u/lambda)^c)
  result = lambda*gamma(1+1/c)*gamma(k-1/c) / gamma(k) * beta_inc(1+1/c, k-1/c, 1-uu) + u * uu^k
  return(result)
}
######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
burr.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations

  c_old = expert_old$get_params()$shape2
  lambda_old = expert_old$get_params()$scale

  logY_e_obs = log(ye)

  pos_idx = (ye!=0)
  term_zkz = z_e_obs[pos_idx]
  term_zkz_logY = z_e_obs[pos_idx] * logY_e_obs[pos_idx]

  logparams_new = stats::optim(par = c(log(c_old), log(lambda_old)),
                             fn = burr_optim_params_exact,
                             ye = ye[pos_idx],
                             z_e_obs = z_e_obs[pos_idx],
                             sum_term_zkzlogy = sum(term_zkz_logY),
                             penalty = penalty,
                             pen_params = pen_params,
                             # lower =  max(log(k_old)-2.0, 0.0),
                             # upper = log(k_old)+2.0,
                             method = "L-BFGS-B",
                             lower = c(max(log(c_old)-2.0, 0.0), log(lambda_old)-2),
                             upper = c(log(c_old)+2, log(lambda_old)+2),
                             control = list(fnscale = +1))$par

  c_new = exp(logparams_new[1])
  lambda_new = exp(logparams_new[2])

  lpow_e_obs = log1p((ye/lambda_new)^c_new)
  term_zkz = z_e_obs[pos_idx]
  term_zkzlpow = z_e_obs[pos_idx] * lpow_e_obs[pos_idx]

  k_new = burr_lpow_to_k_exact(sum(term_zkz), sum(term_zkzlpow),
                               penalty = penalty, pen_params = pen_params)

  return(list(shape1 = k_new, shape2 = c_new, scale = lambda_new))
}

burr.EM_notexact <- function(expert_old,
                              tl, yl, yu, tu,
                              exposure,
                              z_e_obs, z_e_lat, k_e,
                              penalty, pen_params) {

  # Old parameters
  k_old = expert_old$get_params()$shape1
  c_old = expert_old$get_params()$shape2
  lambda_old = expert_old$get_params()$scale

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

  # Uncensored observations
  logY_e_obs[!cencor_idx] = log(yl[!cencor_idx])

  # Unique upper/lower bounds for integration: yl+yu
  y_unique = unique(cbind(yl, yu), MARGIN=1)
  y_unique_length = nrow(y_unique)
  y_unique_match = match(data.frame(t(cbind(yl, yu))), data.frame(t(y_unique)))
  yl_unique = y_unique[,1]
  yu_unique = y_unique[,2]

  # Unique values of integration
  logY_e_obs_unique = rep(0, length(y_unique_length))

  logY_e_obs_unique = intBurrLogYObs(k_old, c_old, lambda_old,
                                     log(yl_unique), log(yu_unique))

  # Match to original observations
  logY_e_obs_match = logY_e_obs_unique[y_unique_match]

  logY_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) * logY_e_obs_match[cencor_idx]


  # Unique upper/lower bounds for integration: tl+tu
  t_unique = unique(cbind(tl, tu), MARGIN=1)
  t_unique_length = nrow(t_unique)
  t_unique_match = match(data.frame(t(cbind(tl, tu))), data.frame(t(t_unique)))
  tl_unique = t_unique[,1]
  tu_unique = t_unique[,2]

  # Unique values of integration
  logY_e_lat_unique = rep(0, length(t_unique_length))

  logY_e_lat_unique = intBurrLogYLat(k_old, c_old, lambda_old,
                                     log(tl_unique), log(tu_unique))

  # Match to original observations
  logY_e_lat_match = logY_e_lat_unique[t_unique_match]

  logY_e_lat = exp(-expert_tn_bar) * logY_e_lat_match

  # Get rid of NaN's
  logY_e_obs[is.na(logY_e_obs)] = 0
  logY_e_lat[is.na(logY_e_lat)] = 0
  logY_e_obs[is.infinite(logY_e_obs)] = 0
  logY_e_lat[is.infinite(logY_e_lat)] = 0


  # Update parameters
  pos_idx = (yu!=0)

  term_zkzlogy = XAPlusYZB(z_e_obs[pos_idx],
                           logY_e_obs[pos_idx],
                           z_e_lat[pos_idx], k_e[pos_idx],
                           logY_e_lat[pos_idx])

  logparams_new = stats::optim(par = c(log(c_old), log(lambda_old)),
                               fn = burr_optim_params,
                               d_old = expert_old,
                               tl = tl[pos_idx], yl = yl[pos_idx],
                               yu = yu[pos_idx], tu = tu[pos_idx],
                               expert_ll = expert_ll[pos_idx], expert_tn_bar = expert_tn_bar[pos_idx],
                               z_e_obs = z_e_obs[pos_idx], z_e_lat = z_e_lat[pos_idx], k_e = k_e[pos_idx],
                               sum_term_zkzlogy = sum(term_zkzlogy),
                               penalty = penalty,
                               pen_params = pen_params,
                               # lower =  max(log(k_old)-2.0, 0.0),
                               # upper = log(k_old)+2.0,
                               method = "L-BFGS-B",
                               lower = c(max(log(c_old)-2.0, 0.0), log(lambda_old)-2),
                               upper = c(log(c_old)+2, log(lambda_old)+2),
                               control = list(fnscale = +1))$par


  c_new = exp(logparams_new[1])
  lambda_new = exp(logparams_new[2])

  # Get new k
  k_old = expert_old$get_params()$shape1
  c_old = expert_old$get_params()$shape2
  lambda_old = expert_old$get_params()$scale

  cencor_idx = (yl!=yu)
  lpowY_e_obs = rep(0, length(yl))
  lpowY_e_lat = rep(0, length(yl))

  # Uncensored observations
  lpowY_e_obs[!cencor_idx] = log1p((yu[!cencor_idx]/lambda_new)^c_new)

  # Unique upper/lower bounds for integration: yl+yu
  y_unique = unique(cbind(yl, yu), MARGIN=1)
  y_unique_length = nrow(y_unique)
  y_unique_match = match(data.frame(t(cbind(yl, yu))), data.frame(t(y_unique)))
  yl_unique = y_unique[,1]
  yu_unique = y_unique[,2]

  # Unique values of integration
  lpowY_e_obs_unique = rep(0, length(y_unique_length))

  lpowY_e_obs_unique = intBurrPolYObs(k_old, c_old, lambda_old,
                                      c_new, lambda_new, log(yl_unique), log(yu_unique))

  # Match to original observations
  lpowY_e_obs_match = lpowY_e_obs_unique[y_unique_match]

  lpowY_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) * lpowY_e_obs_match[cencor_idx]

  # Unique upper/lower bounds for integration: tl+tu
  t_unique = unique(cbind(tl, tu), MARGIN=1)
  t_unique_length = nrow(t_unique)
  t_unique_match = match(data.frame(t(cbind(tl, tu))), data.frame(t(t_unique)))
  tl_unique = t_unique[,1]
  tu_unique = t_unique[,2]

  # Unique values of integration
  lpowY_e_lat_unique = rep(0, length(t_unique_length))

  lpowY_e_lat_unique = intBurrPolYObs(k_old, c_old, lambda_old,
                                      c_new, lambda_new, log(tl_unique), log(tu_unique))

  # Match to original observations
  lpowY_e_lat_match = lpowY_e_lat_unique[t_unique_match]

  lpowY_e_lat = exp(-expert_tn_bar) * lpowY_e_lat_match

  # Get rid of NA's
  lpowY_e_obs[is.na(lpowY_e_obs)] = 0
  lpowY_e_lat[is.na(lpowY_e_lat)] = 0
  lpowY_e_obs[is.infinite(lpowY_e_obs)] = 0
  lpowY_e_lat[is.infinite(lpowY_e_lat)] = 0

  term_zkz = XPlusYZ(z_e_obs, z_e_lat, k_e)
  term_zkzlpow = XAPlusYZB(z_e_obs, lpowY_e_obs, z_e_lat, k_e, lpowY_e_lat)

  k_new = burr_lpow_to_k_exact(sum(term_zkz), sum(term_zkzlpow),
                               penalty = penalty, pen_params = pen_params)


  return(list(shape1 = k_new, shape2 = c_new, scale = lambda_new))
}


burr_lpow_to_k_exact <- function(sum_term_zkz, sum_term_zkzlpow,
                                 penalty = TRUE, pen_params = c(1, Inf, 1, Inf, 1, Inf)) {
  if(penalty){
    return( (sum_term_zkz+(pen_params[1]-1)) / (sum_term_zkzlpow + 1/pen_params[2]) )
  }else{
    return( sum_term_zkz / sum_term_zkzlpow )
  }
}

burr_optim_params_exact <- function(lognew, ye, z_e_obs, sum_term_zkzlogy,
                                    penalty = TRUE, pen_params = c(1, Inf, 1, Inf, 1, Inf)) {
  c_tmp = exp(lognew[1])
  lambda_tmp = exp(lognew[2])

  lpow_e_obs = log1p((ye/lambda_tmp)^c_tmp)

  term_zkz = z_e_obs
  term_zkzlpow = z_e_obs * lpow_e_obs

  k_tmp = burr_lpow_to_k_exact(sum(term_zkz), sum(term_zkzlpow),
                               penalty = penalty, pen_params = pen_params)

  obj = sum(term_zkz)*(log(k_tmp)+lognew[1]-c_tmp*lognew[2]) +
    sum_term_zkzlogy*c_tmp - (k_tmp+1)*sum(term_zkzlpow)
  p = ifelse(penalty, (pen_params[1]-1)*log(k_tmp) - k_tmp/pen_params[2] +
               (pen_params[3]-1)*lognew[1] - c_tmp/pen_params[4] +
               (pen_params[5]-1)*lognew[2] - lambda_tmp/pen_params[6], 0)

  return((obj+p)*(-1))
}



burr_optim_params <- function(lognew, d_old,
                              tl, yl, yu, tu,
                              expert_ll, expert_tn_bar,
                              z_e_obs, z_e_lat, k_e,
                              sum_term_zkzlogy,
                              penalty = TRUE, pen_params = c(1, Inf, 1, Inf, 1, Inf)) {

  k_old = d_old$get_params()$shape1
  c_old = d_old$get_params()$shape2
  lambda_old = d_old$get_params()$scale

  c_tmp = exp(lognew[1])
  lambda_tmp = exp(lognew[2])

  cencor_idx = (yl!=yu)
  lpowY_e_obs = rep(0, length(yl))
  lpowY_e_lat = rep(0, length(yl))

  # Uncensored observations
  lpowY_e_obs[!cencor_idx] = log1p((yu[!cencor_idx]/lambda_tmp)^c_tmp)

  # Unique upper/lower bounds for integration: yl+yu
  y_unique = unique(cbind(yl, yu), MARGIN=1)
  y_unique_length = nrow(y_unique)
  y_unique_match = match(data.frame(t(cbind(yl, yu))), data.frame(t(y_unique)))
  yl_unique = y_unique[,1]
  yu_unique = y_unique[,2]

  # Unique values of integration
  lpowY_e_obs_unique = rep(0, length(y_unique_length))

  lpowY_e_obs_unique = intBurrPolYObs(k_old, c_old, lambda_old,
                                     c_tmp, lambda_tmp, log(yl_unique), log(yu_unique))

  # Match to original observations
  lpowY_e_obs_match = lpowY_e_obs_unique[y_unique_match]

  lpowY_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) * lpowY_e_obs_match[cencor_idx]

  # Unique upper/lower bounds for integration: tl+tu
  t_unique = unique(cbind(tl, tu), MARGIN=1)
  t_unique_length = nrow(t_unique)
  t_unique_match = match(data.frame(t(cbind(tl, tu))), data.frame(t(t_unique)))
  tl_unique = t_unique[,1]
  tu_unique = t_unique[,2]

  # Unique values of integration
  lpowY_e_lat_unique = rep(0, length(t_unique_length))

  lpowY_e_lat_unique = intBurrPolYObs(k_old, c_old, lambda_old,
                                      c_tmp, lambda_tmp, log(tl_unique), log(tu_unique))

  # Match to original observations
  lpowY_e_lat_match = lpowY_e_lat_unique[t_unique_match]

  lpowY_e_lat = exp(-expert_tn_bar) * lpowY_e_lat_match

  # Get rid of NA's
  lpowY_e_obs[is.na(lpowY_e_obs)] = 0
  lpowY_e_lat[is.na(lpowY_e_lat)] = 0
  lpowY_e_obs[is.infinite(lpowY_e_obs)] = 0
  lpowY_e_lat[is.infinite(lpowY_e_lat)] = 0

  term_zkz = XPlusYZ(z_e_obs, z_e_lat, k_e)
  term_zkzlpow = XAPlusYZB(z_e_obs, lpowY_e_obs, z_e_lat, k_e, lpowY_e_lat)

  k_tmp = burr_lpow_to_k_exact(sum(term_zkz), sum(term_zkzlpow),
                               penalty = penalty, pen_params = pen_params)

  obj = sum(term_zkz)*(log(k_tmp)+lognew[1]-c_tmp*lognew[2]) +
    sum_term_zkzlogy*c_tmp - (k_tmp+1)*sum(term_zkzlpow)
  p = ifelse(penalty, (pen_params[1]-1)*log(k_tmp) - k_tmp/pen_params[2] +
               (pen_params[3]-1)*lognew[1] - c_tmp/pen_params[4] +
               (pen_params[5]-1)*lognew[2] - lambda_tmp/pen_params[6], 0)

  return((obj+p)*(-1))

}



######################################################################################
# Register in the ExpertLibrary Object
######################################################################################
