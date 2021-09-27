# gamma Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the gamma
######################################################################################
gamma.params_init <- function(y){
  # Estimate all the parameters that needed for further calculation
  pos_index = which(y > 0)
  mu = mean(y[pos_index])
  variance = var(y[pos_index])
  
  theta_init = variance/mu
  k_init = mu/theta_init
  
  return( list(shape = k_init, scale = theta_init) )
}

gamma.exposurize <- function(params, exposure){
  # Calculate the exposures
  return( params )
}

gamma.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}

######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
gamma.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  return( gamma.logpdf(params = params, x = y) )
}

gamma.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # Check dim tl == yl == tu == yu
  # Check value
  
  # Find indexes of unequal yl & yu: Not exact observations, but censored
  expert.ll = expert.tn = expert.tn.bar = rep(-Inf, length(yu))
  censor.idx = (yl!=yu)
  no.trunc.idx = (tl==tu)
  
  # Expert LL calculation
  ######################################################################################
  prob.log.yu = gamma.logcdf(params, yu[censor.idx])
  prob.log.yl = gamma.logcdf(params, yl[censor.idx])
  
  # Compute loglikelihood for expert j, first for y
  # likelihood of censored interval: some easy algebra
  expert.ll[censor.idx] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl)
  # exact likelihood
  expert.ll[!censor.idx] = gamma.logpdf(params, yu[!censor.idx])
  ######################################################################################
  
  # Expert TN Calculation
  ######################################################################################
  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu = gamma.logcdf(params, tu)
  prob.log.tl = gamma.logcdf(params, tl)
  
  # Normalizing factor for truncation limits, in log
  expert.tn = prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)
  
  # Deal with no truncation case
  expert.tn[no.trunc.idx] = gamma.logpdf(params, tu[no.trunc.idx])
  
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

gamma.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  if(!length(penalty_params)) { penalty_params = c(2.0, 10.0, 2.0, 10.0) }
  p = penalty_params
  d.k = params[["shape"]]
  d.theta = params[["scale"]]
  return( (p[1]-1)*log(d.k) - d.k/p[2] + (p[3]-1)*log(d.theta) - d.theta/p[4] )
}

######################################################################################
# dgamma, pgamma, qgamma and rgamma implementations.
######################################################################################
gamma.simulation <- function(params, n) {
  # simulated n points based on the gamma params
  return( rgamma(n, shape = params[["shape"]], scale = params[["scale"]]) )
}

gamma.mean <- function(params) {
  return( params[["shape"]] * params[["scale"]] )
}

gamma.variance <- function(params) {
  return( params[["shape"]] * params[["scale"]]^2 )
}

gamma.logpdf <- function(params, x) {
  return( dgamma(x, shape = params[["shape"]], scale = params[["scale"]], log = T) )
}

gamma.pdf <- function(params, x) {
  return( dgamma(x, shape = params[["shape"]], scale = params[["scale"]], log = F) )
}

gamma.logcdf <- function(params, q) {
  return( pgamma(q, shape = params[["shape"]], scale = params[["scale"]], log.p = T) )
}

gamma.cdf <- function(params, q) {
  return( pgamma(q, shape = params[["shape"]], scale = params[["scale"]], log.p = F) )
}

gamma.quantile <- function(params, p) {
  return( qgamma(p, shape = params[["shape"]], scale = params[["scale"]]) )
}

gamma.lev <- function(params, u) {
  k = params[["shape"]]
  theta = params[["scale"]]
  result = theta * k * gammainc(k+1, u/theta) + u*(1-gammainc(k, u/theta))
  return(result)
}
######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
gamma.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations

  k_old = expert_old$get_params()$shape
  theta_old = expert_old$get_params()$scale

  Y_e_obs = ye
  logY_e_obs = log(ye)

  pos_idx = (ye!=0)
  term_zkz = z_e_obs[pos_idx]
  term_zkz_Y = z_e_obs[pos_idx] * Y_e_obs[pos_idx]
  term_zkz_logY = z_e_obs[pos_idx] * logY_e_obs[pos_idx]

  logk_new = stats::optimise(f = gamma_optim_k,
                             sum_term_zkz = sum(term_zkz),
                             sum_term_zkzy = sum(term_zkz_Y),
                             sum_term_zkzlogy = sum(term_zkz_logY),
                             penalty = penalty, pen_params = pen_params,
                             lower =  max(log(k_old)-2.0, 0.0),
                             upper = log(k_old)+2.0,
                             tol = .Machine$double.eps^0.05,
                             maximum = FALSE)$minimum
  k_new = exp(logk_new)
  theta_new = gamma_k_to_theta(k_new, sum(term_zkz), sum(term_zkz_Y),
                               penalty = penalty, pen_params = pen_params)

  return(list(shape = k_new, scale = theta_new))
}

gamma.EM_notexact <- function(expert_old, 
                              tl, yl, yu, tu, 
                              exposure, 
                              z_e_obs, z_e_lat, k_e,
                              penalty, pen_params) {
  
  # Old parameters
  k_old = expert_old$get_params()$shape
  theta_old = expert_old$get_params()$scale
  
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
  y_e_obs = rep(0, length(yl))
  y_e_lat = rep(0, length(yl))
  logY_e_obs = rep(0, length(yl))
  logY_e_lat = rep(0, length(yl))
  
  # Uncensored observations
  y_e_obs[!cencor_idx] = yl[!cencor_idx]
  logY_e_obs[!cencor_idx] = log(yl[!cencor_idx])
  
  # Unique upper/lower bounds for integration: yl+yu
  y_unique = unique(cbind(yl, yu), MARGIN=1)
  y_unique_length = nrow(y_unique)
  y_unique_match = match(data.frame(t(cbind(yl, yu))), data.frame(t(y_unique)))
  yl_unique = y_unique[,1]
  yu_unique = y_unique[,2]
  
  # Unique values of integration
  y_e_obs_unique = rep(0, length(y_unique_length))
  logY_e_obs_unique = rep(0, length(y_unique_length))
  
  y_e_obs_unique = intGammaYObs(k_old, theta_old, yl_unique, yu_unique)
  logY_e_obs_unique = intGammaLogYObs(k_old, theta_old, log(yl_unique), log(yu_unique))
  
  # Match to original observations
  y_e_obs_match = y_e_obs_unique[y_unique_match]
  logY_e_obs_match = logY_e_obs_unique[y_unique_match]
  
  y_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) * y_e_obs_match[cencor_idx]
  logY_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) * logY_e_obs_match[cencor_idx]
  
  
  # Unique upper/lower bounds for integration: tl+tu
  t_unique = unique(cbind(tl, tu), MARGIN=1)
  t_unique_length = nrow(t_unique)
  t_unique_match = match(data.frame(t(cbind(tl, tu))), data.frame(t(t_unique)))
  tl_unique = t_unique[,1]
  tu_unique = t_unique[,2]
  
  # Unique values of integration
  y_e_lat_unique = rep(0, length(t_unique_length))
  logY_e_lat_unique = rep(0, length(t_unique_length))
  
  y_e_lat_unique = intGammaYLat(k_old, theta_old, tl_unique, tu_unique)
  logY_e_lat_unique = intGammaLogYLat(k_old, theta_old, log(tl_unique), log(tu_unique))
  
  # Match to original observations
  y_e_lat_match = y_e_lat_unique[t_unique_match]
  logY_e_lat_match = logY_e_lat_unique[t_unique_match]
  
  y_e_lat = exp(-expert_tn_bar) * y_e_lat_match
  logY_e_lat = exp(-expert_tn_bar) * logY_e_lat_match
  
  # Get rid of NaN's
  y_e_obs[is.na(y_e_obs)] = 0
  y_e_lat[is.na(y_e_lat)] = 0
  logY_e_obs[is.na(logY_e_obs)] = 0
  logY_e_lat[is.na(logY_e_lat)] = 0
  
  
  # Update parameters
  pos_idx = (yu!=0)
  
  term_zkz = XPlusYZ(z_e_obs[pos_idx], z_e_lat[pos_idx], k_e[pos_idx])
  term_zkz_Y = XAPlusYZB(z_e_obs[pos_idx], y_e_obs[pos_idx], z_e_lat[pos_idx], k_e[pos_idx], y_e_lat[pos_idx])
  term_zkz_logY = XAPlusYZB(z_e_obs[pos_idx], logY_e_obs[pos_idx], z_e_lat[pos_idx], k_e[pos_idx], logY_e_lat[pos_idx])
  
  logk_new = stats::optimise(f = gamma_optim_k,
                             sum_term_zkz = sum(term_zkz),
                             sum_term_zkzy = sum(term_zkz_Y),
                             sum_term_zkzlogy = sum(term_zkz_logY),
                             penalty = penalty, pen_params = pen_params,
                             lower =  max(log(k_old)-2.0, 0.0),
                             upper = log(k_old)+2.0,
                             tol = .Machine$double.eps^0.05,
                             maximum = FALSE)$minimum
  k_new = exp(logk_new)
  theta_new = gamma_k_to_theta(k_new, sum(term_zkz), sum(term_zkz_Y),
                               penalty = penalty, pen_params = pen_params)
  
  return(list(shape = k_new, scale = theta_new))
}


gamma_k_to_theta <- function(k, sum_term_zkz, sum_term_zkzy,
                              penalty = TRUE, pen_params = c(1, Inf, 1, Inf)) {
  if(penalty){
    hyper_1 = pen_params[3]
    hyper_2 = pen_params[4]
    if((hyper_1 == 1) & (hyper_2 == Inf)){
      return(sum_term_zkzy/(sum_term_zkz*k))
    }
    quad_1 = ( k * sum_term_zkz - (hyper_1 - 1) )^2
    quad_2 = (4/hyper_2) * sum_term_zkzy
    return((hyper_2/2) * ( (hyper_1-1) - k*sum_term_zkz + sqrt(quad_1 + quad_2) ))
  }else{
    return(sum_term_zkzy/(sum_term_zkz*k))
  }
}

gamma_optim_k <- function(logk, sum_term_zkz, sum_term_zkzy, sum_term_zkzlogy,
                          penalty = TRUE, pen_params = c(1, Inf, 1, Inf)) {
  k_tmp = exp(logk)
  theta_tmp = gamma_k_to_theta(k_tmp, sum_term_zkz, sum_term_zkzy,
                               penalty = penalty, pen_params = pen_params)
  obj = (k_tmp-1)*sum_term_zkzlogy - (1/theta_tmp)*sum_term_zkzy - (k_tmp*log(theta_tmp) + lgamma(k_tmp))*sum_term_zkz
  p = ifelse(penalty, (pen_params[1]-1)*log(k_tmp) - k_tmp/pen_params[2] + (pen_params[3]-1)*log(theta_tmp) - theta_tmp/pen_params[4] , 0.0)
  return((obj+p)*(-1))
}

intGammaYObs <- function(m, theta, lower, upper) {
  diff.dist.untrunc = pgamma(upper, shape = m+1, scale = theta, lower.tail = TRUE, log.p = FALSE) - 
    pgamma(lower, shape = m+1, scale = theta, lower.tail = TRUE, log.p = FALSE)
  return(diff.dist.untrunc*m*theta)
}

intGammaYLat <- function(m, theta, lower, upper) {
  diff.dist.untrunc = pgamma(upper, shape = m+1, scale = theta, lower.tail = TRUE, log.p = FALSE) - 
    pgamma(lower, shape = m+1, scale = theta, lower.tail = TRUE, log.p = FALSE)
  return((1-diff.dist.untrunc)*m*theta)
}

######################################################################################
# Register the gamma at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################