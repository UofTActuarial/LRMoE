# weibull Draft
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the distribuion
######################################################################################
weibull.params_init <- function(y){
  # Calculate all the parameters that needed for further calculation
  y = y[which(y > 0)]
  init_obj <- function(logk) {
    n = length(y)
    k_tmp = exp(logk)
    theta_tmp = ( sum(y^k_tmp) / n)^(1/k_tmp)
    return(-1 * (n*logk - n*log(theta_tmp) + k_tmp*sum(log(y)) - n*(k_tmp-1)*log(theta_tmp) - 1/(theta_tmp^k_tmp)*sum(y^k_tmp) ))
  }
  
  logk_init = optim(par = c(1), init_obj, method = "L-BFGS-B")
  k_init = exp(logk_init$par[1])
  theta_init = ( sum(y^k_init) / length(y) )^(1/k_init)
  
  return( list(shape = k_init, scale = theta_init) )
}

weibull.exposurize <- function(params, exposure){
  # Calculate the exposures
  return( params )
}

weibull.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
weibull.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  return( weibull.logpdf(params = params, x = y) )
}

weibull.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # Check dim tl == yl == tu == yu
  # Check value
  
  # Find indexes of unequal yl & yu: Not exact observations, but censored
  expert.ll = expert.tn = expert.tn.bar = rep(-Inf, length(yu))
  censor.idx = (yl!=yu)
  no.trunc.idx = (tl==tu)
  
  # Expert LL calculation
  ######################################################################################
  prob.log.yu = weibull.logcdf(params, yu[censor.idx])
  prob.log.yl = weibull.logcdf(params, yl[censor.idx])
  
  # Compute loglikelihood for expert j, first for y
  # likelihood of censored interval: some easy algebra
  expert.ll[censor.idx] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl)
  # exact likelihood
  expert.ll[!censor.idx] = weibull.logpdf(params, yu[!censor.idx])
  ######################################################################################
  
  ######################################################################################
  # Deal with numerical underflow: prob.log.yu and prob.log.yl can both be -Inf
  expert.ll[which(is.na(expert.ll))] = -Inf
  ######################################################################################
  
  # Expert TN Calculation
  ######################################################################################
  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu = weibull.logcdf(params, tu)
  prob.log.tl = weibull.logcdf(params, tl)
  
  # Normalizing factor for truncation limits, in log
  expert.tn= prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)
  
  # Deal with no truncation case
  expert.tn[no.trunc.idx] = weibull.logpdf(params, tu[no.trunc.idx])
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

weibull.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  if(!length(penalty_params)) { penalty_params = c(2.0, 10.0, 1.0, Inf)}
  p = penalty_params
  k = params[["shape"]]
  theta = params[["scale"]]
  
  return( (p[1] - 1)*log(k) -k/p[2] + (p[3] - 1)*log(theta) - theta/p[4] )
}

######################################################################################
# dweibull, pweibull, qweibull and rweibull implementations.
######################################################################################
weibull.simulation <- function(params, n) {
  # simulated n points based on the weibull params
  return( rweibull(n, shape = params[["shape"]], scale = params[["scale"]]) )
}

weibull.mean <- function(params) {
  # Calculate the mean based on the params
  k = params[["shape"]]
  theta = params[["scale"]]
  return( theta*gamma(1 + 1/k) )
}

weibull.variance <- function(params) {
  # Calculate the variance based on the params
  k = params[["shape"]]
  theta = params[["scale"]]
  return( theta^2 * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2) )
}

weibull.logpdf <- function(params, x) {
  # Return the log pdf based on the input x
  return( dweibull(x, shape = params[["shape"]], scale = params[["scale"]], log = T) )
}

weibull.pdf <- function(params, x) {
  # Return the pdf based on the input x
  return( dweibull(x, shape = params[["shape"]], scale = params[["scale"]], log = F) )
}

weibull.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  return( pweibull(q, shape = params[["shape"]], scale = params[["scale"]], log.p = T) )
}

weibull.cdf <- function(params, q) {
  # return the cdf based on the input x
  return( pweibull(q, shape = params[["shape"]], scale = params[["scale"]], log.p = F) )
}

weibull.quantile <- function(params, p) {
  # return the percentage points based on the value of p
  return( qweibull(p, shape = params[["shape"]], scale = params[["scale"]]) )
}

weibull.lev <- function(params, u) {
  k = params[["shape"]]
  theta = params[["scale"]]
  
  result = theta * gammainc(1/k+1, (u/theta)^k) * gamma(1/k+1) + u*exp(-(u/theta)^k)
  return(result)
}

######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
weibull.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations
  
  k_old = expert_old$get_params()$shape
  theta_old = expert_old$get_params()$scale
  
  Y_e_obs = ye
  logY_e_obs = log(ye)
  
  pos_idx = (ye!=0)
  
  logk_new = stats::optimise(f = weibull_optim_k_exact,
                             ye = ye[pos_idx],
                             z_e_obs = z_e_obs[pos_idx],
                             logY_e_obs = logY_e_obs[pos_idx],
                             penalty = penalty,
                             pen_params = pen_params,
                             lower =  max(log(k_old)-2.0, 0.0),
                             upper = log(k_old)+2.0,
                             tol = .Machine$double.eps^0.05,
                             maximum = FALSE)$minimum
  k_new = exp(logk_new)
  
  powY_e_obs = ye^(k_new)
  
  term_zkz = z_e_obs
  term_zkz_powY = z_e_obs * powY_e_obs
  
  theta_new = weibull_k_to_lambda(k_new, sum(term_zkz), sum(term_zkz_powY),
                               penalty = penalty, pen_params = pen_params)
  
  return(list(shape = k_new, scale = theta_new))
}

weibull.EM_notexact <- function(expert_old, 
                              tl, yl, yu, tu, 
                              exposure, 
                              z_e_obs, z_e_lat, k_e,
                              penalty, pen_params) {
  
  # Old parameters
  k_old = expert_old$get_params()$shape
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
  
  logY_e_obs_unique = intWeibullLogYObs(k_old, lambda_old, (yl_unique), (yu_unique))
  
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
  
  logY_e_lat_unique = intWeibullLogYLat(k_old, lambda_old, (tl_unique), (tu_unique))
  
  # Match to original observations
  logY_e_lat_match = logY_e_lat_unique[t_unique_match]
  
  logY_e_lat = exp(-expert_tn_bar) * logY_e_lat_match
  
  # Get rid of NaN's
  logY_e_obs[is.na(logY_e_obs)] = 0
  logY_e_lat[is.na(logY_e_lat)] = 0
  
  
  # Update parameters
  pos_idx = (yu!=0)
  
  logk_new = stats::optimise(f = weibull_optim_k,
                             d_old = expert_old,
                             tl = tl[pos_idx], yl = yl[pos_idx], yu = yu[pos_idx], tu = tu[pos_idx],
                             expert_ll = expert_ll[pos_idx], expert_tn_bar = expert_tn_bar[pos_idx],
                             z_e_obs = z_e_obs[pos_idx], z_e_lat = z_e_lat[pos_idx], k_e = k_e[pos_idx],
                             logY_e_obs = logY_e_obs[pos_idx], logY_e_lat = logY_e_lat[pos_idx],
                             lower =  max(log(k_old)-2.0, 0.0),
                             upper = log(k_old)+2.0,
                             tol = .Machine$double.eps^0.05,
                             maximum = FALSE)$minimum
  k_new = exp(logk_new)
  
  
  censor_idx = (yl!=yu)
  powY_e_obs = rep(0, length(yl))
  powY_e_lat = rep(0, length(yl))
  
  # Untruncated and uncensored case
  powY_e_obs[!censor_idx] = (yu[!censor_idx])^(k_new)
  
  # Untrencated and uncensored case
  powY_e_obs[censor_idx] = exp(-expert_ll[censor_idx]) *
    lambda_old^(k_new) *
    (gammainc((k_old+k_new)/k_old, (yl[censor_idx]/lambda_old)^(k_old)) -
       gammainc((k_old+k_new)/k_old, (yu[censor_idx]/lambda_old)^(k_old)) )
  
  # Truncated case
  powY_e_lat = exp(-expert_tn_bar) *
    lambda_old^(k_new) *
    (gammainc((k_old+k_new)/k_old, (0 /lambda_old)^(k_old)) -
       gammainc((k_old+k_new)/k_old, (tl/lambda_old)^(k_old)) +
       gammainc((k_old+k_new)/k_old, (tu/lambda_old)^(k_old)) -
       gammainc((k_old+k_new)/k_old, (Inf/lambda_old)^(k_old)))
  
  powY_e_obs[is.na(powY_e_obs)] = 0
  powY_e_lat[is.na(powY_e_lat)] = 0
  
  term_zkz = XPlusYZ(z_e_obs, z_e_lat, k_e)
  term_zkz_logY = XAPlusYZB(z_e_obs, logY_e_obs, z_e_lat, k_e, logY_e_lat)
  term_zkz_powY = XAPlusYZB(z_e_obs, powY_e_obs, z_e_lat, k_e, powY_e_lat)
  
  lambda_new = weibull_k_to_lambda(k_new, sum(term_zkz), sum(term_zkz_powY),
                                  penalty = penalty, pen_params = pen_params)
  
  return(list(shape = k_new, scale = lambda_new))
}

weibull_k_to_lambda <- function(k, sum_term_zkz, sum_term_zkzpowy,
                                penalty = TRUE, pen_params = c(1, Inf, 1, Inf)) {
  return((sum_term_zkzpowy/sum_term_zkz)^(1/k))
}

weibull_optim_k_exact <- function(logk, ye, z_e_obs, logY_e_obs,
                                  penalty = TRUE, pen_params = c(1, Inf, 1, Inf)) {
  k_tmp = exp(logk)
  
  powY_e_obs = ye^(k_tmp)
  
  term_zkz = z_e_obs
  term_zkz_logY = z_e_obs * logY_e_obs
  term_zkz_powY = z_e_obs * powY_e_obs
  
  theta_tmp = weibull_k_to_lambda(k_tmp, sum(term_zkz), sum(term_zkz_powY),
                                  penalty = penalty, pen_params = pen_params)
  
  obj = logk*sum(term_zkz) - k_tmp*log(theta_tmp)*sum(term_zkz) +
    (k_tmp-1)*sum(term_zkz_logY) - (theta_tmp^(-k_tmp))*sum(term_zkz_powY)
  p = ifelse(penalty, (pen_params[1]-1)*logk - k_tmp/pen_params[2] +
               (pen_params[3]-1)*log(theta_tmp) - theta_tmp/pen_params[4], 0)
  
  return((obj+p)*(-1))
}

weibull_optim_k <- function(logk, d_old,
                            tl, yl, yu, tu,
                            expert_ll, expert_tn_bar,
                            z_e_obs, z_e_lat, k_e,
                            logY_e_obs, logY_e_lat,
                            penalty = TRUE, pen_params = c(1, Inf, 1, Inf)) {
  
  k_tmp = exp(logk)
  
  k_old = d_old$get_params()$shape
  lambda_old = d_old$get_params()$scale
  
  censor_idx = (yl!=yu)
  powY_e_obs = rep(0, length(yl))
  powY_e_lat = rep(0, length(yl))
  
  # Untruncated and uncensored case
  powY_e_obs[!censor_idx] = (yu[!censor_idx])^(k_tmp)
  
  # Untrencated and uncensored case
  powY_e_obs[censor_idx] = exp(-expert_ll[censor_idx]) *
    lambda_old^(k_tmp) *
    (gammainc((k_old+k_tmp)/k_old, (yl[censor_idx]/lambda_old)^(k_old)) -
       gammainc((k_old+k_tmp)/k_old, (yu[censor_idx]/lambda_old)^(k_old)) )
  
  # Truncated case
  powY_e_lat = exp(-expert_tn_bar) *
    lambda_old^(k_tmp) *
    (gammainc((k_old+k_tmp)/k_old, (0 /lambda_old)^(k_old)) -
       gammainc((k_old+k_tmp)/k_old, (tl/lambda_old)^(k_old)) +
       gammainc((k_old+k_tmp)/k_old, (tu/lambda_old)^(k_old)) -
       gammainc((k_old+k_tmp)/k_old, (Inf/lambda_old)^(k_old)))
  
  powY_e_obs[is.na(powY_e_obs)] = 0
  powY_e_lat[is.na(powY_e_lat)] = 0
  
  term_zkz = XPlusYZ(z_e_obs, z_e_lat, k_e)
  term_zkz_logY = XAPlusYZB(z_e_obs, logY_e_obs, z_e_lat, k_e, logY_e_lat)
  term_zkz_powY = XAPlusYZB(z_e_obs, powY_e_obs, z_e_lat, k_e, powY_e_lat)
  
  theta_tmp = weibull_k_to_lambda(k_tmp, sum(term_zkz), sum(term_zkz_powY),
                                  penalty = penalty, pen_params = pen_params)
  
  obj = logk*sum(term_zkz) - k_tmp*log(theta_tmp)*sum(term_zkz) +
    (k_tmp-1)*sum(term_zkz_logY) - (theta_tmp^(-k_tmp))*sum(term_zkz_powY)
  p = ifelse(penalty, (pen_params[1]-1)*logk - k_tmp/pen_params[2] +
               (pen_params[3]-1)*log(theta_tmp) - theta_tmp/pen_params[4], 0)
  
  return((obj+p)*(-1))
}
######################################################################################
# Register in the ExpertLibrary Object
######################################################################################