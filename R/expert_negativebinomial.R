# negativebinomial Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the negativebinomial
######################################################################################
negativebinomial.params_init <- function(y){
  # Calculate all the parameters that needed for further calculation
  mu = mean(y)
  variance = var(y)
  p_init = mu / variance
  n_init = mu * p_init/(1-p_init)
  return( list(n = n_init, p = p_init) )
}

negativebinomial.exposurize <- function(params, exposure){
  return( list(n = params[["n"]] * exposure, p = params[["p"]]) )
}

negativebinomial.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
negativebinomial.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  return( negativebinomial.logpdf(params = params, x = y) )
}

negativebinomial.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # Check dim tl == yl == tu == yu
  # Check value
  
  # Find indexes of unequal yl & yu: Not exact observations, but censored
  expert.ll = expert.tn = expert.tn.bar = rep(-Inf, length(yu))
  censor.idx = (yl!=yu)
  no.trunc.idx = (tl==tu)
  
  # Expert LL calculation
  ######################################################################################
  prob.log.yu = negativebinomial.logcdf(params, yu[censor.idx])
  prob.log.yl = negativebinomial.logcdf(params, ceiling(yl[censor.idx])-1)
  
  # Compute loglikelihood for expert j, first for y
  # likelihood of censored interval: some easy algebra
  expert.ll[censor.idx] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl)
  # exact likelihood
  expert.ll[!censor.idx] = negativebinomial.logpdf(params, yu[!censor.idx])
  ######################################################################################
  
  # Expert TN Calculation
  ######################################################################################
  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu = negativebinomial.logcdf(params, tu)
  prob.log.tl = negativebinomial.logcdf(params, ceiling(tl)-1)
  
  # Normalizing factor for truncation limits, in log
  expert.tn = prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)
  
  # Deal with no truncation case
  expert.tn[no.trunc.idx] = negativebinomial.logpdf(params, tu[no.trunc.idx])
  ######################################################################################
  
  # Expert TN Bar Calculation
  ######################################################################################
  expert.tn.bar[!no.trunc.idx] = log1mexp(-expert.tn[!no.trunc.idx])
  expert.tn.bar[no.trunc.idx] = log1mexp(-expert.tn[no.trunc.idx])
  ######################################################################################
  
  # Return values
  return( list(expert_ll = expert.ll, expert_tn = expert.tn, expert_tn_bar = expert.tn.bar) )
}

negativebinomial.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  if(!length(penalty_params)) { penalty_params = c(2.0, 10.0) }
  p = penalty_params
  d.n = params[["n"]]
  d.p = params[["p"]]
  return( (p[1]-1)*log(d.n) - d.n/p[2] )
}

######################################################################################
# dnegativebinomial, pnegativebinomial, qnegativebinomial and rnegativebinomial implementations.
######################################################################################
negativebinomial.simulation <- function(params, n) {
  # simulated n points based on the negativebinomial params
  return( rnbinom(n, size = params[["n"]], prob = params[["p"]]))
}

negativebinomial.mean <- function(params) {
  # Calculate the mean based on the params
  n = params[["n"]]
  p = params[["p"]]
  return( n*p/(1-p) )
}

negativebinomial.variance <- function(params) {
  # Calculate the variance based on the params
  n = params[["n"]]
  p = params[["p"]]
  return( n*p/(1-p)^2 )
}

negativebinomial.logpdf <- function(params, x) {
  # Return the log pdf based on the input x
  return( dnbinom(x, size = params[["n"]], prob = params[["p"]], log = T) )
}

negativebinomial.pdf <- function(params, x) {
  # Return the pdf based on the input x
  return( dnbinom(x, size = params[["n"]], prob = params[["p"]], log = F) )
}

negativebinomial.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  return( pnbinom(q, size = params[["n"]], prob = params[["p"]], log.p = T) )
}

negativebinomial.cdf <- function(params, q) {
  # return the cdf based on the input x
  return( pnbinom(q, size = params[["n"]], prob = params[["p"]], log.p = F) )
}

negativebinomial.quantile <- function(params, p) {
  # return the percentage points based on the value of p
  return( qnbinom(p, size = params[["n"]], prob = params[["p"]]) )
}

######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
# negativebinomial._EStep <- function() {
#   # Perform the E step
#   NULL
# }
# 
# negativebinomial._MStep <- function() {
#   # Perform the M step
#   NULL
# }
# 
# negativebinomial.compute_EM <- function() {
#   # Perform the EM optimization
#   NULL
# }

negativebinomial.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations
  
  n_old = expert_old$get_params()$n
  p_old = expert_old$get_params()$p
  
  logn_new = stats::optimise(f = negativebinomial_optim_n_exact,
                             ye = ye,
                             exposure = exposure,
                             z_e_obs = z_e_obs,
                             penalty = penalty,
                             pen_params = pen_params,
                             lower =  log(n_old)-0.5,
                             upper = log(n_old)+0.5,
                             tol = .Machine$double.eps^0.05,
                             maximum = FALSE)$minimum
  n_new = exp(logn_new)
  
  term_zkz = z_e_obs * n_new * exposure
  term_zkz_Y = z_e_obs * ye
  
  p_new = negativebinomial_n_to_p(sum(term_zkz), sum(term_zkz_Y),
                                  penalty = penalty, pen_params = pen_params)
  
  if((p_new^n_new > 0.999999) | (is.na(p_new^n_new))){
    n_new = n_old
    p_new = p_old
  }
  
  return(list(n = n_new, p = p_new))
}

negativebinomial_n_to_p <- function(sum_term_zkz, sum_term_zkzy,
                                    penalty = TRUE, pen_params = c(2, 10)) {
  return(1 - sum_term_zkzy/(sum_term_zkz + sum_term_zkzy))
}

negativebinomial_optim_n_exact <- function(logn, ye, exposure, z_e_obs,
                                           penalty = TRUE, pen_params = c(2, 10)) {
  
  n_tmp = exp(logn)
  
  Y_e_obs = ye
  logY_e_obs = lgamma(ye + n_tmp*exposure)
  
  term_zkz = z_e_obs * n_tmp * exposure
  term_zkz_Y = z_e_obs * Y_e_obs
  term_zkz_logY = z_e_obs * logY_e_obs
  
  p_tmp = negativebinomial_n_to_p(sum(term_zkz), sum(term_zkz_Y),
                                  penalty = penalty, pen_params = pen_params)
  
  obj = sum(term_zkz_logY) - sum(z_e_obs * lgamma(n_tmp * exposure)) +
    sum(term_zkz)*log(p_tmp) + sum(term_zkz_Y)*log(1-p_tmp)
  p = ifelse(penalty, (pen_params[1]-1)*logn - n_tmp/pen_params[2], 0.0)
  return((obj+p)*(-1))
  
}

negativebinomial.EM_notexact <- function(expert_old, 
                                tl, yl, yu, tu, 
                                exposure, 
                                z_e_obs, z_e_lat, k_e,
                                penalty, pen_params) {
  
  # Old parameters
  n_old = expert_old$get_params()$n
  p_old = expert_old$get_params()$p
  
  # Old loglikelihoods
  expert_ll = rep(-Inf, length(yl))
  expert_tn_bar = rep(-Inf, length(yl))
  # Additional E-step
  y_e_obs = rep(0, length(yl))
  y_e_lat = rep(0, length(yl))
  
  for(i in 1:length(yl)){
    expert_expo = expert_old$exposurize(exposure[i])
    n_expo = expert_expo$get_params()$n
    p_expo = expert_expo$get_params()$p
    result_set = expert_expo$ll_not_exact(tl[i], yl[i], yu[i], tu[i])
    expert_ll[i] = result_set[["expert_ll"]]
    expert_tn_bar[i] = result_set[["expert_tn_bar"]]
    y_e_obs[i] = exp(-expert_ll[i]) * sumNegativeBinomialYObs(n_expo, p_expo, yl[i], yu[i])
    y_e_lat[i] = exp(-expert_tn_bar[i]) * sumNegativeBinomialYLat(n_expo, p_expo, tl[i], tu[i])
  }
  
  y_e_obs[is.na(y_e_obs)] = 0
  y_e_lat[is.na(y_e_lat)] = 0
  y_e_obs[is.infinite(y_e_obs)] = 0
  y_e_lat[is.infinite(y_e_lat)] = 0
  
  logn_new = stats::optimise(f = negativebinomial_optim_n,
                             expert_old = expert_old,
                             tl = tl, yl = yl, yu = yu, tu = tu,
                             exposure = exposure,
                             z_e_obs = z_e_obs, z_e_lat = z_e_lat, k_e = k_e,
                             y_e_obs = y_e_obs, y_e_lat = y_e_lat,
                             penalty = penalty,
                             pen_params = pen_params,
                             lower =  log(n_old)-0.5,
                             upper = log(n_old)+0.5,
                             tol = .Machine$double.eps^0.05,
                             maximum = FALSE)$minimum
  n_new = exp(logn_new)
  
  term_zkz = XPlusYZ(z_e_obs, z_e_lat, k_e) * exposure * n_new
  term_zkz_Y = XAPlusYZB(z_e_obs, y_e_obs, z_e_lat, k_e, y_e_lat)
  
  p_new = negativebinomial_n_to_p(sum(term_zkz), sum(term_zkz_Y),
                                  penalty = penalty, pen_params = pen_params)
  
  if ((p_new^n_new > 0.999999) | (is.na(p_new^n_new))){
    n_new = n_old
    p_new = p_old
  }
 
  return(list(n = n_new, p = p_new))
}

negativebinomial_optim_n <- function(logn, 
                                     expert_old,
                                     tl, yl, yu, tu,
                                     exposure,
                                     z_e_obs, z_e_lat, k_e,
                                     y_e_obs, y_e_lat,
                                     penalty = TRUE, pen_params = c(2, 10)) {
  
  n_tmp = exp(logn)
  
  # Old loglikelihoods
  expert_ll = rep(-Inf, length(yl))
  expert_tn_bar = rep(-Inf, length(yl))
  # Additional E-step
  logY_e_obs = rep(0, length(yl))
  logY_e_lat = rep(0, length(yl))
  
  for(i in 1:length(yl)){
    expert_expo = expert_old$exposurize(exposure[i])
    n_expo = expert_expo$get_params()$n
    p_expo = expert_expo$get_params()$p
    result_set = expert_expo$ll_not_exact(tl[i], yl[i], yu[i], tu[i])
    expert_ll[i] = result_set[["expert_ll"]]
    expert_tn_bar[i] = result_set[["expert_tn_bar"]]
    logY_e_obs[i] = exp(-expert_ll[i]) * 
      sumNegativeBinomialLfacYObs(n_expo, p_expo, n_tmp*exposure[i], yl[i], yu[i])
    logY_e_lat[i] = exp(-expert_tn_bar[i]) * 
      sumNegativeBinomialLfacYLat(n_expo, p_expo, n_tmp*exposure[i], tl[i], tu[i])
  }
  
  logY_e_obs[is.na(logY_e_obs)] = 0
  logY_e_lat[is.na(logY_e_lat)] = 0
  logY_e_obs[is.infinite(logY_e_obs)] = 0
  logY_e_lat[is.infinite(logY_e_lat)] = 0
  
  
  
  term_zkz = XPlusYZ(z_e_obs, z_e_lat, k_e) * exposure * n_tmp
  term_zkz_Y = XAPlusYZB(z_e_obs, y_e_obs, z_e_lat, k_e, y_e_lat)
  term_zkz_logY = XAPlusYZB(z_e_obs, logY_e_obs, z_e_lat, k_e, logY_e_lat)
  
  p_tmp = negativebinomial_n_to_p(sum(term_zkz), sum(term_zkz_Y),
                                  penalty = penalty, pen_params = pen_params)
  
  obj = sum(term_zkz_logY) - sum(z_e_obs * lgamma(n_tmp * exposure)) +
    sum(term_zkz)*log(p_tmp) + sum(term_zkz_Y)*log(1-p_tmp)
  p = ifelse(penalty, (pen_params[1]-1)*logn - n_tmp/pen_params[2], 0.0)
  return((obj+p)*(-1))
  
}

######################################################################################
# Register the negativebinomial at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################