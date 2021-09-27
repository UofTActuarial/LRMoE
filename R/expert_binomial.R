# Binomial Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the binomial
######################################################################################
binomial.params_init <- function(y){
  # Estimate all the parameters that needed for further calculation
  n_init = max(y) + 2
  p_init = mean(y)/length(y)
  return( list(n = n_init, p = p_init) )
}

binomial.exposurize <- function(params, exposure){
  # Calculate the exposures
  return( params )
}

binomial.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
binomial.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  return( binomial.logpdf(params = params, x = y) )
}

binomial.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # Check dim tl == yl == tu == yu
  # Check value
  
  # Find indexes of unequal yl & yu: Not exact observations, but censored
  expert.ll = expert.tn = expert.tn.bar = rep(-Inf, length(yu))
  censor.idx = (yl!=yu)
  no.trunc.idx = (tl==tu)
  
  # Expert LL calculation
  ######################################################################################
  prob.log.yu = binomial.logcdf(params, yu[censor.idx])
  prob.log.yl = binomial.logcdf(params, ceiling(yl[censor.idx])-1)
  
  # Compute loglikelihood for expert j, first for y
  # likelihood of censored interval: some easy algebra
  expert.ll[censor.idx] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl)
  # exact likelihood
  expert.ll[!censor.idx] = binomial.logpdf(params, yu[!censor.idx])
  ######################################################################################
  
  # Expert TN Calculation
  ######################################################################################
  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu=binomial.logcdf(params, tu)
  prob.log.tl=binomial.logcdf(params, ceiling(tl)-1)
  
  # Normalizing factor for truncation limits, in log
  expert.tn = prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)
  
  # Deal with no truncation case
  expert.tn[no.trunc.idx] = binomial.logpdf(params, tu[no.trunc.idx])
  ######################################################################################
  
  # Expert TN Bar Calculation
  ######################################################################################
  expert.tn.bar[!no.trunc.idx] = log1mexp(-expert.tn[!no.trunc.idx])
  expert.tn.bar[no.trunc.idx] = log1mexp(-expert.tn[no.trunc.idx])
  ######################################################################################
  
  # Return values
  return( list(expert_ll = expert.ll, expert_tn = expert.tn, expert_tn_bar = expert.tn.bar) )
}

binomial.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  return(0)
  # stop("There are no approprate penalty parameters applied, please check `binomial.penalty()`")
}

######################################################################################
# dbinomial, pbinomial, qbinomial and rbinomial implementations.
######################################################################################
binomial.simulation <- function(params, n) {
  return( rbinom(n, size = params[["n"]], prob = params[["p"]]) )
}

binomial.mean <- function(params) {
  # Calculate the mean based on the params
  return( params[["n"]] * params[["p"]] )
}

binomial.variance <- function(params) {
  # Calculate the variance based on the params
  return( params[["n"]] * params[["p"]] * (1 - params[["p"]]) )
}

binomial.logpdf <- function(params, x) {
  # Return the log pdf based on the input x
  return( dbinom(x, size = params[["n"]], prob = params[["p"]], log = TRUE) )
}

binomial.pdf <- function(params, x) {
  # Return the pdf based on the input x
  return( dbinom(x, size = params[["n"]], prob = params[["p"]]) )
}

binomial.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  return( pbinom(q, size = params[["n"]], prob = params[["p"]], log.p = TRUE) )
}

binomial.cdf <- function(params, q) {
  # return the cdf based on the input x
  return( pbinom(q, size = params[["n"]], prob = params[["p"]]) )
}

binomial.quantile <- function(params, p) {
  # return the percentage points based on the value of p
  return( qbinom(p, size = params[["n"]], prob = params[["p"]]) ) 
}

######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
# binomial._EStep <- function() {
#   # Perform the E step
#   NULL
# }
# 
# binomial._MStep <- function() {
#   # Perform the M step
#   NULL
# }
# 
# binomial.compute_EM <- function() {
#   # Perform the EM optimization
#   NULL
# }

binomial.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations
  
  n_old = expert_old$get_params()$n
  p_old = expert_old$get_params()$p
  
  term_zkz_Y = z_e_obs * ye
  term_zkz_n_Y = z_e_obs * (n_old - ye)
  
  p_new = sum(term_zkz_Y) / (sum(term_zkz_Y) + sum(term_zkz_n_Y))
  n_new = n_old
  
  if(((1-p_new)^n_new > 0.999999) | (is.na((1-p_new)^n_new))){
    n_new = n_old
    p_new = p_old
  }
  
  return(list(n = n_new, p = p_new))
}

binomial.EM_notexact <- function(expert_old, 
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
    y_e_obs[i] = exp(-expert_ll[i]) * sumBinomialYObs(n_expo, p_expo, yl[i], yu[i])
    y_e_lat[i] = exp(-expert_tn_bar[i]) * (n_expo*p_expo - sumBinomialYObs(n_expo, p_expo, tl[i], tu[i]))
  }
  
  y_e_obs[is.na(y_e_obs)] = 0
  y_e_lat[is.na(y_e_lat)] = 0
  y_e_obs[is.infinite(y_e_obs)] = 0
  y_e_lat[is.infinite(y_e_lat)] = 0
  
  term_zkz_Y = XAPlusYZB(z_e_obs, y_e_obs, z_e_lat, k_e, y_e_lat)
  term_zkz_n_Y = XAPlusYZB(z_e_obs, n_old - y_e_obs, z_e_lat, k_e, n_old - y_e_lat)
  
  p_new = sum(term_zkz_Y) / (sum(term_zkz_Y) + sum(term_zkz_n_Y))
  n_new = n_old
  
  if(((1-p_new)^n_new > 0.999999) | (is.na((1-p_new)^n_new))){
    n_new = n_old
    p_new = p_old
  }
  
  return(list(n = n_new, p = p_new))
}
######################################################################################
# Register the binomial at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################