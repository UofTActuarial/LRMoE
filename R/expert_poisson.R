######################################################################################
# Initialize the parameters of the distribuion
######################################################################################
poisson.params_init <- function(y) {
  lambda = mean(y)
  return( list(lambda = lambda) )
}

poisson.exposurize <- function(params, exposure) {
  lambda = params[["lambda"]]
  return( list(lambda = lambda * exposure) )
}

poisson.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
poisson.expert_ll_exact <- function(y, params) {
  return( poisson.logpdf(params = params, x = y) )
}

poisson.expert_ll_not_exact <- function(tl, yl, yu, tu, params){
  # Check dim tl == yl == tu == yu
  # Check value

  expert.ll = expert.tn = expert.tn.bar = rep(-Inf, length(yu))
  censor.idx = (yl!=yu)
  no.trunc.idx = (tl==tu)

  # Expert LL calculation
  ######################################################################################
  prob.log.yu = poisson.logcdf(params, yu[censor.idx])
  prob.log.yl = poisson.logcdf(params, ceiling(yl[censor.idx])-1)

  expert.ll[censor.idx] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl)
  expert.ll[!censor.idx] = poisson.logpdf(params, yu[!censor.idx])
  ######################################################################################

  # Expert TN Calculation
  ######################################################################################
  prob.log.tu = poisson.logcdf(params, tu)
  prob.log.tl = poisson.logcdf(params, ceiling(tl)-1)
  expert.tn = prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)
  expert.tn[no.trunc.idx] = poisson.logpdf(params, tu[no.trunc.idx])
  ######################################################################################

  # Expert TN Bar Calculation
  ######################################################################################
  expert.tn.bar[!no.trunc.idx] = log1mexp(-expert.tn[!no.trunc.idx])
  expert.tn.bar[no.trunc.idx] = log1mexp(-expert.tn[no.trunc.idx])
  ######################################################################################

  return( list(expert_ll = expert.ll, expert_tn = expert.tn, expert_tn_bar = expert.tn.bar) )
}

poisson.penalty <- function(params, penalty_params) {
  if(!length(penalty_params)) { penalty_params = c(2.0, 1.0) }
  penalty = (penalty_params[1] - 1)*log(params[["lambda"]]) - params[["lambda"]]/penalty_params[2]
  return(penalty)
}

poisson.default_penalty <- function() {
  return(c(2.0, 1.0))
}

######################################################################################
# ddistribution, pdistribution, qdistribution and rdistribution implementations.
######################################################################################
poisson.simulation <- function(params, n) {
  # simulated n points based on the distribution params
  return( rpois(n = n, lambda = params[["lambda"]]) )
}

poisson.mean <- function(params) {
  return( params[["lambda"]] )
}

poisson.variance <- function(params) {
  return( params[["lambda"]] )
}

poisson.logpdf <- function(params, x) {
  return( dpois(x = x, lambda = params[["lambda"]], log = T) )
}

poisson.pdf <- function(params, x) {
  return( dpois(x = x, lambda = params[["lambda"]], log = F) )
}

poisson.logcdf <- function(params, q) {
  return( ppois(q = q, lambda = params[["lambda"]], log.p = T) )
}

poisson.cdf <- function(params, q) {
  return( ppois(q = q, lambda = params[["lambda"]], log.p = F) )
}

poisson.quantile <- function(params, p) {
  return( qpois(p = p, lambda = params[["lambda"]]) )
}

######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
poisson.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  lambda_old = expert_old$get_params()$lambda
  term_zkz = z_e_obs * exposure
  term_zkz_Y = z_e_obs * ye
  lambda_new = ifelse(penalty,
                      (sum(term_zkz_Y) - (pen_params[1]-1)) / (sum(term_zkz) + 1/pen_params[2]),
                      sum(term_zkz_Y) / sum(term_zkz)
                      )

  if(lambda_new < 1e-4){
    lambda_new = lambda_old
  }
  return(list(lambda = lambda_new))
}

poisson.EM_notexact <- function(expert_old,
                              tl, yl, yu, tu,
                              exposure,
                              z_e_obs, z_e_lat, k_e,
                              penalty, pen_params) {

  # Old parameters
  lambda_old = expert_old$get_params()$lambda

  # Old loglikelihoods
  expert_ll = rep(-Inf, length(yl))
  expert_tn_bar = rep(-Inf, length(yl))
  # Additional E-step
  cencor_idx = (yl!=yu)
  y_e_obs = rep(0, length(yl))
  y_e_lat = rep(0, length(yl))

  for(i in 1:length(yl)){
    expert_expo = expert_old$exposurize(exposure[i])
    lambda_expo = expert_expo$get_params()$lambda
    result_set = expert_expo$ll_not_exact(tl[i], yl[i], yu[i], tu[i])
    expert_ll[i] = result_set[["expert_ll"]]
    expert_tn_bar[i] = result_set[["expert_tn_bar"]]
    y_e_obs[i] = ifelse(yl[i]==yu[i],
                        yl[i],
                        exp(-expert_ll[i]) * sumPoissonYObs(lambda_expo, yl[i], yu[i]))
    y_e_lat[i] = ifelse(tl[i]==tu[i],
                        lambda_expo - tl[i],
                        exp(-expert_tn_bar[i]) * sumPoissonYLat(lambda_expo, tl[i], tu[i]))
  }

  y_e_obs[is.na(y_e_obs)] = 0
  y_e_lat[is.na(y_e_lat)] = 0

  term_zkz = XPlusYZ(z_e_obs, z_e_lat, k_e) * exposure
  term_zkz_Y = XAPlusYZB(z_e_obs, y_e_obs, z_e_lat, k_e, y_e_lat)

  lambda_new = ifelse(penalty,
                      (sum(term_zkz_Y) - (pen_params[1]-1)) / (sum(term_zkz) + 1/pen_params[2]),
                      sum(term_zkz_Y) / sum(term_zkz)
                      )

  return(list(lambda = lambda_new))
}
######################################################################################
# Register the distribution at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################
