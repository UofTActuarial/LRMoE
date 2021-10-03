# inversegaussian Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the distribution
######################################################################################
inversegaussian.params_init <- function(y){
  pos_index = which(y > 0)
  mu = mean(y[pos_index])
  variance = var(y[pos_index])
  mu_init = mu
  lambda_init = mu^3/variance

  return( list(mean = mu_init, shape = lambda_init) )
}

inversegaussian.exposurize <- function(params, exposure){
  return( params )
}

inversegaussian.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
inversegaussian.expert_ll_exact <- function(y, params){
  return( inversegaussian.logpdf(params = params, x = y) )
}

inversegaussian.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # Check dim tl == yl == tu == yu
  # Check value

  # Find indexes of unequal yl & yu: Not exact observations, but censored
  expert.ll = expert.tn = expert.tn.bar = rep(-Inf, length(yu))
  censor.idx = (yl!=yu)
  no.trunc.idx = (tl==tu)

  # Expert LL calculation
  ######################################################################################
  prob.log.yu = inversegaussian.logcdf(params, yu[censor.idx])
  prob.log.yl = inversegaussian.logcdf(params, yl[censor.idx])

  # Compute loglikelihood for expert j, first for y
  # likelihood of censored interval: some easy algebra
  expert.ll[censor.idx] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl)
  # exact likelihood
  expert.ll[!censor.idx] = inversegaussian.logpdf(params, yu[!censor.idx])
  ######################################################################################

  # Expert TN Calculation
  ######################################################################################
  # Compute loglikelihood for expert j, then for truncation limits t
  prob.log.tu=inversegaussian.logcdf(params, tu)
  prob.log.tl=inversegaussian.logcdf(params, tl)

  # Normalizing factor for truncation limits, in log
  expert.tn = prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)

  # Deal with no truncation case
  expert.tn[no.trunc.idx] = inversegaussian.logpdf(params, tu[no.trunc.idx])

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

inversegaussian.penalty <- function(params, penalty_params) {
  if(!length(penalty_params)) { penalty_params = c(1.0, Inf, 1.0, Inf) }
  p = penalty_params
  mu = params[["mean"]]
  shape = params[["shape"]]
  return( (p[1] - 1)*log(mu) - mu/p[2] + (p[3] - 1)*log(shape) - shape/p[4] )
}

inversegaussian.default_penalty <- function() {
  return(c(1.0, Inf, 1.0, Inf))
}

######################################################################################
# dinversegaussian, pinversegaussian, qinversegaussian and rinversegaussian implementations.
######################################################################################
inversegaussian.simulation <- function(params, n) {
  return( rinvgauss(n, mean = params[["mean"]], shape = params[["shape"]]) )
}

inversegaussian.mean <- function(params) {
  return( params[["mean"]] )
}

inversegaussian.variance <- function(params) {
  return( params[["mean"]]^3/params[["shape"]] )
}

inversegaussian.logpdf <- function(params, x) {
  return( dinvgauss(x, mean = params[["mean"]], shape = params[["shape"]], log = T) )
}

inversegaussian.pdf <- function(params, x) {
  return( dinvgauss(x, mean = params[["mean"]], shape = params[["shape"]], log = F) )
}

inversegaussian.logcdf <- function(params, q) {
  return( pinvgauss(q, mean = params[["mean"]], shape = params[["shape"]], log.p = T) )
}

inversegaussian.cdf <- function(params, q) {
  return( pinvgauss(q, mean = params[["mean"]], shape = params[["shape"]], log.p = F) )
}

inversegaussian.quantile <- function(params, p) {
  return( qinvgauss(p, mean = params[["mean"]], shape = params[["shape"]]) )
}

inversegaussian.lev <- function(params, u) {
  mu = params[["mean"]]
  lambda = params[["shape"]]

  z = (u - mu)/mu
  y = (u + mu)/mu
  result = u - mu*z*pnorm(q = z*(lambda/u)^0.5) - mu * y * exp(2*lambda/mu) * pnorm(-y*(lambda/u)*0.5)
}
######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
inversegaussian.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations
  Y_e_obs = ye
  invY_e_obs = 1/ye

  pos_idx = (ye!=0)
  term_zkz = z_e_obs[pos_idx]
  term_zkz_Y = z_e_obs[pos_idx] * Y_e_obs[pos_idx]
  term_zkz_invY = z_e_obs[pos_idx] * invY_e_obs[pos_idx]

  mean_new = sum(term_zkz_Y)/sum(term_zkz)
  shape_new = sum(term_zkz) / ( sum(term_zkz_Y)/(mean_new^2) - 2*sum(term_zkz)/mean_new + sum(term_zkz_invY) )

  return(list(mean = mean_new, shape = shape_new))
}

inversegaussian.EM_notexact <- function(expert_old,
                                        tl, yl, yu, tu,
                                        exposure,
                                        z_e_obs, z_e_lat, k_e,
                                        penalty, pen_params) {

  # Old parameters
  mean_old = expert_old$get_params()$mean
  shape_old = expert_old$get_params()$shape

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
  # logY_e_obs = rep(0, length(yl))
  # logY_e_lat = rep(0, length(yl))
  invY_e_obs = rep(0, length(yl))
  invY_e_lat = rep(0, length(yl))

  # Uncensored observations
  y_e_obs[!cencor_idx] = yl[!cencor_idx]
  # logY_e_obs[!cencor_idx] = log(yl[!cencor_idx])
  invY_e_obs[!cencor_idx] = 1/yl[!cencor_idx]

  # Unique upper/lower bounds for integration: yl+yu
  y_unique = unique(cbind(yl, yu), MARGIN=1)
  y_unique_length = nrow(y_unique)
  y_unique_match = match(data.frame(t(cbind(yl, yu))), data.frame(t(y_unique)))
  yl_unique = y_unique[,1]
  yu_unique = y_unique[,2]

  # Unique values of integration
  y_e_obs_unique = rep(0, length(y_unique_length))
  # logY_e_obs_unique = rep(0, length(y_unique_length))
  invY_e_obs_unique = rep(0, length(y_unique_length))

  y_e_obs_unique = intInvGaussYObs(mean_old, shape_old, log(yl_unique), log(yu_unique))
  # logY_e_obs_unique = intInvGaussLogYObs(mean_old, shape_old, log(yl_unique), log(yu_unique))
  invY_e_obs_unique = intInvGaussInvYObs(mean_old, shape_old, log(yl_unique), log(yu_unique))

  # Match to original observations
  y_e_obs_match = y_e_obs_unique[y_unique_match]
  # logY_e_obs_match = logY_e_obs_unique[y_unique_match]
  invY_e_obs_match = invY_e_obs_unique[y_unique_match]

  y_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) * y_e_obs_match[cencor_idx]
  # logY_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) * logY_e_obs_match[cencor_idx]
  invY_e_obs[cencor_idx] = exp(-expert_ll[cencor_idx]) * invY_e_obs_match[cencor_idx]


  # Unique upper/lower bounds for integration: tl+tu
  t_unique = unique(cbind(tl, tu), MARGIN=1)
  t_unique_length = nrow(t_unique)
  t_unique_match = match(data.frame(t(cbind(tl, tu))), data.frame(t(t_unique)))
  tl_unique = t_unique[,1]
  tu_unique = t_unique[,2]

  # Unique values of integration
  y_e_lat_unique = rep(0, length(t_unique_length))
  # logY_e_lat_unique = rep(0, length(t_unique_length))
  invY_e_lat_unique = rep(0, length(t_unique_length))

  y_e_lat_unique = intInvGaussYLat(mean_old, shape_old, log(tl_unique), log(tu_unique))
  # logY_e_lat_unique = intInvGaussLogYLat(mean_old, shape_old, log(tl_unique), log(tu_unique))
  invY_e_lat_unique = intInvGaussInvYLat(mean_old, shape_old, log(tl_unique), log(tu_unique))

  # Match to original observations
  y_e_lat_match = y_e_lat_unique[t_unique_match]
  # logY_e_lat_match = logY_e_lat_unique[t_unique_match]
  invY_e_lat_match = invY_e_lat_unique[t_unique_match]

  y_e_lat = exp(-expert_tn_bar) * y_e_lat_match
  # logY_e_lat = exp(-expert_tn_bar) * logY_e_lat_match
  invY_e_lat = exp(-expert_tn_bar) * invY_e_lat_match

  # Get rid of NaN's
  y_e_obs[is.na(y_e_obs)] = 0
  y_e_lat[is.na(y_e_lat)] = 0
  # logY_e_obs[is.na(logY_e_obs)] = 0
  # logY_e_lat[is.na(logY_e_lat)] = 0
  invY_e_obs[is.na(invY_e_obs)] = 0
  invY_e_lat[is.na(invY_e_lat)] = 0

  # Update parameters
  pos_idx = (yu!=0)

  mean_new = inversegaussian_mean_update(z_e_obs[pos_idx], z_e_lat[pos_idx], k_e[pos_idx],
                                         y_e_obs[pos_idx], y_e_lat[pos_idx])

  shape_new = inversegaussian_mean_to_shape(mean_new,
                                            z_e_obs[pos_idx], z_e_lat[pos_idx], k_e[pos_idx],
                                            y_e_obs[pos_idx], y_e_lat[pos_idx],
                                            invY_e_obs[pos_idx], invY_e_lat[pos_idx])

  return(list(mean = mean_new, shape = shape_new))
}


inversegaussian_mean_update <- function(z.obs.j, z.lat.j, k.e, y.e.obs.j, y.e.lat.j)
{
  sum.one = XPlusYZ(z.obs.j, z.lat.j, k.e)
  sum.two = XAPlusYZB(z.obs.j, y.e.obs.j, z.lat.j, k.e, y.e.lat.j)
  return( sum(sum.two)/sum(sum.one) )
}

inversegaussian_mean_to_shape <- function(mean.mu,
                                          z.obs.j, z.lat.j, k.e,
                                          y.e.obs.j, y.e.lat.j,
                                          y.inv.e.obs.j, y.inv.e.lat.j)
{
  sum.one = XPlusYZ(z.obs.j, z.lat.j, k.e)
  sum.two = XAPlusYZB(z.obs.j, y.e.obs.j, z.lat.j, k.e, y.e.lat.j)
  sum.three = XAPlusYZB(z.obs.j, y.inv.e.obs.j, z.lat.j, k.e, y.inv.e.lat.j)
  shape.lambda.temp = sum(sum.one) / ( (1/(mean.mu^2))*sum(sum.two) - (2/(mean.mu))*sum(sum.one) + sum(sum.three) )

  return(shape.lambda.temp)
}
######################################################################################
# Register in the ExpertLibrary Object at zzz.R
######################################################################################
