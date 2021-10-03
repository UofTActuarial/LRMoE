# zigamma Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the zigamma
######################################################################################
zigamma.params_init <- function(y){
  # Estimate all the parameters that needed for further calculation
  p0_init = sum(y==0)/sum(y>=0)
  y = y[which(y > 0)]
  result = gamma.params_init(y)
  return( c(result, list(p_zero = p0_init)) )
}

zigamma.exposurize <- function(params, exposure){
  # Calculate the exposures
  return( params )
}

zigamma.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}

######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
zigamma.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  p0 = params[["p_zero"]]
  return( ifelse(y == 0, log(p0), log(1-p0) + gamma.logpdf(params, x = y)) )
}

zigamma.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  return( faster_zi_result(tl, tu, yl, yu, params, "gamma") )
}

zigamma.penalty <- function(params, penalty_params) {
  return( gamma.penalty(params, penalty_params) )
}

zigamma.default_penalty <- function() {
  return(gamma.default_penalty())
}

######################################################################################
# dzigamma, pzigamma, qzigamma and rzigamma implementations.
######################################################################################
zigamma.simulation <- function(params, n) {
  # simulated n points based on the zigamma params
  p0 = params[["p_zero"]]
  return( (1 - rbinom(n,1,p0)) * gamma.simulation(params, n) )
}

zigamma.mean <- function(params) {
  p0 = params[["p_zero"]]
  return( (1-p0) * gamma.mean(params) )
}

zigamma.variance <- function(params) {
  p0 = params[["p_zero"]]
  return( (1-p0)*gamma.variance(params) + p0*(1-p0)*gamma.mean(params)^2 )
}

zigamma.logpdf <- function(params, x) {
  return( ifelse(params[["shape"]] < 1 & x <= 0, -Inf, gamma.logpdf(params, x)) )
}

zigamma.pdf <- function(params, x) {
  return( ifelse(params[["shape"]] < 1 & x <= 0, 0, gamma.pdf(params, x)) )
}

zigamma.logcdf <- function(params, q) {
  return( ifelse(params[["shape"]] < 1 & q <= 0, -Inf, gamma.logpdf(params, q)) )
}

zigamma.cdf <- function(params, q) {
  return( ifelse(params[["shape"]] < 1 & q <= 0, 0, gamma.cdf(params, q)) )
}

zigamma.quantile <- function(params, p) {
  p0 = params[["p_zero"]]
  return( ifelse(p0 >= 0, 0, gamma.quantile(params, p - p0)) )
}

zigamma.lev <- function(params, u) {
  p0 = params[["p_zero"]]
  result = (1-p0)*gamma.lev(params, u)
  return(result)
}
######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
zigamma.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  tmp_exp = ExpertFunction$new("gamma",
                               list(shape = expert_old$get_params()$shape,
                                    scale = expert_old$get_params()$scale),
                               pen_params)

  expert_ll_pos = tmp_exp$ll_exact(ye)

  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(ye, p_old, expert_ll_pos)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, 0.0, 0.0, 0.0)

  tmp_update = gamma.EM_exact(tmp_exp, ye, exposure, z_pos_e_obs,
                                penalty, pen_params)

  return(list(p_zero = p_new,
              shape = tmp_update$shape,
              scale = tmp_update$scale))
}

zigamma.EM_notexact <- function(expert_old,
                                tl, yl, yu, tu,
                                exposure,
                                z_e_obs, z_e_lat, k_e,
                                penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  tmp_exp = ExpertFunction$new("gamma",
                               list(shape = expert_old$get_params()$shape,
                                    scale = expert_old$get_params()$scale),
                               pen_params)

  expert_ll = rep(-Inf, length(yl))
  expert_tn = rep(-Inf, length(yl))
  expert_tn_bar = rep(-Inf, length(yl))

  for(i in 1:length(yl)){
    expert_expo = tmp_exp$exposurize(exposure[i])
    result_set = expert_expo$ll_not_exact(tl[i], yl[i], yu[i], tu[i])
    expert_ll[i] = result_set[["expert_ll"]]
    expert_tn[i] = result_set[["expert_tn"]]
    expert_tn_bar[i] = result_set[["expert_tn_bar"]]
  }

  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(yl, p_old, expert_ll)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  z_zero_e_lat = z_e_lat * EM_E_z_zero_lat(tl, p_old, expert_tn_bar)
  z_pos_e_lat = z_e_lat - z_zero_e_lat
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, z_zero_e_lat, z_pos_e_lat, k_e)

  tmp_update = gamma.EM_notexact(expert_old = tmp_exp,
                                 tl = tl, yl = yl, yu = yu, tu = tu,
                                 exposure = exposure,
                                 z_e_obs = z_pos_e_obs, z_e_lat = z_pos_e_lat,
                                 k_e = k_e,
                                 penalty = penalty, pen_params = pen_params)

  return(list(p_zero = p_new,
              shape = tmp_update$shape,
              scale = tmp_update$scale))
}
######################################################################################
# Register the zigamma at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################
