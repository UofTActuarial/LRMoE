# zinegativebinomial Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the zinegativebinomial
######################################################################################
zinegativebinomial.params_init <- function(y){
  # Calculate all the parameters that needed for further calculation
  pos_idx = which(y > 0)
  mu = mean(y[pos_idx])
  variance = var(y[pos_idx])
  p_init = mu / variance
  n_init = mu * p_init/(1-p_init)
  p0_init = 1 - mean(y)/mu
  return( list(n = n_init, p = p_init, p_zero = p0_init) )
}

zinegativebinomial.exposurize <- function(params, exposure){
  result = negativebinomial.exposurize(params, exposure)
  p_zero = params[["p_zero"]]
  return( c(result, list(p_zero = p_zero)) )
}

zinegativebinomial.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
zinegativebinomial.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  p0 = params[["p_zero"]]
  return( ifelse(y == 0,
                 log(p0 + (1 - p0)*negativebinomial.pdf(params, x = y)),
                 log(1 - p0) + negativebinomial.logpdf(params, x = y)
                )
        )
}

zinegativebinomial.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # If some observations are not exact, calculate its corresponding likelihood
  return( faster_zi_result(tl, tu, yl, yu, params, "negativebinomial") )
}

zinegativebinomial.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  return( negativebinomial.penalty(params, penalty_params) )
}

zinegativebinomial.default_penalty <- function() {
  return(negativebinomial.default_penalty())
}

######################################################################################
# dzinegativebinomial, pzinegativebinomial, qzinegativebinomial and rzinegativebinomial implementations.
######################################################################################
zinegativebinomial.simulation <- function(params, n) {
  # simulated n points based on the zinegativebinomial params
  p0 = params[["p_zero"]]
  return( (1 - rbinom(n,1,p0))*negativebinomial.simulation(params, n) )
}

zinegativebinomial.mean <- function(params) {
  # Calculate the mean based on the params
  return( (1 - params[["p_zero"]])*negativebinomial.mean(params) )
}

zinegativebinomial.variance <- function(params) {
  # Calculate the variance based on the params
  p0 = params[["p_zero"]]
  return( (1-p0)*negativebinomial.variance(params) + p0*(1-p0)*negativebinomial.mean(params)^2 )
}

zinegativebinomial.logpdf <- function(params, x) {
  # Return the log pdf based on the input x
  return( ifelse(is.infinite(x), -Inf, negativebinomial.logpdf(params, x)) )
}

zinegativebinomial.pdf <- function(params, x) {
  # Return the pdf based on the input x
  return( ifelse(is.infinite(x), 0, negativebinomial.pdf(params, x)) )
}

zinegativebinomial.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  return( ifelse(is.infinite(q), 0, negativebinomial.logcdf(params, q)) )
}

zinegativebinomial.cdf <- function(params, q) {
  # return the cdf based on the input x
  return( ifelse(is.infinite(q), 1, negativebinomial.cdf(params, q)) )
}

zinegativebinomial.quantile <- function(params, p) {
  # return the percentage points based on the value of p
  p0 = params[["p_zero"]]
  return( ifelse(p <= p0, 0, negativebinomial.quantile(params, p - p0)) )
}

######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
# zinegativebinomial._EStep <- function() {
#   # Perform the E step
#   NULL
# }
#
# zinegativebinomial._MStep <- function() {
#   # Perform the M step
#   NULL
# }
#
# zinegativebinomial.compute_EM <- function() {
#   # Perform the EM optimization
#   NULL
# }

zinegativebinomial.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  tmp_exp = ExpertFunction$new("negativebinomial",
                               list(n = expert_old$get_params()$n,
                                    p = expert_old$get_params()$p),
                               pen_params)

  expert_ll_pos = tmp_exp$ll_exact(ye)

  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(ye, p_old, expert_ll_pos)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, 0.0, 0.0, 0.0)

  tmp_update = negativebinomial.EM_exact(tmp_exp, ye, exposure, z_pos_e_obs,
                              penalty, pen_params)

  return(list(p_zero = p_new,
              n = tmp_update$n,
              p = tmp_update$p))
}

zinegativebinomial.EM_notexact <- function(expert_old,
                                  tl, yl, yu, tu,
                                  exposure,
                                  z_e_obs, z_e_lat, k_e,
                                  penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  if(p_old > 0.999999){
    return(list(p_zero = p_old,
                n = expert_old$n,
                p = expert_old$p))
  }

  tmp_exp = ExpertFunction$new("negativebinomial",
                               list(n = expert_old$get_params()$n,
                                    p = expert_old$get_params()$p),
                               pen_params)

  expert_ll = rep(-Inf, length(yl))
  expert_tn_bar = rep(-Inf, length(yl))

  for(i in 1:length(yl)){
    expert_expo = tmp_exp$exposurize(exposure[i])
    result_set = expert_expo$ll_not_exact(tl[i], yl[i], yu[i], tu[i])
    expert_ll[i] = result_set[["expert_ll"]]
    expert_tn_bar[i] = result_set[["expert_tn_bar"]]
  }

  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(yl, p_old, expert_ll)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  z_zero_e_lat = z_e_lat * EM_E_z_zero_lat(tl, p_old, expert_tn_bar)
  z_pos_e_lat = z_e_lat - z_zero_e_lat
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, z_zero_e_lat, z_pos_e_lat, k_e)

  tmp_update = negativebinomial.EM_notexact(expert_old = tmp_exp,
                                   tl = tl, yl = yl, yu = yu, tu = tu,
                                   exposure = exposure,
                                   z_e_obs = z_pos_e_obs, z_e_lat = z_pos_e_lat,
                                   k_e = k_e,
                                   penalty = penalty, pen_params = pen_params)

  return(list(p_zero = p_new,
              n = tmp_update$n,
              p = tmp_update$p))
}
######################################################################################
# Register the zinegativebinomial at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################
