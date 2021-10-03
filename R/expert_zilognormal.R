# zilognormal Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the zilognormal
######################################################################################
zilognormal.params_init <- function(y){
  # Estimate all the parameters that needed for further calculation
  p0 = sum(y == 0) / sum(y >= 0)
  y = y[which(y > 0)]
  result = lognormal.params_init(y)
  return( c(result, list(p_zero = p0)) )
}

zilognormal.exposurize <- function(params, exposure){
  # Calculate the exposures
  return( params )
}

zilognormal.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
zilognormal.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  return( ifelse(y == 0, log(params[["p_zero"]]), log(1 - params[["p_zero"]]) + lognormal.logpdf(params, x=y)) )
}

zilognormal.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # If some observations are not exact, calculate its corresponding likelihood
  return( faster_zi_result(tl, tu, yl, yu, params, "lognormal") )
}

zilognormal.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  return( lognormal.penalty(params, penalty_params) )
}

zilognormal.default_penalty <- function() {
  return(lognormal.default_penalty())
}

######################################################################################
# dzilognormal, pzilognormal, qzilognormal and rzilognormal implementations.
######################################################################################
zilognormal.simulation <- function(params, n) {
  p0 = params[["p_zero"]]
  return( (1 - rbinom(n,1,p0)) * lognormal.simulation(params, n) )
}

zilognormal.mean <- function(params) {
  # Calculate the mean based on the params
  return( (1 - params[["p_zero"]])*lognormal.mean(params) )
}

zilognormal.variance <- function(params) {
  # Calculate the variance based on the params
  p0 = params[["p_zero"]]
  return( (1-p0)*lognormal.variance(params) + p0*(1-p0)*lognormal.mean(params)^2 )
}

zilognormal.logpdf <- function(params, x) {
  # Return the log pdf based on the input x
  return( lognormal.logpdf(params, x) )
}

zilognormal.pdf <- function(params, x) {
  # Return the pdf based on the input x
  return( lognormal.pdf(params, x) )
}

zilognormal.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  return( lognormal.logcdf(params, q) )
}

zilognormal.cdf <- function(params, q) {
  # return the cdf based on the input x
  return( lognormal.cdf(params, q) )
}

zilognormal.quantile <- function(params, p) {
  # return the percentage points based on the value of p
  p0 = params[["p_zero"]]
  return( ifelse(p0 >= p, 0, lognormal.quantile(params, p - p0)) )
}

zilognormal.lev <- function(params, u) {
  p0 = params[["p_zero"]]
  result = (1-p0)*lognormal.lev(params, u)
  return(result)
}
######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
# zilognormal._EStep <- function() {
#   # Perform the E step
#   NULL
# }
#
# zilognormal._MStep <- function() {
#   # Perform the M step
#   NULL
# }
#
# zilognormal.compute_EM <- function() {
#   # Perform the EM optimization
#   NULL
# }

zilognormal.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  tmp_exp = ExpertFunction$new("lognormal",
                               list(meanlog = expert_old$get_params()$meanlog,
                                    sdlog = expert_old$get_params()$sdlog),
                               pen_params)

  expert_ll_pos = tmp_exp$ll_exact(ye)

  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(ye, p_old, expert_ll_pos)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, 0.0, 0.0, 0.0)

  tmp_update = lognormal.EM_exact(tmp_exp, ye, exposure, z_pos_e_obs,
                              penalty, pen_params)

  return(list(p_zero = p_new,
              meanlog = tmp_update$meanlog,
              sdlog = tmp_update$sdlog))
}

zilognormal.EM_notexact <- function(expert_old,
                                tl, yl, yu, tu,
                                exposure,
                                z_e_obs, z_e_lat, k_e,
                                penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  tmp_exp = ExpertFunction$new("lognormal",
                               list(meanlog = expert_old$get_params()$meanlog,
                                    sdlog = expert_old$get_params()$sdlog),
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

  tmp_update = lognormal.EM_notexact(expert_old = tmp_exp,
                                 tl = tl, yl = yl, yu = yu, tu = tu,
                                 exposure = exposure,
                                 z_e_obs = z_pos_e_obs, z_e_lat = z_pos_e_lat,
                                 k_e = k_e,
                                 penalty = penalty, pen_params = pen_params)

  return(list(p_zero = p_new,
              meanlog = tmp_update$meanlog,
              sdlog = tmp_update$sdlog))
}
######################################################################################
# Register the zilognormal at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################
