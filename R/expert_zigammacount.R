# zigammacount Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the zigammacount
######################################################################################
zigammacount.params_init <- function(y){
  # Calculate all the parameters that needed for further calculation
  p0 = sum(y == 0)/sum(y >= 0)
  y = y[which(y > 0)]
  result = gammacount.params_init(y)
  return( c(result, list(p_zero = p0)) )
}

zigammacount.exposurize <- function(params, exposure){
  # Calculate the exposures
  return( list(p_zero = params[["p_zero"]], m = params[["m"]], s = params[["s"]]/exposure) )
}

zigammacount.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
zigammacount.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  p0 = params[["p_zero"]]
  return( ifelse(y == 0, 
                 log(p0 + (1-p0)*gammacount.pdf(params, x = y)),
                 log(1-p0) + gammacount.logpdf(params, x = y)
                 )
          )
}

zigammacount.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  return( faster_zi_result(tl, tu, yl, yu, params, "gammacount") )
}

zigammacount.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  return( gammacount.penalty(params, penalty_params) )
}

######################################################################################
# dzigammacount, pzigammacount, qzigammacount and rzigammacount implementations.
######################################################################################
zigammacount.simulation <- function(params, n) {
  p0 = params[["p_zero"]]
  return( (1 - rbinom(n, 1, p0))*gammacount.simulation(params, n) )
}

zigammacount.mean <- function(params) {
  return( (1 - params[["p_zero"]])*gammacount.mean(params) )
}

zigammacount.variance <- function(params) {
  # Calculate the variance based on the params
  p0 = params[["p_zero"]]
  return( (1 - p0)*gammacount.variance(params) + p0*(1-p0)*gammacount.mean(params)^2 )
}

zigammacount.logpdf <- function(params, x) {
  return( ifelse(is.infinite(x), 0, gammacount.logpdf(params, x)) )
}

zigammacount.pdf <- function(params, x) {
  # Return the pdf based on the input x
  return( ifelse(is.infinite(x), -Inf, gammacount.pdf(params, x)) )
}

zigammacount.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  return( ifelse(is.infinite(q), 0, gammacount.logcdf(params, q)) )
}

zigammacount.cdf <- function(params, q) {
  # return the cdf based on the input x
  return( ifelse(is.infinite(q), 1, gammacount.cdf(params, q)) )
}

zigammacount.quantile <- function(params, p) {
  return( ifelse(params[["p_zero"]] >= p, 0, gammacount.quantile(params, p - params[["p_zero"]])) )
}

######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
# zigammacount._EStep <- function() {
#   # Perform the E step
#   NULL
# }
# 
# zigammacount._MStep <- function() {
#   # Perform the M step
#   NULL
# }
# 
# zigammacount.compute_EM <- function() {
#   # Perform the EM optimization
#   NULL
# }

zigammacount.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations
  
  p_old = expert_old$get_params()$p_zero
  
  tmp_exp = ExpertFunction$new("gammacount", 
                               list(m = expert_old$get_params()$m,
                                    s = expert_old$get_params()$s),
                               pen_params)
  
  expert_ll_pos = tmp_exp$ll_exact(ye)
  
  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(ye, p_old, expert_ll_pos)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, 0.0, 0.0, 0.0)
  
  tmp_update = gammacount.EM_exact(tmp_exp, ye, exposure, z_pos_e_obs,
                                         penalty, pen_params)
  
  return(list(p_zero = p_new, 
              m = tmp_update$m, 
              s = tmp_update$s))
}

zigammacount.EM_notexact <- function(expert_old, 
                                           tl, yl, yu, tu, 
                                           exposure, 
                                           z_e_obs, z_e_lat, k_e,
                                           penalty, pen_params) {
  # Perform the EM optimization with exact observations
  
  p_old = expert_old$get_params()$p_zero
  
  if(p_old > 0.999999){
    return(list(p_zero = p_old, 
                m = expert_old$m, 
                s = expert_old$s))
  }
  
  tmp_exp = ExpertFunction$new("gammacount", 
                               list(m = expert_old$get_params()$m,
                                    s = expert_old$get_params()$s),
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
  
  tmp_update = gammacount.EM_notexact(expert_old = tmp_exp, 
                                            tl = tl, yl = yl, yu = yu, tu = tu, 
                                            exposure = exposure, 
                                            z_e_obs = z_pos_e_obs, z_e_lat = z_pos_e_lat, 
                                            k_e = k_e,
                                            penalty = penalty, pen_params = pen_params)
  
  return(list(p_zero = p_new, 
              m = tmp_update$m, 
              s = tmp_update$s))
}

######################################################################################
# Register the zigammacount at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################