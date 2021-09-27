# zipoisson Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the zipoisson
######################################################################################
zipoisson.params_init <- function(y){
  # Calculate all the parameters that needed for further calculation
  mu = mean(y)
  variance = var(y)
  lp = (variance - mu)/mu
  tmp = lp/mu
  p0 = tmp/(1+tmp)
  lambda = lp/p0
  return(list(p_zero = p0, lambda = lambda))
}

zipoisson.exposurize <- function(params, exposure){
  # Calculate the exposures
  p0 = params[["p_zero"]]
  result = poisson.exposurize(params, exposure)
  return( c(result, list(p_zero = p0)) )
}

zipoisson.set_params <- function(params){
  # Check the parameters are valid for distribution
  return(params)
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
zipoisson.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  p0 = params[["p_zero"]]
  return( ifelse(y == 0,
                 log(p0 + (1 - p0)*poisson.pdf(params, x = y)),
                 log(1 - p0) + poisson.logpdf(params, x = y)
                 )
          )
}

zipoisson.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # If some observations are not exact, calculate its corresponding likelihood
  return( faster_zi_result(tl, tu, yl, yu, params, "poisson") )
}

zipoisson.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  return(poisson.penalty(params, penalty_params))
}

######################################################################################
# dzipoisson, pzipoisson, qzipoisson and rzipoisson implementations.
######################################################################################
zipoisson.simulation <- function(params, n) {
  # simulated n points based on the zipoisson params
  return( (1-rbinom(n,1,params[["p_zero"]]))*poisson.simulation(params,n) )
}

zipoisson.mean <- function(params) {
  # Calculate the mean based on the params
  p0 = params[["p_zero"]]
  return( (1-p0)*poisson.mean(params) )
}

zipoisson.variance <- function(params) {
  # Calculate the variance based on the params
  p0 = params[["p_zero"]]
  return( (1-p0)*poisson.variance(params) + p0*(1-p0)*poisson.mean(params)^2 )
}

zipoisson.logpdf <- function(params, x) {
  # Return the log pdf based on the input x
  return( ifelse(is.infinite(x), 0, poisson.logpdf(params, x)) )
}

zipoisson.pdf <- function(params, x) {
  # Return the pdf based on the input x
  return( ifelse(is.infinite(x), -Inf, poisson.pdf(params, x)) )
}

zipoisson.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  return( ifelse(is.infinite(q), 0, poisson.logcdf(params, q)) )
}

zipoisson.cdf <- function(params, q) {
  # return the cdf based on the input x
  return( ifelse(is.infinite(q), 1, poisson.cdf(params, q)) )
}

zipoisson.quantile <- function(params, p) {
  # return the percentage points based on the value of p
  p0 = params[["p_zero"]]
  return( ifelse(p0 >= p, 0, poisson.quantile(params, p - p0)) )
}

######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
zipoisson.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations
  
  p_old = expert_old$get_params()$p_zero
  
  tmp_exp = ExpertFunction$new("poisson", 
                               list(lambda = expert_old$get_params()$lambda),
                               pen_params)
  
  expert_ll_pos = tmp_exp$ll_exact(ye)
  
  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(ye, p_old, expert_ll_pos)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, 0.0, 0.0, 0.0)
  
  tmp_update = poisson.EM_exact(tmp_exp, ye, exposure, z_pos_e_obs,
                                         penalty, pen_params)
  
  return(list(p_zero = p_new, 
              lambda = tmp_update$lambda))
}

zipoisson.EM_notexact <- function(expert_old, 
                                tl, yl, yu, tu, 
                                exposure, 
                                z_e_obs, z_e_lat, k_e,
                                penalty, pen_params) {
  # Perform the EM optimization with exact observations
  
  p_old = expert_old$get_params()$p_zero
  
  if(p_old > 0.999999){
    return(list(p_zero = p_old, 
                lambda = expert_old$lambda))
  }
  
  tmp_exp = ExpertFunction$new("poisson", 
                               list(lambda = expert_old$get_params()$lambda),
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
  
  tmp_update = poisson.EM_notexact(expert_old = tmp_exp, 
                                 tl = tl, yl = yl, yu = yu, tu = tu, 
                                 exposure = exposure, 
                                 z_e_obs = z_pos_e_obs, z_e_lat = z_pos_e_lat, 
                                 k_e = k_e,
                                 penalty = penalty, pen_params = pen_params)
  
  return(list(p_zero = p_new, 
              lambda = tmp_update$lambda))
}
######################################################################################
# Register the zipoisson at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################