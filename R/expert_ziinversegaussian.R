# ziinversegaussian Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the distribution
######################################################################################
ziinversegaussian.params_init <- function(y){
  p0 = sum(y == 0) / sum(y >= 0)
  y = y[which(y > 0)]
  result = inversegaussian.params_init(y)
  return( c(result, list(p_zero = p0)) )
}

ziinversegaussian.exposurize <- function(params, exposure){
  return( params )
}

ziinversegaussian.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
ziinversegaussian.expert_ll_exact <- function(y, params){
  return( ifelse(y == 0, log(params[["p_zero"]]), log(1-params[["p_zero"]]) + inversegaussian.logpdf(params, x = y)) )
}

ziinversegaussian.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  return( return( faster_zi_result(tl, tu, yl, yu, params, "inversegaussian") ) )
}

ziinversegaussian.penalty <- function(params, penalty_params) {
  return( inversegaussian.penalty(params, penalty_params) )
}

######################################################################################
# dziinversegaussian, pziinversegaussian, qziinversegaussian and rziinversegaussian implementations.
######################################################################################
ziinversegaussian.simulation <- function(params, n) {
  return( (1-rbinom(n,1,params[["p_zero"]]))*inversegaussian.simulation(params, n) )
}

ziinversegaussian.mean <- function(params) {
  return( (1-params[["p_zero"]])*inversegaussian.mean(params) )
}

ziinversegaussian.variance <- function(params) {
  p0 = params[["p_zero"]]
  return( (1-p0)*inversegaussian.variance(params) + p0*(1-p0)*inversegaussian.mean(params)^2 )
}

ziinversegaussian.logpdf <- function(params, x) {
  return( ifelse(is.infinite(x), -Inf, inversegaussian.logpdf(params, x)) )
}

ziinversegaussian.pdf <- function(params, x) {
  return( ifelse(is.infinite(x), 0, inversegaussian.pdf(params, x)) )
}

ziinversegaussian.logcdf <- function(params, q) {
  return( ifelse(is.infinite(q), 0, inversegaussian.logcdf(params, q)) )
}

ziinversegaussian.cdf <- function(params, q) {
  return( ifelse(is.infinite(q), 1, inversegaussian.cdf(params, q)) )
}

ziinversegaussian.quantile <- function(params, p) {
  p0 = params[["p_zero"]]
  return( ifelse(p0 >= p, 0, inversegaussian.quantile(params, p - p0)) )
}

ziinversegaussian.lev <- function(params, u) {
  p0 = params[["p_zero"]]
  result = (1-p0)*inversegaussian.lev(params,u)
  return(result)
}
######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
ziinversegaussian.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations
  
  p_old = expert_old$get_params()$p_zero
  
  tmp_exp = ExpertFunction$new("inversegaussian", 
                               list(mean = expert_old$get_params()$mean,
                                    shape = expert_old$get_params()$shape),
                               pen_params)
  
  expert_ll_pos = tmp_exp$ll_exact(ye)
  
  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(ye, p_old, expert_ll_pos)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, 0.0, 0.0, 0.0)
  
  tmp_update = inversegaussian.EM_exact(tmp_exp, ye, exposure, z_pos_e_obs,
                                        penalty, pen_params)
  
  return(list(p_zero = p_new, 
              mean = tmp_update$mean, 
              shape = tmp_update$shape))
}

ziinversegaussian.EM_notexact <- function(expert_old, 
                                    tl, yl, yu, tu, 
                                    exposure, 
                                    z_e_obs, z_e_lat, k_e,
                                    penalty, pen_params) {
  # Perform the EM optimization with exact observations
  
  p_old = expert_old$get_params()$p_zero
  
  tmp_exp = ExpertFunction$new("inversegaussian", 
                               list(mean = expert_old$get_params()$mean,
                                    shape = expert_old$get_params()$shape),
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
  
  tmp_update = inversegaussian.EM_notexact(expert_old = tmp_exp, 
                                     tl = tl, yl = yl, yu = yu, tu = tu, 
                                     exposure = exposure, 
                                     z_e_obs = z_pos_e_obs, z_e_lat = z_pos_e_lat, 
                                     k_e = k_e,
                                     penalty = penalty, pen_params = pen_params)
  
  return(list(p_zero = p_new, 
              mean = tmp_update$mean, 
              shape = tmp_update$shape))
}
######################################################################################
# Register in the ExpertLibrary Object at zzz.R
######################################################################################