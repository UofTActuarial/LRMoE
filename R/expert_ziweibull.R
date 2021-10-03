# ziweibull Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the ziweibull
######################################################################################
ziweibull.params_init <- function(y){
  # Estimate all the parameters that needed for further calculation
  p0 = sum(y == 0)/sum(y >= 0)
  y = y[which(y>0)]
  result = weibull.params_init(y)
  return( c(result, list(p_zero = p0)) )
}

ziweibull.exposurize <- function(params, exposure){
  # Calculate the exposures
  return( params )
}

ziweibull.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}
######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
ziweibull.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  p0 = params[["p_zero"]]
  return( ifelse(y == 0, log(p0), log(1-p0) + weibull.logpdf(params, x = y)) )
}

ziweibull.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # If some observations are not exact, calculate its corresponding likelihood
  return( faster_zi_result(tl, tu, yl, yu, params, "weibull") )
}

ziweibull.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters
  return( weibull.penalty(params, penalty_params) )
}

ziweibull.default_penalty <- function() {
  return(weibull.default_penalty())
}

######################################################################################
# dziweibull, pziweibull, qziweibull and rziweibull implementations.
######################################################################################
ziweibull.simulation <- function(params, n) {
  # simulated n points based on the ziweibull params
  return( (1 - rbinom(n,1,params[["p_zero"]]))*weibull.simulation(params, n) )
}

ziweibull.mean <- function(params) {
  # Calculate the mean based on the params
  p0 = params[["p_zero"]]
  return( (1-p0)*weibull.mean(params) )
}

ziweibull.variance <- function(params) {
  # Calculate the variance based on the params
  p0 = params[["p_zero"]]
  return( (1-p0)*weibull.variance(params) + p0*(1-p0)*weibull.mean(params)^2 )
}

ziweibull.logpdf <- function(params, x) {
  # Return the log pdf based on the input x
  return(weibull.logpdf(params, x))
}

ziweibull.pdf <- function(params, x) {
  # Return the pdf based on the input x
  return(weibull.pdf(params, x))
}

ziweibull.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  return(weibull.logcdf(params,q))
}

ziweibull.cdf <- function(params, q) {
  # return the cdf based on the input x
  return(weibull.cdf(params, q))
}

ziweibull.quantile <- function(params, p) {
  # return the percentage points based on the value of p
  p0 = params[["p_zero"]]
  return( ifelse(p0 >= p, 0, weibull.quantile(params, p-p0)) )
}

ziweibull.lev <- function(params, u) {
  p0 = params[["p_zero"]]
  result = (1-p0)*weibull.lev(params, u)
  return(result)
}
######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
ziweibull.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  tmp_exp = ExpertFunction$new("weibull",
                               list(shape = expert_old$get_params()$shape,
                                    scale = expert_old$get_params()$scale),
                               pen_params)

  expert_ll_pos = tmp_exp$ll_exact(ye)

  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(ye, p_old, expert_ll_pos)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, 0.0, 0.0, 0.0)

  tmp_update = weibull.EM_exact(tmp_exp, ye, exposure, z_pos_e_obs,
                                penalty, pen_params)

  return(list(p_zero = p_new,
              shape = tmp_update$shape,
              scale = tmp_update$scale))
}

ziweibull.EM_notexact <- function(expert_old,
                                    tl, yl, yu, tu,
                                    exposure,
                                    z_e_obs, z_e_lat, k_e,
                                    penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  tmp_exp = ExpertFunction$new("weibull",
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

  tmp_update = weibull.EM_notexact(expert_old = tmp_exp,
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
# Register the ziweibull at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################
