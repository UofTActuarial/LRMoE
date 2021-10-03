# ziburr Expert Function
######################################################################################
# Initialize the parameters of the distribuion
######################################################################################
ziburr.params_init <- function(y){
  # There is no clear way we could estimate ziburr distribution params from observations, so we would use MLE to estimate
  # Calculate all the parameters that needed for further calculation
  p_init = sum(y == 0) / sum(y>=0)
  y = y[which(y>0)]
  result = burr.params_init(y)
  return( c(result, list(p_zero = p_init)) )
}

ziburr.exposurize <- function(params, exposure){
  return( params )
}

ziburr.set_params <- function(params){
  # Check the parameters are valid for distribution
  return( params )
}

######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
ziburr.expert_ll_exact <- function(y, params){
  return(ifelse(y == 0,
                 log(params[["p_zero"]]),
                 log(1 - params[["p_zero"]]) + burr.logpdf(params, x = y))
        )
}

ziburr.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  return( faster_zi_result(tl, tu, yl, yu, params, "burr") )
}

ziburr.penalty <- function(params, penalty_params) {
  return( burr.penalty(params, penalty_params) )
}

ziburr.default_penalty <- function() {
  return(burr.default_penalty())
}

######################################################################################
# ddistribution, pdistribution, qdistribution and rdistribution implementations.
######################################################################################
ziburr.simulation <- function(params, n) {
  # simulated n points based on the distribution params
  p0 = params[["p_zero"]]
  return( (1 - rbinom(n,1,p0)) * burr.simulation(params, n) )
}

ziburr.mean <- function(params) {
  return( (1 - params[["p_zero"]])*burr.mean(params) )
}

ziburr.variance <- function(params) {
  p0 = params[["p_zero"]]
  return( (1-p0)*burr.variance(params) + p0*(1-p0)*burr.mean(params)^2 )
}

ziburr.logpdf <- function(params, x) {
  return( burr.logpdf(params, x) )
}

ziburr.pdf <- function(params, x) {
  return( burr.pdf(params, x) )
}

ziburr.logcdf <- function(params, q) {
  return( burr.logcdf(params, q) )
}

ziburr.cdf <- function(params, q) {
  return( burr.cdf(params, q) )
}

ziburr.quantile <- function(params, p) {
  p0 = params[["p_zero"]]
  return( ifelse(p0 >= p, 0, burr.quantile(params, p)) )
}

ziburr.lev <- function(params, u) {
  p0 = params[["p_zero"]]
  result = (1-p0) * burr.lev(params, u)
  return(result)
}
######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################

ziburr.EM_exact <- function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  tmp_exp = ExpertFunction$new("burr",
                               list(shape1 = expert_old$get_params()$shape1,
                                    shape2 = expert_old$get_params()$shape2,
                                    scale = expert_old$get_params()$scale),
                               pen_params)

  expert_ll_pos = tmp_exp$ll_exact(ye)

  z_zero_e_obs = z_e_obs * EM_E_z_zero_obs(ye, p_old, expert_ll_pos)
  z_pos_e_obs = z_e_obs - z_zero_e_obs
  p_new = EM_M_zero(z_zero_e_obs, z_pos_e_obs, 0.0, 0.0, 0.0)

  tmp_update = burr.EM_exact(tmp_exp, ye, exposure, z_pos_e_obs,
                              penalty, pen_params)

  return(list(p_zero = p_new,
              shape1 = tmp_update$shape1,
              shape2 = tmp_update$shape2,
              scale = tmp_update$scale))
}

ziburr.EM_notexact <- function(expert_old,
                                          tl, yl, yu, tu,
                                          exposure,
                                          z_e_obs, z_e_lat, k_e,
                                          penalty, pen_params) {
  # Perform the EM optimization with exact observations

  p_old = expert_old$get_params()$p_zero

  tmp_exp = ExpertFunction$new("burr",
                               list(shape1 = expert_old$get_params()$shape1,
                                    shape2 = expert_old$get_params()$shape2,
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

  tmp_update = burr.EM_notexact(expert_old = tmp_exp,
                                           tl = tl, yl = yl, yu = yu, tu = tu,
                                           exposure = exposure,
                                           z_e_obs = z_pos_e_obs, z_e_lat = z_pos_e_lat,
                                           k_e = k_e,
                                           penalty = penalty, pen_params = pen_params)

  return(list(p_zero = p_new,
              shape1 = tmp_update$shape1,
              shape2 = tmp_update$shape2,
              scale = tmp_update$scale))
}
######################################################################################
# Register in the ExpertLibrary Object
######################################################################################
