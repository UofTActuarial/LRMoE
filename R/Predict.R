#' @title predict_class_prior
#' @description: Predicts the latent class probabilities, given covariates `X` and logit regression coefficients `alpha`.
#' 
#' @param X
#' A matrix of covariates.
#' @param alpha
#' A matrix of logit regression coefficients.
#' 
#' @return prob
#' A matrix of latent class probabilities.
#' @return max_prob_idx
#' A matrix of the most likely latent class for each observation.
predict_class_prior <- function(X, alpha) {
  prob = exp(GateLogit(X, alpha))
  return( list(prob = prob, max_prob_index = apply(prob, 1, which.max)) )
}

#' @title predict_class_posterior
#' @description: Predicts the latent class probabilities, given observations `Y`, covariates `X`, 
#' logit regression coefficients `alpha` and a specified `model` of expert functions.
#' 
#' @param Y
#' A matrix of responses.
#' @param X
#' A matrix of covariates.
#' @param alpha
#' A matrix of logit regression coefficients.
#' @param model
#' A matrix specifying the expert functions.
#' 
#' @param exact_Y
#' `true` or `false` (default), indicating if `Y` is observed exactly or with censoring and truncation.
#' @param exposure_past
#' A vector indicating the time exposure (past) of each observation. If nothing is supplied, it is set to 1.0 by default.
#' 
#' @return prob:
#' A matrix of latent class probabilities.
#' @return max_prob_idx:
#' A matrix of the most likely latent class for each observation.
predict_class_posterior <- function(Y, X, alpha, model, exposure_past = list(), exact_Y = F) {
  if(!length(exposure_past)) { exposure_past = rep(1, dim(X)[1])}
  gate = GateLogit(X, alpha)
  
  if(exact_Y) {ll_list = LogLikelihoodExact(Y, gate, model, exposure_past)}
  else{ll_list = LogLikelihoodNotExact(Y, gate, model, exposure_past)}
  
  z_e_obs = EM_E_z_obs(ll_list$gate_expert_ll_comp, ll_list$gate_expert_ll)
  return( list(prob = z_e_obs, max_prob_index = apply(z_e_obs, 1, which.max)) )
}

#' @title 
#' predict_mean_prior
#' @description
#' Predicts the mean values of response, given covariates `X`, logit regression coefficients `alpha` and a specified `model` of expert functions.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model An [ExpertMatrix]
#' @param exposure_future A vector indicating the time exposure (future) of each observation. If nothing is supplied,
#' it is set to 1.0 by default.
#' @return result:
#' A matrix of predicted mean values of response, based on prior probabilities.
predict_mean_prior <- function(X, alpha, model, exposure_future = list()) {
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1])}
  weights = predict_class_prior(X,alpha)$prob
  result = matrix(0, dim(X)[1], model$nrow)
  
  for(i in c(1:dim(X)[1])){
    tem_expert_matrix = model$exposurize(exposure_future[i])
    tem_mean = tem_expert_matrix$get_mean()
    result[i,] = tem_mean %*% weights[i,]
  }
  
  return(result)
}

#' @title 
#' predict_mean_posterior
#' @description: Predicts the mean value of the responses, given observations `Y`, covariates `X`, 
#' logit regression coefficients `alpha` and a specified `model` of expert functions.
#' 
#' @param Y A matrix of responses.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model A matrix specifying the expert functions.
#' 
#' @param exact_Y Bool variable. indicating if `Y` is observed exactly or with censoring and truncation. Default set to be False
#' @param exposure_past A vector indicating the time exposure (past) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' @param exposure_future A vector indicating the time exposure (future) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' @return result:
#' A matrix of predicted mean values of response, based on posterior probabilities.
predict_mean_posterior <- function(Y, X, alpha, model, exposure_past = list(), exposure_future = list(), exact_Y = F) {
  if(!length(exposure_past)) { exposure_past = rep(1, dim(X)[1])}
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1])}
  
  weights = weights = predict_class_posterior(Y, X, alpha, model, exposure_past, exact_Y)$prob
  result = matrix(0, dim(X)[1], model$nrow)
  
  for(i in c(1:dim(X)[1])){
    tem_expert_matrix = model$exposurize(exposure_future[i])
    tem_mean = tem_expert_matrix$get_mean()
    result[i,] = tem_mean %*% weights[i,]
  }
  
  return(result)
}

#' @title predict_var_prior
#' @description 
#' Predicts the variance of response, given covariates `X`,
#' logit regression coefficients `alpha` and a specified `model` of expert functions.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model A matrix specifying the expert functions.
#' @param exposure_future A vector indicating the time exposure of each observation. If nothing is supplied,
#' it is set to 1.0 by default.
#' @return result A matrix of predicted variance of response, based on prior probabilities.
predict_var_prior <- function(X, alpha, model, exposure_future = list()) {
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1])}
  weights = predict_class_prior(X,alpha)$prob
  g_mean = predict_mean_prior(X, alpha, model, exposure_future)
  
  var_c_mean = matrix(0, dim(X)[1], model$nrow)
  mean_c_var = matrix(0, dim(X)[1], model$nrow)
  
  for(i in c(1:dim(X)[1])) {
    tem_expert_matrix = model$exposurize(exposure_future[i])
    c_mean = tem_expert_matrix$get_mean()
    c_var = tem_expert_matrix$get_variance()
    
    var_c_mean[i,] = (c_mean - g_mean[i,])^2 %*% weights[i,]
    mean_c_var[i,] = c_var %*% weights[i,]
  }
  
  return(var_c_mean + mean_c_var)
}

#' @title predict_var_posterior
#' @description 
#' Predicts the variance of response, given observations `Y`, covariates `X`, 
#' logit regression coefficients `alpha` and a specified `model` of expert functions.
#' 
#' @param Y A matrix of responses.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model A matrix specifying the expert functions.
#' 
#' @param exact_Y Bool variable. indicating if `Y` is observed exactly or with censoring and truncation. Default set to be False
#' @param exposure_past A vector indicating the time exposure (past) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' @param exposure_future A vector indicating the time exposure (future) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' @return result A matrix of predicted variance of response, based on posterior probabilities.
predict_var_posterior <- function(Y, X, alpha, model, exposure_past = list(), exposure_future = list(), exact_Y = F) {
  if(!length(exposure_past)) { exposure_past = rep(1, dim(X)[1])}
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1])}
  weights = predict_class_posterior(Y, X, alpha, model, exposure_past, exact_Y = exact_Y)$prob
  g_mean = predict_mean_posterior(Y, X, alpha, model, exposure_past, exposure_future, exact_Y)
  
  var_c_mean = matrix(0, dim(X)[1], model$nrow)
  mean_c_var = matrix(0, dim(X)[1], model$nrow)
  
  for(i in c(1:dim(X)[1])) {
    tem_expert_matrix = model$exposurize(exposure_future[i])
    c_mean = tem_expert_matrix$get_mean()
    c_var = tem_expert_matrix$get_variance()
    
    var_c_mean[i,] = (c_mean - g_mean[i,])^2 %*% weights[i,]
    mean_c_var[i,] = c_var %*% weights[i,]
  }
  
  return(var_c_mean + mean_c_var)
}

#' @title predict_limit_prior
#' @description Predicts the limit expected value (LEV) of response, given covariates `X`, 
#' logit regression coefficients `alpha` and a specified `model` of expert functions.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model A matrix specifying the expert functions.
#' @param limit A matrix specifying the cutoff point.
#' @param exposure_future A vector indicating the time exposure (future) of each observation. If nothing is supplied, it is set to 1.0 by default.
#' @return result: A matrix of predicted limit expected value of response, based on prior probabilities.
predict_limit_prior <- function(X, alpha, model, limit, exposure_future = list()) {
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1])}
  weights = predict_class_prior(X,alpha)$prob
  result = matrix(0, dim(X)[1], model$nrow)
  
  for(i in c(1:dim(X)[1])) {
    tem_expert_matrix = model$exposurize(exposure_future[i])
    tem_means_matrix = matrix(NA, nrow = tem_expert_matrix$nrow, ncol = tem_expert_matrix$ncol)
    for(row in c(1, tem_expert_matrix$nrow)){
      for(col in c(1, tem_expert_matrix$col)){
        tem_means_matrix[row, col] = tem_expert_matrix$select(row,col)$get_lev(limit[i,row])
      }
    }
    result[i,] = tem_means_matrix %*% weights[i,]
  }
  return(result)
}

#' @title predict_limit_posterior
#' @description Predicts the limit expected value (LEV) of response, given observations `Y`, covariates `X`, 
#' logit regression coefficients `alpha` and a specified `model` of expert functions.
#' @param Y A matrix of responses.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model A matrix specifying the expert functions.
#' @param limit A matrix specifying the cutoff point.
#' 
#' @param exact_Y Bool variable. indicating if `Y` is observed exactly or with censoring and truncation. Default set to be False
#' @param exposure_past A vector indicating the time exposure (past) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' @param exposure_future A vector indicating the time exposure (future) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' 
#' @return result: A matrix of predicted limit expected value of response, based on posterior probabilities.
predict_limit_posterior <- function(Y, X, alpha, model, limit, exposure_past = list(), exposure_future = list(), exact_Y = F) {
  if(!length(exposure_past)) { exposure_past = rep(1, dim(X)[1])}
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1])}
  weights = predict_class_posterior(Y, X, alpha, model, exposure_past, exact_Y = exact_Y)$prob
  result = matrix(0, dim(X)[1], model$nrow)
  
  for(i in c(1:dim(X)[1])) {
    tem_expert_matrix = model$exposurize(exposure_future[i])
    tem_means_matrix = matrix(NA, nrow = tem_expert_matrix$nrow, ncol = tem_expert_matrix$ncol)
    for(row in c(1, tem_expert_matrix$nrow)){
      for(col in c(1, tem_expert_matrix$col)){
        tem_means_matrix[row, col] = tem_expert_matrix$select(row,col)$get_lev(limit[i,row])
      }
    }
    result[i,] = tem_means_matrix %*% weights[i,]
  }
  return(result)
}

#' @title predict_excess_prior
#' @description Predicts the excess expectation of response, given covariates `X`, 
#' logit regression coefficients `alpha` and a specified `model` of expert functions.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model A matrix specifying the expert functions.
#' @param limit A matrix specifying the cutoff point.
#' @param exposure_future A vector indicating the time exposure (future) of each observation. If nothing is supplied, it is set to 1.0 by default.
#' @return result: A matrix of predicted excess expectation of response, based on prior probabilities.
predict_excess_prior <- function(X, alpha, model, limit, exposure_future = list()) {
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1])}
  weights = predict_class_prior(X,alpha)$prob
  result = matrix(0, dim(X)[1], model$nrow)
  
  for(i in c(1:dim(X)[1])) {
    tem_expert_matrix = model$exposurize(exposure_future[i])
    tem_means_matrix = matrix(NA, nrow = tem_expert_matrix$nrow, ncol = tem_expert_matrix$ncol)
    for(row in c(1, tem_expert_matrix$nrow)){
      for(col in c(1, tem_expert_matrix$col)){
        tem_means_matrix[row, col] = tem_expert_matrix$select(row,col)$get_excess(limit[i,row])
      }
    }
    result[i,] = tem_means_matrix %*% weights[i,]
  }
  return(result)
}

#' @title predict_excess_posterior
#' @description Predicts the excess expectation of response, given covariates `X`, 
#' logit regression coefficients `alpha` and a specified `model` of expert functions.
#' @param Y A matrix of responses.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model A matrix specifying the expert functions.
#' @param limit A matrix specifying the cutoff point.
#' 
#' @param exact_Y Bool variable. indicating if `Y` is observed exactly or with censoring and truncation. Default set to be False
#' @param exposure_past A vector indicating the time exposure (past) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' @param exposure_future A vector indicating the time exposure (future) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' 
#' @return result: A matrix of predicted excess expectation of response, based on posterior probabilities.
predict_excess_posterior <- function(Y, X, alpha, model, limit, exposure_past = list(), exposure_future = list(), exact_Y = F) {
  if(!length(exposure_past)) { exposure_past = rep(1, dim(X)[1])}
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1])}
  weights = predict_class_posterior(Y, X, alpha, model, exposure_past, exact_Y = exact_Y)$prob
  result = matrix(0, dim(X)[1], model$nrow)
  
  for(i in c(1:dim(X)[1])) {
    tem_expert_matrix = model$exposurize(exposure_future[i])
    tem_means_matrix = matrix(NA, nrow = tem_expert_matrix$nrow, ncol = tem_expert_matrix$ncol)
    for(row in c(1, tem_expert_matrix$nrow)){
      for(col in c(1, tem_expert_matrix$col)){
        tem_means_matrix[row, col] = tem_expert_matrix$select(row,col)$get_excess(limit[i,row])
      }
    }
    result[i,] = tem_means_matrix %*% weights[i,]
  }
  return(result)
}

# Helper Function For Internal Usage
# Solve a quantile of a mixture model
# Bisection method seems to give the most stable results
solve_continuous_mix_quantile <- function(weights, experts, p) {
  experts_ll = rep(0, length(experts))
  for(i in c(1:length(experts))){ experts_ll[i] = experts[[i]]$get_cdf(0) }
  p0 = sum(weights * experts_ll)
  
  if(p <= p0) {
    return(0)
  } else {
    experts_quantile = rep(0, length(experts))
    init_guess = 0
    for(counter in c(1:1000)){
      for(i in c(1:length(experts))){ experts_quantile[i] = experts[[i]]$get_quantile(p) }
      init_guess_tem = max(experts_quantile)
      if(init_guess_tem > init_guess){init_guess = init_guess_tem}
    }
    tryCatch(
    {
      target_function <- function(y) {
        experts_ll = rep(0, length(experts))
        for(i in c(1:length(experts))){ experts_ll[i] = experts[[i]]$get_cdf(y) }
        return(sum(weights * experts_ll) - p)
      }
      VaR = uniroot(target_function, lower = 0, upper = init_guess + 500)
    },
    error=function(cond) {
      message("solve_continuous_mix_quantile has an optimization error")
      message(cond)
      return(NA)
    })
  }
}

# Helper Function For Internal Usage
calc_continuous_CTE <- function(weights, experts, p, VaR) {
  expert_mean = c()
  for(expert in experts) { expert_mean = c(expert_mean, expert$get_mean()) }
  m = sum(weights * expert_mean)
  
  lev_list = c() # Get the lev for all the experts
  for(expert in experts) { lev_list = c(lev_list, expert$get_lev(VaR)) }
  lim_ev = sum(weights * lev_list)
  
  return(VaR + (m - lim_ev)/(1-p))
}

#' @title predict_VaRCTE_prior
#' @description Predicts the `p`-th value-at-risk (VaR) and conditional tail expectation (CTE) of response, 
#' given observations `Y`, covariates `X`, logit regression coefficients `alpha` and a specified `model` of expert functions.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model A matrix specifying the expert functions.
#' @param p A matrix of probabilities.
#' 
#' @param exposure_future A vector indicating the time exposure (future) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' 
#' @return result: list(`VaR`, `CTE`)
#' `VaR`: A matrix of predicted VaR of response, based on prior probabilities.
#' `CTE`: A matrix of predicted CTE of response, based on prior probabilities.
predict_VaRCTE_prior <- function(X, alpha, model, p, exposure_future = list()) {
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1]) }
  weights = predict_class_prior(X, alpha)$prob
  VaR = matrix(NA, nrow = dim(X)[1], ncol = model$nrow)
  CTE = matrix(NA, nrow = dim(X)[1], ncol = model$nrow)
  
  for(i in c(1:dim(X)[1])) {
    tem_expert_matrix = model$exposurize(exposure_future[i])
    for(i in c(1, dim(X)[1])){
      for(k in c(1, model$nrow)){
        VaR[i, k] = solve_continuous_mix_quantile(weights[i,],
                                                  tem_expert_matrix$select(row_index = k), p[i,k])
        CTE[i, k] = calc_continuous_CTE(weights[i,],
                                        tem_expert_matrix$select(row_index = k), p[i,k], VaR[i,k])
      }
    }
  }
  
  return(list(VaR = VaR, CTE = CTE))
}

#' @title predict_VaRCTE_posterior
#' @description Predicts the `p`-th value-at-risk (VaR) and conditional tail expectation (CTE) of response, 
#' given observations `Y`, covariates `X`, logit regression coefficients `alpha` and a specified `model` of expert functions.
#' 
#' @param Y A matrix of responses.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param model A matrix specifying the expert functions.
#' @param p A matrix of probabilities.
#' 
#' @param exact_Y Bool variable. indicating if `Y` is observed exactly or with censoring and truncation. Default set to be False
#' @param exposure_past A vector indicating the time exposure (past) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' @param exposure_future A vector indicating the time exposure (future) of each observation. If nothing is supplied,it is set to 1.0 by default.
#' 
#' @return result: list(`VaR`, `CTE`)
#' `VaR`: A matrix of predicted VaR of response, based on posterior probabilities.
#' `CTE`: A matrix of predicted CTE of response, based on posterior probabilities.
predict_VaRCTE_posterior <- function(Y, X, alpha, model, p, exposure_past = list(), exposure_future = list(), exact_Y = F) {
  if(!length(exposure_future)) { exposure_future = rep(1, dim(X)[1]) }
  if(!length(exposure_past)) { exposure_past = rep(1, dim(X)[1]) }
  
  weights = predict_class_posterior(Y, X, alpha, model, exposure_past = exposure_past, exact_Y = F)$prob
  VaR = matrix(NA, nrow = dim(X)[1], ncol = model$nrow)
  CTE = matrix(NA, nrow = dim(X)[1], ncol = model$nrow)
  
  for(i in c(1:dim(X)[1])) {
    tem_expert_matrix = model$exposurize(exposure_future[i])
    for(i in c(1, dim(X)[1])){
      for(k in c(1, model$nrow)){
        VaR[i, k] = solve_continuous_mix_quantile(weights[i,],
                                                  tem_expert_matrix$select(row_index = k), p[i,k])
        CTE[i, k] = calc_continuous_CTE(weights[i,],
                                        tem_expert_matrix$select(row_index = k), p[i,k], VaR[i,k])
      }
    }
  }
  
  return(list(VaR = VaR, CTE = CTE))
}