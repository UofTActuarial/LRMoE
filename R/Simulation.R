#' Simulation based on Logit Gating Function
#' @name sim_logit_gating
# 
#' @description Given alpha and x, generate a probability matrix. Give the simulation matrix with each row
#' of the probability matrix to be multinomial probability.
#' 
#' @param probs The probability matrix with N x P dimentions. 
#'
#' @return result The result matrix with N X P results
sim_logit_gating <- function(probs) {
  result = matrix(0, nrow = dim(probs)[1], ncol = dim(probs)[2])
  for(index in c(1:dim(probs)[1])) {
    result[index,] = rmultinom(1, 1, probs[index,])
  }
  return(result)
}

#' Simulations of Exposurized Models
#' @name sim_exposurize_model
#' 
#' @description Give one simulation for each expert in the expert_matrix after they got exposurized.
#' 
#' @param expert_matrix (`ExpertMatrix`)\cr
#' The Expert Matrix that will be exposurized and simulated.
#' @param exposure_list (`numeric`)\cr
#' The list of exposure values. Make sure this is a row vector with length to be `expert_matrix$ncol`
#' @param selected_index (`numeric`)\cr
#' The list that indicate which latent class is simulated.
#'
#' @return result (`matrix`)\cr
#' The result matrix with dimension to be `expert_matrix$nrow` rows and `expert_matrix$ncol` columns
sim_exposurize_model <- function(expert_matrix, exposure_list, selected_index) {
  result = matrix(0, nrow = length(selected_index), ncol = expert_matrix$nrow)
  for(i in c(1:length(selected_index))){
    expert_row = expert_matrix$select(col_index = selected_index[i])
    for(d in c(1:length(expert_row))){
      exposurized_expert = expert_row[[d]]$exposurize(exposure_list[i])
      result[i,d] = exposurized_expert$simulate(1)
    }
  }
  return(result)
}

#' Generate the simulations of the expert_matrix
#'
#' @param alpha (`matrix`)
#' A g * P matrix. Logit regression coefficients.
#' @param x (`matrix`)
#' An N * P covariate matrix, where N is sample size. The first column MUST be 1.
#' @param expert_matrix (`ExpertMatrix`)
#' A D*g expert matrix. 
#' @param exposure
#' A N*1 vector that contain the exposure value for each column in expert_matrix
#'
#' @return A N*D simulation matrix of expert_matrix.
#' 
#' @export
#'
#' @examples
#' #alpha = matrix(runif(20,-1,1), nrow = 4)
#' #x = matrix(runif(35, -1, 1), nrow = 7)
#' #params = matrix(list( list(meanlog = 100, sdlog = 2),  list(meanlog = 100, sdlog  = 2),
#' #                      list(meanlog = 333, sdlog = 1),  list(meanlog = 533, sdlog = 1),
#' #                      list(meanlog = 1, sdlog = 2),  list(meanlog = 1, sdlog  = 2),
#' #                      list(meanlog = 3, sdlog = 1),  list(meanlog = 5, sdlog = 1)),
#' #                      nrow = 2, byrow = T)
#' # expert_names = matrix( c("lognormal", "lognormal", "lognormal", "lognormal",
#' #                        "lognormal", "lognormal", "lognormal", "lognormal"), nrow = 2, byrow = T)
#' #expert_matrix = ExpertMatrix$new(expert_matrix = expert_names, expert_params_matrix = params)
#' #for(expert in expert_matrix$expert_matrix) {
#' # expert$set_penalty_params(c(1,1,1))
#' #  expert$initialize_penalty()
#' #}
#' #exposure_list = c(1,2,1,2,2,2,1)
#' #sim_dataset(alpha, x, expert_matrix, exposure = exposure_list)

sim_dataset <- function(alpha, x, expert_matrix, exposure = rep(1.0, dim(x)[1])) {
  probs = exp(GateLogit(x, alpha))
  gating_sim = sim_logit_gating(probs)
  selected_index = apply(gating_sim, 1, which.max)
  simulations = sim_exposurize_model(expert_matrix, exposure, selected_index)
  return( simulations ) # <- what is going on here
}

