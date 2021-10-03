#' Calculate the log likelihood
#'
#' @description
#' Calculate the log likelihood based on the observations, gate matrix, model matrix and exposure matrix.
#'
#' @param Y (`matrix`)\cr
#' The numerical observation matrix with IxJ dimentions
#' @param gate (`matrix`)\cr
#' The numerical gate matrix with IxD dimentions
#' @param model (`ExpertMatrix`)\cr
#' The JxD ExpertMatrix S6 Class which contained ExpertFunctions
#' @param exposure (`matrix`)\cr
#' The Ix1 numerical vector which specify the exposure for each observation.
#' @keywords internal
LogLikelihoodExact <- function(Y, gate, model, exposure){
  expert_ll_comp = matrix(0, nrow = nrow(Y), ncol = model$ncol)
  for(i in 1:nrow(Y)){ # Loop through observation 1:n
    expert_ll_dim_comp = matrix(-Inf, model$nrow, model$ncol)
    for(j in 1:model$ncol){ # Loop through component
      for(d in 1:model$nrow){ # Loop through dimension
        expert_expo = model$select(row_index = d, col_index = j)$exposurize(exposure[i])
        expert_ll_dim_comp[d,j] = expert_expo$ll_exact(Y[i,d])
      }
    }
    expert_ll_comp[i,] = colSums(expert_ll_dim_comp)
  }

  gate_expert_ll_comp = gate + expert_ll_comp
  gate_expert_ll = rowLogSumExps(gate_expert_ll_comp)
  norm_gate_expert_ll = gate_expert_ll
  ll = sum(norm_gate_expert_ll)

  return(list(
    expert_ll_comp = expert_ll_comp,
    gate_expert_ll_comp = gate_expert_ll_comp,
    gate_expert_ll = gate_expert_ll,
    norm_gate_expert_ll = norm_gate_expert_ll,
    ll = ll
))
}

#' Calculate the log likelihood
#'
#' @description
#' Calculate the log likelihood based on the observations, gate matrix, model matrix and exposure matrix.
#' Keep in mind our observation matrix now have a different format
#'
#' @param Y (`array`)\cr
#' The observation matrix with Ix4J dimentions, where each entry is a tuple of this format (tl, yl, tu, yu).
#' Keep in mind it means there will be three dimention for this matrix, and the last dim need to be 4.
#' @param gate (`matrix`)\cr
#' The numerical gate matrix with IxD dimentions
#' @param model (`ExpertMatrix`)\cr
#' The JxD ExpertMatrix S6 Class which contained ExpertFunctions
#' @param exposure (`matrix`)\cr
#' The Ix1 numerical vector which specify the exposure for each observation.
#' @keywords internal
LogLikelihoodNotExact <- function(Y, gate, model, exposure){

  expert_ll_comp = matrix(0, nrow = nrow(Y), ncol = model$ncol)
  expert_tn_comp = matrix(0, nrow = nrow(Y), ncol = model$ncol)
  expert_tn_bar_comp = matrix(0, nrow = nrow(Y), ncol = model$ncol)

  for(i in 1:nrow(Y)){ # Loop through observation 1:n
    expert_ll_dim_comp = matrix(-Inf, model$nrow, model$ncol)
    expert_tn_dim_comp = matrix(-Inf, model$nrow, model$ncol)
    expert_tn_bar_dim_comp = matrix(-Inf, model$nrow, model$ncol)

    for(j in 1:model$ncol){ # Loop through component
      for(d in 1:model$nrow){ # Loop through dimension
        expert_expo = model$select(row_index = d, col_index = j)$exposurize(exposure[i])
        result_set = expert_expo$ll_not_exact(Y[i,4*(d-1) + 1], Y[i,4*(d-1) + 2], Y[i,4*(d-1) + 3], Y[i,4*(d-1) + 4])
        expert_ll_dim_comp[d,j] = result_set[["expert_ll"]]
        expert_tn_dim_comp[d,j] = result_set[["expert_tn"]]
        expert_tn_bar_dim_comp[d,j] = result_set[["expert_tn_bar"]]
      }
    }
    expert_ll_comp[i,] = colSums(expert_ll_dim_comp)
    expert_tn_comp[i,] = colSums(expert_tn_dim_comp)
    expert_tn_bar_comp[i,] = colSums(expert_tn_bar_dim_comp)
  }

  gate_expert_ll_comp = gate + expert_ll_comp
  gate_expert_tn_comp = gate + expert_tn_comp
  gate_expert_tn_bar_comp = gate + log1mexp(-expert_tn_comp)
  gate_expert_tn_bar_comp_k = gate + log1mexp(-expert_tn_bar_comp)
  gate_expert_tn_bar_comp_z_lat = gate + expert_tn_bar_comp

  gate_expert_ll = rowLogSumExps(gate_expert_ll_comp)
  gate_expert_tn = rowLogSumExps(gate_expert_tn_comp)
  gate_expert_tn_bar = rowLogSumExps(gate_expert_tn_bar_comp)
  gate_expert_tn_bar_k = rowLogSumExps(gate_expert_tn_bar_comp_k)
  gate_expert_tn_bar_z_lat = rowLogSumExps(gate_expert_tn_bar_comp_z_lat)

  norm_gate_expert_ll = gate_expert_ll - gate_expert_tn_bar_k
  ll = sum(norm_gate_expert_ll)

  return(list(
    expert_ll_comp = expert_ll_comp,
    expert_tn_comp = expert_tn_comp,
    expert_tn_bar_comp = expert_tn_bar_comp,

    gate_expert_ll_comp = gate_expert_ll_comp,
    gate_expert_tn_comp = gate_expert_tn_comp,
    gate_expert_tn_bar_comp = gate_expert_tn_bar_comp,
    gate_expert_tn_bar_comp_k = gate_expert_tn_bar_comp_k,
    gate_expert_tn_bar_comp_z_lat = gate_expert_tn_bar_comp_z_lat,

    gate_expert_ll = gate_expert_ll,
    gate_expert_tn = gate_expert_tn,
    gate_expert_tn_bar = gate_expert_tn_bar,
    gate_expert_tn_bar_k = gate_expert_tn_bar_k,
    gate_expert_tn_bar_z_lat = gate_expert_tn_bar_z_lat,

    norm_gate_expert_ll = norm_gate_expert_ll,

    ll = ll
  ))
}
