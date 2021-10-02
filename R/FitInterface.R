
FitLRMoE = function(Y, X, alpha_init,
                    comp_dist, params_list,
                    exposure = NULL,
                    exact_Y = FALSE,
                    penalty = TRUE, pen_alpha = 5.0, pen_params = NULL,
                    eps = 1e-3,
                    alpha_iter_max = 3, ecm_iter_max = 200,
                    grad_jump = TRUE, grad_seq = NULL,
                    print_steps = TRUE)
{

  Y = as.matrix(Y)
  X = as.matrix(X)

  if(is.null(exposure)){
    exposure= rep(1, nrow(X))
    warning("No exposure provided. The default value exposure=1 is used for all observations.")
  }

  model_guess = ExpertMatrix$new(comp_dist, params_list)

  if(penalty == FALSE){
    pen_alpha = Inf
  }else{
    if(missing(pen_alpha)){
      pen_alpha = 5.0
      warning("No pen_alpha provided. The default value pen_alpha=5.0 is used.")
    }
    if(is.null(pen_params)){
      pen_params = matrix(list(NULL),
                          nrow = model_guess$nrow, ncol = model_guess$ncol)
      warning("No pen_params provided. The default value pen_params is used.")
    }
  }

  model_guess$set_penalty_params(pen_params)

  if(exact_Y == TRUE){
    result = FitExact(Y = Y, X = X, alpha = alpha_init, model = model_guess,
                       exposure = exposure,
                       penalty = penalty, pen_alpha = pen_alpha, pen_params = pen_params,
                       eps = eps,
                       alpha_iter_max = alpha_iter_max, ecm_iter_max = ecm_iter_max,
                       grad_jump = grad_jump, grad_seq = grad_seq,
                       print_steps = print_steps)
  }else{
    result = FitNotExact(Y = Y, X = X, alpha = alpha_init, model = model_guess,
                      exposure = exposure,
                      penalty = penalty, pen_alpha = pen_alpha, pen_params = pen_params,
                      eps = eps,
                      alpha_iter_max = alpha_iter_max, ecm_iter_max = ecm_iter_max,
                      grad_jump = grad_jump, grad_seq = grad_seq,
                      print_steps = print_steps)
  }
  return(result)
}
