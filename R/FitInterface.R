#' Main fitting function of LRMoE
#'
#' @param Y A N by d (\code{exact_Y=T}) or N by 4d (\code{exact_Y=F}) matrix of numerics,
#'          where N is sample size and d is the dimension of each obsevation.
#'          If the size is N by 4d, Each block of four columns should be organized as \code{(tl, yl, yu, tu)}, representing the
#'          truncation lower bound, censoring lower bound, censoring upper bound and truncation upper bound.
#' @param X A N*P matrix of numerics, where P is the number of covariates.
#'          The first column of \code{X} should be 1, which is the intercept.
#' @param alpha_init A g*P matrix of numerics, which contains initial guess of the logit regression coefficients.
#'                   The last row should all be zero, representing the default latent class.
#'                   If no initialization is provided, all coefficients are set to zero.
#' @param comp_dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param params_list A d * g matrix of list with paramster names and values,
#'                    which is the initial parameter guess for the corresponding \code{comp_dist}.
#' @param exposure A vector of length N, representing the exposure of the observations (how long it has been observed).
#' @param exact_Y TRUE/FALSE: whether \code{Y} is observed exactly, or with censoring and/or truncation.
#' @param penalty TRUE/FALSE: whether the parameters are penalized for their magnitude.
#'                Default (and recommended) is TRUE.
#' @param pen_alpha A numeric, which contains penalties for \code{alpha_init}.
#'                    If \code{penalty=T} but no \code{pen_alpha} is provided, a constant is used.
#' @param pen_params A list of length d, where each element is a sublist of length g.
#'                     Each sublist contains one numeric vector, which is the corresponding penalty for \code{params.init}.
#' @param eps Stopping criteria for loglikelihood convergence. Default is \code{1e-03}.
#' @param alpha_iter_max Maximum number of iterations for updating alpha. Defauls is 5.
#' @param ecm_iter_max Maximum number of iterations for ECM. Default is 200.
#' @param grad_jump TRUE/FALSE: whether to use an approximated gradient jump to speed up convergence.
#' @param grad_seq How are the gradient sequence selected. Default is \code{2^(seq(8)-1)-1}.
#' @param print_steps TRUE/FALSE: whether paramater updates are printed on screen. Default is TRUE.
#'
#'
#' @export
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
      pen_params = matrix(list(list(NULL)),
                          nrow = model_guess$nrow, ncol = model_guess$ncol)
      for(j in 1:nrow(pen_params)){
        for(k in 1:ncol(pen_params)){
          pen_params[j,k][[1]] = model_guess$select(j,k)$default_penalty()
        }
      }
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
