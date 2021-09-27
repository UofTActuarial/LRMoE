

FitNotExact = function(Y, X, alpha, model,
                    exposure = NULL,
                    penalty = TRUE, pen_alpha = 5.0, pen_params = NULL,
                    eps = 1e-3,
                    alpha_iter_max = 3, ecm_iter_max = 200,
                    grad_jump = TRUE, grad_seq = NULL,
                    print_steps = TRUE)
{
  gate_init = GateLogit(X, alpha)
  ll_np_list = LogLikelihoodNotExact(Y, gate_init, model, exposure)
  ll_init_np = ll_np_list$ll
  ll_penalty = model$get_penalty_value()
  ll_init = ll_init_np + ll_penalty
  
  if(print_steps){
    print(paste("Initial loglik: ", ll_init_np, " (no penalty), ",
                ll_init, " (with penalty)", sep = ""))
  }
  
  alpha_em = alpha
  # model_em = model$clone()
  model_em = model
  gate_em = GateLogit(X, alpha_em)
  ll_em_list = LogLikelihoodNotExact(Y, gate_em, model_em, exposure)
  ll_em_np = ll_em_list$ll
  ll_em = ll_init
  ll_em_old = -Inf
  iter = 0
  
  while((ll_em - ll_em_old > eps) & (iter < ecm_iter_max)){
    
    iter = iter + 1
    ll_em_np_old = ll_em_np
    ll_em_old = ll_em
    
    # E-step
    z_e_obs = EM_E_z_obs(ll_em_list$gate_expert_ll_comp, ll_em_list$gate_expert_ll)
      # exp( sweep(ll_em_list$gate_expert_ll_comp, 1, ll_em_list$gate_expert_ll, FUN = "-", check.margin = FALSE) )
      # exp(XColMinusY(ll_em_list$gate_expert_ll_comp, ll_em_list$gate_expert_ll))
    z_e_lat = EM_E_z_lat(ll_em_list$gate_expert_tn_bar_comp_z_lat, ll_em_list$gate_expert_tn_bar_z_lat)
    k_e = EM_E_k(ll_em_list$gate_expert_tn_bar_k)
    
    # M-step: alpha
    ll_em_temp = ll_em
    alpha_em = EMMalpha(X, alpha_em, list(z.e.obs=z_e_obs, z.e.lat = z_e_lat, k.e = k_e), 
                        alpha_iter_max, penalty, pen_alpha)
    gate_em = GateLogit(X, alpha_em)
    ll_em_list = LogLikelihoodNotExact(Y, gate_em, model_em, exposure)
    ll_em_np = ll_em_list$ll
    ll_em_penalty = model_em$get_penalty_value()
    ll_em = ll_em_np + ll_em_penalty
    
    diff = ifelse(ll_em - ll_em_temp>0, "+", "-")
    pct = abs(ll_em - ll_em_temp)/abs(ll_em_old) * 100
    if(print_steps){
      print(paste("Iteration: ", iter, " ,",
                  " updating alpha: ", ll_em_temp, " -> ", ll_em,
                  " (", diff, pct, "%)",
                  sep = ""))
    }
    
    ll_em_temp = ll_em
    
    # M-step: Expert functions
    for(d in c(1:model_em$nrow)){
      for(j in c(1:model_em$ncol)){
        params_old = model_em$select(d,j)$get_params()
        params_new = model_em$select(d,j)$EM_notexact(
          expert_old = model_em$select(d,j), 
          tl = Y[,4*(d-1)+1], yl = Y[,4*(d-1)+2], yu = Y[,4*(d-1)+3], tu = Y[,4*(d-1)+4],
          exposure = exposure, 
          z_e_obs = z_e_obs[,j], z_e_lat = z_e_lat[,j], k_e = k_e,
          penalty = penalty, pen_params = pen_params[d,j][[1]]
        )
        model_em$select(d,j)$set_params(params_new)
        
        print(params_old)
        print(model_em$select(d,j)$get_params())
        
        ll_em_list = LogLikelihoodNotExact(Y, gate_em, model_em, exposure)
        ll_em_np = ll_em_list$ll
        ll_em_penalty = model_em$get_penalty_value()
        ll_em = ll_em_np + ll_em_penalty
        
        diff = ifelse(ll_em - ll_em_temp>0, "+", "-")
        pct = abs(ll_em - ll_em_temp)/abs(ll_em_old) * 100
        if(print_steps){
          print(paste("Iteration: ", iter, " ,",
                      " updating expert[", d, ", ", j, "]: ", ll_em_temp, " -> ", ll_em,
                      " (", diff, pct, "%)",
                      sep = ""))
        }
        ll_em_temp = ll_em
        
      }
    }
    
    alpha_em = alpha_em
    # model_em = model_em$clone()
    gate_em = GateLogit(X, alpha_em)
    ll_em_list = LogLikelihoodNotExact(Y, gate_em, model_em, exposure)
    ll_em_np = ll_em_list$ll
    ll_em_penalty = model_em$get_penalty_value()
    ll_em = ll_em_np + ll_em_penalty
    
  }
  
  
  
  converge = ifelse(ll_em - ll_em_old > eps, FALSE, TRUE)
  AIC = NaN # To Do
  BIC = NaN # To Do
  
  return(list(alpha_fit = alpha_em, model_fit = model_em,
              converge = converge, iter = iter,
              ll_np = ll_em_np, ll = ll_em,
              AIC = AIC, BIC = BIC))
}