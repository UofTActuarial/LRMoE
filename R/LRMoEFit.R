## Main fitting function
#' Fit an LRMoE model.
#'
#' @param Y A N*4d matrix of numerics, where N is sample size and d is the dimension of each obsevation.
#'          Each block of four columns should be organized as \code{(tl, yl, yu, tu)}, representing the
#'          truncation lower bound, censoring lower bound, censoring upper bound and truncation upper bound.
#' @param X A N*P matrix of numerics, where P is the number of covariates.
#'          The first column of \code{X} should be 1, which is the intercept.
#' @param n.comp A numeric which indicates the number of experts desired to fit the data. Default value is 2.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#'                  See below for more details.
#' @param alpha.init A g*P matrix of numerics, which contains initial guess of the logit regression coefficients.
#'                   The last row should all be zero, representing the default latent class.
#'                   If no initialization is provided, all coefficients are set to zero.
#' @param zero.init A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#'                  If the corresponding entry in \code{comp.dist} is not zero-inflated, zero value must by supplied.
#' @param params.init A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the initial parameter guess for the corresponding \code{comp.dist}.
#' @param penalty TRUE/FALSE: whether the parameters are penalized for their magnitude.
#'                Default (and recommended) is TRUE.
#' @param hyper.alpha A numeric, which contains penalties for \code{alpha.init}.
#'                    If \code{penalty=T} but no \code{hyper.alpha} is provided, a constant is used.
#' @param hyper.params A list of length d, where each element is a sublist of length g.
#'                     Each sublist contains one numeric vector, which is the corresponding penalty for \code{params.init}.
#' @param eps Stopping criteria for loglikelihood convergence. Default is \code{1e-03}.
#' @param alpha.iter.max Maximum number of iterations for updating alpha. Defauls is 5.
#' @param ecm.iter.max Maximum number of iterations for ECM. Default is 200.
#' @param grad.jump TRUE/FALSE: whether to use an approximated gradient jump to speed up convergence.
#' @param grad.period How often should \code{grad.jump} occur. Default is every 5 iterations.
#' @param grad.seq How are the gradient sequence selected. Default is \code{2^(seq(8)-1)-1}.
#' @param print TRUE/FALSE: whether paramater updates are printed on screen. Default is TRUE.
#'
#'
#' @export LRMoEFit
#'
LRMoEFit = function(Y, X, n.comp = 2, comp.dist = NULL,
                     alpha.init = NULL,
                     zero.init = NULL, params.init = NULL,
                     penalty = TRUE, hyper.alpha = NULL, hyper.params = NULL,
                     eps=1e-03, alpha.iter.max=5, ecm.iter.max=200,
                     grad.jump = TRUE,
                     grad.period = 5, grad.seq = 2^(seq(8)-1)-1,
                     print = TRUE)
{
  # Initialization: Constants
  dim.m = ncol(Y)/4 # no. of dimensions of observation
  sample.size.n = nrow(Y) # sample size
  n.covar.p = ncol(X) # no. of covariates
  g = n.comp # no. of component distributions

  # Initialization: Convert data into matrix
  Y = data.matrix(Y)
  X = data.matrix(X)

  # Initialization: Some parameter penalizations are specially treated.
  # See penalty switch for more detail.
  # invgauss, lnorm: NO penalty at all, for numerical speed.
  # weibull: scale parameter has no penalty, for numerical speed

  # Initialization: NULL hyper parameters
  if(is.null(hyper.alpha)){hyper.alpha = 1}

  if(is.null(hyper.params)){
    hyper.params = params.init # Tempopary place holder, with same dimension as params.init
    for(k in 1:dim.m){ for(j in 1:n.comp) { hyper.params[[k]][[j]] = DimCompExpertPenaltyInit(comp.dist[k,j]) }  }
  }


  # Initialization: NULL arguments not provided by user
  #############################################################
  #############################################################

  # Initialize Gating function, Expert functions and Loglik
  if(is.null(alpha.init)){alpha.init = matrix(0, nrow = n.comp, ncol = n.covar.p)}
  gate.ll.init = GateLogit(X, alpha.init)
  expert.list.init = DimCompExpertLL(Y, comp.dist, zero.init, params.init)
  ll.list.init = GateExpertLL(alpha.init, gate.ll.init, expert.list.init, penalty, hyper.alpha, hyper.params)


  # Initilization for ECM algorithm
  iter = 0
  alpha.em = alpha.init # This copy goes into EM
  zero.em = zero.init # This copy goes into EM
  for(k in 1:dim.m){zero.em[k,] = expert.list.init[[k]]$ zero.prob}
  params.em = params.init # This copy goes into EM
  gate.ll.em = gate.ll.init # This copy goes into EM
  expert.list.em = DimCompExpertLL(Y, comp.dist, zero.em, params.em) # This copy goes into EM
  ll.list.em = GateExpertLL(alpha.em, gate.ll.em, expert.list.em, penalty, hyper.alpha, hyper.params) # This copy goes into EM
  loglik.em = ll.list.em$ll

  if(print) {print(paste("Initial")); print(loglik.em)}
  loglik.em.old = -Inf

  # Iterations of EM algorithm
  # Stopping criteria: Insignificant improvement of loglik.em, or reaching iter.max
  while( (loglik.em - loglik.em.old > eps) & iter<ecm.iter.max )
  {
    # Take results from the previous iteration, and update the iteration count
    iter = iter + 1
    loglik.em.old = loglik.em

    # E-Step: z.obs, z.lat, k.e
    comp.zkz.e.list = EMEzkz(gate.ll.em, expert.list.em, ll.list.em)

    # M-Step: gating weights, i.e. alpha
    alpha.old = alpha.em
    alpha.em = EMMalpha(X, alpha.old, comp.zkz.e.list, alpha.iter.max, penalty, hyper.alpha)

    # Check that loglik increases after M-step of alpha
    gate.ll.em = GateLogit(X, alpha.em)
    expert.list.em = DimCompExpertLL(Y, comp.dist, zero.em, params.em)
    ll.list.em = GateExpertLL(alpha.em, gate.ll.em, expert.list.em, penalty, hyper.alpha, hyper.params) # This copy goes into EM
    loglik.em = ll.list.em$ll
    if(print) {print(paste("Iter:", iter, "(Update of alpha). Loglik = ", round(loglik.em, digits=6))); print(alpha.em); print("Updating Components: ");}
    # if(print) {print(paste("Iteration:", iter, "(Update of alpha). Loglik =", loglik.em))}


    # M-Step: by dimension and by component, all other parameters
    for(k in 1:dim.m)
    {
      for(j in 1:n.comp)
      {
        # Read information of component distribution
        # Strings and scalar/vector of parameters
        comp.kj.dist            = expert.list.em[[k]]$dim.comp.dist[j]
        comp.kj.zero.inflation  = expert.list.em[[k]]$zero.inflation[j]
        comp.kj.zero.prob.old   = zero.em[k,j]
        comp.kj.params.old      = params.em[[k]][[j]]
        # Vectors. Dim: sample.size.n * 1
        comp.kj.pos.expert.ll     = expert.list.em[[k]]$pos.expert.ll[,j]
        comp.kj.pos.expert.tn     = expert.list.em[[k]]$pos.expert.tn[,j]
        comp.kj.pos.expert.tn.bar = expert.list.em[[k]]$pos.expert.tn.bar[,j]
        comp.kj.expert.ll         = expert.list.em[[k]]$pos.expert.ll[,j]
        comp.kj.expert.tn         = expert.list.em[[k]]$pos.expert.tn[,j]
        comp.kj.expert.tn.bar     = expert.list.em[[k]]$pos.expert.tn.bar[,j]

        # Observations needed for loglik maximization
        tl.k = Y[,(4*(k-1)+1)]
        yl.k = Y[,(4*(k-1)+2)]
        yu.k = Y[,(4*(k-1)+3)]
        tu.k = Y[,(4*(k-1)+4)]

        # E-Step: First find conditional expectation: z.zero...., z.pos.... . Dimension: sample.size.n * 1
        z.zero.e.obs  = array(0, dim = c(sample.size.n, 1))
        z.zero.e.lat  = array(0, dim = c(sample.size.n, 1))
        z.pos.e.obs   = array(0, dim = c(sample.size.n, 1))
        z.pos.e.obs   = array(0, dim = c(sample.size.n, 1))

        z.zero.e.obs = comp.zkz.e.list$z.e.obs[,j] * EMEzzeroobs(yl.k, comp.kj.zero.inflation, comp.kj.zero.prob.old, comp.kj.pos.expert.ll)
        z.pos.e.obs = comp.zkz.e.list$z.e.obs[,j] - z.zero.e.obs

        z.zero.e.lat = comp.zkz.e.list$z.e.lat[,j] * EMEzzerolat(tl.k, comp.kj.zero.inflation, comp.kj.zero.prob.old, comp.kj.pos.expert.tn.bar)
        z.pos.e.lat = comp.zkz.e.list$z.e.lat[,j] - z.zero.e.lat

        # M-Step: loglik maximization wrt zero probabilities
        if(comp.kj.zero.inflation==TRUE) # Only do this for zero-inflated distributions
        {
          gate.ll.em = GateLogit(X, alpha.em) # alpha.em has been updated
          expert.list.em = DimCompExpertLL(Y, comp.dist, zero.em, params.em)
          ll.list.em = GateExpertLL(alpha.em, gate.ll.em, expert.list.em, penalty, hyper.alpha, hyper.params)
          loglik.em.before.kj = ll.list.em$ll

          comp.kj.zero.prob.old = zero.em[k,j]
          comp.kj.zero.prob.new = EMMzero(z.zero.e.obs, z.pos.e.obs, z.zero.e.lat, z.pos.e.lat, comp.zkz.e.list$k.e)
          zero.em[k,j] = comp.kj.zero.prob.new

          gate.ll.em = GateLogit(X, alpha.em) # alpha.em has been updated
          expert.list.em = DimCompExpertLL(Y, comp.dist, zero.em, params.em) # zero.em[k,j] has been updated
          ll.list.em = GateExpertLL(alpha.em, gate.ll.em, expert.list.em, penalty, hyper.alpha, hyper.params)
          loglik.em.after.kj = ll.list.em$ll

          if(print){ print(paste("Dim ", k, " Comp ", j, "z: ", loglik.em.before.kj, "  ", loglik.em.after.kj)); }

          if(loglik.em.after.kj < loglik.em.before.kj){
            if(print){
              print(paste("No update on component zero: Dim ", k, " Comp ", j, " !" ));
              print(paste("Decrease in loglik = ", loglik.em.before.kj - loglik.em.after.kj, " or ", (loglik.em.before.kj - loglik.em.after.kj)/abs(loglik.em.before.kj)*100, "%"));
              print(paste("Intended Update: ", comp.kj.zero.prob.new));# ; print (cat(params.em[[k]][[j]], sep = " "))
            }
            zero.em[k,j] = comp.kj.zero.prob.old
          } # No update, if the loglik decreases due to numerical imprecision.
        }

        # M-Step: loglik maximization wrt to component distribution parameters
        # I should only use those y's that are POSITIVE for severity distributions
        # All y's should be used for frequency distributions
        # This subtlety is dealt with at a lower level of optimization, i.e. component specific.

        # There may also be numerical imprecisions when calling optim, etc.
        # Add a insurance step!

        gate.ll.em = GateLogit(X, alpha.em) # alpha.em has been updated
        expert.list.em = DimCompExpertLL(Y, comp.dist, zero.em, params.em) # zero.em has been updated
        ll.list.em = GateExpertLL(alpha.em, gate.ll.em, expert.list.em, penalty, hyper.alpha, hyper.params)
        loglik.em.before.kj = ll.list.em$ll

        params.em.kj.old = params.em[[k]][[j]] # Before update

        params.em[[k]][[j]] = comp.kj.params.m.recur(comp.kj.dist, comp.kj.params.old,
                                                     tl.k, yl.k, yu.k, tu.k,
                                                     comp.kj.pos.expert.ll, comp.kj.pos.expert.tn, comp.kj.pos.expert.tn.bar,
                                                     z.pos.e.obs, z.pos.e.lat, comp.zkz.e.list$k.e,
                                                     penalty, hyper.params[[k]][[j]])

        params.em.kj.new = params.em[[k]][[j]] # After update

        # gate.ll.em = gate.logit(X, alpha.em)
        # expert.list.em = expert.loglik.dim.comp(Y, comp.dist, zero.em, params.em) # params.em[[k]][[j]] has been updated
        # ll.list.em.temp = gate.expert.loglik(alpha.em, gate.ll.em, expert.list.em, penalty, hyper.alpha, hyper.params)
        # loglik.em.temp = ll.list.em.temp$ll
        #
        # if(loglik.em.temp < loglik.em){ # No update, if the loglik decreases due to numerical imprecision.
        #   params.em[[k]][[j]] = params.em.kj.old
        # }

        # no.jump.list = c("gamma", "ZI-gamma", "burr", "ZI-burr", "nbinom", "ZI-nbinom")
        no.jump.list = c("nbinom", "ZI-nbinom", "binom", "ZI-binom")

        if( grad.jump==TRUE & (iter%%grad.period)==0 & !(comp.kj.dist %in% no.jump.list) ){
          diff.log.params = log(params.em.kj.new) - log(params.em.kj.old)
          grad.length = length(grad.seq)
          ll.seq = rep(loglik.em, grad.length)

          for(t in 1:grad.length)
          {
            # Focus only on k-j-th parameter vector
            params.temp = params.em
            params.temp[[k]][[j]] = params.em.kj.old * exp(grad.seq[t]*diff.log.params)

            gate.ll.em = GateLogit(X, alpha.em)
            expert.list.em = DimCompExpertLL(Y, comp.dist, zero.em, params.temp) # params.em[[k]][[j]] has been updated with gradient
            ll.list.em.temp = GateExpertLL(alpha.em, gate.ll.em, expert.list.em, penalty, hyper.alpha, hyper.params)
            ll.seq[t] = ifelse(ll.list.em.temp$ll>0, NaN, ll.list.em.temp$ll) # Prevent spurious loglik
          }

          grad.rate = grad.seq[which.max(ll.seq)]

          params.em[[k]][[j]]  = params.em.kj.old * exp(grad.rate*diff.log.params)
        }

        gate.ll.em = GateLogit(X, alpha.em)
        expert.list.em = DimCompExpertLL(Y, comp.dist, zero.em, params.em) # params.em[[k]][[j]] has been updated
        ll.list.em.temp = GateExpertLL(alpha.em, gate.ll.em, expert.list.em, penalty, hyper.alpha, hyper.params)
        loglik.em.after.kj = ll.list.em.temp$ll

        if(print){ print(paste("Dim ", k, " Comp ", j, " : ", loglik.em.before.kj, "  ", loglik.em.after.kj)); }

        if(loglik.em.after.kj < loglik.em.before.kj){ # No update, if the loglik decreases due to numerical imprecision.
          if(print){
            print(paste("No update on component: Dim ", k, " Comp ", j, " !" ));
            print(paste("Decrease in loglik = ", loglik.em.before.kj - loglik.em.after.kj, " or ", (loglik.em.before.kj - loglik.em.after.kj)/abs(loglik.em.before.kj)*100, "%"));
            print(paste("Intended Update: ", params.em[[k]][[j]]))# ; print (cat(params.em[[k]][[j]], sep = " "))
          }
          params.em[[k]][[j]] = params.em.kj.old
        }

        if(is.na(loglik.em.after.kj)){ # No update, if the loglik becomes NA
          if(print){
            print(paste("No update on component: Dim ", k, " Comp ", j, " !" ));
            print("Loglik becomes NA!");
            print(paste("Intended Update: ", params.em[[k]][[j]]))# ; print (cat(params.em[[k]][[j]], sep = " "))
          }
          params.em[[k]][[j]] = params.em.kj.old
        }

        if(loglik.em.after.kj>0 | is.infinite(loglik.em.after.kj)){ # No update, if the loglik becomes Inf. (Haven't figured out why)
          if(print){
            print(paste("No update on component: Dim ", k, " Comp ", j, " !" ));
            print(paste("Loglik becomes invalid! ", loglik.em.after.kj));
            print(paste("Intended Update: ", params.em[[k]][[j]]))# ; print (cat(params.em[[k]][[j]], sep = " "))
          }
          params.em[[k]][[j]] = params.em.kj.old
        }

      }
    }

    # Update loglik after one run of ECM
    gate.ll.em = GateLogit(X, alpha.em)
    expert.list.em = DimCompExpertLL(Y, comp.dist, zero.em, params.em)
    ll.list.em = GateExpertLL(alpha.em, gate.ll.em, expert.list.em, penalty, hyper.alpha, hyper.params) # This copy goes into EM
    loglik.em = ll.list.em$ll
    loglik.em.np = ll.list.em$ll.np
    if(print) {if(grad.jump==TRUE & (iter%%grad.period)==0){print("[[[ ! GRADIENT JUMP ! ]]]")};
      print(paste("Iter:", iter, "(Update of component params). Loglik = ", round(loglik.em, digits=6)));
      print(zero.em); print(params.em)}

  }

  # Information criteria
  AIC = -2*loglik.em.np + 2*(count.alpha(alpha.init) + count.zero(comp.dist) + count.pos.params(comp.dist))
  BIC = -2*loglik.em.np + log(sample.size.n)*(count.alpha(alpha.init) + count.zero(comp.dist) + count.pos.params(comp.dist))

  # Naming all output
  # rownames(zero.init) = rownames(zero.em) = paste("dim", c(1:dim.m), sep = " ")
  # rownames(zero.init) = paste("dim", c(1:dim.m), sep = " ")
  rownames(comp.dist) = paste("dim", c(1:dim.m), sep = " ")
  # colnames(zero.init) = colnames(zero.em) = paste("comp", c(1:n.comp), sep = " ")
  # colnames(zero.init) =  paste("comp", c(1:n.comp), sep = " ")
  rownames(alpha.init) = rownames(alpha.em) = colnames(comp.dist) = paste("comp", c(1:n.comp), sep = " ")
  if(!is.null(colnames(X))){
    colnames(alpha.init) = colnames(alpha.em) = colnames(X)
  }
  for(k in 1:dim.m){
    names(params.init)[k] = names(params.em)[k] = paste("dim", k, sep = " ")
    for(j in 1:n.comp){
      names(params.init[[k]])[j] = names(params.em[[k]])[j] = paste("comp", j, sep = " ")
      names(params.init[[k]][[j]]) = names(params.em[[k]][[j]]) = PosName(comp.dist[k,j])
    }
  }


  # Return results
  return(list(
    n.comp = n.comp, comp.dist = comp.dist,
    alpha.init = alpha.init, zero.init = zero.init, params.init = params.init,
    alpha.fit = alpha.em, zero.fit = zero.em, params.fit = params.em,
    ll = loglik.em, ll.np = loglik.em.np,
    AIC = AIC, BIC = BIC,
    iter = iter
  ))

}
