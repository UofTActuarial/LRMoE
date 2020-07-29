## M-Step for alpha
#' ECM: M-Step for logit regression coefficients \code{alpha}.
#'
#' @param X A N*P matrix of numerical covariates.
#' @param alpha A g*P matrix of old logit regression coefficients.
#' @param comp.zkz.e.list An object returned by \code{\link{EMEzkz}}.
#' @param alpha.iter.max Numeric: maximum number of iterations.
#' @param penalty TRUE/FALSE, which indicates whether penalty is applied.
#' @param hyper.alpha A numeric of penalty applied to \code{alpha}.
#'
#' @return \code{alpha.new} Updated logit regression coefficients.
#'
#' @importFrom matrixStats rowLogSumExps
#'
#' @keywords internal
#'
#' @export EMMalpha
EMMalpha = function(X, alpha, comp.zkz.e.list, alpha.iter.max,
                         penalty, hyper.alpha)
{
  comp.zpzk = XPlusYColTimesZ(comp.zkz.e.list$z.e.obs, comp.zkz.e.list$z.e.lat, comp.zkz.e.list$k.e)
    # comp.zkz.e.list$z.e.obs + sweep(comp.zkz.e.list$z.e.lat, 1, comp.zkz.e.list$k.e, FUN = "*", check.margin = FALSE)

  sample.size.n = nrow(X)
  n.covar.p = ncol(X)
  n.comp = nrow(alpha)
  iter=array(0, dim = c(n.comp, 1))
  alpha.new = alpha
  alpha.old = alpha - Inf
  comp.zpzk.marg = apply(comp.zpzk, 1, sum)

  for (j in 1:(n.comp-1)) # The last component's alpha's are always kept at zero (reference category).
  {
    while ((iter[j]<=alpha.iter.max)&(sum((alpha.old[j,]-alpha.new[j,])^2)>10^(-8))) # Stopping criteria: (alpha.iter.max) iterations, or small difference
    {
      alpha.old[j,]=alpha.new[j,]
      gate.body=tcrossprod(X,alpha.new)
      pp = exp(gate.body-rowLogSumExps(gate.body))
      qqj = exp(rowLogSumExps(array(gate.body[,-j],dim=c(sample.size.n,n.comp-1)))-rowLogSumExps(gate.body))
      dQ = EMalphadQ(X, comp.zpzk[,j], comp.zpzk.marg, pp[,j]) - if(penalty){alpha.new[j,]/hyper.alpha^2} else{0}
        # apply(sweep(X,1,comp.zpzk[,j]-comp.zpzk.marg*exp(gate.body[,j]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE),2,sum)-if(penalty){alpha.new[j,]/hyper.alpha^2} else{0}
      dQ2 = EMalphadQ2(X, comp.zpzk.marg, pp[,j], qqj) - if(penalty){diag(1/hyper.alpha^2,nrow = n.covar.p, ncol = n.covar.p)} else{diag(10^(-7),nrow = n.covar.p, ncol = n.covar.p)}
        # -crossprod(sweep(X,1,comp.zpzk.marg*exp(rowLogSumExps(array(gate.body[,-j],dim=c(sample.size.n,n.comp-1)))+gate.body[,j]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE),X)-if(penalty){diag(1/hyper.alpha^2,nrow = n.covar.p, ncol = n.covar.p)} else{diag(10^(-7),nrow = n.covar.p, ncol = n.covar.p)}

      alpha.new[j,]=alpha.new[j,]-crossprod(dQ,solve(dQ2))
      iter[j] = iter[j]+1
    }
  }

  return(alpha.new)
}
