## Gating Function: Multiple Logit
#' Computes the logit regression weights in log.
#'
#' @param x An N * P covariate matrix, where N is sample size. The first column MUST be 1.
#' @param alpha A g * P matrix. Logit regression coefficients.
#' @return The log of \code{tcrossprod(x,alpha)}, with each row normalized by \code{rowLogSumExps}
#'
#' @importFrom matrixStats rowLogSumExps
#`
#' @keywords internal
#'
#' @export GateLogit
GateLogit = function(x, alpha)
{
  gate.body=tcrossprod(x,alpha) # = x %*% t(alpha)
  # return(sweep(gate.body, 1, rowLogSumExps(gate.body), FUN = "-", check.margin = FALSE))
  return(XColMinusY(gate.body, rowLogSumExps(gate.body)))
}
