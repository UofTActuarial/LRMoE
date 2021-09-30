#' @title Compute f(a) = log(1 +/- exp(-a)) Numerically Optimally
#' @name log1mexp
#'
#' @description
#' Compute f(a) = log(1 +/- exp(-a)) Numerically Optimally
#' This helper function was a directly copy from `copula` package. It was initialized here intend to maintain the minimum package complexity
#'
#' @param a \cr
#' numeric vector of positive values
#'
#' @param cutoff \cr
#' log(2) is “optimal”, but the exact value is unimportant, and anything in (0.5, 1) is fine.
log1mexp <- function(a, cutoff = log(2))
{
  if (has.na <- any(ina <- is.na(a))) {
    y <- a
    a <- a[ok <- !ina]
  }
  if (any(a < 0))
    warning("'a' >= 0 needed")
  tst <- a <= cutoff
  r <- a
  r[tst] <- log(-expm1(-a[tst]))
  r[!tst] <- log1p(-exp(-a[!tst]))
  if (has.na) {
    y[ok] <- r
    y
  }
  else r
}

# The helper function to avoid duplicate code in faster_zi_result
faster_zi_process <- function(condition, p0, dist_result) {
  result = ifelse(condition, log(p0 + (1-p0)*exp(dist_result)), log(0.0 + (1-p0)*exp(dist_result)))
  return( result )
}

# The helper function to avoid duplicate code in ZI distribution Scirpts. (ll_not_exact)
faster_zi_result <- function(tl, tu, yl, yu, params, distribution) {
  # Call the None ZI Distribution.expert_ll_not_exact
  p0 = params[["p_zero"]]
  result = do.call(paste0(distribution, ".expert_ll_not_exact"), list(tl = tl, tu = tu, yl = yl, yu = yu, params = params))
  # Zero Inflation Processing
  ######################################################################################
  expert_ll = faster_zi_process(yl == 0, p0, result[["expert_ll"]])
  expert_tn = faster_zi_process(tl == 0, p0, result[["expert_tn"]])
  expert_tn_bar = faster_zi_process(tl > 0, p0, result[["expert_tn_bar"]])
  ######################################################################################

  # Return values
  return( list(expert_ll = expert_ll, expert_tn = expert_tn, expert_tn_bar = expert_tn_bar) )
}

# The helper function that calculate the incomplete beta function
beta_inc <- function(x, a, b) {
  pbeta(x,a,b)*beta(a,b)
}

# The print helper function for named list
print_nl <- function(nl) {
  paste(names(nl),nl,sep="=",collapse=";  " )
}

# The print helper function for vectors
print_vc <- function(vc) {
  vec_str = paste(vc,collapse=", ")
  paste("c(", vec_str, ")")
}

# The helper function to count the number of parameters in alpha
count_alpha <- function(alpha) {
  return((nrow(alpha)-1)*ncol(alpha))
}
