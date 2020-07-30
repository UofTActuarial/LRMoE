context("EMMzero")
library(Rcpp)
# library(RcppEigen)
# library(inline)
library(matrixStats)

zobs = matrix(runif(50000, 0, 1), ncol = 1)
ke = matrix(runif(50000, 0, 1), ncol = 1)
zlat = matrix(runif(50000, 0, 1), ncol = 1)

oldver = function(z.zero.e.obs, z.pos.e.obs, z.zero.e.lat, z.pos.e.lat, k.e)
{
  term.zero = sweep(matrix(z.zero.e.obs), 1, sweep(matrix(z.zero.e.lat), 1, matrix(k.e), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE)
  # z.zero.e.obs  + sweep(z.zero.e.lat, 1, k.e, FUN = "*", check.margin = FALSE)
  term.pos  = sweep(matrix(z.pos.e.obs), 1, sweep(matrix(z.pos.e.lat), 1, matrix(k.e), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE)
  # z.pos.e.obs   + sweep(z.pos.e.lat, 1, k.e, FUN = "*", check.margin = FALSE)

  numerator   = apply(term.zero, 2, sum)
  denominator = numerator + apply(term.pos, 2, sum) # apply(term.zero, 2, sum) + apply(term.pos, 2, sum)
  return( numerator / denominator )
}


cppFunction("
SEXP XPlusYColTimesZ(SEXP x, SEXP y, SEXP z) {
  NumericMatrix xx(x);
  NumericMatrix yy(y);
  NumericVector zz(z);
  NumericMatrix result(xx.nrow(), xx.ncol());

  for(int j=0; j<xx.ncol(); j++){
    result(_,j) = xx(_,j) + yy(_,j)*zz;
  }
  return result;
}")

test_that("EMMzero", {
  expect_equal(XPlusYColTimesZ(zobs, zlat, ke),
               sweep(matrix(zobs), 1, sweep(matrix(zlat), 1, matrix(ke), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE))
})

# library(microbenchmark)
# bench = microbenchmark(
#   XPlusYColTimesZ(zobs, zlat, ke),
#   sweep(matrix(zobs), 1, sweep(matrix(zlat), 1, matrix(ke), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE)
#
# )
#
# library(ggplot2)
# autoplot(bench)

temp1 = XPlusYColTimesZ(zobs, zlat, ke)
temp2 = sweep(matrix(zobs), 1, sweep(matrix(zlat), 1, matrix(ke), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE)

test_that("EMMzero", {
  expect_equal(sum(temp1),
               apply(temp1, 2, sum))
  expect_equal(sum(temp2),
               apply(temp2, 2, sum))
})

# library(microbenchmark)
# bench = microbenchmark(
#   sum(temp1),
#   apply(temp1, 2, sum)
# )
#
# library(ggplot2)
# autoplot(bench)


# Very significant improvement!
