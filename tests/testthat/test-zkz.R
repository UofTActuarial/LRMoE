context("EMEzkz")
library(Rcpp)

ll.ind = matrix(runif(250000, -10, 10), ncol = 5)
comp.aggre.ll.ind = matrix(runif(50000, -10, 10), ncol = 1)
ll.tn.bar = matrix(runif(250000, -10, 10), ncol = 5)
comp.aggre.ll.tn.bar = matrix(runif(50000, -10, 10), ncol = 1)

cppFunction("
SEXP tempf(SEXP x, SEXP y) {
  NumericMatrix xx(x);
  NumericVector yy(y);
  NumericMatrix result(xx.nrow(), xx.ncol());

  for(int j=0; j<xx.ncol(); j++){
    result(_,j) = xx(_,j) - yy;
  }
  return result;
}")

test_that("EMEzkz", {
  expect_equal(exp(tempf(ll.ind, comp.aggre.ll.ind)), exp( sweep(ll.ind, 1, comp.aggre.ll.ind, FUN = "-", check.margin = FALSE) ))
  expect_equal(exp(tempf(ll.tn.bar, comp.aggre.ll.tn.bar)), exp( sweep(ll.tn.bar, 1, comp.aggre.ll.tn.bar, FUN = "-", check.margin = FALSE) ))
  expect_equal(exp(-comp.aggre.ll.tn.bar)-1, expm1(-comp.aggre.ll.tn.bar))
})

# library(microbenchmark)
# bench = microbenchmark(
#   exp(tempf(ll.ind, comp.aggre.ll.ind)),
#   exp( sweep(ll.ind, 1, comp.aggre.ll.ind, FUN = "-", check.margin = FALSE) )
# )
#
# library(ggplot2)
# autoplot(bench)

# Very significant improvement!
