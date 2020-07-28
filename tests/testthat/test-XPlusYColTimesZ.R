context("XPlusYColTimesZ")
library(Rcpp)

zobs = matrix(runif(250000, 0, 1), ncol = 5)
ke = matrix(runif(50000, 0, 1), ncol = 1)
zlat = matrix(runif(250000, 0, 1), ncol = 5)

# zobs + sweep(zlat, 1, ke, FUN = "*", check.margin = FALSE)

cppFunction("
SEXP tempf(SEXP x, SEXP y, SEXP z) {
  NumericMatrix xx(x);
  NumericMatrix yy(y);
  NumericVector zz(z);
  NumericMatrix result(xx.nrow(), xx.ncol());

  for(int j=0; j<xx.ncol(); j++){
    result(_,j) = xx(_,j) + yy(_,j)*zz;
  }
  return result;
}")

test_that("XPlusYColTimesZ", {
  expect_equal(tempf(zobs, zlat, ke), zobs + sweep(zlat, 1, ke, FUN = "*", check.margin = FALSE))
})

# library(microbenchmark)
# bench = microbenchmark(
#   tempf(zobs, zlat, ke),
#   zobs + sweep(zlat, 1, ke, FUN = "*", check.margin = FALSE)
# )
#
# library(ggplot2)
# autoplot(bench)

# Very significant improvement!
