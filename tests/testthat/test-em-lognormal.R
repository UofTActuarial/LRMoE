context("EMMLognormal")
library(Rcpp)


zobs = matrix(runif(50000, 0, 1), ncol = 1)
ke = matrix(runif(50000, 0, 1), ncol = 1)
zlat = matrix(runif(50000, 0, 1), ncol = 1)

Yobs = matrix(runif(50000, 0, 1), ncol = 1)
Ylat = matrix(runif(50000, 0, 1), ncol = 1)

# z.e.obs[pos.idx] * y.log.e.obs[pos.idx] + k.e[pos.idx] * z.e.lat[pos.idx] * y.log.e.lat[pos.idx]

cppFunction("
SEXP tempf(SEXP x, SEXP a, SEXP y, SEXP z, SEXP b) {
  NumericVector xx(x);
  NumericVector yy(y);
  NumericVector zz(z);
  NumericVector aa(a);
  NumericVector bb(b);
  NumericVector result(xx.length());

  // for(int j=0; j<xx.ncol(); j++){
  //   result(_,j) = xx(_,j)*aa + yy(_,j)*zz*bb;
  // }

  result = xx*aa + yy*zz*bb;
  return result;
}")

temp1 = tempf(zobs, Yobs, zlat, ke, Ylat)
temp2 = zobs * Yobs + zlat * ke * Ylat
test_that("EMMLognormal", {
  expect_equal(tempf(zobs, Yobs, zlat, ke, Ylat), c(zobs * Yobs + zlat * ke * Ylat))
})

# library(microbenchmark)
# bench = microbenchmark(
#   tempf(zobs, Yobs, zlat, ke, Ylat),
#   zobs * Yobs + zlat * ke * Ylat
# )
#
# library(ggplot2)
# autoplot(bench)

# Slight improvement!



# data("LRMoEDemoData")
#
# head(Y.obs)
# summary(Y.obs)
# tempY = Y.obs[,5:8]
#
# meanlog = 2
# sdlog = 1
#
# censor.idx = (tempY[,2]!=tempY[,3])

