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

data("LRMoEDemoData")
head(X)

alpha = matrix(c(0.5, 0.25, -0.05, 0.3, -0.2,
                 0, 0, 0, 0, 0),
               nrow = 2, byrow = T)

head(exp(GateLogit(X, alpha)))
hist(exp(GateLogit(X, alpha))[,1])

comp.dist = matrix(c("lnorm", "ZI-lnorm"),
                   nrow = 1, byrow = T)

zero.prob = matrix(c(0, 0.10),
                   nrow = 1, byrow = T)

params.list = list(list(c(3, 0.25), c(1, 0.5)))

simy = SimYSet(X, alpha, comp.dist, zero.prob, params.list)

hist(simy, breaks = 200, xlim = c(0, 50))

YY = cbind(rep(0, nrow(X)), simy, simy, rep(Inf, nrow(X)))

YY[c(251:500),2] = simy[c(250:500)] * 0.75
YY[c(501:750),3] = simy[c(250:500)] * 1.50

YY[c(751:1000),1] = simy[c(751:1000)] * 0.25

YY[c(1001:1250),4] = simy[c(1001:1250)] * 2

YY[c(1251:1500),1] = simy[c(1251:1500)] * 0.25
YY[c(1251:1500),2] = simy[c(1251:1500)] * 0.75
YY[c(1251:1500),3] = simy[c(1251:1500)] * 1.25
YY[c(1251:1500),4] = simy[c(1251:1500)] * 2


alpha.init = alpha = matrix(c(0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0),
                            nrow = 2, byrow = T)
comp.dist = comp.dist
zero.init = matrix(c(0, 0.50),
                   nrow = 1, byrow = T)
params.init = list(list(c(2.5, 0.1), c(2, 0.8)))

modelfit = LRMoEFit(YY, X, 2, comp.dist, alpha.init, zero.init, params.init)
