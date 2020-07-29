context("EMMalpha")
library(Rcpp)
# library(RcppEigen)
# library(inline)
library(matrixStats)

zobs = matrix(runif(300000, 0, 1), ncol = 6)
ke = matrix(runif(50000, 0, 1), ncol = 1)
zlat = matrix(runif(300000, 0, 1), ncol = 6)

X = matrix(runif(300000, 0, 1), ncol = 6)
alpha = matrix(runif(30, 0, 1), ncol = 6)

comp.zpzk = zobs + sweep(zlat, 1, ke, FUN = "*", check.margin = FALSE)
comp.zpzk.sum = apply(comp.zpzk, 1, sum)
gate.body=tcrossprod(X, alpha)
pp = exp(gate.body-rowLogSumExps(gate.body))

# Old dQ
# apply(sweep(X,1,comp.zpzk[,j]-comp.zpzk.marg*exp(gate.body[,j]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE),
#       2,sum)
# apply(sweep(X,1,comp.zpzk[,1]-comp.zpzk.sum*exp(gate.body[,1]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), 2,sum)

cppFunction("
SEXP tempdQ(SEXP x, SEXP zj, SEXP z, SEXP p) {
  NumericMatrix xx(x);
  NumericVector zzj(zj);
  NumericVector zz(z);
  NumericVector pp(p);
  NumericVector result(xx.ncol());

  for(int j=0; j<xx.ncol(); j++){
    result(j) = sum(xx(_,j)*(zzj - zz*pp));
  }
  return result;

}")

test_that("EMMalpha", {
  expect_equal(tempdQ(X, comp.zpzk[,1], comp.zpzk.sum, pp[,1]), apply(sweep(X,1,comp.zpzk[,1]-comp.zpzk.sum*exp(gate.body[,1]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), 2,sum))
  expect_equal(tempdQ(X, comp.zpzk[,2], comp.zpzk.sum, pp[,2]), apply(sweep(X,1,comp.zpzk[,2]-comp.zpzk.sum*exp(gate.body[,2]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), 2,sum))
  expect_equal(tempdQ(X, comp.zpzk[,3], comp.zpzk.sum, pp[,3]), apply(sweep(X,1,comp.zpzk[,3]-comp.zpzk.sum*exp(gate.body[,3]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), 2,sum))
  expect_equal(tempdQ(X, comp.zpzk[,4], comp.zpzk.sum, pp[,4]), apply(sweep(X,1,comp.zpzk[,4]-comp.zpzk.sum*exp(gate.body[,4]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), 2,sum))
  expect_equal(tempdQ(X, comp.zpzk[,5], comp.zpzk.sum, pp[,5]), apply(sweep(X,1,comp.zpzk[,5]-comp.zpzk.sum*exp(gate.body[,5]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), 2,sum))
})

# library(microbenchmark)
# bench = microbenchmark(
#   tempdQ(X, comp.zpzk[,1], comp.zpzk.sum, pp[,1]),
#   apply(sweep(X,1,comp.zpzk[,1]-comp.zpzk.sum*exp(gate.body[,1]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), 2,sum)
# )
#
# library(ggplot2)
# autoplot(bench)

# Very significant improvement!

# Old dQ2
# -crossprod(sweep(X,1,comp.zpzk.marg*exp(rowLogSumExps(array(gate.body[,-j],dim=c(sample.size.n,n.comp-1)))+gate.body[,j]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE),X)
# -crossprod(sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-j],dim=c(sample.size.n,n.comp-1)))+gate.body[,j]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), X)

cppFunction("
SEXP tempdQ2left(SEXP x, SEXP z, SEXP p, SEXP q) {
  NumericMatrix xx(x);
  NumericVector zz(z);
  NumericVector pp(p);
  NumericVector qq(q);
  NumericMatrix result(xx.nrow(), xx.ncol());

  for(int j=0; j<xx.ncol(); j++){
    result(_,j) = xx(_,j) * (zz*pp*qq);
  }
  return result;

}")

test_that("EMMalpha", {
  expect_equal(tempdQ2left(X, comp.zpzk.sum, pp[,1], exp(rowLogSumExps(array(gate.body[,-1],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-1],dim=c(50000,5-1)))+gate.body[,1]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE))

  expect_equal(tempdQ2left(X, comp.zpzk.sum, pp[,2], exp(rowLogSumExps(array(gate.body[,-2],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-2],dim=c(50000,5-1)))+gate.body[,2]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE))

  expect_equal(tempdQ2left(X, comp.zpzk.sum, pp[,3], exp(rowLogSumExps(array(gate.body[,-3],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-3],dim=c(50000,5-1)))+gate.body[,3]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE))

  expect_equal(tempdQ2left(X, comp.zpzk.sum, pp[,4], exp(rowLogSumExps(array(gate.body[,-4],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-4],dim=c(50000,5-1)))+gate.body[,4]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE))

  expect_equal(tempdQ2left(X, comp.zpzk.sum, pp[,5], exp(rowLogSumExps(array(gate.body[,-5],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-5],dim=c(50000,5-1)))+gate.body[,5]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE))
})

# library(microbenchmark)
# qq1 = exp(rowLogSumExps(array(gate.body[,-1],dim=c(50000,5-1)))-rowLogSumExps(gate.body))
# qqq1 = exp(rowLogSumExps(array(gate.body[,-1],dim=c(50000,5-1)))+gate.body[,1]-2*rowLogSumExps(gate.body))
# bench = microbenchmark(
#   tempdQ2left(X, comp.zpzk.sum, pp[,1], qq1),
#   sweep(X,1,comp.zpzk.sum*qqq1,FUN="*",check.margin=FALSE)
# )
#
# library(ggplot2)
# autoplot(bench)

# Very significant improvement!

cppFunction(
"
SEXP tempdQ2(SEXP x, SEXP z, SEXP p, SEXP q) {
  NumericMatrix xx(x);
  NumericVector zz(z);
  NumericVector pp(p);
  NumericVector qq(q);
  NumericMatrix temp(xx.nrow(), xx.ncol());
  NumericMatrix result(xx.ncol(), xx.ncol());

  // for(int j=0; j<xx.ncol(); j++){
  //   temp(_,j) = xx(_,j) * (zz*pp*qq);
  // }

  for(int i=0; i<result.nrow(); i++){
    temp(_,i) = xx(_,i) * (zz*pp*qq);
    for(int j=0; j<result.ncol(); j++){
      result(i,j) = -1.0* sum(temp(_,i)*xx(_,j));
    }
  }

  return result;

}")

test_that("EMMalpha", {
  expect_equal(tempdQ2(X, comp.zpzk.sum, pp[,1], exp(rowLogSumExps(array(gate.body[,-1],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               -crossprod(sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-1],dim=c(50000,5-1)))+gate.body[,1]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), X))

  expect_equal(tempdQ2(X, comp.zpzk.sum, pp[,2], exp(rowLogSumExps(array(gate.body[,-2],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               -crossprod(sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-2],dim=c(50000,5-1)))+gate.body[,2]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), X))

  expect_equal(tempdQ2(X, comp.zpzk.sum, pp[,3], exp(rowLogSumExps(array(gate.body[,-3],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               -crossprod(sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-3],dim=c(50000,5-1)))+gate.body[,3]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), X))

  expect_equal(tempdQ2(X, comp.zpzk.sum, pp[,4], exp(rowLogSumExps(array(gate.body[,-4],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               -crossprod(sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-4],dim=c(50000,5-1)))+gate.body[,4]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), X))

  expect_equal(tempdQ2(X, comp.zpzk.sum, pp[,5], exp(rowLogSumExps(array(gate.body[,-5],dim=c(50000,5-1)))-rowLogSumExps(gate.body))),
               -crossprod(sweep(X,1,comp.zpzk.sum*exp(rowLogSumExps(array(gate.body[,-5],dim=c(50000,5-1)))+gate.body[,5]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE), X))
})

# library(microbenchmark)
# qq1 = exp(rowLogSumExps(array(gate.body[,-1],dim=c(50000,5-1)))-rowLogSumExps(gate.body))
# qqq1 = exp(rowLogSumExps(array(gate.body[,-1],dim=c(50000,5-1)))+gate.body[,1]-2*rowLogSumExps(gate.body))
# bench = microbenchmark(
#   tempdQ2(X, comp.zpzk.sum, pp[,1], qq1),
#   -crossprod(sweep(X,1,comp.zpzk.sum*qqq1,FUN="*",check.margin=FALSE), X)
# )
#
# library(ggplot2)
# autoplot(bench)

# Very significant improvement!
