# context("EMEzzeroobs")
# library(Rcpp)
# library(copula)
#
# # ll.ind = matrix(runif(250000, -10, 10), ncol = 5)
# # comp.aggre.ll.ind = matrix(runif(50000, -10, 10), ncol = 1)
# # ll.tn.bar = matrix(runif(250000, -10, 10), ncol = 5)
# # comp.aggre.ll.tn.bar = matrix(runif(50000, -10, 10), ncol = 1)
#
# yl.k = rpois(100000, 5)
# comp.kj.zero.inflation = T
# comp.kj.zero.prob.old = 0.30
# comp.kj.pos.expert.ll = runif(100000, -10, 10)
#
# oldver = function(yl.k, comp.kj.zero.inflation, comp.kj.zero.prob.old, comp.kj.pos.expert.ll)
# {
#   sample.size.n = length(yl.k)
#   temp = array(0, dim = c(sample.size.n, 1))
#
#   if(comp.kj.zero.inflation == TRUE) # Otherwise, there is no possibility of zero component.
#   {
#     obs.y.zero.idx = (yl.k==0) # yl=0: possible for observation to be from zero component
#     # Update z.zero only for those indices. Otherwise, they stay zero.
#     temp[obs.y.zero.idx] = exp( - log1pexp( log(1/comp.kj.zero.prob.old-1) + comp.kj.pos.expert.ll[obs.y.zero.idx]) )
#       # comp.kj.zero.prob.old / (comp.kj.zero.prob.old + (1-comp.kj.zero.prob.old)*exp(comp.kj.pos.expert.ll[obs.y.zero.idx]))
#   }
#
#   return(temp)
# }
#
# cppFunction("
# SEXP tempf(SEXP yl, SEXP ll, SEXP pr) {
#   NumericVector yyl(yl);
#   NumericVector posl(ll);
#   NumericVector ppr(pr);
#   double pp = ppr(0);
#  NumericVector result = ifelse(yyl == 0.0, pp/(pp+(1-pp)*exp(posl)), 0.0);
#
#   return result;
# }")
#
# test_that("EMEzzeroobs", {
#   expect_equal(tempf(yl.k, comp.kj.pos.expert.ll, comp.kj.zero.prob.old), c(oldver(yl.k, comp.kj.zero.inflation, comp.kj.zero.prob.old, comp.kj.pos.expert.ll)) )
# })
#
# library(microbenchmark)
# bench = microbenchmark(
#   tempf(yl.k, comp.kj.pos.expert.ll, comp.kj.zero.prob.old),
#   oldver(yl.k, comp.kj.zero.inflation, comp.kj.zero.prob.old, comp.kj.pos.expert.ll),
# times = 500
# )
#
# library(ggplot2)
# autoplot(bench)
#
# # No significant improvement!
