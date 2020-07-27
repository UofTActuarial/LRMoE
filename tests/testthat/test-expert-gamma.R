context("ExpertGamma")
library(copula)

y = rgamma(100000, shape = 5, scale = 2)
tl = y*0.50
yl = y*0.75
yu = y*1.25
tu = y*1.50

expert.gamma = function(tl, yl, yu, tu, g = 1, m, theta)
{
  # Initialization: return value are N * g matrices
  expert.gamma.ll=expert.gamma.tn=expert.gamma.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=pgamma(yu[censor.idx],shape=m[j],scale=theta[j],log.p=TRUE)
    prob.log.yl=pgamma(yl[censor.idx],shape=m[j],scale=theta[j],log.p=TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.gamma.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.gamma.ll[!censor.idx,j]=dgamma(yu[!censor.idx],shape=m[j],scale=theta[j],log=TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=pgamma(tu,shape=m[j],scale=theta[j],log.p=TRUE)
    prob.log.tl=pgamma(tl,shape=m[j],scale=theta[j],log.p=TRUE)

    # Normalizing factor for truncation limits, in log
    expert.gamma.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.gamma.tn[no.trunc.idx,j] = dgamma(tu[no.trunc.idx],shape=m[j],scale=theta[j],log=TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    zero.idx = (tu==0)
    expert.gamma.ll[zero.idx,j]=(-Inf)
    expert.gamma.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.gamma.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.gamma.tn[!no.trunc.idx,j])
    expert.gamma.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.gamma.ll=expert.gamma.ll, expert.gamma.tn=expert.gamma.tn, expert.gamma.tn.bar=expert.gamma.tn.bar)
}

test_that("Gamma Expert", {
  expect_equal(ExpertGamma(tl, yl, yu, tu, 5, 2)[[1]], expert.gamma(tl, yl, yu, tu, 1, 5, 2)[[1]])
  expect_equal(ExpertGamma(tl, yl, yu, tu, 5, 2)[[2]], expert.gamma(tl, yl, yu, tu, 1, 5, 2)[[2]])
  expect_equal(ExpertGamma(tl, yl, yu, tu, 5, 2)[[3]], expert.gamma(tl, yl, yu, tu, 1, 5, 2)[[3]])
})

# library(microbenchmark)
# bench = microbenchmark(
#   ExpertGamma(tl, yl, yu, tu, 5, 2),
#   expert.gamma(tl, yl, yu, tu, 1, 5, 2)
# )
#
# library(ggplot2)
# autoplot(bench)
