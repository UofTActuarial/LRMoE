context("ExpertLognormal")
library(copula)

y = rlnorm(100000, 0, 1)
tl = y*0.50
yl = y*0.75
yu = y*1.25
tu = rep(Inf, 100000)

expert.lognormal = function(tl, yl, yu, tu, g = 1, meanlog, sdlog)
{
  # Initialization: return value are N * g matrices
  expert.lognormal.ll=expert.lognormal.tn=expert.lognormal.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=plnorm(yu[censor.idx],meanlog = meanlog[j],sdlog = sdlog[j],log.p=TRUE)
    prob.log.yl=plnorm(yl[censor.idx],meanlog = meanlog[j],sdlog = sdlog[j],log.p=TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.lognormal.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.lognormal.ll[!censor.idx,j]=dlnorm(yu[!censor.idx],meanlog = meanlog[j],sdlog = sdlog[j],log=TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=plnorm(tu,meanlog = meanlog[j],sdlog = sdlog[j],log.p=TRUE)
    prob.log.tl=plnorm(tl,meanlog = meanlog[j],sdlog = sdlog[j],log.p=TRUE)

    # Normalizing factor for truncation limits, in log
    expert.lognormal.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.lognormal.tn[no.trunc.idx,j] = dlnorm(tu[no.trunc.idx],meanlog = meanlog[j],sdlog = sdlog[j],log=TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    zero.idx = (tu==0)
    expert.lognormal.ll[zero.idx,j]=(-Inf)
    expert.lognormal.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.lognormal.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.lognormal.tn[!no.trunc.idx,j])
    expert.lognormal.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.lognormal.ll=expert.lognormal.ll, expert.lognormal.tn=expert.lognormal.tn, expert.lognormal.tn.bar=expert.lognormal.tn.bar)
}

test_that("Lornormal Expert", {
  expect_equal(ExpertLognormal(tl, yl, yu, tu, 0, 1)[[1]], expert.lognormal(tl, yl, yu, tu, 1, 0, 1)[[1]])
  expect_equal(ExpertLognormal(tl, yl, yu, tu, 0, 1)[[2]], expert.lognormal(tl, yl, yu, tu, 1, 0, 1)[[2]])
  expect_equal(ExpertLognormal(tl, yl, yu, tu, 0, 1)[[3]], expert.lognormal(tl, yl, yu, tu, 1, 0, 1)[[3]])
})

# library(microbenchmark)
# bench = microbenchmark(
#   ExpertLognormal(tl, yl, yu, tu, 0, 1),
#   expert.lognormal(tl, yl, yu, tu, 1, 0, 1)
# )
#
# library(ggplot2)
# autoplot(bench)
