context("ExpertInvgauss")
library(copula)
library(statmod)

y = statmod::rinvgauss(100000, 10, 10)
tl = y*0.50
yl = y*0.75
yu = y*1.25
tu = rep(Inf, 100000)

expert.invgauss = function(tl, yl, yu, tu, g = 1, mean.mu, shape.lamba)
{
  # Initialization: return value are N * g matrices
  expert.invgauss.ll=expert.invgauss.tn=expert.invgauss.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=statmod::pinvgauss(yu[censor.idx],mean = mean.mu[j],shape = shape.lamba[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.yl=statmod::pinvgauss(yl[censor.idx],mean = mean.mu[j],shape = shape.lamba[j], lower.tail = TRUE, log.p = TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.invgauss.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.invgauss.ll[!censor.idx,j]=statmod::dinvgauss(yu[!censor.idx],mean = mean.mu[j],shape = shape.lamba[j], log = TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=statmod::pinvgauss(tu,mean = mean.mu[j],shape = shape.lamba[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.tl=statmod::pinvgauss(tl,mean = mean.mu[j],shape = shape.lamba[j], lower.tail = TRUE, log.p = TRUE)

    # Normalizing factor for truncation limits, in log
    expert.invgauss.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.invgauss.tn[no.trunc.idx,j] = statmod::dinvgauss(tu[no.trunc.idx],mean = mean.mu[j],shape = shape.lamba[j], log = TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    zero.idx = (tu==0)
    expert.invgauss.ll[zero.idx,j]=(-Inf)
    expert.invgauss.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.invgauss.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.invgauss.tn[!no.trunc.idx,j])
    expert.invgauss.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.invgauss.ll=expert.invgauss.ll, expert.invgauss.tn=expert.invgauss.tn, expert.invgauss.tn.bar=expert.invgauss.tn.bar)
}

test_that("Invgauss Expert", {
  expect_equal(ExpertInvgauss(tl, yl, yu, tu, 10, 10)[[1]], expert.invgauss(tl, yl, yu, tu, 1, 10, 10)[[1]])
  expect_equal(ExpertInvgauss(tl, yl, yu, tu, 10, 10)[[2]], expert.invgauss(tl, yl, yu, tu, 1, 10, 10)[[2]])
  expect_equal(ExpertInvgauss(tl, yl, yu, tu, 10, 10)[[3]], expert.invgauss(tl, yl, yu, tu, 1, 10, 10)[[3]])
})
#
# library(microbenchmark)
# bench = microbenchmark(
#   ExpertInvgauss(tl, yl, yu, tu, 10, 10),
#   expert.invgauss(tl, yl, yu, tu, 1, 10, 10)
# )
#
# library(ggplot2)
# autoplot(bench)
