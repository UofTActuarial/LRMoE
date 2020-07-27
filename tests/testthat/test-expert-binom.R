context("ExpertBinomial")
library(copula)

y = rbinom(100000, 20, 0.5)
tl = y # floor(y*0.50)
yl = y # floor(y*0.75)
yu = y # ceiling(y*1.25)
tu = y # rep(Inf, 100000)

expert.binom = function(tl, yl, yu, tu, g = 1, size.n, prob.p)
{
  # Initialization: return value are N * g matrices
  expert.binom.ll=expert.binom.tn=expert.binom.tn.bar=array(-Inf, dim=c(length(yu),g))

  for (j in 1:g) # Loop trough each expert component
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=pbinom(yu[censor.idx], size = size.n[j], prob = prob.p[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.yl=pbinom(ceiling(yl[censor.idx])-1, size = size.n[j], prob = prob.p[j], lower.tail = TRUE, log.p = TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.binom.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.binom.ll[!censor.idx,j]=dbinom(yu[!censor.idx], size = size.n[j], prob = prob.p[j], log = TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=pbinom(tu, size = size.n[j], prob = prob.p[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.tl=pbinom(ceiling(tl)-1, size = size.n[j], prob = prob.p[j], lower.tail = TRUE, log.p = TRUE)

    # Normalizing factor for truncation limits, in log
    expert.binom.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.binom.tn[no.trunc.idx,j] = dbinom(tu[no.trunc.idx], size = size.n[j], prob = prob.p[j], log = TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case: The following code treating zeros is NOT applicable for frequency distributions!!!
    # zero.idx = (tu==0)
    # expert.binom.ll[zero.idx,j]=(-Inf)
    # expert.binom.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.binom.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.binom.tn[!no.trunc.idx,j])
    expert.binom.tn.bar[no.trunc.idx,j] = log1mexp(-expert.binom.tn[no.trunc.idx,j])
    # expert.nbinom.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.binom.ll=expert.binom.ll, expert.binom.tn=expert.binom.tn, expert.binom.tn.bar=expert.binom.tn.bar)
}

test_that("Binomial Expert", {
  expect_equal(ExpertBinom(tl, yl, yu, tu, 20, 0.5)[[1]], expert.binom(tl, yl, yu, tu, 1, 20, 0.5)[[1]])
  expect_equal(ExpertBinom(tl, yl, yu, tu, 20, 0.5)[[2]], expert.binom(tl, yl, yu, tu, 1, 20, 0.5)[[2]])
  expect_equal(ExpertBinom(tl, yl, yu, tu, 20, 0.5)[[3]], expert.binom(tl, yl, yu, tu, 1, 20, 0.5)[[3]])
})

# library(microbenchmark)
# bench = microbenchmark(
#   ExpertBinom(tl, yl, yu, tu, 20, 0.5),
#   expert.binom(tl, yl, yu, tu, 1, 20, 0.5)
# )
#
# library(ggplot2)
# autoplot(bench)
