context("ExpertBurr")
library(copula)
library(actuar)

y = actuar::rburr(100000, shape1 = 1, shape2 = 2, scale = 10)
tl = y*0.50
yl = y*0.75
yu = y*1.25
tu = rep(Inf, 100000)

expert.burr = function(tl, yl, yu, tu, g = 1, shape1.k, shape2.c, scale.lambda)
{
  # Initialization: return value are N * g matrices
  expert.burr.ll=expert.burr.tn=expert.burr.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=actuar::pburr(yu[censor.idx],shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.yl=actuar::pburr(yl[censor.idx],shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    # prob.log.yu = log1mexp( shape1.k[j] * log(1+(yu[censor.idx]/scale.lambda[j])^(shape2.c[j])) )
    # prob.log.yl = log1mexp( shape1.k[j] * log(1+(yl[censor.idx]/scale.lambda[j])^(shape2.c[j])) )

    # Compute loglikelihood for expert j, first for y
    expert.burr.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra

    ###################################################################
    # Deal with numerical underflow: prob.log.yu and prob.log.yl can both be -Inf
    NA.idx = which(is.na(expert.burr.ll[,j]))
    expert.burr.ll[NA.idx, j] = -Inf
    ###################################################################

    expert.burr.ll[!censor.idx,j]=actuar::dburr(yu[!censor.idx],shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], log = TRUE) # exact likelihood
    # expert.burr.ll[!censor.idx,j] = log(shape1.k[j]*shape2.c[j]) + (shape2.c[j]-1) * log(yu[!censor.idx]) - shape2.c[j] * log(scale.lambda[j]) +
    #                                 (-shape1.k[j]-1) * log1pexp(shape2.c[j] * log(yu[!censor.idx]/scale.lambda[j]) )

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=actuar::pburr(tu,shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.tl=actuar::pburr(tl,shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    # prob.log.tu = log1mexp( shape1.k[j] * log(1+(tu/scale.lambda[j])^(shape2.c[j])) )
    # prob.log.tl = log1mexp( shape1.k[j] * log(1+(tl/scale.lambda[j])^(shape2.c[j])) )

    # Normalizing factor for truncation limits, in log
    expert.burr.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with numerical underflow: prob.log.tu and prob.log.tl can both be -Inf
    NA.idx = which(is.na(expert.burr.tn[,j]))
    expert.burr.tn[NA.idx, j] = -Inf
    ###################################################################

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    # expert.burr.tn[no.trunc.idx,j] = log(shape1.k[j]*shape2.c[j]) + (shape2.c[j]-1) * log(tu[no.trunc.idx]) - shape2.c[j] * log(scale.lambda[j]) +
    #   (-shape1.k[j]-1) * log1pexp(shape2.c[j] * log(tu[no.trunc.idx]/scale.lambda[j]) )
    expert.burr.tn[no.trunc.idx,j] = actuar::dburr(tu[no.trunc.idx],shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], log = TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    zero.idx = (tu==0)
    expert.burr.ll[zero.idx,j]=(-Inf)
    expert.burr.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.burr.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.burr.tn[!no.trunc.idx,j])
    expert.burr.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.burr.ll=expert.burr.ll, expert.burr.tn=expert.burr.tn, expert.burr.tn.bar=expert.burr.tn.bar)
}

test_that("Burr Expert", {
  expect_equal(ExpertBurr(tl, yl, yu, tu, 1, 2, 10)[[1]], expert.burr(tl, yl, yu, tu, 1, 1, 2, 10)[[1]])
  expect_equal(ExpertBurr(tl, yl, yu, tu, 1, 2, 10)[[2]], expert.burr(tl, yl, yu, tu, 1, 1, 2, 10)[[2]])
  expect_equal(ExpertBurr(tl, yl, yu, tu, 1, 2, 10)[[3]], expert.burr(tl, yl, yu, tu, 1, 1, 2, 10)[[3]])
})

# library(microbenchmark)
# bench = microbenchmark(
#   ExpertBurr(tl, yl, yu, tu, 1, 2, 10),
#   expert.burr(tl, yl, yu, tu, 1, 1, 2, 10)
# )
#
# library(ggplot2)
# autoplot(bench)
