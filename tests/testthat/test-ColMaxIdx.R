# context("ColMaxIdx")
# library(Rcpp)
#
# w = matrix(runif(1000000, 0, 1), ncol = 5)
#
# temp1 = ColMaxIdx(w)
# temp2 = matrix(apply(w, 1, FUN = "which.max"))
#
#
# test_that("ColMaxIdx", {
#   expect_equal(temp1, temp2)
# })
#
# library(microbenchmark)
# bench = microbenchmark(
#   ColMaxIdx(w),
#   matrix(apply(w, 1, FUN = "which.max"))
# )
#
# library(ggplot2)
# autoplot(bench)

# Very significant improvement!
