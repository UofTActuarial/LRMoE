context("Gating")
library(matrixStats)

X = matrix(runif(500000, -1, 1), ncol = 5)
alpha = matrix(runif(50, -1, 1), ncol = 5)

gate.body=tcrossprod(X, alpha)

test_that("Gating function", {
  expect_equal(GateLogit(X, alpha), sweep(gate.body, 1, rowLogSumExps(gate.body), FUN = "-", check.margin = FALSE))
})

gate.old = function(x, alpha){
  gate.body=tcrossprod(X, alpha)
  return(sweep(gate.body, 1, rowLogSumExps(gate.body), FUN = "-", check.margin = FALSE))
}
#
# library(microbenchmark)
# bench = microbenchmark(
#   GateLogit(X, alpha),
#   gate.old(X, alpha)
# )
#
# library(ggplot2)
# autoplot(bench)
