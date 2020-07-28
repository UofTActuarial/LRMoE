context("GammaCount")
library(copula)
library(rmutil)

y = c(0:999)

test_that("GammaCount", {
  expect_equal(dgammacount(y, 2, 2, F), rmutil::dgammacount(y, 2, 2, F))
  # expect_equal(dgammacount(y, 2, 2, T), rmutil::dgammacount(y, 2, 2, T)) # Should fail this one
  expect_equal(pgammacount(y, 2, 2, F), rmutil::pgammacount(y, 2, 2))
})
