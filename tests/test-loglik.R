# data(LRMoEDemoData)
#
# head(X)
# head(Y)
#
# alpha.init = matrix(runif(15, -1, 1), nrow = 3, byrow = T)
# comp.dist = matrix(c("poisson", "nbinom", "ZI-gammacount",
#                      "lnorm", "ZI-burr", "invgauss"),
#                    nrow = 2, byrow = T)
# zero.init = matrix(c(0, 0, 0.5,
#                      0, 0.5, 0),
#                    nrow = 2, byrow = T)
# params.init = list(
#   list(c(10), c(10, 0.5), c(5, 5)),
#   list(c(2, 1), c(1, 2, 10), c(10, 10))
# )
# params.pen = list(
#   list(c(2, 1), c(2, 1), c(2, 1, 2, 1)),
#   list(c(Inf, 1, Inf), c(2, 1, 2, 1, 2, 1), c(1, Inf, 1, Inf))
# )

# gate.ll = GateLogit(X.obs, alpha.init)
# expert.list = DimCompExpertLL(Y.obs, comp.dist, zero.init, params.init)
# loglik.list = GateExpertLL(alpha.init, gate.ll, expert.list, T, 5, params.pen)
#
# zkz.list = EMEzkz(gate.ll, expert.list, loglik.list)

# temp1 = PredictClassPrior(X, alpha.init, T)
# temp2 = PredictClassPrior(X, alpha.init, F)
# temp3 = PredictClassPosterior(Y, X, alpha.init, comp.dist, zero.init, params.init, T)
# temp4 = PredictClassPosterior(Y, X, alpha.init, comp.dist, zero.init, params.init, F)

# temp1 = PredictMeanPrior(X, alpha.init, comp.dist, zero.init, params.init)
# temp2 = PredictMeanPosterior(Y, X, alpha.init, comp.dist, zero.init, params.init)

# temp1 = PredictVariancePrior(X, alpha.init, comp.dist, zero.init, params.init)
# temp2 = PredictVariancePosterior(Y, X, alpha.init, comp.dist, zero.init, params.init)

# Returns numbers!




