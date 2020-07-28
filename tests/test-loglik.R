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
#
# gate.ll = GateLogit(X.obs, alpha.init)
# expert.list = DimCompExpertLL(Y.obs, comp.dist, zero.init, params.init)
# loglik.list = GateExpertLL(alpha.init, gate.ll, expert.list, T, 5, params.pen)

# Returns numbers!
