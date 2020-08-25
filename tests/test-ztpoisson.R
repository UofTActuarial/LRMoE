# # Get covariates X
# data("LRMoEDemoData")
#
# alpha = matrix(c(0.5, 1, -0.04, 0.05, 0.02,
#                  0, 0, 0, 0, 0),
#                nrow = 2, byrow = T)
#
# comp.dist = matrix(c("ztpoisson", "ZI-ztpoisson"),
#                    nrow = 1, byrow = T)
#
# zero.prob = matrix(c(0, 0.25),
#                    nrow = 1, byrow = T)
#
# params.list = list(list(c(10), c(25)))
#
#
#
# raw_y = SimYSet(X, alpha, comp.dist, zero.prob, params.list)
#
# hist(raw_y)
#
#
# # Exact observations
# tl = rep(0, length(raw_y))
# yl = raw_y
# yu = raw_y
# tu = rep(Inf, length(raw_y))
#
# Y = matrix(c(tl, yl, yu, tu), byrow = F, ncol = 4)
#
# head(Y)
#
#
# alpha.init = matrix(c(0, 0, 0, 0, 0,
#                       0, 0, 0, 0, 0),
#                     nrow = 2, byrow = T)
#
# comp.dist = comp.dist
#
# zero.init = matrix(c(0, 0.5),
#                    nrow = 1, byrow = T)
#
# params.init = list(list(c(8), c(15)))
#
# model.fit = LRMoEFit(Y, X, 2, comp.dist, alpha.init, zero.init, params.init, print = T)
#
#
# # With censoring and truncation
# tl = floor(raw_y * 0.50) # rep(0, length(raw_y))
# yl = floor(raw_y * 0.80)
# yu = ceiling(raw_y * 1.20)
# tu = ceiling(raw_y * 2.0) # rep(Inf, length(raw_y))
#
# Y = matrix(c(tl, yl, yu, tu), byrow = F, ncol = 4)
#
# head(Y)
#
# model.fit = LRMoEFit(Y, X, 2, comp.dist, alpha.init, zero.init, params.init, print = T)
