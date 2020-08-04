context("EMMGammaCount")
# library(Rcpp)

# n = 20
# p = 0.30
#
# int.y.fcn = function(lower, upper, size.n.j, prob.p.j)
# {
#   lower.bound = ceiling(lower)
#   upper.bound = floor(upper)
#
#   if (upper!=Inf) {
#     y.series = c((lower.bound):(upper.bound))
#     dens.series = dnbinom(y.series, size = size.n.j, prob = prob.p.j, log = FALSE)
#     result = sum(y.series * dens.series)
#   }else{
#     if(lower.bound<=1) {
#       y.series = c(0)
#     }else{
#       y.series = c((0):(lower.bound-1))
#     }
#     dens.series = dnbinom(y.series, size = size.n.j, prob = prob.p.j, log = FALSE)
#     result = size.n.j*(1-prob.p.j)/prob.p.j - sum(y.series * dens.series)
#   }
#
#   return(sum(result))
# }
#
# int.y.log.fac.fcn = function(size.n.new, lower, upper, size.n.j, prob.p.j)
# {
#   lower.bound = ceiling(lower)
#   upper.bound = floor(upper)
#
#   if (upper!=Inf) {
#     y.series = c((lower.bound):(upper.bound))
#     lfac.series = lfactorial(y.series + size.n.new - 1)
#     dens.series = dnbinom(y.series, size = size.n.j, prob = prob.p.j, log = FALSE)
#     result = sum(lfac.series * dens.series)
#   }else{
#     # I cannot find a closed-form expression for lfactorial(y+size.n.new-1).
#     # Throw away all tail values with very small probabilities.
#     # The lfactorial function grows very slow, so this should be fine.
#     upper.bound.finite = qnbinom(1-1e-10, size = size.n.j, prob = prob.p.j, lower.tail = TRUE)
#     y.series = c((lower.bound):(upper.bound.finite))
#     lfac.series = lfactorial(y.series + size.n.new - 1)
#     dens.series = dnbinom(y.series, size.n.j, prob.p.j, log = FALSE)
#     result = sum(lfac.series * dens.series)
#   }
#   return(result)
# }
#
# yl = rnbinom(10000, size= n, prob = p)
# yu = yl + rnbinom(10000, size= n, prob = p)
#
# nn = 25
#
# yl[1:1000] = 0
# yu[1001:2000] = Inf
# yl[2001:3000] = yu[2001:3000]
#
# temp1 = mapply(function(x, y) ifelse(x!=y,
#                                      int.y.fcn(x, y, n, p),
#                                      x*dnbinom(x, size=n, prob=p)),
#                yl, yu)
# temp2 = sumNegativeBinomialYObs(n, p, (yl), (yu))
#
# temp3 = mapply(function(x, y) ifelse(x!=y,
#                                      int.y.log.fac.fcn(size.n.new = nn, lower = x, upper = y, size.n.j = n, prob.p.j = p),
#                                      lfactorial(x+nn-1)*dnbinom(x, size=n, prob=p)),
#                yl, yu)
#
# temp4 = sumNegativeBinomialLfacYObs(n, p, nn, (yl), (yu))
#
# temp5 = mapply(function(x, y) ifelse(x!=y,
#                                      int.y.log.fac.fcn(size.n.new = nn, lower = x, upper = y, size.n.j = n, prob.p.j = p),
#                                      0),
#                # rep(0,tn.unique.length), tl.unique) +
#                rep(0,length(yl)), (ceiling(yl)-1)*(ceiling(yl)-1 >=0) ) +
#   mapply(function(x, y) ifelse(x!=y,
#                                int.y.log.fac.fcn(size.n.new = nn, lower = x, upper = y, size.n.j = n, prob.p.j = p),
#                                0),
#          # tu.unique, rep(Inf,tn.unique.length))
#          floor(yu)+1, rep(Inf,length(yu)))
#
# temp6 = sumNegativeBinomialLfacYLat(n, p, nn, (yl), (yu))

# test_that("EMMGammaCount", {
#   expect_equal(temp1, temp2)
# })

# library(microbenchmark)
# bench = microbenchmark(
#
#   sumNegativeBinomialLfacYLat(n, p, nn, (yl), (yu)),
#
#   mapply(function(x, y) ifelse(x!=y,
#                                int.y.log.fac.fcn(size.n.new = nn, lower = x, upper = y, size.n.j = n, prob.p.j = p),
#                                0),
#          # rep(0,tn.unique.length), tl.unique) +
#          rep(0,length(yl)), (ceiling(yl)-1)*(ceiling(yl)-1 >=0) ) +
#     mapply(function(x, y) ifelse(x!=y,
#                                  int.y.log.fac.fcn(size.n.new = nn, lower = x, upper = y, size.n.j = n, prob.p.j = p),
#                                  0),
#            # tu.unique, rep(Inf,tn.unique.length))
#            floor(yu)+1, rep(Inf,length(yu)))
# )

# library(ggplot2)
# autoplot(bench)


# Very good improvement!




# data("LRMoEDemoData")
# head(X)
#
# alpha = matrix(c(0.5, 0.25, -0.05, 0.3, -0.2,
#                  0, 0, 0, 0, 0),
#                nrow = 2, byrow = T)
#
# head(exp(GateLogit(X, alpha)))
# hist(exp(GateLogit(X, alpha))[,1])
#
# comp.dist = matrix(c("gammacount", "ZI-gammacount"),
#                    nrow = 1, byrow = T)
#
# zero.prob = matrix(c(0, 0.10),
#                    nrow = 1, byrow = T)
#
# params.list = list(list(c(20, 0.3), c(35, 1)))
#
# simy = SimYSet(X, alpha, comp.dist, zero.prob, params.list)
#
# hist(simy, breaks = 300, xlim = c(0, 100))
#
# YY = cbind(rep(0, nrow(X)), simy, simy, rep(Inf, nrow(X)))
#
# YY[c(1000:2000),2] = floor(simy[c(1000:2000)] * 0.75)
# YY[c(2001:3000),3] = ceiling(simy[c(2001:3000)] * 1.50)
#
# YY[c(3001:4000),1] = floor(simy[c(3001:4000)] * 0.25)
#
# # YY[c(4001:5000),4] = ceiling(simy[c(4001:5000)] * 1.5)
#
# YY[c(5001:7000),1] = floor(simy[c(5001:7000)] * 0.25)
# YY[c(5001:7000),2] = floor(simy[c(5001:7000)] * 0.75)
# YY[c(5001:7000),3] = ceiling(simy[c(5001:7000)] * 1.25)
# # YY[c(5001:7000),4] = ceiling(simy[c(5001:7000)] * 2)
#
#
# alpha.init = alpha = matrix(c(0, 0, 0, 0, 0,
#                               0, 0, 0, 0, 0),
#                             nrow = 2, byrow = T)
# comp.dist = comp.dist
# zero.init = matrix(c(0, 0.50),
#                    nrow = 1, byrow = T)
# params.init = list(list(c(20, 0.2), c(25, 0.8)))
#
#
# modelfit = LRMoEFit(YY, X, 2, comp.dist, alpha.init, zero.init, params.init)

