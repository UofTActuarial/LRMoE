context("EMMInverseGaussian")
# library(Rcpp)

# mu = 15
# lambda = 10
# int.y.fcn = function(u, mean.mu, shape.lambda) # u = log(y) for numerical stability
# {
#   u.dens.log = statmod::dinvgauss(exp(u), mean = mean.mu, shape = shape.lambda, log = TRUE) + u
#   temp = exp(u) * exp(u.dens.log)
#   temp[which(u==Inf)] = 0
#   temp[which(u==-Inf)] = 0
#   temp[which(is.na(temp))] = 0
#   return(temp)
# }
#
# int.y.log.fcn = function(u, mean.mu, shape.lambda) # u = log(y) for numerical stability
# {
#   u.dens.log = statmod::dinvgauss(exp(u), mean = mean.mu, shape = shape.lambda, log = TRUE) + u
#   temp = u * exp(u.dens.log)
#   temp[which(u==Inf)] = 0
#   temp[which(u==-Inf)] = 0
#   temp[which(is.na(temp))] = 0
#   return(temp)
# }
#
# int.y.inv.fcn = function(u, mean.mu, shape.lambda) # u = log(y) for numerical stability
# {
#   u.dens.log = statmod::dinvgauss(exp(u), mean = mean.mu, shape = shape.lambda, log = TRUE) + u
#   temp = exp(-u) * exp(u.dens.log)
#   temp[which(u==Inf)] = 0
#   temp[which(u==-Inf)] = 0
#   temp[which(is.na(temp))] = 0
#   return(temp)
# }
#
# yl = rinvgauss(10000, mean = mu, shape = lambda)
# yu = yl + rinvgauss(10000, mean = mu, shape = lambda)
#
# yl[1:1000] = 0
# yu[1001:2000] = Inf
# yl[2001:3000] = yu[2001:3000]
#
# temp1 = mapply(function(x, y) ifelse(x!=y,
#                                      integrate(int.y.log.fcn, log(x), log(y),
#                                                mean.mu=mu, shape.lambda=lambda,
#                                                rel.tol=.Machine$double.eps^0.5)$value,
#                                      0),
#                yl, yu)
# temp2 = intInvGaussLogYObs(mu, lambda, log(yl), log(yu))
#
# temp3 = sapply(yl,
#                function(x) ifelse(x!=0,
#                                   integrate(int.y.log.fcn, -Inf, log(x),
#                                             mean.mu=mu, shape.lambda=lambda,
#                                             rel.tol=.Machine$double.eps^0.5)$value,
#                                   0))+
#         sapply(yu,
#                function(x) ifelse(x!=Inf,
#                                   integrate(int.y.log.fcn, log(x), Inf,
#                                             mean.mu=mu, shape.lambda=lambda,
#                                             rel.tol=.Machine$double.eps^0.5)$value,
#                                   0))
# temp4 = intInvGaussLogYLat(mu, lambda, log(yl), log(yu))
#
# temp5 = mapply(function(x, y) ifelse(x!=y,
#                                     integrate(int.y.fcn, log(x), log(y),
#                                               mean.mu=mu, shape.lambda=lambda,
#                                               rel.tol=.Machine$double.eps^0.5)$value,
#                                     0),
#               yl, yu)
# temp6 = intInvGaussYObs(mu, lambda, log(yl), log(yu))
#
# temp7 = sapply(yl,
#                function(x) ifelse(x!=0,
#                                   integrate(int.y.fcn, -Inf, log(x),
#                                             mean.mu=mu, shape.lambda=lambda,
#                                             rel.tol=.Machine$double.eps^0.5)$value,
#                                   0))+
#         sapply(yu,
#                function(x) ifelse(x!=Inf,
#                                   integrate(int.y.fcn, log(x), Inf,
#                                             mean.mu=mu, shape.lambda=lambda,
#                                             rel.tol=.Machine$double.eps^0.5)$value,
#                                   0))
# temp8 = intInvGaussYLat(mu, lambda, log(yl), log(yu))
#
# temp9 = mapply(function(x, y) ifelse(x!=y,
#                                      integrate(int.y.inv.fcn, log(x), log(y),
#                                                mean.mu=mu, shape.lambda=lambda,
#                                                rel.tol=.Machine$double.eps^0.5)$value,
#                                      0),
#                yl, yu)
# temp10 = intInvGaussInvYObs(mu, lambda, log(yl), log(yu))
#
# temp11 = sapply(yl,
#                 function(x) ifelse(x!=0,
#                                    integrate(int.y.inv.fcn, -Inf, log(x),
#                                              mean.mu=mu, shape.lambda=lambda,
#                                              rel.tol=.Machine$double.eps^0.5)$value,
#                                    0))+
#         sapply(yu,
#                function(x) ifelse(x!=Inf,
#                                   integrate(int.y.inv.fcn, log(x), Inf,
#                                             mean.mu=mu, shape.lambda=lambda,
#                                             rel.tol=.Machine$double.eps^0.5)$value,
#                                   0))
# temp12 = intInvGaussInvYLat(mu, lambda, log(yl), log(yu))
#
# test_that("EMMInverseGaussian", {
#   expect_equal(temp1, temp2)
#   expect_equal(temp3, temp4)
#   expect_equal(temp5, temp6)
#   expect_equal(temp7, temp8)
#   expect_equal(temp9, temp10)
#   expect_equal(temp11, temp12)
# })

# library(microbenchmark)
# bench = microbenchmark(
#
#   intInvGaussLogYObs(mu, lambda, log(yl), log(yu)),
#
#   mapply(function(x, y) ifelse(x!=y,
#                                integrate(int.y.log.fcn, log(x), log(y),
#                                          mean.mu=mu, shape.lambda=lambda,
#                                          rel.tol=.Machine$double.eps^0.5)$value,
#                                0),
#          yl, yu)
# )
#
# library(ggplot2)
# autoplot(bench)

# library(microbenchmark)
# bench = microbenchmark(
#
#   intInvGaussLogYLat(mu, lambda, log(yl), log(yu)),
#
#   sapply(yl,
#          function(x) ifelse(x!=0,
#                             integrate(int.y.log.fcn, -Inf, log(x),
#                                       mean.mu=mu, shape.lambda=lambda,
#                                       rel.tol=.Machine$double.eps^0.5)$value,
#                             0))+
#     sapply(yu,
#            function(x) ifelse(x!=Inf,
#                               integrate(int.y.log.fcn, log(x), Inf,
#                                         mean.mu=mu, shape.lambda=lambda,
#                                         rel.tol=.Machine$double.eps^0.5)$value,
#                               0)),
#   times = 10
# )
#
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
# comp.dist = matrix(c("invgauss", "ZI-invgauss"),
#                    nrow = 1, byrow = T)
#
# zero.prob = matrix(c(0, 0.10),
#                    nrow = 1, byrow = T)
#
# params.list = list(list(c(10, 20), c(25, 10)))
#
# simy = SimYSet(X, alpha, comp.dist, zero.prob, params.list)
#
# hist(simy, breaks = 300, xlim = c(0, 100))
#
# YY = cbind(rep(0, nrow(X)), simy, simy, rep(Inf, nrow(X)))
#
# YY[c(1000:2000),2] = simy[c(1000:2000)] * 0.75
# YY[c(2001:3000),3] = simy[c(2001:3000)] * 1.50
#
# YY[c(3001:4000),1] = simy[c(3001:4000)] * 0.25
#
# YY[c(4001:5000),4] = simy[c(4001:5000)] * 2
#
# YY[c(5001:7000),1] = simy[c(5001:7000)] * 0.25
# YY[c(5001:7000),2] = simy[c(5001:7000)] * 0.75
# YY[c(5001:7000),3] = simy[c(5001:7000)] * 1.25
# YY[c(5001:7000),4] = simy[c(5001:7000)] * 2
#
#
# alpha.init = alpha = matrix(c(0, 0, 0, 0, 0,
#                               0, 0, 0, 0, 0),
#                             nrow = 2, byrow = T)
# comp.dist = comp.dist
# zero.init = matrix(c(0, 0.50),
#                    nrow = 1, byrow = T)
# params.init = list(list(c(15, 15), c(20, 8)))
#
#
# modelfit = LRMoEFit(YY, X, 2, comp.dist, alpha.init, zero.init, params.init)
#
