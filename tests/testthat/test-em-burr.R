context("EMMBurr")
# library(Rcpp)
# library(expint)
# library(actuar)
#
# k = 3
# c = 1
# lambda = 10
# int.y.log.fcn = function(u, shape1.k.j, shape2.c.j, scale.lambda.j) # u = log(y) for numerical stability
# {
#   u.density.log = u + actuar::dburr(exp(u), shape1 = shape1.k.j, shape2 = shape2.c.j, scale = scale.lambda.j, log = TRUE)
#   temp = u * exp(u.density.log)
#   temp[which(u==Inf)] = 0
#   temp[which(u==-Inf)] = 0
#   temp[which(is.na(temp))] = 0
#
#   temp[which(temp==Inf)] = 0
#
#   return(temp)
# }
#
# int.y.pol.fcn = function(u, shape2.c.new, scale.lambda.new,
#                          shape1.k.j, shape2.c.j, scale.lambda.j) # u = log(y) for numerical stability
# {
#   # Use existing package so that I don't have to manually type in the density!
#   u.density.log = u + actuar::dburr(exp(u), shape1 = shape1.k.j, shape2 = shape2.c.j, scale = scale.lambda.j, log = TRUE)
#   temp = log(1+(exp(u)/scale.lambda.new)^shape2.c.new) * exp(u.density.log)
#   temp[which(u==Inf)] = 0
#   temp[which(u==-Inf)] = 0
#   temp[which(is.na(temp))] = 0
#
#   temp[which(temp==Inf)] = 0
#   return(temp)
# }
#
# yl = rburr(10000, shape1 = k, shape2 = c, scale = lambda)
# yu = yl + rburr(10000, shape1 = k, shape2 = c, scale = lambda)
#
# yl[1:1000] = 0
# yu[1001:2000] = Inf
# yl[2001:3000] = yu[2001:3000]
#
# temp1 = mapply(function(x, y) ifelse(x!=y,
#                                      integrate(int.y.log.fcn, log(x), log(y),
#                                                shape1.k.j = k, shape2.c.j = c, scale.lambda.j = lambda,
#                                                rel.tol=.Machine$double.eps^0.5)$value,
#                                      0),
#                yl, yu)
# temp2 = intBurrLogYObs(k, c, lambda, log(yl), log(yu))
#
#
# temp3 = sapply(yl,
#                function(x) ifelse(x!=0,
#                                   integrate(int.y.log.fcn, -Inf, log(x),
#                                             shape1.k.j = k, shape2.c.j = c, scale.lambda.j = lambda,
#                                             rel.tol=.Machine$double.eps^0.5)$value,
#                                   0))+
#   sapply(yu,
#          function(x) ifelse(x!=Inf,
#                             integrate(int.y.log.fcn, log(x), Inf, # actuar::qburr(1-1e-09, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = T, log.p = F), # Inf causes a problem
#                                       shape1.k.j = k, shape2.c.j = c, scale.lambda.j = lambda,
#                                       rel.tol=.Machine$double.eps^0.5)$value,
#                             0))
#
# temp4 = intBurrLogYLat(k, c, lambda, log(yl), log(yu))
#
# summary(temp1-temp2)
# summary(temp3-temp4)
#
# cc = 2
# ll = 15
#
# temp5 = mapply(function(x, y) ifelse(x!=y,
#                                     integrate(int.y.pol.fcn, log(x), log(y),
#                                               shape2.c.new = cc, scale.lambda.new = ll,
#                                               shape1.k.j = k, shape2.c.j = c, scale.lambda.j = lambda,
#                                               rel.tol=.Machine$double.eps^0.5)$value,
#                                     0),
#               yl, yu)
#
# temp6 = intBurrPolYObs(k, c, lambda, cc, ll, log(yl), log(yu))
#
# temp7 = sapply(yl,
#                function(x) ifelse(x!=0,
#                                   integrate(int.y.pol.fcn, -Inf, log(x),
#                                             shape2.c.new = cc, scale.lambda.new = ll,
#                                             shape1.k.j = k, shape2.c.j = c, scale.lambda.j = lambda,
#                                             rel.tol=.Machine$double.eps^0.5)$value,
#                                   0))+
#   sapply(yu,
#          function(x) ifelse(x!=Inf,
#                             integrate(int.y.pol.fcn, log(x), Inf, # actuar::qburr(1-1e-09, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = T, log.p = F), # Inf causes a problem
#                                       shape2.c.new = cc, scale.lambda.new = ll,
#                                       shape1.k.j = k, shape2.c.j = c, scale.lambda.j = lambda,
#                                       rel.tol=.Machine$double.eps^0.5)$value,
#                             0))
#
# temp8 = intBurrPolYLat(k, c, lambda, cc, ll, log(yl), log(yu))
#
#
# test_that("EMMBurr", {
#   expect_equal(temp1, temp2)
#   expect_equal(temp3, temp4)
#   expect_equal(temp5, temp6)
#   expect_equal(temp7, temp8)
# })

# library(microbenchmark)
# bench = microbenchmark(
#
#   intBurrPolYLat(k, c, lambda, cc, ll, log(yl), log(yu)),
#
#   sapply(yl,
#          function(x) ifelse(x!=0,
#                             integrate(int.y.pol.fcn, -Inf, log(x),
#                                       shape2.c.new = cc, scale.lambda.new = ll,
#                                       shape1.k.j = k, shape2.c.j = c, scale.lambda.j = lambda,
#                                       rel.tol=.Machine$double.eps^0.5)$value,
#                             0))+
#     sapply(yu,
#            function(x) ifelse(x!=Inf,
#                               integrate(int.y.pol.fcn, log(x), Inf, # actuar::qburr(1-1e-09, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = T, log.p = F), # Inf causes a problem
#                                         shape2.c.new = cc, scale.lambda.new = ll,
#                                         shape1.k.j = k, shape2.c.j = c, scale.lambda.j = lambda,
#                                         rel.tol=.Machine$double.eps^0.5)$value,
#                               0))
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
# comp.dist = matrix(c("burr", "ZI-burr"),
#                    nrow = 1, byrow = T)
#
# zero.prob = matrix(c(0, 0.10),
#                    nrow = 1, byrow = T)
#
# params.list = list(list(c(1, 3, 10), c(2, 5, 30)))
#
# simy = SimYSet(X, alpha, comp.dist, zero.prob, params.list)
#
# hist(simy, breaks = 100, xlim = c(0, 50))
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
# params.init = list(list(c(1.5, 2.5, 15), c(3, 4, 20)))
#
#
# modelfit = LRMoEFit(YY, X, 2, comp.dist, alpha.init, zero.init, params.init)

