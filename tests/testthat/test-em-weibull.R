context("EMMWeibull")
# library(Rcpp)

k = 3
lambda = 10
int.y.log.fcn = function(u, shape.k.j, scale.lambda.j) # u = log(y) for numerical stability # NO
{
  u.density.log = stats::dweibull(u, shape = shape.k.j, scale = scale.lambda.j, log = TRUE)
    # stats::dweibull(exp(u), shape = shape.k.j, scale = scale.lambda.j, log = TRUE) + u
  # Numerical underflow of dweibull for small u, and certain set of parameters. e.g. (4, 35)
  temp = log(u)*exp(u.density.log)
  # temp[which(u==Inf)] = 0
  # temp[which(u==-Inf)] = 0
  temp[which(u==0)] = 0
  temp[which(u==Inf)] = 0
  temp[which(is.na(temp))] = 0
  return(temp)
}


yl = rweibull(10000, k, lambda)
yu = yl + rweibull(10000, k, lambda)

yl[1:1000] = 0
yu[1001:2000] = Inf
yl[2001:3000] = yu[2001:3000]

temp1 = mapply(function(x, y) ifelse(x!=y,
                                     integrate(int.y.log.fcn, x, y,
                                               shape.k.j=k, scale.lambda.j = lambda,
                                               rel.tol=.Machine$double.eps^0.5)$value,
                                     0),
               yl, yu)
temp2 = intWeibullLogYObs(k, lambda, (yl), (yu))

# temp3 = sapply(yl,
#                function(x) ifelse(x!=0,
#                                   integrate(int.y.log.fcn, -Inf, (x),
#                                             shape.k.j=k, scale.lambda.j = lambda,
#                                             rel.tol=.Machine$double.eps^0.5)$value,
#                                   0))+
#   sapply(yu,
#          function(x) ifelse(x!=Inf,
#                             integrate(int.y.log.fcn, (x), stats::qweibull(1-1e-09, shape = k, scale = lambda, log.p = F), # Inf is a problem
#                                       shape.k.j=k, scale.lambda.j = lambda,
#                                       rel.tol=.Machine$double.eps^0.5)$value,
#                             0))
# temp4 = intWeibullLogYLat(k, lambda, (yl), (yu))


test_that("EMMWeibull", {
  expect_equal(temp1, temp2)
  # expect_equal(temp3, temp4)
})

# library(microbenchmark)
# bench = microbenchmark(
#
#   intWeibullLogYObs(k, lambda, (yl), (yu)),
#
#   mapply(function(x, y) ifelse(x!=y,
#                                integrate(int.y.log.fcn, x, y,
#                                          shape.k.j=k, scale.lambda.j = lambda,
#                                          rel.tol=.Machine$double.eps^0.5)$value,
#                                0),
#          yl, yu)
# )
#
# library(ggplot2)
# autoplot(bench)


# Very good improvement!




data("LRMoEDemoData")
head(X)

alpha = matrix(c(0.5, 0.25, -0.05, 0.3, -0.2,
                 0, 0, 0, 0, 0),
               nrow = 2, byrow = T)

head(exp(GateLogit(X, alpha)))
hist(exp(GateLogit(X, alpha))[,1])

comp.dist = matrix(c("weibull", "ZI-weibull"),
                   nrow = 1, byrow = T)

zero.prob = matrix(c(0, 0.10),
                   nrow = 1, byrow = T)

params.list = list(list(c(3, 10), c(5, 30)))

simy = SimYSet(X, alpha, comp.dist, zero.prob, params.list)

hist(simy, breaks = 100, xlim = c(0, 50))

YY = cbind(rep(0, nrow(X)), simy, simy, rep(Inf, nrow(X)))

YY[c(1000:2000),2] = simy[c(1000:2000)] * 0.75
YY[c(2001:3000),3] = simy[c(2001:3000)] * 1.50

YY[c(3001:4000),1] = simy[c(3001:4000)] * 0.25

YY[c(4001:5000),4] = simy[c(4001:5000)] * 2

YY[c(5001:7000),1] = simy[c(5001:7000)] * 0.25
YY[c(5001:7000),2] = simy[c(5001:7000)] * 0.75
YY[c(5001:7000),3] = simy[c(5001:7000)] * 1.25
YY[c(5001:7000),4] = simy[c(5001:7000)] * 2


alpha.init = alpha = matrix(c(0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0),
                            nrow = 2, byrow = T)
comp.dist = comp.dist
zero.init = matrix(c(0, 0.50),
                   nrow = 1, byrow = T)
params.init = list(list(c(2, 15), c(4, 20)))


modelfit = LRMoEFit(YY, X, 2, comp.dist, alpha.init, zero.init, params.init)

