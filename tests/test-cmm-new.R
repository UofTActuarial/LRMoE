# data("LRMoEDemoData")
#
# library(fitdistrplus)
# library(actuar)
#
#
# d = 2
#
# Yd = pmin((Y.obs[,4*(d-1)+2] + Y.obs[,4*(d-1)+3])/2, Y.obs[,4*(d-1)+2])
#
# temp = mledist(Yd[which(Yd>0)], "weibull")
#
# mean.pos = mean(Yd[which(Yd>0)])
# var.pos = var(Yd[which(Yd>0)])
#
# mean.pos^2/var.pos
# var.pos/mean.pos
#
#
# tempf = function(params, y){
#   return(-sum(dburr(y, shape1 = exp(params[1]), shape2 = exp(params[2]), scale = exp(params[3]), log = T)))
# }
#
#
# temp = optim(log(c(2, 1, 10)), tempf, y = Yd[which(Yd>0)], method = "Nelder-Mead")
#
#
# tt = CMMInit(Y.obs[,1:4], X.obs, 2, 'F')
#
# tt = CMMInit(Y.obs, X.obs, 2, c("F", "S"))
#
# tt = CMMSeverity(Y.obs, X.obs, 2)
