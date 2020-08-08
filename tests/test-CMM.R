# data("LRMoEDemoData")
#
# X.norm = X[,2:5]
# X.norm[,1] = (X.norm[,1]-mean(X.norm[,1]))/sd(X.norm[,1])
# X.norm[,2] = (X.norm[,2]-mean(X.norm[,2]))/sd(X.norm[,2])
# X.norm[,3] = (X.norm[,3]-mean(X.norm[,3]))/sd(X.norm[,3])
# X.norm[,4] = (X.norm[,4]-mean(X.norm[,4]))/sd(X.norm[,4])
# summary(X.norm)

# init.list = CMMSeverity(Y[,6], X.norm, 3)
# init.list = CMMSeverity(Y[,6], X.norm, 4)
#
# init.list = CMMFrequency(Y[,2], X.norm, 3)
# init.list = CMMFrequency(Y[,2], X.norm, 4)



# temp1 = CMMSeverity(Y[, 5:8], X.norm, 4)
# temp2 = CMMFrequency(Y[, 1:4], X.norm, 4)
#
#
# temp3 = CMMInit(Y, X.norm, 4, c('F', 'S'))
