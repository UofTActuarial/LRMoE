## Severity Initialization
#' Initializes parameter for severity distributions using CMM.
#'
#' @param Y A vector of response variables.
#' @param X A matrix of covariates (normalized for running \code{\link[stats]{kmeans}}).
#' @param n.comp Number of components
#' @return A list of parameter initialization.
#'
#' @importFrom EnvStats skewness kurtosis
#' @importFrom stats var kmeans
#'
#' @export CMMSeverity
CMMSeverity = function(Y, X, n.comp)
{
  km = kmeans(data.matrix(X), n.comp)
  cluster = km$cluster

  result = NULL

  for(j in 1:n.comp){
    subset = Y[which(cluster==j)]
    subset.pos = subset[which(subset>0)]

    cluster.prop = length(subset) / length(Y)

    zero.prop = 1 - length(subset.pos)/length(subset)

    mean.pos = mean(subset.pos)
    var.pos = var(subset.pos)
    cv.pos = sqrt(var.pos)/mean.pos
    skew.pos = skewness(subset.pos)
    kurt.pos = kurtosis(subset.pos)

    gamma.init = c(shape = mean.pos^2/var.pos, scale = var.pos/mean.pos)
    lnorm.init = c(meanlog = log(mean.pos)- 0.5*(log(var.pos/(mean.pos^2) + 1)), sdlog = sqrt(log(var.pos/(mean.pos^2) + 1)))
    invgauss.init = c(mean = mean.pos, scale = (mean.pos)^3/var.pos )
    weibull.init = c(shape = 3, scale = mean.pos * 1.5 )# Ad-hoc
    burr.init = c(shape1 = 2, shape2 = 5, scale = mean.pos/0.136) # Ad-hoc. 0.136 = Beta(4.8, 1.2)

    result[[length(result)+1]] = list(cluster.prop = cluster.prop, zero.prop = zero.prop, mean.pos = mean.pos, var.pos = var.pos, cv.pos = cv.pos, skew.pos = skew.pos, kurt.pos = kurt.pos,
                                      gamma.init = gamma.init, lnorm.init = lnorm.init, invgauss.init = invgauss.init,
                                      weibull.init = weibull.init, burr.init = burr.init)
  }

  return(result)
}


## Frequency Initialization
#' Initializes parameter for frequency distributions using CMM.
#'
#' @param Y A vector of response variables.
#' @param X A matrix of covariates (normalized for running \code{\link[stats]{kmeans}}).
#' @param n.comp Number of components
#' @return A list of parameter initialization.
#'
#' @importFrom EnvStats skewness kurtosis
#' @importFrom stats var kmeans
#'
#' @export CMMFrequency
CMMFrequency = function(Y, X, n.comp)
{
  km = kmeans(data.matrix(X), n.comp)
  cluster = km$cluster

  result = NULL

  for(j in 1:n.comp){
    subset = Y[which(cluster==j)]
    subset.pos = subset[which(subset>0)]

    cluster.prop = length(subset) / length(Y)

    zero.prop = 1 - length(subset.pos)/length(subset)

    mean.pos = mean(subset.pos)
    var.pos = var(subset.pos)
    cv.pos = sqrt(var.pos)/mean.pos
    skew.pos = skewness(subset.pos, method = "moment")
    kurt.pos = kurtosis(subset.pos, method = "moment")

    poisson.init = c(lambda = mean.pos)
    ztpoisson.init = c(lambda = mean.pos)
    nbinom.init = c(size.n = mean.pos * (mean.pos/var.pos/(1-mean.pos/var.pos)), prob.p = mean.pos/var.pos)
    binom.init = c(size.n = mean.pos/(1-var.pos/mean.pos), prob.p = 1-var.pos/mean.pos)
    gammacount.init = c(shape = mean.pos * (mean.pos/var.pos/(1-mean.pos/var.pos)), scale = mean.pos/var.pos ) # ad-hoc, same as nbinom


    result[[length(result)+1]] = list(cluster.prop = cluster.prop, zero.prop = zero.prop, mean.pos = mean.pos, var.pos = var.pos, cv.pos = cv.pos, skew.pos = skew.pos, kurt.pos = kurt.pos,
                                      poisson.init = poisson.init, ztpoisson.init = ztpoisson.init,
                                      nbinom.init = nbinom.init, binom.init = binom.init, gammacount.init = gammacount.init)
  }

  return(result)
}
