#' Modified \code{\link[rmutil]{GammaCount}} cdf for better numerical accuracy.
#'
#' @param q Vector of quantiles.
#' @param m,s Paramaters of Gamma count distribution.
#' @param log.p TRUE/FALSE: whether log.p should be returned.
#'
#' @seealso \code{\link[rmutil]{GammaCount}}.
#'
#' @importFrom stats pgamma
#' @importFrom copula log1mexp
#'
#' @export pgammacount
#'
pgammacount = function(q, m, s, log.p=FALSE)
{
  if(any(m<=0))stop("m must be positive")
  if(any(s<=0))stop("s must be positive")

  neg.idx = which(q<0)
  q = floor(q) # convert to integers
  if(log.p==TRUE){
    temp = log1mexp(-pgamma(m*s,(q+1)*s,1,log.p=TRUE))
    temp[neg.idx] = -Inf
    return(temp)
  }else{
    temp = 1-pgamma(m*s,(q+1)*s,1)
    temp[neg.idx] = 0
    return(temp)
  }
}

#' Modified \code{\link[rmutil]{GammaCount}} pmf for better numerical accuracy.
#'
#' @param y Vector of gamma count values.
#' @param m,s Paramaters of Gamma count distribution.
#' @param log TRUE/FALSE: whether log density should be returned.
#'
#' @seealso \code{\link[rmutil]{GammaCount}}.
#'
#' @importFrom stats pgamma
#' @importFrom copula log1mexp
#'
#' @export dgammacount
#'
dgammacount = function(y, m, s, log=FALSE)
{
  if(any(m<=0))stop("m must be positive")
  if(any(s<=0))stop("s must be positive")

  neg.idx = which(y<0)
  nonint.idx = which(y!=floor(y))
  tmp = ifelse(y==0,pgamma(m*s,(y+1)*s,1,log.p=TRUE,lower.tail=FALSE),
               # pgamma(m*s,y*s+(y==0),1,log=TRUE) + log1mexp(pgamma(m*s,y*s+(y==0),1,log.p=TRUE)-pgamma(m*s,(y+1)*s,1,log.p=TRUE)) )
               pgamma(m*s,y*s,1,log.p=TRUE) + log1mexp( pgamma(m*s,y*s,1,log.p=TRUE)-pgamma(m*s,(y+1)*s,1,log.p=TRUE)) )

  tmp[neg.idx] = -Inf
  tmp[nonint.idx] = -Inf

  if(log==TRUE){
    return(tmp)
  }else{
    return(exp(tmp))
  }
}

#' Calculates moments of \code{\link[rmutil]{GammaCount}} using finite approximation.
#'
#' @param order A vector of positive power indices.
#' @param m,s Paramaters of Gamma count distribution.
#' @param tol Cut-off probability threshold. Values above (1-tol) are discarded.
#'
#' @seealso \code{\link[rmutil]{GammaCount}}.
#'
#' @return A vector of Gamma count distribution moments.
#'
#' @importFrom rmutil qgammacount
#'
#' @export mgammacount
mgammacount = function(order, m, s, tol = 1e-10)
{
  upper = rmutil::qgammacount(p = 1 - tol, m = m, s= s)

  y.series = c(1:upper)
  prob.series = dgammacount(y.series, m = m, s= s, log = FALSE)

  means = array(0, length(order))
  for(j in 1:length(order))
  {
    means[j] = sum((y.series^order[j])*prob.series)
  }

  return(means)
}
