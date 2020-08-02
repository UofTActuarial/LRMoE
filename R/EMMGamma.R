## ECM algorithm of Gamma expert
#' ECM: M-Step for Gamma expert.
#'
#'
#' @importFrom stats pgamma dgamma optimise integrate
#'
#' @keywords internal
#'
#' @export EMMGamma
EMMGamma = function(params.old,
                                tl, yl, yu, tu,
                                expert.ll, expert.tn, expert.tn.bar,
                                z.e.obs, z.e.lat, k.e,
                                penalty, hyper.params)
{
  # Constants
  sample.size.n = length(tl)
  # Take old parameters
  m = params.old[1]
  theta = params.old[2]
  # Take hyper parameters
  nu.m = hyper.params[1]
  lambda.m = hyper.params[2]
  nu.theta = hyper.params[3]
  lambda.theta = hyper.params[4]
  # Value to return
  params.new = params.old


  # E-Step: conditional expectations for y, log(y)
  censor.idx = (yl!=yu)
  y.e.obs = rep(0, sample.size.n) # array(0, dim = c(sample.size.n, 1))
  y.e.lat = rep(0, sample.size.n) # array(0, dim = c(sample.size.n, 1))
  y.log.e.obs = rep(0, sample.size.n) # array(0, dim = c(sample.size.n, 1))
  y.log.e.lat = rep(0, sample.size.n) # array(0, dim = c(sample.size.n, 1))

  # Conditional expectation of y: untruncated and uncensored case.
  y.e.obs[!censor.idx] = yu[!censor.idx]
  # Conditional expectation of y: untruncated but censored case.
  diff.dist.untrunc = pgamma(yu[censor.idx], shape = m+1, scale = theta, lower.tail = TRUE, log.p = FALSE) - pgamma(yl[censor.idx], shape = m+1, scale = theta, lower.tail = TRUE, log.p = FALSE)
  y.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * (diff.dist.untrunc * m * theta)
  # Conditional expecration of y: truncated case.
  diff.dist.trunc = pgamma(tu, shape = m+1, scale = theta, lower.tail = TRUE, log.p = FALSE) - pgamma(tl, shape = m+1, scale = theta, lower.tail = TRUE, log.p = FALSE)
  y.e.lat = exp(-expert.tn.bar) * ( (1-diff.dist.trunc) * m*theta)
  y.e.lat[is.na(y.e.lat)] = 0 # Hardcode zeros to prevent NaN's

  # Conditional expectation of log(y): untruncated but censored case.
  # A function to numerically integrate for log(y)
  # Q.y.log = function(u, m, theta) # Change of variable: u = log(y) for numerical stability
  # {
  #   dens.u.log = dgamma(exp(u), shape = m, scale = theta, log = TRUE) + u
  #   return( u *exp(dens.u.log))
  # }
  # Find unique integration limits, for numerical speed
  y.unique = unique(cbind(yl,yu),MARGIN=1)
  y.unique.length = nrow(y.unique)
  y.unique.match = match(data.frame(t(cbind(yl,yu))),data.frame(t(y.unique)))
  yl.unique = y.unique[,1]
  yu.unique = y.unique[,2]
  tn.unique = unique(cbind(tl,tu),MARGIN=1)
  tn.unique.length = nrow(tn.unique)
  tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
  tl.unique = tn.unique[,1]
  tu.unique = tn.unique[,2]
  y.log.e.obs.unique=array(0, dim = c(y.unique.length, 1))
  y.log.e.lat.unique=array(0, dim = c(tn.unique.length, 1))

  y.log.e.obs.unique[,1]= intGammaLogYObs(m, theta, log(yl.unique), log(yu.unique))
    # mapply(function(x, y) ifelse(x!=y,
    #                              integrate(Q.y.log, log(x), log(y),
    #                                        m = m, theta = theta,
    #                                        rel.tol=.Machine$double.eps^0.5)$value,
    #                              0),
    #        yl.unique, yu.unique)
  y.log.e.obs = y.log.e.obs.unique[y.unique.match,]
  y.log.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * y.log.e.obs[censor.idx]

  # Conditional expectation of log(y): untruncated and uncensored case.
  y.log.e.obs[!censor.idx] = log(yu[!censor.idx])
  y.log.e.obs[is.na(y.log.e.obs)] = 0 # Hardcode zeros to prevent NaN's

  # Conditional expectation of log(y): truncated case.
  y.log.e.lat.unique[,1]= intGammaLogYLat(m, theta, log(tl.unique), log(tu.unique))
    # sapply(tl.unique,
    #        function(x) ifelse(x!=0,
    #                           integrate(Q.y.log, -Inf, log(x),
    #                                     m = m, theta = theta,
    #                                     rel.tol=.Machine$double.eps^0.5)$value,
    #                           0))+
    # sapply(tu.unique,
    #        function(x) ifelse(x!=Inf,
    #                           integrate(Q.y.log, log(x), Inf,
    #                                     m = m, theta = theta,
    #                                     rel.tol=.Machine$double.eps^0.5)$value,
    #                           0))

  y.log.e.lat = y.log.e.lat.unique[tn.unique.match,]
  y.log.e.lat = exp(-expert.tn.bar) * y.log.e.lat
  y.log.e.lat[is.na(y.log.e.lat)] = 0 # Hardcode zeros to prevent NaN's


  # M-Step: maximization of loglik
  # Numerical optimization is needed. Optimal theta is a function of optimal m.
  theta.m = function(m, z.e.obs, z.e.lat, k.e, y.e.obs, y.e.lat, nu.m, lambda.m, nu.theta, lambda.theta)
  {
    term.zkz = XPlusYZ(z.e.obs, z.e.lat, k.e)
    # sweep(matrix(z.e.obs), 1,
    #                  sweep(matrix(k.e), 1, matrix(z.e.lat), FUN = "*", check.margin = FALSE),
    #                  FUN = "+", check.margin = FALSE)
    # z.e.obs + k.e * z.e.lat
    term.zkz.y = XAPlusYZB(z.e.obs, y.e.obs, z.e.lat, k.e, y.e.lat)
    # sweep(
    #   sweep(matrix(z.e.obs), 1, matrix(y.e.obs), FUN = "*", check.margin = FALSE), 1,
    #   sweep(matrix(k.e), 1,
    #         sweep(matrix(z.e.lat), 1, matrix(y.e.lat), FUN = "*", check.margin = FALSE),
    #         FUN = "*", check.margin = FALSE),
    #   FUN = "+", check.margin = FALSE)
    # z.e.obs * y.e.obs + k.e * z.e.lat * y.e.lat

    quad.first = ( m*sum(term.zkz) - (nu.theta-1) )^2
      # ( m*apply(term.zkz, 2, FUN = "sum") - (nu.theta-1) )^2
    quad.second = (4/lambda.theta) * sum(term.zkz.y)
      # (4/lambda.theta) * apply(term.zkz.y, 2, FUN = "sum")

    theta = (lambda.theta/2) * ( (nu.theta-1) - m*sum(term.zkz) + sqrt(quad.first + quad.second)  )
      # (lambda.theta/2) * ( (nu.theta-1) - m*apply(term.zkz, 2, FUN = "sum") + sqrt(quad.first + quad.second)  )

    return(theta)
  }
  Q.m = function(m, z.e.obs, z.e.lat, k.e, y.e.obs, y.e.lat, y.log.e.obs, y.log.e.lat, nu.m, lambda.m, nu.theta, lambda.theta)
  {
    term.zkz = XPlusYZ(z.e.obs, z.e.lat, k.e)
      # sweep(matrix(z.e.obs), 1,
      #                sweep(matrix(k.e), 1, matrix(z.e.lat), FUN = "*", check.margin = FALSE),
      #                FUN = "+", check.margin = FALSE)
      # z.e.obs + k.e * z.e.lat
    term.zkz.y = XAPlusYZB(z.e.obs, y.e.obs, z.e.lat, k.e, y.e.lat)
      # sweep(
      # sweep(matrix(z.e.obs), 1, matrix(y.e.obs), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.e.lat), 1, matrix(y.e.lat), FUN = "*", check.margin = FALSE),
      #       FUN = "*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
      # z.e.obs * y.e.obs + k.e * z.e.lat * y.e.lat
    term.zkz.logy = XAPlusYZB(z.e.obs, y.log.e.obs, z.e.lat, k.e, y.log.e.lat)
      # sweep(
      # sweep(matrix(z.e.obs), 1, matrix(y.log.e.obs), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.e.lat), 1, matrix(y.log.e.lat), FUN = "*", check.margin = FALSE),
      #       FUN = "*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
      # z.e.obs * y.log.e.obs + k.e * z.e.lat * y.log.e.lat

    theta = theta.m(m, z.e.obs, z.e.lat, k.e, y.e.obs, y.e.lat, nu.m, lambda.m, nu.theta, lambda.theta)

    temp = (m-1)*sum(term.zkz.logy) - (1/theta)*sum(term.zkz.y) - (m*log(theta) + lgamma(m))*sum(term.zkz) + ((nu.m-1)*log(m) - m/lambda.m) + ((nu.theta-1)*log(theta) - theta/lambda.theta)
      # (m-1)*apply(term.zkz.logy, 2, FUN = "sum") - (1/theta)*apply(term.zkz.y, 2, FUN = "sum") - (m*log(theta) + lgamma(m))*apply(term.zkz, 2, FUN = "sum") + ((nu.m-1)*log(m) - m/lambda.m) + ((nu.theta-1)*log(theta) - theta/lambda.theta)

    return(temp)
  }

  # I should only optimize using "positive" observations.
  pos.idx = (yu!=0)
  m.new = optimise(f = Q.m, interval = c(0.5*m, 10*m),
                   z.e.obs = z.e.obs[pos.idx], z.e.lat = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                   y.e.obs = y.e.obs[pos.idx], y.e.lat = y.e.lat[pos.idx], y.log.e.obs = y.log.e.obs[pos.idx], y.log.e.lat = y.log.e.lat[pos.idx],
                   nu.m = nu.m, lambda.m = lambda.m, nu.theta = nu.theta, lambda.theta = lambda.theta,
                   tol = .Machine$double.eps^0.05,
                   maximum = TRUE)$maximum
  theta.new = theta.m(m.new, z.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx], y.e.obs[pos.idx], y.e.lat[pos.idx],
                      nu.m = nu.m, lambda.m = lambda.m, nu.theta = nu.theta, lambda.theta = lambda.theta)

  # Update the parameters and return the result
  params.new[1] = m.new
  params.new[2] = theta.new
  return(params.new)

}
