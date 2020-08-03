## ECM algorithm of Weibull expert
#' ECM: M-Step for Weibull expert.
#'
#'
#' @importFrom stats dweibull optimize integrate
#' @importFrom expint gammainc
#'
#' @keywords internal
#'
#' @export EMMWeibull
EMMWeibull = function(params.old,
                          tl, yl, yu, tu,
                          expert.ll, expert.tn, expert.tn.bar,
                          z.e.obs, z.e.lat, k.e,
                          penalty, hyper.params)
{
  # Constants
  sample.size.n = length(tl)
  # Take old parameters
  shape.k = params.old[1]
  scale.lambda = params.old[2]
  # Take hyper parameters
  hyper.k.1 = hyper.params[1]
  hyper.k.2 = hyper.params[2]
  hyper.lambda.1   = hyper.params[3]
  hyper.lambda.2   = hyper.params[4]
  # Value to return
  params.new = params.old


  # E-Step: Conditional expectation of log(y)
  censor.idx = (yl!=yu)
  y.log.e.obs = array(0, dim = c(sample.size.n, 1))
  y.log.e.lat = array(0, dim = c(sample.size.n, 1))

  # Untruncated and uncensored case
  y.log.e.obs[!censor.idx, 1] = log(yl[!censor.idx])

  # Function to numerically integrate for censored and truncated case
  int.y.log.fcn = function(u, shape.k.j, scale.lambda.j) # u = log(y) for numerical stability
  {
    u.density.log = stats::dweibull(exp(u), shape = shape.k.j, scale = scale.lambda.j, log = TRUE) + u
    # Numerical underflow of dweibull for small u, and certain set of parameters. e.g. (4, 35)
    temp = u*exp(u.density.log)
    temp[which(u==Inf)] = 0
    temp[which(u==-Inf)] = 0
    temp[which(is.na(temp))] = 0
    return(temp)
  }

  # First find unique upper and lower bounds of integration: untruncated case
  y.unique = unique(cbind(yl,yu),MARGIN=1)
  y.unique.length = nrow(y.unique)
  y.unique.match = match(data.frame(t(cbind(yl,yu))),data.frame(t(y.unique)))
  yl.unique = y.unique[,1]
  yu.unique = y.unique[,2]

  # Integration: untruncated case
  y.log.e.obs.unique = array(0,dim=c(y.unique.length, 1))

  y.log.e.obs.unique[,1]=
    mapply(function(x, y) ifelse(x!=y,
                                 integrate(int.y.log.fcn, log(x), log(y),
                                           shape.k.j=shape.k, scale.lambda.j = scale.lambda,
                                           rel.tol=.Machine$double.eps^0.5)$value,
                                 0),
           yl.unique, yu.unique)
  # Match to all observations of y
  temp.y.log.e.obs = array(0, dim = c(sample.size.n, 1))
  temp.y.log.e.obs = y.log.e.obs.unique[y.unique.match,]
  # Only take out those values of censored cases
  y.log.e.obs[censor.idx,1] = exp(-expert.ll[censor.idx]) * temp.y.log.e.obs[censor.idx]

  # Similarly, find unique upper and lower bounds of integration: truncated case
  tn.unique = unique(cbind(tl,tu),MARGIN=1)
  tn.unique.length = nrow(tn.unique)
  tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
  tl.unique = tn.unique[,1]
  tu.unique = tn.unique[,2]

  # Integration
  y.log.e.lat.unique = array(0,dim=c(tn.unique.length, 1))

  y.log.e.lat.unique[,1]=
    sapply(tl.unique,
           function(x) ifelse(x!=0,
                              integrate(int.y.log.fcn, -Inf, log(x),
                                        shape.k.j=shape.k, scale.lambda.j = scale.lambda,
                                        rel.tol=.Machine$double.eps^0.5)$value,
                              0))+
    sapply(tu.unique,
           function(x) ifelse(x!=Inf,
                              integrate(int.y.log.fcn, log(x), stats::qweibull(1-1e-09, shape = shape.k, scale = scale.lambda, log.p = F), # Inf is a problem
                                        shape.k.j=shape.k, scale.lambda.j = scale.lambda,
                                        rel.tol=.Machine$double.eps^0.5)$value,
                              0))
  # Match to all observations of y
  temp.y.log.e.lat = array(0, dim = c(sample.size.n, 1))
  temp.y.log.e.lat = y.log.e.lat.unique[tn.unique.match,]
  # Conditional expectation of log(y): truncated
  y.log.e.lat = exp(-expert.tn.bar) * temp.y.log.e.lat

  # Finally, get rid of NaN's
  y.log.e.obs[is.nan(y.log.e.obs)] = 0
  y.log.e.lat[is.nan(y.log.e.lat)] = 0


  # M-Step: Maximization of loglik
  # Optimal scale.lambda is a function of shape.k
  scale.shape = function(shape.k, z.e.obs, z.e.lat, k.e, y.pow.e.obs, y.pow.e.lat)
  {
    term.zkz = XPlusYZ(z.e.obs, z.e.lat, k.e)
      # sweep(matrix(z.e.obs), 1,
      #                sweep(matrix(k.e), 1, matrix(z.e.lat), FUN = "*", check.margin = FALSE),
      #                FUN = "+", check.margin = FALSE)
    # z.e.obs + k.e * z.e.lat
    term.zkz.ypow = XAPlusYZB(z.e.obs, y.pow.e.obs, z.e.lat, k.e, y.pow.e.lat)
      # sweep(
      # sweep(matrix(z.e.obs), 1, matrix(y.pow.e.obs), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.e.lat), 1, matrix(y.pow.e.lat), FUN = "*", check.margin = FALSE),
      #       FUN = "*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.e.obs * y.pow.e.obs + k.e * z.e.lat * y.pow.e.lat

    temp.scale = ( sum(term.zkz.ypow)/sum(term.zkz) )^(1/shape.k)
      # ( apply(term.zkz.ypow, 2, FUN = "sum")/apply(term.zkz, 2, FUN = "sum") )^(1/shape.k)
    return(temp.scale)
  }

  Q.shape = function(shape.k.new, shape.k.old, scale.lambda.old,
                     z.e.obs, z.e.lat, k.e,
                     y.log.e.obs, y.log.e.lat,
                     tl, yl, yu, tu, expert.ll, expert.tn, expert.tn.bar,
                     penalty, nu.shape, theta.shape)
  {
    # E-Step: conditional expectation of y.pow
    sample.size.n = length(yu)
    censor.idx = (yl!=yu)
    y.pow.e.obs = array(0, dim = c(sample.size.n, 1))
    y.pow.e.lat = array(0, dim = c(sample.size.n, 1))

    # Untruncated and uncensored case
    y.pow.e.obs[!censor.idx] = (yu[!censor.idx])^(shape.k.new)

    # Untrencated and uncensored case
    y.pow.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * scale.lambda.old^(shape.k.new) *
      (gammainc((shape.k.old+shape.k.new)/shape.k.old, (yl[censor.idx]/scale.lambda.old)^(shape.k.old)) -
         gammainc((shape.k.old+shape.k.new)/shape.k.old, (yu[censor.idx]/scale.lambda.old)^(shape.k.old)) )

    # Truncated case
    y.pow.e.lat = exp(-expert.tn.bar) * scale.lambda.old^(shape.k.new) *
      (gammainc((shape.k.old+shape.k.new)/shape.k.old, (0 /scale.lambda.old)^(shape.k.old)) -
         gammainc((shape.k.old+shape.k.new)/shape.k.old, (tl/scale.lambda.old)^(shape.k.old)) +
         gammainc((shape.k.old+shape.k.new)/shape.k.old, (tu/scale.lambda.old)^(shape.k.old)) -
         gammainc((shape.k.old+shape.k.new)/shape.k.old, (Inf/scale.lambda.old)^(shape.k.old)))

    # Finally, get rid of NaN's
    y.pow.e.obs[is.nan(y.pow.e.obs)] = 0
    y.pow.e.lat[is.nan(y.pow.e.lat)] = 0

    scale.shape.new = scale.shape(shape.k.new, z.e.obs, z.e.lat, k.e, y.pow.e.obs, y.pow.e.lat)

    # Compute the function calue to optimize
    term.zkz = XPlusYZ(z.e.obs, z.e.lat, k.e)
      # sweep(matrix(z.e.obs), 1,
      #                sweep(matrix(k.e), 1, matrix(z.e.lat), FUN = "*", check.margin = FALSE),
      #                FUN = "+", check.margin = FALSE)
    # z.e.obs + k.e * z.e.lat
    term.zkz.ylog = XAPlusYZB(z.e.obs, y.log.e.obs, z.e.lat, k.e, y.log.e.lat)
      # sweep(
      # sweep(matrix(z.e.obs), 1, matrix(y.log.e.obs), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.e.lat), 1, matrix(y.log.e.lat), FUN = "*", check.margin = FALSE),
      #       FUN = "*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.e.obs * y.log.e.obs + k.e * z.e.lat * y.log.e.lat
    term.zkz.ypow = XAPlusYZB(z.e.obs, y.pow.e.obs, z.e.lat, k.e, y.pow.e.lat)
      # sweep(
      # sweep(matrix(z.e.obs), 1, matrix(y.pow.e.obs), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.e.lat), 1, matrix(y.pow.e.lat), FUN = "*", check.margin = FALSE),
      #       FUN = "*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.e.obs * y.pow.e.obs + k.e * z.e.lat * y.pow.e.lat

    result = log(shape.k.new)*sum(term.zkz) - shape.k.new*log(scale.shape.new)*sum(term.zkz) + (shape.k.new-1)*sum(term.zkz.ylog) - scale.shape.new^(-shape.k.new)*sum(term.zkz.ypow)
      # log(shape.k.new)*apply(term.zkz, 2, FUN = "sum") - shape.k.new*log(scale.shape.new)*apply(term.zkz, 2, FUN = "sum") + (shape.k.new-1)*apply(term.zkz.ylog, 2, FUN = "sum") - scale.shape.new^(-shape.k.new)*apply(term.zkz.ypow, 2, FUN = "sum")

    if(penalty==TRUE)
    {
      result = result + (nu.shape-1)*log(shape.k.new) - shape.k.new/theta.shape
    }

    return(result)

  }

  # I should only use "positive" observations.
  pos.idx = (yu!=0)

  # Update of shape.k
  shape.k.new = optimize(f = Q.shape, lower = max(0.5*shape.k, 1.0000), upper = 2*shape.k,
                         # lower = max(0.5*shape.k, 1.0000), upper = 5*shape.k, # interval = c(0.5*shape.k, 5*shape.k),
                         shape.k.old = shape.k, scale.lambda.old = scale.lambda,
                         z.e.obs = z.e.obs[pos.idx], z.e.lat = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                         y.log.e.obs = y.log.e.obs[pos.idx], y.log.e.lat = y.log.e.lat[pos.idx],
                         tl = tl[pos.idx], yl = yl[pos.idx], yu = yu[pos.idx], tu = tu[pos.idx],
                         expert.ll = expert.ll[pos.idx], expert.tn = expert.tn[pos.idx], expert.tn.bar = expert.tn.bar[pos.idx],
                         penalty = penalty, nu.shape = hyper.k.1, theta.shape = hyper.k.2,
                         maximum = TRUE)$maximum

  # Update of scale.lambda
  shape.k.old = shape.k
  scale.lambda.old = scale.lambda

  censor.idx = (yl!=yu)
  y.pow.e.obs = array(0, dim = c(sample.size.n, 1))
  y.pow.e.lat = array(0, dim = c(sample.size.n, 1))

  # Untruncated and uncensored case
  y.pow.e.obs[!censor.idx] = (yu[!censor.idx])^(shape.k.new)

  # Untrencated and uncensored case
  y.pow.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * scale.lambda.old^(shape.k.new) *
    (gammainc((shape.k.old+shape.k.new)/shape.k.old, (yl[censor.idx]/scale.lambda.old)^(shape.k.old)) -
       gammainc((shape.k.old+shape.k.new)/shape.k.old, (yu[censor.idx]/scale.lambda.old)^(shape.k.old)) )

  # Truncated case
  y.pow.e.lat = exp(-expert.tn.bar) * scale.lambda.old^(shape.k.new) *
    (gammainc((shape.k.old+shape.k.new)/shape.k.old, (0 /scale.lambda.old)^(shape.k.old)) -
       gammainc((shape.k.old+shape.k.new)/shape.k.old, (tl/scale.lambda.old)^(shape.k.old)) +
       gammainc((shape.k.old+shape.k.new)/shape.k.old, (tu/scale.lambda.old)^(shape.k.old)) -
       gammainc((shape.k.old+shape.k.new)/shape.k.old, (Inf/scale.lambda.old)^(shape.k.old)))

  # Finally, get rid of NaN's
  y.pow.e.obs[is.nan(y.pow.e.obs)] = 0
  y.pow.e.lat[is.nan(y.pow.e.lat)] = 0

  scale.lambda.new = scale.shape(shape.k.new, z.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx], y.pow.e.obs[pos.idx], y.pow.e.lat[pos.idx])

  # Return results
  params.new[1] = shape.k.new
  params.new[2] = scale.lambda.new
  return(params.new)

}
