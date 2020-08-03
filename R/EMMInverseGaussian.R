## ECM algorithm of Inverse Gaussian expert
#' ECM: M-Step for Inverse Gaussian expert.
#'
#'
#' @importFrom stats optimise integrate
#' @importFrom statmod dinvgauss
#'
#' @keywords internal
#'
#' @export EMMInverseGaussian
EMMInverseGaussian = function(params.old,
                                   tl, yl, yu, tu,
                                   expert.ll, expert.tn, expert.tn.bar,
                                   z.e.obs, z.e.lat, k.e,
                                   penalty, hyper.params)
{
  # Constants
  sample.size.n = length(tl)
  # Take old parameters
  mean.mu = params.old[1]
  shape.lambda = params.old[2]
  # Take hyper parameters
  hyper.mean.mu.1 = hyper.params[1]
  hyper.mean.mu.2 = hyper.params[2]
  hyper.nu.lambda.1 = hyper.params[3]
  hyper.nu.lambda.2 = hyper.params[4]
  # Value to return
  params.new = params.old

  # E-Step: conditional expectations for y, log(y), 1/y
  censor.idx = (yl!=yu)
  y.e.obs = array(0, dim = c(sample.size.n, 1))
  y.e.lat = array(0, dim = c(sample.size.n, 1))
  y.log.e.obs = array(0, dim = c(sample.size.n, 1))
  y.log.e.lat = array(0, dim = c(sample.size.n, 1))
  y.inv.e.obs = array(0, dim = c(sample.size.n, 1))
  y.inv.e.lat = array(0, dim = c(sample.size.n, 1))

  # Conditional expectation of y: untruncated and uncensored case
  y.e.obs[!censor.idx, 1] = yl[!censor.idx]
  # Conditional expectation of log(y): untruncated and uncensored case
  y.log.e.obs[!censor.idx, 1] = log(yl[!censor.idx])
  # Conditional expectation of inv(y): untruncated and uncensored case
  y.inv.e.obs[!censor.idx, 1] = 1/yl[!censor.idx]

  # Conditional expectation of y, log(y) and inv(y): untruncated but cencored case
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

  # First find unique upper and lower bounds of integration
  y.unique = unique(cbind(yl,yu),MARGIN=1)
  y.unique.length = nrow(y.unique)
  y.unique.match = match(data.frame(t(cbind(yl,yu))),data.frame(t(y.unique)))
  yl.unique = y.unique[,1]
  yu.unique = y.unique[,2]

  # Integration
  y.e.obs.unique = array(0,dim=c(y.unique.length,1))
  y.log.e.obs.unique = array(0,dim=c(y.unique.length,1))
  y.inv.e.obs.unique = array(0,dim=c(y.unique.length,1))

  y.e.obs.unique[,1]= intInvGaussYObs(mean.mu, shape.lambda, log(yl.unique), log(yu.unique))
    # mapply(function(x, y) ifelse(x!=y,
    #                              integrate(int.y.fcn, log(x), log(y),
    #                                        mean.mu=mean.mu, shape.lambda=shape.lambda,
    #                                        rel.tol=.Machine$double.eps^0.5)$value,
    #                              0),
    #        yl.unique, yu.unique)

  y.log.e.obs.unique[,1]= intInvGaussLogYObs(mean.mu, shape.lambda, log(yl.unique), log(yu.unique))
    # mapply(function(x, y) ifelse(x!=y,
    #                              integrate(int.y.log.fcn, log(x), log(y),
    #                                        mean.mu=mean.mu, shape.lambda=shape.lambda,
    #                                        rel.tol=.Machine$double.eps^0.5)$value,
    #                              0),
    #        yl.unique, yu.unique)

  y.inv.e.obs.unique[,1]= intInvGaussInvYObs(mean.mu, shape.lambda, log(yl.unique), log(yu.unique))
    # mapply(function(x, y) ifelse(x!=y,
    #                              integrate(int.y.inv.fcn, log(x), log(y),
    #                                        mean.mu=mean.mu, shape.lambda=shape.lambda,
    #                                        rel.tol=.Machine$double.eps^0.5)$value,
    #                              0),
    #        yl.unique, yu.unique)

  # Match to all observations of y
  temp.y.e.obs = temp.y.log.e.obs = temp.y.inv.e.obs = array(0, dim = c(sample.size.n, 1))
  temp.y.e.obs = y.e.obs.unique[y.unique.match,]
  temp.y.log.e.obs = y.log.e.obs.unique[y.unique.match,]
  temp.y.inv.e.obs = y.inv.e.obs.unique[y.unique.match,]

  # Conditional expectation of y: untruncated and censored case
  y.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * temp.y.e.obs[censor.idx]
  # Conditional expectation of log(y): untruncated and censored case
  y.log.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * temp.y.log.e.obs[censor.idx]
  # Conditional expectation of inv(y): untruncated and censored case
  y.inv.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * temp.y.inv.e.obs[censor.idx]


  # Conditional expectation of y, log(y) and inv(y): truncated
  # Reuse the integration functions above
  tn.unique = unique(cbind(tl,tu),MARGIN=1)
  tn.unique.length = nrow(tn.unique)
  tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
  tl.unique = tn.unique[,1]
  tu.unique = tn.unique[,2]

  # Integration
  y.e.lat.unique = array(0,dim=c(tn.unique.length,1))
  y.log.e.lat.unique = array(0,dim=c(tn.unique.length,1))
  y.inv.e.lat.unique = array(0,dim=c(tn.unique.length,1))

  y.e.lat.unique[,1]= intInvGaussYLat(mean.mu, shape.lambda, log(tl.unique), log(tu.unique))
    # sapply(tl.unique,
    #        function(x) ifelse(x!=0,
    #                           integrate(int.y.fcn, -Inf, log(x),
    #                                     mean.mu=mean.mu, shape.lambda=shape.lambda,
    #                                     rel.tol=.Machine$double.eps^0.5)$value,
    #                           0))+
    # sapply(tu.unique,
    #        function(x) ifelse(x!=Inf,
    #                           integrate(int.y.fcn, log(x), Inf,
    #                                     mean.mu=mean.mu, shape.lambda=shape.lambda,
    #                                     rel.tol=.Machine$double.eps^0.5)$value,
    #                           0))

  y.log.e.lat.unique[,1]= intInvGaussLogYLat(mean.mu, shape.lambda, log(tl.unique), log(tu.unique))
    # sapply(tl.unique,
    #        function(x) ifelse(x!=0,
    #                           integrate(int.y.log.fcn, -Inf, log(x),
    #                                     mean.mu=mean.mu, shape.lambda=shape.lambda,
    #                                     rel.tol=.Machine$double.eps^0.5)$value,
    #                           0))+
    # sapply(tu.unique,
    #        function(x) ifelse(x!=Inf,
    #                           integrate(int.y.log.fcn, log(x), Inf,
    #                                     mean.mu=mean.mu, shape.lambda=shape.lambda,
    #                                     rel.tol=.Machine$double.eps^0.5)$value,
    #                           0))

  y.inv.e.lat.unique[,1]= intInvGaussInvYLat(mean.mu, shape.lambda, log(tl.unique), log(tu.unique))
    # sapply(tl.unique,
    #        function(x) ifelse(x!=0,
    #                           integrate(int.y.inv.fcn, -Inf, log(x),
    #                                     mean.mu=mean.mu, shape.lambda=shape.lambda,
    #                                     rel.tol=.Machine$double.eps^0.5)$value,
    #                           0))+
    # sapply(tu.unique,
    #        function(x) ifelse(x!=Inf,
    #                           integrate(int.y.inv.fcn, log(x), Inf,
    #                                     mean.mu=mean.mu, shape.lambda=shape.lambda,
    #                                     rel.tol=.Machine$double.eps^0.5)$value,
    #                           0))

  # Match to all observations of y
  temp.y.e.lat = temp.y.log.e.lat = temp.y.inv.e.lat = array(0, dim = c(sample.size.n, 1))
  temp.y.e.lat = y.e.lat.unique[tn.unique.match,]
  temp.y.log.e.lat = y.log.e.lat.unique[tn.unique.match,]
  temp.y.inv.e.lat = y.inv.e.lat.unique[tn.unique.match,]

  # Conditional expectation of y: truncated
  y.e.lat = exp(-expert.tn.bar) * temp.y.e.lat
  # Conditional expectation of log(y): truncated
  y.log.e.lat = exp(-expert.tn.bar) * temp.y.log.e.lat
  # Conditional expectation of inv(y): truncated
  y.inv.e.lat = exp(-expert.tn.bar) * temp.y.inv.e.lat

  # Finally, get rid of NaN's
  y.e.obs[is.na(y.e.obs)] = 0
  y.e.lat[is.na(y.e.lat)] = 0
  y.log.e.obs[is.na(y.log.e.obs)] = 0
  y.log.e.lat[is.na(y.log.e.lat)] = 0
  y.inv.e.obs[is.na(y.inv.e.obs)] = 0
  y.inv.e.lat[is.na(y.inv.e.lat)] = 0


  ## M-step: mean.mu and shape.lambda
  # Closed-form solution exists, if there is no penalty
  mu.update = function(z.obs.j, z.lat.j, k.e, y.e.obs.j, y.e.lat.j)
  {
    sum.one = XPlusYZ(z.obs.j, z.lat.j, k.e)
      # sweep(matrix(z.obs.j), 1,
      #               sweep(matrix(k.e), 1, matrix(z.lat.j), FUN = "*", check.margin = FALSE),
      #               FUN = "+", check.margin = FALSE)
    # z.obs.j + k.e * z.lat.j
    sum.two = XAPlusYZB(z.obs.j, y.e.obs.j, z.lat.j, k.e, y.e.lat.j)
      # sweep(
      # sweep(matrix(z.obs.j), 1, matrix(y.e.obs.j), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.lat.j), 1, matrix(y.e.lat.j), FUN = "*", check.margin = FALSE),
      #       FUN = "*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.obs.j * y.e.obs.j + k.e * z.lat.j * y.e.lat.j

    return( sum(sum.two)/sum(sum.one) )
    # return( apply(sum.two, 2, FUN = "sum")/apply(sum.one, 2, FUN = "sum") )
  }

  lambda.mu = function(mean.mu, z.obs.j, z.lat.j, k.e, y.e.obs.j, y.e.lat.j, y.inv.e.obs.j, y.inv.e.lat.j)
  {
    sum.one = XPlusYZ(z.obs.j, z.lat.j, k.e)
      # sweep(matrix(z.obs.j), 1,
      #               sweep(matrix(k.e), 1, matrix(z.lat.j), FUN = "*", check.margin = FALSE),
      #               FUN = "+", check.margin = FALSE)
    # z.obs.j + k.e * z.lat.j
    sum.two = XAPlusYZB(z.obs.j, y.e.obs.j, z.lat.j, k.e, y.e.lat.j)
      # sweep(
      # sweep(matrix(z.obs.j), 1, matrix(y.e.obs.j), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.lat.j), 1, matrix(y.e.lat.j), FUN = "*", check.margin = FALSE),
      #       FUN = "*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.obs.j * y.e.obs.j + k.e * z.lat.j * y.e.lat.j
    sum.three = XAPlusYZB(z.obs.j, y.inv.e.obs.j, z.lat.j, k.e, y.inv.e.lat.j)
      # sweep(
      # sweep(matrix(z.obs.j), 1, matrix(y.inv.e.obs.j), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.lat.j), 1, matrix(y.inv.e.lat.j), FUN = "*", check.margin = FALSE),
      #       FUN = "*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.obs.j * y.inv.e.obs.j + k.e * z.lat.j * y.inv.e.lat.j
    shape.lambda.temp = sum(sum.one) / ( (1/(mean.mu^2))*sum(sum.two) - (2/(mean.mu))*sum(sum.one) + sum(sum.three) )
      # apply(sum.one, 2, FUN = "sum") / ( (1/(mean.mu^2))*apply(sum.two, 2, FUN = "sum") - (2/(mean.mu))*apply(sum.one, 2, FUN = "sum") + apply(sum.three, 2, FUN = "sum") )

    return(shape.lambda.temp)
  }

  # I should only maximize using the "positive" observations
  pos.idx = (yu!=0)

  mean.mu.new = mu.update(z.obs.j = z.e.obs[pos.idx], z.lat.j = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                          y.e.obs.j = y.e.obs[pos.idx], y.e.lat.j = y.e.lat[pos.idx])

  shape.lambda.new = lambda.mu(mean.mu.new,
                               z.obs.j = z.e.obs[pos.idx], z.lat.j = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                               y.e.obs.j = y.e.obs[pos.idx], y.e.lat.j = y.e.lat[pos.idx], y.inv.e.obs.j = y.inv.e.obs[pos.idx], y.inv.e.lat.j = y.inv.e.lat[pos.idx])


  # Update the parameters and return the result
  params.new[1] = mean.mu.new
  params.new[2] = shape.lambda.new
  return(params.new)

}
