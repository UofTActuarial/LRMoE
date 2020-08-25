## ECM algorithm of Inverse Gaussian expert
#' ECM: M-Step for Inverse Gaussian expert.
#'
#' @importFrom stats optim integrate
#' @importFrom actuar dburr qburr pburr
#'
#' @keywords internal
#'
#' @export EMMBurr
EMMBurr = function(params.old,
                     tl, yl, yu, tu,
                     expert.ll, expert.tn, expert.tn.bar,
                     z.e.obs, z.e.lat, k.e,
                     penalty, hyper.params)
{
  # Constants
  sample.size.n = length(tl)
  # Take old parameters
  shape1.k = params.old[1]
  shape2.c = params.old[2]
  scale.lambda = params.old[3]
  # Take hyper parameters
  hyper.k.1 = hyper.params[1]
  hyper.k.2 = hyper.params[2]
  hyper.c.1 = hyper.params[3]
  hyper.c.2 = hyper.params[4]
  hyper.lambda.1 = hyper.params[5]
  hyper.lambda.2 = hyper.params[6]
  # Value to return
  params.new = params.old


  # E-Step: conditional expectation of log(y)
  censor.idx = (yl!=yu)
  y.log.e.obs = array(0, dim = c(sample.size.n ,1))
  y.log.e.lat = array(0, dim = c(sample.size.n, 1))

  # Untruncated and uncensored case
  y.log.e.obs[!censor.idx, 1] = log(yu[!censor.idx])

  # Function to numerically integrate for censored and truncated case
  int.y.log.fcn = function(u, shape1.k.j, shape2.c.j, scale.lambda.j) # u = log(y) for numerical stability
  {
    u.density.log = u + actuar::dburr(exp(u), shape1 = shape1.k.j, shape2 = shape2.c.j, scale = scale.lambda.j, log = TRUE)
    temp = u * exp(u.density.log)
    temp[which(u==Inf)] = 0
    temp[which(u==-Inf)] = 0
    temp[which(is.na(temp))] = 0

    temp[which(temp==Inf)] = 0

    return(temp)
  }

  # First find unique upper and lower bounds of integration: untruncated case
  y.unique = unique(cbind(yl,yu),MARGIN=1)
  y.unique.length = nrow(y.unique)
  y.unique.match = match(data.frame(t(cbind(yl,yu))),data.frame(t(y.unique)))
  yl.unique = y.unique[,1]
  yu.unique = y.unique[,2]

  # Integration
  y.log.e.obs.unique = array(0,dim=c(y.unique.length,1))

  y.log.e.obs.unique[,1]= # intBurrLogYObs(shape1.k, shape2.c, scale.lambda, log(yl.unique), log(yu.unique))
    mapply(function(x, y) ifelse(x!=y,
                                 integrate(int.y.log.fcn, log(x), log(y),
                                           shape1.k.j = shape1.k, shape2.c.j = shape2.c, scale.lambda.j = scale.lambda,
                                           rel.tol=.Machine$double.eps^0.5)$value,
                                 0),
           yl.unique, yu.unique)
  # Match to all observations of y
  temp.y.log.e.obs = array(0, dim = c(sample.size.n, 1))
  temp.y.log.e.obs = y.log.e.obs.unique[y.unique.match,]

  y.log.e.obs[censor.idx, 1] = exp(-expert.ll[censor.idx]) * temp.y.log.e.obs[censor.idx]

  # Similarly, find unique upper and lower bounds of integration: truncated case
  tn.unique = unique(cbind(tl,tu),MARGIN=1)
  tn.unique.length = nrow(tn.unique)
  tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
  tl.unique = tn.unique[,1]
  tu.unique = tn.unique[,2]

  # Integration
  y.log.e.lat.unique = array(0,dim=c(tn.unique.length, 1))

  y.log.e.lat.unique[,1]= # intBurrLogYObs(shape1.k, shape2.c, scale.lambda, log(tl.unique), log(tu.unique))
    sapply(tl.unique,
           function(x) ifelse(x!=0,
                              integrate(int.y.log.fcn, -Inf, log(x),
                                        shape1.k.j = shape1.k, shape2.c.j = shape2.c, scale.lambda.j = scale.lambda,
                                        rel.tol=.Machine$double.eps^0.5)$value,
                              0))+
    sapply(tu.unique,
           function(x) ifelse(x!=Inf,
                              integrate(int.y.log.fcn, log(x), Inf, # actuar::qburr(1-1e-09, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = T, log.p = F), # Inf causes a problem
                                        shape1.k.j = shape1.k, shape2.c.j = shape2.c, scale.lambda.j = scale.lambda,
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

  # M-Step: maximization of loglik
  # Optimal shape1.k is a function of y.pol
  shape1.logpow = function(z.e.obs, z.e.lat, k.e, y.pol.e.obs, y.pol.e.lat, penalty, hyper.k.1, hyper.k.2)
  {
    Inf.idx = which(y.pol.e.obs==-Inf)
    y.pol.e.obs[Inf.idx] = 0
    Inf.idx = which(y.pol.e.lat==-Inf)
    y.pol.e.lat[Inf.idx] = 0

    Inf.idx = which(y.pol.e.obs==Inf)
    y.pol.e.obs[Inf.idx] = 0
    Inf.idx = which(y.pol.e.lat==Inf)
    y.pol.e.lat[Inf.idx] = 0

    term.zkz = XPlusYZ(z.e.obs, z.e.lat, k.e)
      # sweep(matrix(z.e.obs), 1,
      #                sweep(matrix(k.e), 1, matrix(z.e.lat), FUN = "*", check.margin = FALSE),
      #                FUN = "+", check.margin = FALSE)
    # z.e.obs + k.e * z.e.lat
    term.zkz.ypol = XAPlusYZB(z.e.obs, y.pol.e.obs, z.e.lat, k.e, y.pol.e.lat)
      # sweep(
      # sweep(matrix(z.e.obs), 1, matrix(y.pol.e.obs), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.e.lat), 1, matrix(y.pol.e.lat), FUN = "*", check.margin = FALSE),
      #       FUN ="*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.e.obs * y.pol.e.obs + k.e * z.e.lat * y.pol.e.lat

    result = sum(term.zkz)/sum(term.zkz.ypol)
      # apply(term.zkz, 2, FUN = "sum")/apply(term.zkz.ypol, 2, FUN = "sum")

    if(penalty==TRUE){
      result = (sum(term.zkz) +(hyper.k.1-1)) / (sum(term.zkz.ypol) + 1/hyper.k.2 )
        # (apply(term.zkz, 2, FUN = "sum") +(hyper.k.1-1)) / (apply(term.zkz.ypol, 2, FUN = "sum") + 1/hyper.k.2 )
    }
    return(result)
  }
  # A function to numerically integrate for y.pol
  int.y.pol.fcn = function(u, shape2.c.new, scale.lambda.new,
                           shape1.k.j, shape2.c.j, scale.lambda.j) # u = log(y) for numerical stability
  {
    # Use existing package so that I don't have to manually type in the density!
    u.density.log = u + actuar::dburr(exp(u), shape1 = shape1.k.j, shape2 = shape2.c.j, scale = scale.lambda.j, log = TRUE)
    temp = log(1+(exp(u)/scale.lambda.new)^shape2.c.new) * exp(u.density.log)
    temp[which(u==Inf)] = 0
    temp[which(u==-Inf)] = 0
    temp[which(is.na(temp))] = 0

    temp[which(temp==Inf)] = 0
    return(temp)
  }
  # A function to numerically optimize
  Q.T = function(params.new, shape1.k.old, shape2.c.old, scale.lambda.old,
                 z.e.obs, z.e.lat, k.e,
                 y.log.e.obs, y.log.e.lat,
                 tl, yl, yu, tu, expert.ll, expert.tn, expert.tn.bar,
                 penalty,
                 hyper.k.1, hyper.k.2, hyper.c.1, hyper.c.2, hyper.lambda.1, hyper.lambda.2)
  {
    # Constants
    sample.size.n = length(tl)
    censor.idx = (yl!=yu)

    # shape2.c.new = params.new[1]
    # scale.lambda.new = params.new[2]

    shape2.c.new = exp(params.new[1])
    scale.lambda.new = exp(params.new[2])

    # E-Step: conditional expectation of y.pol
    y.pol.e.obs = array(0, dim = c(sample.size.n ,1))
    y.pol.e.lat = array(0, dim = c(sample.size.n, 1))

    # Untruncated and uncensored case
    y.pol.e.obs[!censor.idx] = log1p((yu[!censor.idx]/scale.lambda.new)^shape2.c.new)
      # log(1+(yu[!censor.idx]/scale.lambda.new)^shape2.c.new)

    # First find unique upper and lower bounds of integration: untruncated case
    y.unique = unique(cbind(yl,yu),MARGIN=1)
    y.unique.length = nrow(y.unique)
    y.unique.match = match(data.frame(t(cbind(yl,yu))),data.frame(t(y.unique)))
    yl.unique = y.unique[,1]
    yu.unique = y.unique[,2]

    # Integration
    y.pol.e.obs.unique = array(0,dim=c(y.unique.length,1))

    y.pol.e.obs.unique[,1]= # intBurrPolYObs(shape1.k.old, shape2.c.old, scale.lambda.old,
      #                                      shape2.c.new, scale.lambda.new, log(yl.unique), log(yu.unique))
      mapply(function(x, y) ifelse(x!=y,
                                   integrate(int.y.pol.fcn, log(x), log(y),
                                             shape2.c.new = shape2.c.new, scale.lambda.new = scale.lambda.new,
                                             shape1.k.j = shape1.k.old, shape2.c.j = shape2.c.old, scale.lambda.j = scale.lambda.old,
                                             rel.tol=.Machine$double.eps^0.5)$value,
                                   0),
             yl.unique, yu.unique)
    # Match to all observations of y
    temp.y.pol.e.obs = array(0, dim = c(sample.size.n, 1))
    temp.y.pol.e.obs = y.pol.e.obs.unique[y.unique.match,]

    y.pol.e.obs[censor.idx, 1] = exp(-expert.ll[censor.idx]) * temp.y.pol.e.obs[censor.idx]

    # Similarly, find unique upper and lower bounds of integration: truncated case
    tn.unique = unique(cbind(tl,tu),MARGIN=1)
    tn.unique.length = nrow(tn.unique)
    tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
    tl.unique = tn.unique[,1]
    tu.unique = tn.unique[,2]

    # Integration
    y.pol.e.lat.unique = array(0,dim=c(tn.unique.length, 1))

    y.pol.e.lat.unique[,1]=# intBurrPolYLat(shape1.k.old, shape2.c.old, scale.lambda.old,
                           #                shape2.c.new, scale.lambda.new, log(tl.unique), log(tu.unique))
      sapply(tl.unique,
             function(x) ifelse(x!=0,
                                integrate(int.y.pol.fcn, -Inf, log(x),
                                          shape2.c.new = shape2.c.new, scale.lambda.new = scale.lambda.new,
                                          shape1.k.j = shape1.k.old, shape2.c.j = shape2.c.old, scale.lambda.j = scale.lambda.old,
                                          rel.tol=.Machine$double.eps^0.5)$value,
                                0))+
      sapply(tu.unique,
             function(x) ifelse(x!=Inf,
                                integrate(int.y.pol.fcn, log(x), Inf, # actuar::qburr(1-1e-09, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = T, log.p = F), # Inf causes a problem
                                          shape2.c.new = shape2.c.new, scale.lambda.new = scale.lambda.new,
                                          shape1.k.j = shape1.k.old, shape2.c.j = shape2.c.old, scale.lambda.j = scale.lambda.old,
                                          rel.tol=.Machine$double.eps^0.5)$value,
                                0))
    # Match to all observations of y
    temp.y.pol.e.lat = array(0, dim = c(sample.size.n, 1))
    temp.y.pol.e.lat = y.pol.e.lat.unique[tn.unique.match,]
    # Conditional expectation of log(y): truncated
    y.pol.e.lat = exp(-expert.tn.bar) * temp.y.pol.e.lat

    # Finally, get rid of NaN's
    y.pol.e.obs[is.nan(y.pol.e.obs)] = 0
    y.pol.e.lat[is.nan(y.pol.e.lat)] = 0

    # Prevent z.e.obs * y.log.e.obs  = 0*(-Inf) case
    Inf.idx = which(y.log.e.obs==-Inf)
    y.log.e.obs[Inf.idx] = 0
    Inf.idx = which(y.log.e.lat==-Inf)
    y.log.e.lat[Inf.idx] = 0
    Inf.idx = which(y.pol.e.obs==-Inf)
    y.pol.e.obs[Inf.idx] = 0
    Inf.idx = which(y.pol.e.lat==-Inf)
    y.pol.e.lat[Inf.idx] = 0

    Inf.idx = which(y.log.e.obs==Inf)
    y.log.e.obs[Inf.idx] = 0
    Inf.idx = which(y.log.e.lat==Inf)
    y.log.e.lat[Inf.idx] = 0
    Inf.idx = which(y.pol.e.obs==Inf)
    y.pol.e.obs[Inf.idx] = 0
    Inf.idx = which(y.pol.e.lat==Inf)
    y.pol.e.lat[Inf.idx] = 0

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
      #       FUN ="*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.e.obs * y.log.e.obs + k.e * z.e.lat * y.log.e.lat
    term.zkz.ypol = XAPlusYZB(z.e.obs, y.pol.e.obs, z.e.lat, k.e, y.pol.e.lat)
      # sweep(
      # sweep(matrix(z.e.obs), 1, matrix(y.pol.e.obs), FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.e.lat), 1, matrix(y.pol.e.lat), FUN = "*", check.margin = FALSE),
      #       FUN ="*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.e.obs * y.pol.e.obs + k.e * z.e.lat * y.pol.e.lat

    shape1.k.new = shape1.logpow(z.e.obs, z.e.lat, k.e, y.pol.e.obs, y.pol.e.lat, penalty, hyper.k.1, hyper.k.2)

    result = (log(shape1.k.new) + log(shape2.c.new) - shape2.c.new*log(scale.lambda.new)) * sum(term.zkz) +
      shape2.c.new*sum(term.zkz.ylog) - (shape1.k.new+1)*sum(term.zkz.ypol)
      # (log(shape1.k.new) + log(shape2.c.new) - shape2.c.new*log(scale.lambda.new)) * apply(term.zkz, 2, FUN = "sum") +
      # shape2.c.new*apply(term.zkz.ylog, 2, FUN = "sum") - (shape1.k.new+1)*apply(term.zkz.ypol, 2, FUN = "sum")

    if(penalty==TRUE)
    {
      result = result + (hyper.k.1-1)*log(shape1.k.new) - shape1.k.new/hyper.k.2
      + (hyper.c.1-1)*log(shape2.c.new) - shape2.c.new/hyper.c.2
      + (hyper.lambda.1-1)*log(scale.lambda.new) - scale.lambda.new/hyper.lambda.2
    }

    return(result * (-1)) # optim is minimization
  }

  # I should only use "positive" observations
  pos.idx = (yu!=0)
  shape1.k.old = shape1.k
  shape2.c.old = shape2.c
  scale.lambda.old = scale.lambda
  params.init = log(c(shape2.c.old, scale.lambda.old)) # c(shape2.c.old, scale.lambda.old)
  temp.update = optim(par = params.init, fn = Q.T,
                      shape1.k.old = shape1.k.old, shape2.c.old = shape2.c.old, scale.lambda.old = scale.lambda.old,
                      z.e.obs = z.e.obs[pos.idx], z.e.lat = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                      y.log.e.obs = y.log.e.obs[pos.idx], y.log.e.lat = y.log.e.lat[pos.idx],
                      tl = tl[pos.idx], yl = yl[pos.idx], yu = yu[pos.idx], tu = tu[pos.idx],
                      expert.ll = expert.ll[pos.idx], expert.tn = expert.tn[pos.idx], expert.tn.bar = expert.tn.bar[pos.idx],
                      penalty = penalty,
                      hyper.k.1 = hyper.k.1, hyper.k.2 = hyper.k.2,
                      hyper.c.1 = hyper.c.1, hyper.c.2 = hyper.c.2,
                      hyper.lambda.1 = hyper.lambda.1, hyper.lambda.2 = hyper.lambda.2,
                      lower = pmax(0.5*params.init, 0), upper = pmax(2*params.init, 0),
                      # method = "L-BFGS-B")$par
                      method = "Nelder-Mead")$par

  # Update of c and lambda
  # shape2.c.new = temp.update[1]
  # scale.lambda.new = temp.update[2]
  shape2.c.new = exp(temp.update[1])
  scale.lambda.new = exp(temp.update[2])

  # Update of k
  # E-Step: conditional expectation of y.pol
  y.pol.e.obs = array(0, dim = c(sample.size.n ,1))
  y.pol.e.lat = array(0, dim = c(sample.size.n, 1))

  # Untruncated and uncensored case
  y.pol.e.obs[!censor.idx] = log1p((yu[!censor.idx]/scale.lambda.new)^shape2.c.new)
    # log(1+(yu[!censor.idx]/scale.lambda.new)^shape2.c.new)

  # First find unique upper and lower bounds of integration: untruncated case
  y.unique = unique(cbind(yl,yu),MARGIN=1)
  y.unique.length = nrow(y.unique)
  y.unique.match = match(data.frame(t(cbind(yl,yu))),data.frame(t(y.unique)))
  yl.unique = y.unique[,1]
  yu.unique = y.unique[,2]

  # Integration
  y.pol.e.obs.unique = array(0,dim=c(y.unique.length,1))

  y.pol.e.obs.unique[,1]= # intBurrPolYObs(shape1.k.old, shape2.c.old, scale.lambda.old,
    #                                      shape2.c.new, scale.lambda.new, log(yl.unique), log(yu.unique))
    mapply(function(x, y) ifelse(x!=y,
                                 integrate(int.y.pol.fcn, log(x), log(y),
                                           shape2.c.new = shape2.c.new, scale.lambda.new = scale.lambda.new,
                                           shape1.k.j = shape1.k.old, shape2.c.j = shape2.c.old, scale.lambda.j = scale.lambda.old,
                                           rel.tol=.Machine$double.eps^0.5)$value,
                                 0),
           yl.unique, yu.unique)
  # Match to all observations of y
  temp.y.pol.e.obs = array(0, dim = c(sample.size.n, 1))
  temp.y.pol.e.obs = y.pol.e.obs.unique[y.unique.match,]

  y.pol.e.obs[censor.idx, 1] = exp(-expert.ll[censor.idx]) * temp.y.pol.e.obs[censor.idx]

  # Similarly, find unique upper and lower bounds of integration: truncated case
  tn.unique = unique(cbind(tl,tu),MARGIN=1)
  tn.unique.length = nrow(tn.unique)
  tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
  tl.unique = tn.unique[,1]
  tu.unique = tn.unique[,2]

  # Integration
  y.pol.e.lat.unique = array(0,dim=c(tn.unique.length, 1))

  y.pol.e.lat.unique[,1]=# intBurrPolYLat(shape1.k.old, shape2.c.old, scale.lambda.old,
                        #                 shape2.c.new, scale.lambda.new, log(tl.unique), log(tu.unique))
    sapply(tl.unique,
           function(x) ifelse(x!=0,
                              integrate(int.y.pol.fcn, -Inf, log(x),
                                        shape2.c.new = shape2.c.new, scale.lambda.new = scale.lambda.new,
                                        shape1.k.j = shape1.k.old, shape2.c.j = shape2.c.old, scale.lambda.j = scale.lambda.old,
                                        rel.tol=.Machine$double.eps^0.5)$value,
                              0))+
    sapply(tu.unique,
           function(x) ifelse(x!=Inf,
                              integrate(int.y.pol.fcn, log(x), Inf, # actuar::qburr(1-1e-09, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = T, log.p = F), # Inf causes a problem
                                        shape2.c.new = shape2.c.new, scale.lambda.new = scale.lambda.new,
                                        shape1.k.j = shape1.k.old, shape2.c.j = shape2.c.old, scale.lambda.j = scale.lambda.old,
                                        rel.tol=.Machine$double.eps^0.5)$value,
                              0))
  # Match to all observations of y
  temp.y.pol.e.lat = array(0, dim = c(sample.size.n, 1))
  temp.y.pol.e.lat = y.pol.e.lat.unique[tn.unique.match,]
  # Conditional expectation of log(y): truncated
  y.pol.e.lat = exp(-expert.tn.bar) * temp.y.pol.e.lat

  # Finally, get rid of NaN's
  y.pol.e.obs[is.nan(y.pol.e.obs)] = 0
  y.pol.e.lat[is.nan(y.pol.e.lat)] = 0

  shape1.k.new = shape1.logpow(z.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx], y.pol.e.obs[pos.idx], y.pol.e.lat[pos.idx],
                               penalty, hyper.k.1, hyper.k.2)

  # Update and return
  params.new[1] = shape1.k.new
  params.new[2] = shape2.c.new
  params.new[3] = scale.lambda.new

  return(params.new)

}


#' ## ECM algorithm of Inverse Gaussian expert
#' #' ECM: M-Step for Inverse Gaussian expert.
#' #'
#' #' @importFrom stats optim integrate
#' #' @importFrom actuar dburr qburr pburr
#' #'
#' #' @keywords internal
#' #'
#' #' @export EMMBurrFast
#' EMMBurrFast = function(params.old,
#'                    tl, yl, yu, tu,
#'                    expert.ll, expert.tn, expert.tn.bar,
#'                    z.e.obs, z.e.lat, k.e,
#'                    penalty, hyper.params)
#' {
#'   # Constants
#'   sample.size.n = length(tl)
#'   # Take old parameters
#'   shape1.k = params.old[1]
#'   shape2.c = params.old[2]
#'   scale.lambda = params.old[3]
#'   # Take hyper parameters
#'   hyper.k.1 = hyper.params[1]
#'   hyper.k.2 = hyper.params[2]
#'   hyper.c.1 = hyper.params[3]
#'   hyper.c.2 = hyper.params[4]
#'   hyper.lambda.1 = hyper.params[5]
#'   hyper.lambda.2 = hyper.params[6]
#'   # Value to return
#'   params.new = params.old
#'
#'   # Use brute-force optimization, due to complex form of pdf/cdf of gammacount
#'   Q.params = function(params.new,
#'                       tl, yl, yu, tu,
#'                       z.e.obs, z.e.lat, k.e,
#'                       penalty,
#'                       hyper.k.1, hyper.k.2, hyper.c.1, hyper.c.2, hyper.lambda.1, hyper.lambda.2)
#'   {
#'     shape1.k = params.new[1]
#'     shape2.c = params.new[2]
#'     scale.lambda = params.new[3]
#'
#'     # Initialization: return value are N * 1 matrices
#'     expert.ll = expert.tn = expert.tn.bar = array(-Inf, dim=c(length(yu),1))
#'
#'     # Find indexes of unequal yl & yu: Not exact observations, but censored
#'     censor.idx = (yl!=yu)
#'     prob.log.yu = actuar::pburr(yu[censor.idx], shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)
#'     prob.log.yl = actuar::pburr(yl[censor.idx], shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)
#'
#'     # Compute loglikelihood for expert j, first for y
#'     expert.ll[censor.idx,1] = prob.log.yu + log1mexp(prob.log.yu - prob.log.yl) # likelihood of censored interval: some easy algebra
#'     expert.ll[!censor.idx,1] = actuar::dburr(yu[!censor.idx], shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, log = TRUE)  # exact likelihood
#'
#'     ###################################################################
#'     # Deal with numerical underflow: prob.log.yu and prob.log.yl can both be -Inf
#'     # NA.idx = which(is.na(expert.ll[,j]))
#'     expert.ll[which(is.na(expert.ll[,1])), 1] = -Inf
#'     ###################################################################
#'
#'     # Compute loglikelihood for expert j, then for truncation limits t
#'     prob.log.tu = actuar::pburr(tu, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)
#'     prob.log.tl = actuar::pburr(tl, shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, lower.tail = TRUE, log.p = TRUE)
#'
#'     # Normalizing factor for truncation limits, in log
#'     expert.tn[,1]= prob.log.tu + log1mexp(prob.log.tu - prob.log.tl)
#'
#'     ###################################################################
#'     # Deal with numerical underflow: prob.log.tu and prob.log.tl can both be -Inf
#'     # NA.idx = which(is.na(expert.tn[,j]))
#'     expert.tn[which(is.na(expert.tn[,1])), 1] = -Inf
#'     ###################################################################
#'
#'     ###################################################################
#'     # Deal with no truncation case
#'     no.trunc.idx = (tl==tu)
#'     expert.tn[no.trunc.idx,1] = actuar::dburr(tu[no.trunc.idx], shape1 = shape1.k, shape2 = shape2.c, scale = scale.lambda, log = TRUE)
#'     ###################################################################
#'
#'     ###################################################################
#'     # Deal with exact zero case
#'     zero.idx = (tu==0)
#'     expert.ll[zero.idx,1]=(-Inf)
#'     expert.tn[zero.idx,1]=(-Inf)
#'     ###################################################################
#'
#'     # Log of Pr(outside of truncation interval)
#'     expert.tn.bar[!no.trunc.idx,1] = log1mexp(-expert.tn[!no.trunc.idx,1])
#'     expert.tn.bar[no.trunc.idx,1] = 0
#'
#'     result = sum(XAPlusYZB(z.e.obs, expert.ll[,1], z.e.lat, k.e, expert.tn[,1]))
#'     if(penalty==TRUE){
#'       result = result + (hyper.k.1-1)*log(shape1.k) - shape1.k/hyper.k.2
#'       + (hyper.c.1-1)*log(shape2.c) - shape2.c/hyper.c.2
#'       + (hyper.lambda.1-1)*log(scale.lambda) - scale.lambda/hyper.lambda.2
#'     }
#'
#'     return(result * (-1))
#'   }
#'
#'   pos.idx = (yu!=0)
#'
#'   temp.params = optim(par = params.old, fn = Q.params,
#'                       tl = tl[pos.idx], yl = yl[pos.idx], yu = yu[pos.idx], tu = tu[pos.idx],
#'                       z.e.obs = z.e.obs[pos.idx], z.e.lat = z.e.lat[pos.idx], k.e = k.e[pos.idx],
#'                       penalty = penalty,
#'                       hyper.k.1 = hyper.k.1, hyper.k.2 = hyper.k.2, hyper.c.1 = hyper.c.1, hyper.c.2 = hyper.c.2,
#'                       hyper.lambda.1 = hyper.lambda.1, hyper.lambda.2 = hyper.lambda.2,
#'                       method = "L-BFGS-B", lower = 0.5*params.old, upper = 2*params.old)$par
#'
#'   params.new = temp.params
#'
#'   return(params.new)
#'
#' }
