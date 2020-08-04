## ECM algorithm of Poisson expert
#' ECM: M-Step for Poisson expert.
#'
#' @importFrom stats dpois optimise integrate
#'
#' @keywords internal
#'
#' @export EMMPoisson
EMMPoisson = function(params.old,
                                  tl, yl, yu, tu,
                                  expert.ll, expert.tn, expert.tn.bar,
                                  z.e.obs, z.e.lat, k.e,
                                  penalty, hyper.params)
{
  # Constants
  sample.size.n = length(tl)
  # Take old parameters
  mean.theta = params.old[1]
  # Take hyper parameters
  hyper.mean.1 = hyper.params[1]
  hyper.mean.2 = hyper.params[2]
  # Value to return
  params.new = params.old


  # E-Step: Conditional expectation of y
  censor.idx = (yl!=yu)
  y.e.obs = array(0, dim = c(sample.size.n, 1))
  y.e.lat = array(0, dim = c(sample.size.n, 1))

  # Conditional expectation of y: untruncated and uncensored case
  y.e.obs[!censor.idx, 1] = yl[!censor.idx]

  # Conditional expectation of y: untruncated but cencored case
  int.y.fcn = function(lower, upper, mean.theta.j)
  {
    lower.bound = ceiling(lower)
    upper.bound = floor(upper)

    if (upper!=Inf) {
      y.series = c((lower.bound):(upper.bound))
      dens.series = dpois(y.series, mean.theta.j, log = FALSE)
      result = sum(y.series * dens.series)
    }else{
      if(lower.bound<=1)
      {
        y.series = c(0)
      }else
      {
        y.series = c((0):(lower.bound-1))
      }
      dens.series = dpois(y.series, mean.theta.j, log = FALSE)
      result = mean.theta.j - sum(y.series * dens.series)
    }

    return(sum(result))
  }

  # First find unique upper and lower bounds of integration
  y.unique = unique(cbind(yl,yu),MARGIN=1)
  y.unique.length = nrow(y.unique)
  y.unique.match = match(data.frame(t(cbind(yl,yu))),data.frame(t(y.unique)))
  yl.unique = y.unique[,1]
  yu.unique = y.unique[,2]

  # Integration
  y.e.obs.unique = array(0,dim=c(y.unique.length,1))
  # y.e.obs.unique = matrix(0, nrow = y.unique.length, ncol = 1)

  y.e.obs.unique[,1] = sumPoissonYObs(mean.theta, (yl.unique), (yu.unique))
    # mapply(function(x, y) ifelse(x!=y,
    #                              int.y.fcn(x, y, mean.theta),
    #                              # 0),
    #                              y*dpois(y, mean.theta) ),
    #        yl.unique, yu.unique)
  # y.e.obs.unique[which(yl.unique!=yu.unique)] =
  #   # mapply(function(x, y){int.y.fcn(x, y, mean.theta)}, yl.unique[yl.unique!=yu.unique], yu.unique[yl.unique!=yu.unique], MoreArgs = list(mean.theta = mean.theta))
  #   sapply(X = cbind(yl.unique[yl.unique!=yu.unique], yu.unique[yl.unique!=yu.unique]),
  #          FUN = function(x){int.y.fcn(x[,1], x[,2], mean.theta)},
  #          mean.theta = mean.theta,
  #          simplify = TRUE)
  # y.e.obs.unique[which(yl.unique==yu.unique)] =
  #   yl.unique[which(yl.unique==yu.unique)] * dpois(yl.unique[which(yl.unique==yu.unique)], mean.theta)
  # y.e.obs.unique = matrix(y.e.obs.unique)

  # Match to all observations of y
  temp.y.e.obs = array(0, dim = c(sample.size.n, 1))
  temp.y.e.obs = y.e.obs.unique[y.unique.match,]
  # Conditional expectation of y: untruncated and censored case
  y.e.obs[censor.idx, 1] = exp(-expert.ll[censor.idx]) * temp.y.e.obs[censor.idx]

  # Conditional expectation of y, log(y) and inv(y): truncated
  # Reuse the integration functions above
  tn.unique = unique(cbind(tl,tu),MARGIN=1)
  tn.unique.length = nrow(tn.unique)
  tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
  tl.unique = tn.unique[,1]
  tu.unique = tn.unique[,2]

  # Integration
  y.e.lat.unique = array(0,dim=c(tn.unique.length,1))

  y.e.lat.unique[,1] = mean.theta - sumPoissonYObs(mean.theta, (tl.unique), (tu.unique))
    # sumPoissonYLat(mean.theta, (tl.unique), (tu.unique))
    # mean.theta - sumPoissonYObs(mean.theta, (tl.unique), (tu.unique))
    # mean.theta -
    # mapply(function(x, y) ifelse(x!=y,
    #                              int.y.fcn(x, y, mean.theta),
    #                              # 0),
    #                              y*dpois(y, mean.theta) ),
    #        tl.unique, tu.unique)
  # y.e.lat.unique[(tl.unique==tu.unique),1] = mean.theta - tl.unique[tl.unique==tu.unique]*dpois(tl.unique[tl.unique==tu.unique], mean.theta)
  # y.e.lat.unique[(tl.unique!=tu.unique),1] = mean.theta - mapply(function(x, y){int.y.fcn(x, y, mean.theta)}, tl.unique[tl.unique!=tu.unique], tu.unique[tl.unique!=tu.unique])

  # Match to all observations of y
  temp.y.e.lat = array(0, dim = c(sample.size.n, 1))
  temp.y.e.lat = y.e.lat.unique[tn.unique.match,]

  # Conditional expectation of y: truncated
  y.e.lat[,1] = exp(-expert.tn.bar) * temp.y.e.lat

  # Finally, get rid of NaN's
  y.e.obs[is.nan(y.e.obs)] = 0
  y.e.lat[is.nan(y.e.lat)] = 0


  # M-Step
  # pos.idx = (yu!=0)
  pos.idx = rep(TRUE, length(yu))

  sum.one = XPlusYZ(z.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx])
    # sweep(matrix(z.e.obs[pos.idx]), 1, sweep(matrix(z.e.lat[pos.idx]), 1, matrix(k.e[pos.idx]), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE)
  # z.e.obs[pos.idx] + k.e[pos.idx] * z.e.lat[pos.idx]
  sum.two = XAPlusYZB(z.e.obs[pos.idx], y.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx], y.e.lat[pos.idx])
    # sweep(sweep(matrix(z.e.obs[pos.idx]), 1, matrix(y.e.obs[pos.idx]), FUN = "*", check.margin = FALSE), 1,
    #               sweep(sweep(matrix(z.e.lat[pos.idx]), 1, matrix(k.e[pos.idx]), FUN = "*", check.margin = FALSE), 1, matrix(y.e.lat[pos.idx]), FUN = "*", check.margin = FALSE),
    #               FUN = "+", check.margin = FALSE)
  # z.e.obs[pos.idx] * y.e.obs[pos.idx] + k.e[pos.idx] * z.e.lat[pos.idx] * y.e.lat[pos.idx]
  mean.theta.new = (sum(sum.two) - (hyper.mean.1-1)) / (sum(sum.one) + 1/hyper.mean.2)
    # (apply(sum.two, 2, FUN = "sum") - (hyper.mean.1-1)) / (apply(sum.one, 2, FUN = "sum") + 1/hyper.mean.2)

  params.new[1] = mean.theta.new

  return(params.new)

}
