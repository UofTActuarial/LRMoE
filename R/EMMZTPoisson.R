## ECM algorithm of Zero-Truncated Poisson expert
#' ECM: M-Step for Zero-Truncated Poisson expert.
#'
#' @importFrom stats dpois optimise integrate
#' @importFrom countreg pztpois dztpois
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
#' @export EMMZTPoisson
EMMZTPoisson = function(params.old,
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
  # int.y.fcn = function(lower, upper, mean.theta.j)
  # {
  #   lower.bound = ceiling(lower)
  #   upper.bound = floor(upper)
  #
  #   if (upper!=Inf) {
  #     y.series = c((lower.bound):(upper.bound))
  #     dens.series = dpois(y.series, mean.theta.j, log = FALSE)
  #     result = sum(y.series * dens.series)
  #   }else{
  #     if(lower.bound<=1)
  #     {
  #       y.series = c(0)
  #     }else
  #     {
  #       y.series = c((0):(lower.bound-1))
  #     }
  #     dens.series = dpois(y.series, mean.theta.j, log = FALSE)
  #     result = mean.theta.j - sum(y.series * dens.series)
  #   }
  #
  #   return(sum(result))
  # }

  # First find unique upper and lower bounds of integration
  y.unique = unique(cbind(yl,yu),MARGIN=1)
  y.unique.length = nrow(y.unique)
  y.unique.match = match(data.frame(t(cbind(yl,yu))),data.frame(t(y.unique)))
  yl.unique = y.unique[,1]
  yu.unique = y.unique[,2]

  # Integration
  y.e.obs.unique = array(0,dim=c(y.unique.length,1))
  # y.e.obs.unique = matrix(0, nrow = y.unique.length, ncol = 1)

  y.e.obs.unique[,1] = sumPoissonYObs(mean.theta, (yl.unique), (yu.unique)) / (1-exp(-mean.theta))

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

  y.e.lat.unique[,1] = ( mean.theta - sumPoissonYObs(mean.theta, (tl.unique), (tu.unique)) ) / (1-exp(-mean.theta))

  # Match to all observations of y
  temp.y.e.lat = array(0, dim = c(sample.size.n, 1))
  temp.y.e.lat = y.e.lat.unique[tn.unique.match,]

  # Conditional expectation of y: truncated
  y.e.lat[,1] = exp(-expert.tn.bar) * temp.y.e.lat

  # Finally, get rid of NaN's
  y.e.obs[is.nan(y.e.obs)] = 0
  y.e.lat[is.nan(y.e.lat)] = 0

  # A function to numerically optimize
  Q.T = function(params.new,
                 z.e.obs, z.e.lat, k.e,
                 y.e.obs, y.e.lat,
                 penalty,
                 hyper.mean.1, hyper.mean.2)
  {

    sum.one = XPlusYZ(z.e.obs, z.e.lat, k.e)
    sum.two = XAPlusYZB(z.e.obs, y.e.obs, z.e.lat, k.e, y.e.lat)

    result = -sum(sum.one)*params.new + sum(sum.two)*log(params.new) - sum(sum.one)*log1mexp(params.new)

    if(penalty==TRUE)
    {
      result = result + (hyper.mean.1-1)*log(params.new) - params.new/hyper.mean.2
    }

    return(result * (-1)) # optim is minimization
  }


  # M-Step
  # pos.idx = (yu!=0)
  pos.idx = rep(TRUE, length(yu))
  temp.update = optim(par = mean.theta, fn = Q.T,
                      z.e.obs = z.e.obs[pos.idx], z.e.lat = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                      y.e.obs = y.e.obs[pos.idx], y.e.lat = y.e.lat[pos.idx],
                      penalty = penalty,
                      hyper.mean.1 = hyper.mean.1, hyper.mean.2 = hyper.mean.2,
                      lower = 0.5*mean.theta, upper = 5*mean.theta,
                      method = "L-BFGS-B")$par

  params.new[1] = temp.update

  return(params.new)

}
