## ECM algorithm of Binomial expert
#' ECM: M-Step for Binomial expert.
#'
#' @importFrom stats dbinom optimise integrate
#'
#' @keywords internal
#'
#' @export EMMBinomial
EMMBinomial = function(params.old,
                                tl, yl, yu, tu,
                                expert.ll, expert.tn, expert.tn.bar,
                                z.e.obs, z.e.lat, k.e,
                                penalty, hyper.params)
{
  # Constants
  sample.size.n = length(tl)
  # Take old parameters
  size.n = params.old[1]
  prob.p = params.old[2]
  # Take hyper parameters
  # hyper.mean.1 = hyper.params[1]
  # hyper.mean.2 = hyper.params[2]
  # Value to return
  params.new = params.old


  # E-Step: Conditional expectation of y
  censor.idx = (yl!=yu)
  y.e.obs = array(0, dim = c(sample.size.n, 1))
  y.e.lat = array(0, dim = c(sample.size.n, 1))

  # Conditional expectation of y: untruncated and uncensored case
  y.e.obs[!censor.idx, 1] = yl[!censor.idx]

  # Conditional expectation of y: untruncated but cencored case
  int.y.fcn = function(lower, upper, size.n.j, prob.p.j)
  {
    lower.bound = ceiling(lower)
    upper.bound = floor(upper)

    if (upper!=Inf) {
      y.series = c((lower.bound):(upper.bound))
      dens.series = dbinom(y.series, size.n.j, prob.p.j, log = FALSE)
      result = sum(y.series * dens.series)
    }else{
      if(lower.bound<=1)
      {
        y.series = c(0)
      }else
      {
        y.series = c((0):(lower.bound-1))
      }
      dens.series = dbinom(y.series, size.n.j, prob.p.j, log = FALSE)
      result = size.n.j*prob.p.j - sum(y.series * dens.series)
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

  y.e.obs.unique[,1]=
    mapply(function(x, y) ifelse(x!=y,
                                 int.y.fcn(x, y, size.n, prob.p),
                                 y*dbinom(y, size.n, prob.p) ),
           yl.unique, yu.unique)

  # Match to all observations of y
  temp.y.e.obs = array(0, dim = c(sample.size.n, 1))
  temp.y.e.obs = y.e.obs.unique[y.unique.match,]
  # Conditional expectation of y: untruncated and censored case
  y.e.obs[censor.idx, 1] = exp(-expert.ll[censor.idx]) * temp.y.e.obs[censor.idx]

  # Conditional expectation of y: truncated
  # Reuse the integration functions above
  tn.unique = unique(cbind(tl,tu),MARGIN=1)
  tn.unique.length = nrow(tn.unique)
  tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
  tl.unique = tn.unique[,1]
  tu.unique = tn.unique[,2]

  # Integration
  y.e.lat.unique = array(0,dim=c(tn.unique.length,1))

  y.e.lat.unique[,1] = size.n*prob.p -
    mapply(function(x, y) ifelse(x!=y,
                                 int.y.fcn(x, y, size.n, prob.p),
                                 y*dbinom(y, size.n, prob.p) ),
           tl.unique, tu.unique)

  # Match to all observations of y
  temp.y.e.lat = array(0, dim = c(sample.size.n, 1))
  temp.y.e.lat = y.e.lat.unique[tn.unique.match,]

  # Conditional expectation of y: truncated
  y.e.lat[,1] = exp(-expert.tn.bar) * temp.y.e.lat

  # Finally, get rid of NaN's
  y.e.obs[is.nan(y.e.obs)] = 0
  y.e.lat[is.nan(y.e.lat)] = 0
  y.e.obs[is.infinite(y.e.obs)] = 0
  y.e.lat[is.infinite(y.e.lat)] = 0


  # M-Step
  # pos.idx = (yu!=0)
  pos.idx = rep(TRUE, length(yu))

  sum.one = XAPlusYZB(z.e.obs[pos.idx], y.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx], y.e.lat[pos.idx])
    # z.e.obs[pos.idx] * y.e.obs[pos.idx] + k.e[pos.idx] * z.e.lat[pos.idx] * y.e.lat[pos.idx]
  sum.two = XAPlusYZB(z.e.obs[pos.idx], (size.n-y.e.obs[pos.idx]), z.e.lat[pos.idx], k.e[pos.idx], (size.n-y.e.lat[pos.idx]))
    # z.e.obs[pos.idx] * (size.n-y.e.obs[pos.idx]) + k.e[pos.idx] * z.e.lat[pos.idx] * (size.n-y.e.lat[pos.idx])
  prob.p.new = (sum(sum.one)) / (sum(sum.one) + sum(sum.two))

  params.new[2] = prob.p.new

  return(params.new)

}
