## ECM algorithm of Gamma Count expert
#' ECM: M-Step for Gamma Count expert.
#'
#' @importFrom rmutil qgammacount
#'
#' @keywords internal
#'
#' @export EMMGammaCount
EMMGammaCount = function(params.old,
                           tl, yl, yu, tu,
                           expert.ll, expert.tn, expert.tn.bar,
                           z.e.obs, z.e.lat, k.e,
                           penalty, hyper.params)
{
  # Constants
  sample.size.n = length(tl)
  # Take old parameters
  shape.m = params.old[1]
  disp.s = params.old[2]
  # Take hyper parameters
  hyper.m.1 = hyper.params[1]
  hyper.m.2 = hyper.params[2]
  hyper.s.1 = hyper.params[3]
  hyper.s.2 = hyper.params[4]
  # Value to return
  params.new = params.old

  # Use brute-force optimization, due to complex form of pdf/cdf of gammacount
  # Q.params = function(params.new,
  #                     tl, yl, yu, tu,
  #                     z.e.obs, z.e.lat, k.e,
  #                     penalty, hyper.m.1, hyper.m.2, hyper.s.1, hyper.s.2)
  # {
  #   censor.idx = (yl!=yu)
  #   sample.size.n = length(tl)
  #
  #   expert.gammacount.ll=expert.gammacount.tn=expert.gammacount.tn.bar=array(-Inf, dim=c(sample.size.n,1))
  #
  #   shape.m.new = params.new[1]
  #   disp.s.new = params.new[2]
  #
  #   prob.log.yu = ifelse(yu[censor.idx]==Inf, 0, pgammacount(yu[censor.idx], m = shape.m.new, s = disp.s.new, log.p = TRUE))
  #     # ifelse(yu[censor.idx]==Inf, 0, pgammacount.new(yu[censor.idx], m = shape.m.new, s = disp.s.new, log.p = TRUE))
  #   prob.log.yl = pgammacount(ceiling(yl[censor.idx])-1, m = shape.m.new, s = disp.s.new, log.p = TRUE)
  #     # ifelse(yl[censor.idx]==0, -Inf, pgammacount.new(ceiling(yl[censor.idx])-1, m = shape.m.new, s = disp.s.new, log.p = TRUE))
  #
  #   expert.gammacount.ll[censor.idx,1]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
  #   expert.gammacount.ll[!censor.idx,1]= dgammacount(yl[!censor.idx], m = shape.m.new, s = disp.s.new, log = TRUE)
  #     # dgammacount.new(yl[!censor.idx], m = shape.m.new, s = disp.s.new, log = TRUE)
  #
  #   prob.log.tu = ifelse(tu==Inf, 0, pgammacount(tu, m = shape.m.new, s = disp.s.new, log.p = TRUE))
  #     # ifelse(tu==Inf, 0, pgammacount.new(tu, m = shape.m.new, s = disp.s.new, log.p = TRUE))
  #   prob.log.tl = pgammacount(ceiling(tl)-1, m = shape.m.new, s = disp.s.new, log.p = TRUE)
  #     # ifelse(tl==0, -Inf, pgammacount.new(ceiling(tl)-1, m = shape.m.new, s = disp.s.new, log.p = TRUE))
  #
  #   expert.gammacount.tn[,1]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)
  #   # Deal with no truncation case
  #   no.trunc.idx = (tl==tu)
  #
  #   # expert.gammacount.tn[no.trunc.idx,1] = rmutil::dgammacount(tu[no.trunc.idx], m = shape.m.new, s = disp.s.new, log = TRUE)
  #   expert.gammacount.tn[no.trunc.idx,1] = dgammacount(tu[no.trunc.idx], m = shape.m.new, s = disp.s.new, log = TRUE)
  #     # dgammacount.new(tu[no.trunc.idx], m = shape.m.new, s = disp.s.new, log = TRUE)
  #
  #   expert.gammacount.tn.bar[!no.trunc.idx,1] = log1mexp(-expert.gammacount.tn[!no.trunc.idx,1])
  #   expert.gammacount.tn.bar[no.trunc.idx,1] = log1mexp(-expert.gammacount.tn[no.trunc.idx,1])
  #
  #   result = sum(XAPlusYZB(z.e.obs, expert.gammacount.ll[,1], z.e.lat, k.e, expert.gammacount.tn[,1]))
  #     # sum(z.e.obs*expert.gammacount.ll) + sum(k.e*z.e.lat*expert.gammacount.tn)
  #   if(penalty==TRUE){
  #     result = result + (hyper.m.1-1)*log(shape.m.new) - shape.m.new/hyper.m.2 + (hyper.s.1-1)*log(disp.s.new) - disp.s.new/hyper.s.2
  #   }
  #
  #   return(result * (-1))
  # }

  int.logf.func = function(lower, upper, params.new, params.old)
  {
    shape.m.new = params.new[1]
    disp.s.new = params.new[2]
    shape.m.old = params.old[1]
    disp.s.old = params.old[2]

    lower.bound = max(ceiling(lower), 0)
    upper.bound = min(floor(upper), qgammacount(1-1e-8, shape.m.old, disp.s.old))

    series = c((lower.bound):(upper.bound))
    logf.series = dgammacount(series, shape.m.new, disp.s.new, log = TRUE)
    prob.series = dgammacount(series, shape.m.old, disp.s.old, log = FALSE)

    return(sum(logf.series*prob.series))
  }

  Q.params = function(params.new, params.old,
                      tl, yl, yu, tu,
                      expert.ll, expert.tn, expert.tn.bar,
                      z.e.obs, z.e.lat, k.e,
                      penalty, hyper.m.1, hyper.m.2, hyper.s.1, hyper.s.2)
  {
    censor.idx = (yl!=yu)
    sample.size.n = length(tl)

    shape.m.new = params.new[1]
    disp.s.new = params.new[2]
    shape.m.old = params.old[1]
    disp.s.old = params.old[2]

    # E-Step: conditional expectations for y, log(y), 1/y
    censor.idx = (yl!=yu)
    y.e.obs = array(0, dim = c(sample.size.n, 1))
    y.e.lat = array(0, dim = c(sample.size.n, 1))

    # Conditional expectation of y: untruncated and uncensored case
    y.e.obs[!censor.idx, 1] = dgammacount(yl[!censor.idx], shape.m.new, disp.s.new, log = TRUE)

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
                                   int.logf.func((x), (y), params.new, params.old),
                                   0),
             yl.unique, yu.unique)

    # Match to all observations of y
    temp.y.e.obs = array(0, dim = c(sample.size.n, 1))
    temp.y.e.obs = y.e.obs.unique[y.unique.match,]

    # Conditional expectation of y: untruncated and censored case
    y.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * temp.y.e.obs[censor.idx]

    # Reuse the integration functions above
    tn.unique = unique(cbind(tl,tu),MARGIN=1)
    tn.unique.length = nrow(tn.unique)
    tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
    tl.unique = tn.unique[,1]
    tu.unique = tn.unique[,2]

    # Integration
    y.e.lat.unique = array(0,dim=c(tn.unique.length,1))

    y.e.lat.unique[,1]=
      sapply(tl.unique,
             function(x) ifelse(x!=0,
                                int.logf.func(0, (x), params.new, params.old),
                                0))+
      sapply(tu.unique,
             function(x) ifelse(x!=Inf,
                                int.logf.func((x), (Inf), params.new, params.old),
                                0))

    # Match to all observations of y
    temp.y.e.lat = array(0, dim = c(sample.size.n, 1))
    temp.y.e.lat = y.e.lat.unique[tn.unique.match,]

    # Conditional expectation of y: truncated
    y.e.lat = exp(-expert.tn.bar) * temp.y.e.lat

    # Finally, get rid of NaN's
    y.e.obs[is.na(y.e.obs)] = 0
    y.e.lat[is.na(y.e.lat)] = 0

    result = sum(XAPlusYZB(z.e.obs, y.e.obs, z.e.lat, k.e, y.e.lat))
      # sum(XAPlusYZB(z.e.obs, expert.gammacount.ll[,1], z.e.lat, k.e, expert.gammacount.tn[,1]))
    # sum(z.e.obs*expert.gammacount.ll) + sum(k.e*z.e.lat*expert.gammacount.tn)

    if(penalty==TRUE){
      result = result + (hyper.m.1-1)*log(shape.m.new) - shape.m.new/hyper.m.2 + (hyper.s.1-1)*log(disp.s.new) - disp.s.new/hyper.s.2
    }

    return(result * (-1))
  }

  # pos.idx = (yu!=0)
  pos.idx = rep(TRUE, length(yu))

  temp.params = optim(par = params.old, fn = Q.params,
                      params.old = params.old,
                      tl = tl[pos.idx], yl = yl[pos.idx], yu = yu[pos.idx], tu = tu[pos.idx],
                      expert.ll = expert.ll[pos.idx], expert.tn = expert.tn[pos.idx], expert.tn.bar = expert.tn.bar[pos.idx],
                      z.e.obs = z.e.obs[pos.idx], z.e.lat = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                      penalty = penalty, hyper.m.1 = hyper.m.1, hyper.m.2 = hyper.m.2, hyper.s.1 = hyper.s.1, hyper.s.2 = hyper.s.2,
                      # method = "L-BFGS-B", lower = 0.5*params.old, upper = 2*params.old)$par
                      method = "Nelder-Mead")$par

  params.new = temp.params


  return(params.new)
}
