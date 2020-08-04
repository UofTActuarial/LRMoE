## ECM algorithm of Gamma Count expert
#' ECM: M-Step for Gamma Count expert.
#'
#' @importFrom stats dnbinom
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
  Q.params = function(params.new,
                      tl, yl, yu, tu,
                      z.e.obs, z.e.lat, k.e,
                      penalty, hyper.m.1, hyper.m.2, hyper.s.1, hyper.s.2)
  {
    censor.idx = (yl!=yu)
    sample.size.n = length(tl)

    expert.gammacount.ll=expert.gammacount.tn=expert.gammacount.tn.bar=array(-Inf, dim=c(sample.size.n,1))

    shape.m.new = params.new[1]
    disp.s.new = params.new[2]

    prob.log.yu = ifelse(yu[censor.idx]==Inf, 0, pgammacount(yu[censor.idx], m = shape.m.new, s = disp.s.new, log.p = TRUE))
      # ifelse(yu[censor.idx]==Inf, 0, pgammacount.new(yu[censor.idx], m = shape.m.new, s = disp.s.new, log.p = TRUE))
    prob.log.yl = pgammacount(ceiling(yl[censor.idx])-1, m = shape.m.new, s = disp.s.new, log.p = TRUE)
      # ifelse(yl[censor.idx]==0, -Inf, pgammacount.new(ceiling(yl[censor.idx])-1, m = shape.m.new, s = disp.s.new, log.p = TRUE))

    expert.gammacount.ll[censor.idx,1]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.gammacount.ll[!censor.idx,1]= dgammacount(yl[!censor.idx], m = shape.m.new, s = disp.s.new, log = TRUE)
      # dgammacount.new(yl[!censor.idx], m = shape.m.new, s = disp.s.new, log = TRUE)

    prob.log.tu = ifelse(tu==Inf, 0, pgammacount(tu, m = shape.m.new, s = disp.s.new, log.p = TRUE))
      # ifelse(tu==Inf, 0, pgammacount.new(tu, m = shape.m.new, s = disp.s.new, log.p = TRUE))
    prob.log.tl = pgammacount(ceiling(tl)-1, m = shape.m.new, s = disp.s.new, log.p = TRUE)
      # ifelse(tl==0, -Inf, pgammacount.new(ceiling(tl)-1, m = shape.m.new, s = disp.s.new, log.p = TRUE))

    expert.gammacount.tn[,1]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)

    # expert.gammacount.tn[no.trunc.idx,1] = rmutil::dgammacount(tu[no.trunc.idx], m = shape.m.new, s = disp.s.new, log = TRUE)
    expert.gammacount.tn[no.trunc.idx,1] = dgammacount(tu[no.trunc.idx], m = shape.m.new, s = disp.s.new, log = TRUE)
      # dgammacount.new(tu[no.trunc.idx], m = shape.m.new, s = disp.s.new, log = TRUE)

    expert.gammacount.tn.bar[!no.trunc.idx,1] = log1mexp(-expert.gammacount.tn[!no.trunc.idx,1])
    expert.gammacount.tn.bar[no.trunc.idx,1] = log1mexp(-expert.gammacount.tn[no.trunc.idx,1])

    result = sum(XAPlusYZB(z.e.obs, expert.gammacount.ll[,1], z.e.lat, k.e, expert.gammacount.tn[,1]))
      # sum(z.e.obs*expert.gammacount.ll) + sum(k.e*z.e.lat*expert.gammacount.tn)
    if(penalty==TRUE){
      result = result + (hyper.m.1-1)*log(shape.m.new) - shape.m.new/hyper.m.2 + (hyper.s.1-1)*log(disp.s.new) - disp.s.new/hyper.s.2
    }

    return(result * (-1))
  }

  # pos.idx = (yu!=0)
  pos.idx = rep(TRUE, length(yu))

  temp.params = optim(par = params.old, fn = Q.params,
                      tl = tl[pos.idx], yl = yl[pos.idx], yu = yu[pos.idx], tu = tu[pos.idx],
                      z.e.obs = z.e.obs[pos.idx], z.e.lat = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                      penalty = penalty, hyper.m.1 = hyper.m.1, hyper.m.2 = hyper.m.2, hyper.s.1 = hyper.s.1, hyper.s.2 = hyper.s.2,
                      method = "L-BFGS-B", lower = 0.5*params.old, upper = 2*params.old)$par

  params.new = temp.params


  return(params.new)
}
