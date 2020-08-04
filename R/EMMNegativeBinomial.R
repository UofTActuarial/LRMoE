## ECM algorithm of Negative Binomial expert
#' ECM: M-Step for Negative Binomial expert.
#'
#' @importFrom stats dnbinom pnbinom optim qnbinom
#' @importFrom NMOF gridSearch
#'
#' @keywords internal
#'
#' @export EMMNegativeBinomial
EMMNegativeBinomial = function(params.old,
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
  hyper.n.1 = hyper.params[1]
  hyper.n.2 = hyper.params[2]
  # Value to return
  params.new = params.old


  # E-Step: conditional expectation of y
  censor.idx = (yl!=yu)
  y.e.obs = array(0, dim = c(sample.size.n, 1))
  y.e.lat = array(0, dim = c(sample.size.n, 1))

  # Conditional expectation of y: untruncated and uncensored case
  y.e.obs[!censor.idx, 1] = yl[!censor.idx]

  # Conditional expectation of y: untruncated but cencored case
  # int.y.fcn = function(lower, upper, size.n.j, prob.p.j)
  # {
  #   lower.bound = ceiling(lower)
  #   upper.bound = floor(upper)
  #
  #   if (upper!=Inf) {
  #     y.series = c((lower.bound):(upper.bound))
  #     dens.series = dnbinom(y.series, size = size.n.j, prob = prob.p.j, log = FALSE)
  #     result = sum(y.series * dens.series)
  #   }else{
  #     if(lower.bound<=1) {
  #       y.series = c(0)
  #     }else{
  #       y.series = c((0):(lower.bound-1))
  #     }
  #     dens.series = dnbinom(y.series, size = size.n.j, prob = prob.p.j, log = FALSE)
  #     result = size.n.j*(1-prob.p.j)/prob.p.j - sum(y.series * dens.series)
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

  y.e.obs.unique[,1]= sumNegativeBinomialYObs(size.n, prob.p, (yl.unique), (yu.unique))
    # mapply(function(x, y) ifelse(x!=y,
    #                              int.y.fcn(x, y, size.n, prob.p),
    #                              0),
    #        yl.unique, yu.unique)

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

  y.e.lat.unique[,1] = size.n*(1-prob.p)/prob.p - sumNegativeBinomialYObs(size.n, prob.p, (tl.unique), (tu.unique))
    # size.n*(1-prob.p)/prob.p -
    # mapply(function(x, y) ifelse(x!=y,
    #                              int.y.fcn(x, y, size.n, prob.p),
    #                              0),
    #        tl.unique, tu.unique)

  # Match to all observations of y
  temp.y.e.lat = array(0, dim = c(sample.size.n, 1))
  temp.y.e.lat = y.e.lat.unique[tn.unique.match,]

  # Conditional expectation of y: truncated
  y.e.lat = exp(-expert.tn.bar) * temp.y.e.lat

  # Finally, get rid of NaN's
  y.e.obs[is.nan(y.e.obs)] = 0
  y.e.lat[is.nan(y.e.lat)] = 0


  # M-Step: Maximization of loglik
  #         It involves maximizing a conditional expectation of lfactorial, which has NO closed-form expression.
  int.y.log.fac.fcn = function(size.n.new, lower, upper, size.n.j, prob.p.j)
  {
    lower.bound = ceiling(lower)
    upper.bound = floor(upper)

    if (upper!=Inf) {
      y.series = c((lower.bound):(upper.bound))
      lfac.series = lfactorial(y.series + size.n.new - 1)
      dens.series = dnbinom(y.series, size = size.n.j, prob = prob.p.j, log = FALSE)
      result = sum(lfac.series * dens.series)
    }else{
      # I cannot find a closed-form expression for lfactorial(y+size.n.new-1).
      # Throw away all tail values with very small probabilities.
      # The lfactorial function grows very slow, so this should be fine.
      upper.bound.finite = qnbinom(1-1e-10, size = size.n.j, prob = prob.p.j, lower.tail = TRUE)
      y.series = c((lower.bound):(upper.bound.finite))
      lfac.series = lfactorial(y.series + size.n.new - 1)
      dens.series = dnbinom(y.series, size.n.j, prob.p.j, log = FALSE)
      result = sum(lfac.series * dens.series)
    }
    return(result)
  }
  # The optimal prob.p is a function of size.n
  prob.size = function(size.n, z.e.obs, z.e.lat, k.e, y.e.obs, y.e.lat)
  {
    term.zkz = XPlusYZ(z.e.obs, z.e.lat, k.e)
      # sweep(matrix(z.e.obs), 1,
      #                sweep(matrix(k.e), 1, matrix(z.e.lat), FUN = "*", check.margin = FALSE),
      #                FUN = "+", check.margin = FALSE)
    # z.e.obs + k.e * z.e.lat
    term.zkz.y = XAPlusYZB(z.e.obs, y.e.obs, z.e.lat, k.e, y.e.lat)
      # sweep(
      # sweep(matrix(z.e.obs), 1, y.e.obs, FUN = "*", check.margin = FALSE), 1,
      # sweep(matrix(k.e), 1,
      #       sweep(matrix(z.e.lat), 1, matrix(y.e.lat),
      #             FUN = "*", check.margin = FALSE),
      #       FUN = "*", check.margin = FALSE),
      # FUN = "+", check.margin = FALSE)
    # z.e.obs * y.e.obs + k.e * z.e.lat * y.e.lat
    temp.p = 1 - sum(term.zkz.y) / (size.n*sum(term.zkz) + sum(term.zkz.y))
      # 1 - apply(term.zkz.y, 2, FUN = "sum") / ( size.n*apply(term.zkz, 2, FUN = "sum") + apply(term.zkz.y, 2, FUN = "sum") )
    return(temp.p)
  }



  # Function to optimize
  Q.size = function(size.n.new, size.n.old, prob.p.old,
                    z.e.obs, z.e.lat, k.e, y.e.obs, y.e.lat,
                    tl, yl, yu, tu, expert.ll, expert.tn, expert.tn.bar,
                    penalty, hyper.n.1, hyper.n.2)
  {
    sample.size.n = length(tl)

    y.lfac.e.obs = array(0, dim = c(sample.size.n, 1))
    y.lfac.e.lat = array(0, dim = c(sample.size.n, 1))

    # Uncensored case
    censor.idx = (yl!=yu)
    y.lfac.e.obs[!censor.idx, 1] = lfactorial(yu[!censor.idx] + size.n.new - 1)

    # First find unique upper and lower bounds of integration
    y.unique = unique(cbind(yl,yu),MARGIN=1)
    y.unique.length = nrow(y.unique)
    y.unique.match = match(data.frame(t(cbind(yl,yu))),data.frame(t(y.unique)))
    yl.unique = y.unique[,1]
    yu.unique = y.unique[,2]

    # Integration
    y.lfac.e.obs.unique = array(0,dim=c(y.unique.length,1))

    y.lfac.e.obs.unique[,1]= sumNegativeBinomialLfacYObs(size.n.old, prob.p.old, size.n.new, (yl.unique), (yu.unique))
      # mapply(function(x, y) ifelse(x!=y,
      #                              int.y.log.fac.fcn(size.n.new = size.n.new, lower = x, upper = y, size.n.j = size.n.old, prob.p.j = prob.p.old),
      #                              0),
      #        yl.unique, yu.unique)

    # Match to all observations of y
    temp.y.lfac.e.obs = array(0, dim = c(sample.size.n, 1))
    temp.y.lfac.e.obs = y.lfac.e.obs.unique[y.unique.match,]

    # Conditional expectation of y: untruncated and censored case
    y.lfac.e.obs[censor.idx] = exp(-expert.ll[censor.idx]) * temp.y.lfac.e.obs[censor.idx]

    # Conditional expectation of y: truncated
    # Reuse the integration functions above
    tn.unique = unique(cbind(tl,tu),MARGIN=1)
    tn.unique.length = nrow(tn.unique)
    tn.unique.match = match(data.frame(t(cbind(tl,tu))),data.frame(t(tn.unique)))
    tl.unique = tn.unique[,1]
    tu.unique = tn.unique[,2]

    # Integration
    y.lfac.e.lat.unique = array(0,dim=c(tn.unique.length,1))

    y.lfac.e.lat.unique[,1]= sumNegativeBinomialLfacYLat(size.n.old, prob.p.old, size.n.new, (tl.unique), (tu.unique))
      # mapply(function(x, y) ifelse(x!=y,
      #                              int.y.log.fac.fcn(size.n.new = size.n.new, lower = x, upper = y, size.n.j = size.n.old, prob.p.j = prob.p.old),
      #                              0),
      #        # rep(0,tn.unique.length), tl.unique) +
      #        rep(0,tn.unique.length), (ceiling(tl.unique)-1)*(ceiling(tl.unique)-1 >=0) ) +
      # mapply(function(x, y) ifelse(x!=y,
      #                              int.y.log.fac.fcn(size.n.new = size.n.new, lower = x, upper = y, size.n.j = size.n.old, prob.p.j = prob.p.old),
      #                              0),
      #        # tu.unique, rep(Inf,tn.unique.length))
      #        floor(tu.unique)+1, rep(Inf,tn.unique.length))

    # Match to all observations of y
    temp.y.lfac.e.lat = array(0, dim = c(sample.size.n, 1))
    temp.y.lfac.e.lat = y.lfac.e.lat.unique[tn.unique.match,]

    # Conditional expectation of y: truncated
    y.lfac.e.lat = exp(-expert.tn.bar) * temp.y.lfac.e.lat

    # Finally, get rid of NaN's
    y.lfac.e.obs[is.nan(y.lfac.e.obs)] = 0
    y.lfac.e.lat[is.nan(y.lfac.e.lat)] = 0


    # Now, construct the function to optimize
    term.zkz = XPlusYZ(z.e.obs, z.e.lat, k.e)
      # sweep(matrix(z.e.obs), 1,
      #                sweep(matrix(k.e), 1, matrix(z.e.lat), FUN = "*", check.margin = FALSE),
      #                FUN = "+", check.margin = FALSE)
    # z.e.obs + k.e * z.e.lat
    term.zkz.y = XAPlusYZB(z.e.obs, y.e.obs, z.e.lat, k.e, y.e.lat)
    #   sweep(
    #   sweep(matrix(z.e.obs), 1, matrix(y.e.obs), FUN = "*", check.margin = FALSE), 1,
    #   sweep(matrix(k.e), 1,
    #         sweep(matrix(z.e.lat), 1, matrix(y.e.lat), FUN = "*", check.margin = FALSE),
    #         FUN = "*", check.margin = FALSE),
    #   FUN = "+", check.margin = FALSE
    # )
    # z.e.obs * y.e.obs + k.e * z.e.lat * y.e.lat
    term.zkz.ylfac = XAPlusYZB(z.e.obs, y.lfac.e.obs, z.e.lat, k.e, y.lfac.e.lat)
    #   sweep(
    #   sweep(matrix(z.e.obs), 1, matrix(y.lfac.e.obs), FUN = "*", check.margin = FALSE), 1,
    #   sweep(matrix(k.e), 1,
    #         sweep(matrix(z.e.lat), 1, matrix(y.lfac.e.lat), FUN = "*", check.margin = FALSE),
    #         FUN = "*", check.margin = FALSE),
    #   FUN = "+", check.margin = FALSE
    # )
    # z.e.obs * y.lfac.e.obs + k.e * z.e.lat * y.lfac.e.lat

    temp.p = prob.size(size.n.new, z.e.obs, z.e.lat, k.e, y.e.obs, y.e.lat)

    result = sum(term.zkz.ylfac) - sum(term.zkz)*lfactorial(size.n.new-1) + sum(term.zkz)*size.n.new*log(temp.p) + sum(term.zkz.y)*log(1-temp.p)
      # apply(term.zkz.ylfac, 2, FUN = "sum") - apply(term.zkz, 2, FUN = "sum")*lfactorial(size.n.new-1) + apply(term.zkz, 2, FUN = "sum")*size.n.new*log(temp.p) + apply(term.zkz.y, 2, FUN = "sum")*log(1-temp.p)
    if(penalty==TRUE)
    {
      result = result + (hyper.n.1 - 1)*log(size.n.new) - size.n.new/hyper.n.2
    }

    return(result*(-1))
  }

  # Optimization:
  # pos.idx = (yu!=0)
  pos.idx = rep(TRUE, length(yu))
  size.n.new = NMOF::gridSearch(Q.size,
                                size.n.old = size.n, prob.p.old = prob.p,
                                z.e.obs = z.e.obs[pos.idx], z.e.lat = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                                y.e.obs = y.e.obs[pos.idx], y.e.lat = y.e.lat[pos.idx],
                                tl = tl[pos.idx], yl = yl[pos.idx], yu = yu[pos.idx], tu = tu[pos.idx],
                                expert.ll = expert.ll[pos.idx], expert.tn = expert.tn[pos.idx], expert.tn.bar = expert.tn.bar[pos.idx],
                                penalty = penalty, hyper.n.1 = hyper.n.1, hyper.n.2 = hyper.n.2,
                                # lower = max(0, size.n-5), upper = size.n+5, npar = 1,
                                levels = list(size.n.new = seq(max(1, size.n-10),size.n+10)),
                                printDetail = FALSE)$minlevels

  prob.p.new = prob.size(size.n.new, z.e.obs[pos.idx], z.e.lat[pos.idx], k.e[pos.idx], y.e.obs[pos.idx], y.e.lat[pos.idx])


  params.new[1] = size.n.new
  params.new[2] = prob.p.new

  return(params.new)

}
