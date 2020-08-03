## ECM: M-Step for updating parameters of positive part of component experts.
#' ECM: M-Step for updating parameters of positive part of component experts.
#'
#' @param comp.kj.dist A string which indicates component distribution.
#' @param comp.kj.params.old A vector of numerics: old parameter values for \code{comp.kj.dist}.
#' @param tl.k,yl.k,yu.k,tu.k Numerical vectors of length N.
#' @param comp.kj.pos.expert.ll,comp.kj.pos.expert.tn,comp.kj.pos.expert.tn.bar See \code{\link{PosExpertLL}}.
#' @param z.e.obs,z.e.lat,k.e A numerical vector returned by \code{\link{EMEzkz}}.
#' @param penalty,hyper.params.kj See \code{\link{GateExpertLL}}.
#'
#' @return Updated parameter values.
#'
#' @keywords internal
#'
#' @export EMMCompParams
EMMCompParams = function(comp.kj.dist, comp.kj.params.old,
                          tl.k, yl.k, yu.k, tu.k,
                          comp.kj.pos.expert.ll, comp.kj.pos.expert.tn, comp.kj.pos.expert.tn.bar,
                          z.e.obs, z.e.lat, k.e,
                          penalty, hyper.params.kj)
{
  temp = NULL
  switch (comp.kj.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = EMMGamma(params.old = comp.kj.params.old,
                                                 tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                 expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                 z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                 penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-gamma"    = {temp = EMMGamma(params.old = comp.kj.params.old,
                                                 tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                 expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                 z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                 penalty = penalty, hyper.params = hyper.params.kj) },
          "invgauss"    = {temp = EMMInverseGaussian(params.old = comp.kj.params.old,
                                                tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-invgauss" = {temp = EMMInverseGaussian(params.old = comp.kj.params.old,
                                                tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                penalty = penalty, hyper.params = hyper.params.kj) },
          "lnorm"       = {temp = EMMLognormal(params.old = comp.kj.params.old,
                                                     tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                     expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                     z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                     penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-lnorm"    = {temp = EMMLognormal(params.old = comp.kj.params.old,
                                                     tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                     expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                     z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                     penalty = penalty, hyper.params = hyper.params.kj) },
          "weibull"     = {temp = EMMWeibull(params.old = comp.kj.params.old,
                                                     tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                     expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                     z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                     penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-weibull"  = {temp = EMMWeibull(params.old = comp.kj.params.old,
                                                     tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                     expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                     z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                     penalty = penalty, hyper.params = hyper.params.kj) },
          "burr"        = {temp = EMMBurr(params.old = comp.kj.params.old,
                                                      tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                      expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                      z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                      penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-burr"     = {temp = EMMBurr(params.old = comp.kj.params.old,
                                                      tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                      expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                      z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                      penalty = penalty, hyper.params = hyper.params.kj) },
          # Frequency distributions & their zero-inflation
          # "poisson"     = {temp = poisson.params.m.recur(poisson.params.old = comp.kj.params.old,
          #                                                tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
          #                                                expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
          #                                                z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
          #                                                penalty = penalty, hyper.params = hyper.params.kj) },
          # "ZI-poisson"  = {temp = poisson.params.m.recur(poisson.params.old = comp.kj.params.old,
          #                                                tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
          #                                                expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
          #                                                z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
          #                                                penalty = penalty, hyper.params = hyper.params.kj) },
          # "nbinom"      = {temp = nbinom.params.m.recur(nbinom.params.old = comp.kj.params.old,
          #                                               tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
          #                                               expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
          #                                               z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
          #                                               penalty = penalty, hyper.params = hyper.params.kj) },
          # "ZI-nbinom"   = {temp = nbinom.params.m.recur(nbinom.params.old = comp.kj.params.old,
          #                                               tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
          #                                               expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
          #                                               z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
          #                                               penalty = penalty, hyper.params = hyper.params.kj) },
          # "binom"       = {temp = binom.params.m.recur(binom.params.old = comp.kj.params.old,
          #                                              tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
          #                                              expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
          #                                              z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
          #                                              penalty = penalty, hyper.params = hyper.params.kj) },
          # "ZI-binom"    = {temp = binom.params.m.recur(binom.params.old = comp.kj.params.old,
          #                                              tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
          #                                              expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
          #                                              z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
          #                                              penalty = penalty, hyper.params = hyper.params.kj) },
          # "gammacount"  = {temp = gammacount.params.m.recur(gammacount.params.old = comp.kj.params.old,
          #                                                   tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
          #                                                   expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
          #                                                   z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
          #                                                   penalty = penalty, hyper.params = hyper.params.kj) },
          # "ZI-gammacount"  = {temp = gammacount.params.m.recur(gammacount.params.old = comp.kj.params.old,
          #                                                      tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
          #                                                      expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
          #                                                      z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
          #                                                      penalty = penalty, hyper.params = hyper.params.kj) },
          # Error
          stop("Invalid distribution!")
  )
  return(temp)
}
