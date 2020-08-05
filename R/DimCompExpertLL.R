## Computation of loglikelihood, by dimension and by component
#' Computes loglikelihood with zero inflation, by dimension and by component.
#'
#' @param Y An \code{N*4d} matrix of numerics, where N is the sample size and d is the dimension of observation.
#'          For each dimension, the four columns shoud be organized as \code{(tl, yl, yu, tu)}.
#' @param comp.dist A \code{d*g} matrix of strings, where d is the dimension of observation and g is the number of components.
#' @param zero.init A \code{d*g} matrix of numerics between 0 and 1, representing zero inflation by dimension and by component.
#' @param params.init A list of length d. Each element in this list is again a sublist of length g.
#'                    Each sublist contains one vector of numerics, representing the parameters of experts by dimension and by component.
#'
#' @return \code{dim.comp.dist}: Same as input \code{comp.dist}, for passing to other functions calling this function.
#' @return \code{zero.inflation}: A d*g matrix of TRUE/FALSE, where TRUE indicates the expert is zero-inflated.
#' @return \code{zero.prob}: Same as input \code{zero.init}.
#' @return \code{pos.params}: Same as input \code{params.init}.
#' @return \code{pos.expert.ll, pos.expert.tn, pos.expert.tn.bar, expert.ll, expert.tn, expert.tn.bar} See \code{\link{PosExpertLL}}.
#'         Here, \code{pos.} indicates the values are for positive parts only, while other values are for zero-inflated versions.
#'
#'
#' @keywords internal
#'
#' @export DimCompExpertLL
DimCompExpertLL = function(Y, comp.dist, zero.init, params.init)
{
  sample.size.n = nrow(Y) # sample size
  n.comp = ncol(comp.dist) # no. of experts
  dim.m = nrow(comp.dist) # no. of dimensions of y
  expert.list = rep(
    list(list(dim.comp.dist = rep("", n.comp),
              zero.inflation = rep(FALSE, n.comp),
              zero.prob = rep(0, n.comp),
              pos.params = list(),
              pos.expert.ll = array(-Inf, dim = c(sample.size.n, n.comp)),
              pos.expert.tn = array(-Inf, dim = c(sample.size.n, n.comp)),
              pos.expert.tn.bar = array(-Inf, dim = c(sample.size.n, n.comp)),
              expert.ll = array(-Inf, dim = c(sample.size.n, n.comp)),
              expert.tn = array(-Inf, dim = c(sample.size.n, n.comp)),
              expert.tn.bar = array(-Inf, dim = c(sample.size.n, n.comp))
    )
    ),
    dim.m) # An object to return. The list is by the dimension of observation y.

  for(k in 1:dim.m) # Loop through dim(y), i.e. the rows
  {
    expert.list[[k]]$dim.comp.dist = comp.dist[k,]
    expert.list[[k]]$pos.params = params.init[[k]]
    for(j in 1:n.comp)
    {
      # Identify zero-inflation components, and initialize zero.prob
      # The initialization is no longer necessary, as it is done in/before the main fitting function
      if(substr(expert.list[[k]]$dim.comp.dist[j], 1, 3)=="ZI-")
      {
        expert.list[[k]]$zero.inflation[j] = TRUE
        if (is.null(zero.init))
        {
          expert.list[[k]]$zero.prob[j] = 0.5
        }else
        {
          expert.list[[k]]$zero.prob[j] = zero.init[k,j]
        }

      }else
      {
        expert.list[[k]]$zero.inflation[j] = FALSE
        expert.list[[k]]$zero.prob[j] = 0
      }

      # Compute the log-likelihood of non-zero part
      # tl = Y[,(4*(k-1)+1)]
      # yl = Y[,(4*(k-1)+2)]
      # yu = Y[,(4*(k-1)+3)]
      # tu = Y[,(4*(k-1)+4)]
      temp.dim.pos.expert.list = PosExpertLL(comp.dist[k,j],
                                             tl = Y[,(4*(k-1)+1)], yl = Y[,(4*(k-1)+2)], yu = Y[,(4*(k-1)+3)], tu = Y[,(4*(k-1)+4)],
                                             params = params.init[[k]][[j]])
      expert.list[[k]]$pos.expert.ll[,j] = temp.dim.pos.expert.list[[1]]
      expert.list[[k]]$pos.expert.tn[,j] = temp.dim.pos.expert.list[[2]]
      expert.list[[k]]$pos.expert.tn.bar[,j] = temp.dim.pos.expert.list[[3]]

      # Compute the log-likelihood of zero and non-zero parts combined
      # The computation is very delicate: I tried multiple times. There are subtle difference between severity and frequency.
      # Frequence is easier, while severity requires more care.
      obs.y.zero.idx = (Y[,(4*(k-1)+2)]==0) # (yl==0) # yl=0: possible for observation to be from zero component
      obs.tn.zero.idx = (Y[,(4*(k-1)+1)]==0)# (tl==0) # tl=0: possible for the truncation interval to include zero
      obs.exact.zero.idx = (Y[,(4*(k-1)+4)]==0) # (tu==0) # tu=0: the ovservation is exactly zero, subset of preceding two cases

      # If modifying this code in future, make sure to change this freq.list, if new frequency distributions are to be added!
      freq.dist = c("poisson", "ZI-poisson", "nbinom", "ZI-nbinom", "binom", "ZI-binom", "gammacount", "ZI-gammacount", "ztpoisson", "ZI-ztpoisson")
      # ll.ind: only possibility is the non-zero component
      expert.list[[k]]$expert.ll[!obs.y.zero.idx,j] = log(0 + (1-expert.list[[k]]$zero.prob[j])*exp(expert.list[[k]]$pos.expert.ll[!obs.y.zero.idx,j]))
      # ll.ind: can come from zero component. Need to separate severity and frequency
      expert.list[[k]]$expert.ll[obs.y.zero.idx,j] = log(expert.list[[k]]$zero.prob[j] + (1-expert.list[[k]]$zero.prob[j])*exp(expert.list[[k]]$pos.expert.ll[obs.y.zero.idx,j]))
      if(expert.list[[k]]$dim.comp.dist[j] %in% freq.dist)
      {
        expert.list[[k]]$expert.ll[obs.exact.zero.idx,j] = log(expert.list[[k]]$zero.prob[j] + (1-expert.list[[k]]$zero.prob[j])*exp(expert.list[[k]]$pos.expert.ll[obs.exact.zero.idx,j]))
      }else
      {
        expert.list[[k]]$expert.ll[obs.exact.zero.idx,j] = log(expert.list[[k]]$zero.prob[j] + 0 )
      }

      # ll.ind.tn: only possibility is non-zero component
      expert.list[[k]]$expert.tn[!obs.tn.zero.idx,j] = log(0 + (1-expert.list[[k]]$zero.prob[j])*exp(expert.list[[k]]$pos.expert.tn[!obs.tn.zero.idx,j]))
      # ll.ind.tn: can come from zero component. Need to separate severity and frequency
      expert.list[[k]]$expert.tn[obs.tn.zero.idx,j] = log(expert.list[[k]]$zero.prob[j] + (1-expert.list[[k]]$zero.prob[j])*exp(expert.list[[k]]$pos.expert.tn[obs.tn.zero.idx,j]))
      if(expert.list[[k]]$dim.comp.dist[j] %in% freq.dist)
      {
        expert.list[[k]]$expert.tn[obs.exact.zero.idx,j] = log(expert.list[[k]]$zero.prob[j] + (1-expert.list[[k]]$zero.prob[j])*exp(expert.list[[k]]$pos.expert.tn[obs.exact.zero.idx,j]))
      }else
      {
        expert.list[[k]]$expert.tn[obs.exact.zero.idx,j] = log(expert.list[[k]]$zero.prob[j] + 0 )
      }

      # ll.ind.tn.bar: (for latent) only possibility is non-zero component
      expert.list[[k]]$expert.tn.bar[obs.tn.zero.idx,j] = log(0 + (1-expert.list[[k]]$zero.prob[j])*exp(expert.list[[k]]$pos.expert.tn.bar[obs.tn.zero.idx,j]))
      # ll.ind.tn.bar: (for latent) can come from zero component. Need to separate severity and frequency
      expert.list[[k]]$expert.tn.bar[!obs.tn.zero.idx,j] = log(expert.list[[k]]$zero.prob[j] + (1-expert.list[[k]]$zero.prob[j])*exp(expert.list[[k]]$pos.expert.tn.bar[!obs.tn.zero.idx,j]))

    }
  }

  return(expert.list)
}
