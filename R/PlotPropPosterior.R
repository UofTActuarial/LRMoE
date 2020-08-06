#' Plots two stacked bar charts of most likely latent class proportions, which contrasts prior and posterior latent class proportions, given a covariate matrix.
#'
#' @param Y A matrix of observed responses for \code{X}.
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param title A text string for plot title.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{PredictClassPosterior}}
#'
#' @import ggplot2
#' @importFrom stats aggregate
#'
#'
#' @export PlotPropPosterior
PlotPropPosterior = function(Y, X, alpha, comp.dist, zero.prob, params.list, title = "Proportion of Latent Classes")
{
  # Old assignments
  alpha.old = alpha
  assign.class.old = PredictClassPrior(X, alpha.old, FALSE)
  df.old = data.frame(class = assign.class.old, count = rep(1, length(assign.class.old)))
  df.old$class = factor(df.old$class, levels = c(1:nrow(alpha.old)))

  df.old.aggre = aggregate(df.old$count, by = list(class = df.old$class), FUN = "sum")
  df.old.plot = data.frame(cbind(case = rep("Prior", nrow(df.old.aggre)), df.old.aggre) )
  df.old.plot$probability = df.old.plot[,3]/sum(df.old.plot[,3])

  # New assignments
  assign.class.new = PredictClassPosterior(Y, X, alpha, comp.dist, zero.prob, params.list, FALSE) # predict.class(X, alpha.new)
  df.new = data.frame(class = assign.class.new, count = rep(1, length(assign.class.new)))
  df.new$class = factor(df.new$class, levels = c(1:nrow(alpha.old)))

  df.new.aggre = aggregate(df.new$count, by = list(class = df.new$class), FUN = "sum")
  df.new.plot = data.frame(cbind(case = rep("Posterior", nrow(df.new.aggre)), df.new.aggre) )
  df.new.plot$probability = df.new.plot[,3]/sum(df.new.plot[,3])

  # df to plot
  df.plot = data.frame(rbind(df.old.plot, df.new.plot))
  df.plot$case = factor(df.plot$case, levels = c("Prior", "Posterior"))

  return(
    ggplot(df.plot, aes_string(fill='class', y='probability', x = 'case')) +
      geom_bar(position="fill", stat="identity") +
      # geom_text(aes(label = round(probability, 2)),
      #           position = position_stack(vjust = 0.5)) +
      xlab("") + ylab("Latent Class Proportion") +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
