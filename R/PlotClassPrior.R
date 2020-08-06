#' Plots a stacked bar chart of latent class probabilities, given a covariate vector.
#'
#' @param X A vector of covariates for one policyholder.
#' @param alpha A matrix of logit regression coeffiecients.
#' @param title A text string for plot title.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{PredictClassPrior}}
#'
#' @import ggplot2
#'
#'
#' @export PlotClassPrior
#'
PlotClassPrior = function(X, alpha, title = "Latent Class Probabilities")
{
  weighting = PredictClassPrior(X, alpha, TRUE)

  df = data.frame(covariate = rep("X", nrow(alpha)), class = as.factor(c(1:nrow(alpha))), probability = t(weighting))

  return(
    ggplot(df, aes(fill=class, y=probability, x = covariate)) +
      geom_bar(position="fill", stat="identity") +
      # geom_text(aes(label = round(probability, 2)),
      #           position = position_stack(vjust = 0.5)) +
      xlab("") + ylab("Latent Class Probabilities") +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
