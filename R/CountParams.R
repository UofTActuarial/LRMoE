#' Count the number of logit regression coefficients.
#'
#' @param alpha A g*P matrix of logit regression coefficients.
#'
#' @return The number of logit regression coefficients.
#'
#'
#' @keywords internal
#'
#' @export CountAlpha
CountAlpha = function(alpha)
{
  return( ncol(alpha)*(nrow(alpha)-1) )
}



#' Count the number of zero probability masses.
#'
#' @param comp.dist A d*g matrix of strings describing the component distributions.
#'
#' @return The number of zero probability masses.
#'
#'
#' @keywords internal
#'
#' @export CountZero
CountZero = function(comp.dist)
{
  result = 0
  for(k in 1:(nrow(comp.dist)))
  {
    for(j in 1:(ncol(comp.dist)))
    {
      if(substr(comp.dist[k,j], 1, 3)=="ZI-")
      {
        result = result + 1
      }
    }
  }
  return(result)
}
