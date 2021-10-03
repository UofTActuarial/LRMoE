#' @import R6
#' @import stats
#' @importFrom matrixStats rowLogSumExps
#' @importFrom rmutil rgammacount
#' @import actuar
#' @useDynLib LRMoE, .registration = TRUE
#' @import RcppEigen
#' @import RcppNumerical
#' @import checkmate
#' @importFrom Rcpp sourceCpp
#' @importFrom expint gammainc
#' @importFrom EnvStats skewness kurtosis
#' @importFrom fitdistrplus mledist
"_PACKAGE"

.onLoad = function(libname, pkgname) {
  op <- options()

  # Register all the expert functions
  .expertlib <<- ExpertLibrary$new()

  # discrete
  .expertlib$register(expert_name = "binomial", continuous = FALSE)
  .expertlib$register(expert_name = "gammacount", continuous = FALSE)
  .expertlib$register(expert_name = "negativebinomial", continuous = FALSE)
  .expertlib$register(expert_name = "poisson", continuous = FALSE)
  .expertlib$register(expert_name = "zibinomial", continuous = FALSE)
  .expertlib$register(expert_name = "zigammacount", continuous = FALSE)
  .expertlib$register(expert_name = "zinegativebinomial", continuous = FALSE)
  .expertlib$register(expert_name = "zipoisson", continuous = FALSE)

  #continuous
  .expertlib$register(expert_name = "lognormal", continuous = TRUE)
  .expertlib$register(expert_name = "gamma", continuous = TRUE)
  .expertlib$register(expert_name = "burr", continuous = TRUE)
  .expertlib$register(expert_name = "inversegaussian", continuous = TRUE)
  .expertlib$register(expert_name = "weibull", continuous = TRUE)
  .expertlib$register(expert_name = "zilognormal", continuous = TRUE)
  .expertlib$register(expert_name = "zigamma", continuous = TRUE)
  .expertlib$register(expert_name = "ziburr", continuous = TRUE)
  .expertlib$register(expert_name = "ziinversegaussian", continuous = TRUE)
  .expertlib$register(expert_name = "ziweibull", continuous = TRUE)
}


#' Demo data for \code{LRMoE}.
#'
#' \itemize{
#'   \item \code{X}: A matrix of covariates required by the LRMoE fitting function (complete).
#'   \item \code{Y}: A matrix of response required by the LRMoE fitting function (complete).
#'   \item \code{X.obs}: A matrix of covariates required by the LRMoE fitting function (after truncation and censoring).
#'   \item \code{Y.obs}: A matrix of response required by the LRMoE fitting function (after truncation and censoring).
#' }
#'
#' @docType data
#' @keywords data
#' @name LRMoEDemoData
#'
'X'

#' @rdname LRMoEDemoData
'Y'

#' @rdname LRMoEDemoData
'X.obs'

#' @rdname LRMoEDemoData
'Y.obs'
