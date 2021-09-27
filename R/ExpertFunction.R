#' @title Expert Function Classes
#' @name ExpertFunction
#'
#' @description 
#' This is the base class for Expert Functions. Users are expected to interact with `ExpertFunction` to perform all the computations.
#'
#' @export
#' @examples
#' # p = ExpertFunction$new("poisson")
#' # p$continuous
#' # p$distribution
#' # Estimate the start parameters and get them
#' # p$initialize_params(c(1,12,2,3,4,6,7,7))
#' # p$get_params()
#' # Initialize the penalty using default penalty parameters and get them
#' # p$initialize_penalty()
#' # p$get_penalty()
#' # Without setiing parameters, the object will use preset parameters
#' # p$get_penalty_params()
#' # Exposurize function, this will not change the parameters of the ExpertFunction
#' # p$exposurize(2)
#' # Calculate the loglikelihood by two methods
#' # p$ll_exact(5)
#' # p$ll_not_exact(1,1,2,2)
ExpertFunction = R6Class("ExpertFunction", cloneable = TRUE,
    public = list(
      #' @field distribution  (`character`)\cr
      #' The distribution name specified by the user.
      distribution = NULL,
      
      #' @field continuous (`logical`)\cr
      #' Given distribution is continuous or not, default set to be `TRUE`.
      continuous = NULL,
      
      #' @description
      #' Create a new ExpertFunction Class
      #'
      #' @param distribution  (`character`)\cr
      #' The distribution name specified by the user.
      #' @param params (`list`)\cr
      #' The parameters for the given distribution.
      #' Could be preset by users or postset by running `ExpertFunction$params_init(observations)`
      #' @param penalty_params (`list`)\cr
      #' The penalty parameters used to determine the penalty applied on likelihood, will be further used in (`initialize_penalty`) methods
      initialize = function( distribution, params = list(), penalty_params = list()) {
        self$distribution = .expertlib$is_valid_expert(tolower(distribution))
        self$continuous = .expertlib$is_continous(self$distribution)
        self$set_params(params)
        self$set_penalty_params(penalty_params)
      },

      #' Calculate the penalty value of the Expert Distribution
      #' 
      #' @param default (`logical`)\cr
      #' Use default penalty parameters if default is TRUE, else use user defined penalty parameters. Default set to be TRUE.
      #' 
      #' @description 
      #' Prerequisite: self$params and self$penalty_params need to be set. Otherwise will return the default penalty value.
      initialize_penalty = function() {
        penalty = ifelse(!length(private$penalty_params),
                         do.call(paste0(self$distribution, ".penalty"),
                                 list(params = private$params, penalty_params = list())), # change to penalty_params
                         do.call(paste0(self$distribution, ".penalty"),
                                 list(params = private$params, penalty_params = private$penalty_params))
                         )
        private$penalty = penalty
      },
      
      #' @description
      #' Estimate all the parameters for this expert function based on observation
      #'
      #' @param y (`numeric`)\cr
      #' The observations from the preset expert distribution
      initialize_params = function(y) {
        params = do.call( paste0(self$distribution, ".params_init"), list(y=y) )
        private$params = self$set_params(params)
      },
      
      #' Set the parameters for the expert function
      #' 
      #' @param expert_params (`logical`)\cr
      #' The parameters for the given distribution.
      set_params = function(expert_params) {
        private$params = do.call( paste0(self$distribution, ".set_params"), list(params = expert_params) )
      },
      
      #' Get the parameters for the expert function
      #'
      #' @param params (`logical`)\cr
      #' The parameters for the given distribution.
      get_params = function() {
        return( private$params )
      },
      
      #' Set the penalty parameters
      #' If there are no expert_penalty_params given, use default penalty parameters.
      #'
      #' @param expert_penalty_params (`logical`)\cr
      #' The parameters for the given distribution. Could be preset by users or postset by running `ExpertFunction$initialize_params(observations)` 
      set_penalty_params = function(expert_penalty_params) {
        private$penalty_params = expert_penalty_params
        if(length(expert_penalty_params)){ self$initialize_penalty() }
      },
      
      #' @description 
      #' Get the penalty parameters
      get_penalty_params = function() {
        return(private$penalty_params)
      },
      
      #' @description 
      #' Get the penalty value, penalty value cannot be set. It can only been initialized.
      get_penalty = function() {
        return(private$penalty)
      },
      
      #' @description 
      #' Get mean of this expert function
      get_mean = function() {
        do.call( paste0(self$distribution, ".mean"), list(params = private$params) )
      },
      
      #' @description 
      #' Get variance of this expert function
      get_variance = function() {
        do.call( paste0(self$distribution, ".variance"), list(params = private$params) )
      },
      
      #' @description 
      #' Get the cdf of this expert function given RV
      #' @param q The value of random variable
      get_cdf = function(q) {
        do.call( paste0(self$distribution, ".cdf"), list(params = private$params, q = q) )
      },
      
      #' @description 
      #' Get the logcdf of this expert function given RV
      #' @param q The value of random variable
      get_logcdf = function(q) {
        do.call( paste0(self$distribution, ".logcdf"), list(params = private$params, q = q) )
      },
      
      #' @description 
      #' Get the pdf of this expert function given RV
      #' @param x The value of random variable
      get_pdf = function(x) {
        do.call( paste0(self$distribution, ".pdf"), list(params = private$params, x = x) )
      },
      
      #' @description 
      #' Get the logpdf of this expert function given RV
      #' @param x The value of random variable
      get_logpdf = function(x) {
        do.call( paste0(self$distribution, ".logpdf"), list(params = private$params, x = x) )
      },
      
      #' @description 
      #' Get the quantile of this expert function given probability
      #' @param p The probability
      get_quantile = function(p) {
        do.call( paste0(self$distribution, ".quantile"), list(params = private$params, p = p) )
      },
      
      #' @description 
      #' Get the Limited Expected Value of the expert function.
      #' @param u (`numeric`)
      get_lev = function(u) {
        if(self$continuous) {
          result = do.call( paste0(self$distribution, ".lev"), list(params = private$params, u = u) )
        } else {
          warning("Operation Warning: Limited Expected Value for discrete function is meaningless")
          result = NA
        }
        return(result)
      },
      
      #' @description 
      #' Get the Excess of the function.
      #' @param u (`numeric`)
      get_excess = function(u) {
        return(self$get_mean() - self$get_lev(u))
      },
      
      #' @description
      #' Exposurize the parameters, return a copy of the object with the updated parameters
      #' Note this function will not change the parameters of the original object
      #' 
      #' @param exposure (`numeric`)\cr
      #' The exposure applied to the parameters, default set to be 1
      exposurize = function(exposure = 1) {
        exposure = do.call( paste0(self$distribution, ".exposurize"), list(params=private$params, exposure = exposure) )
        new_expert = self$clone()
        new_expert$set_params(exposure)
        return( new_expert )
      },
      
      #' @description
      #' Simulate value based on distribution parameters
      #' 
      #' @param n (`numeric`)\cr
      #' The number of simulated values you want
      simulate = function(n = 1) {
        simulations = do.call( paste0(self$distribution, ".simulation"), list(params=private$params, n = n) )
        return( simulations )
      },
      
      #' Calculate the exact log likelihood of the Expert Distribution
      #'
      #' @param y (`numeric`)\cr
      #' The observations from the preset expert distribution
      ll_exact = function(y) {
        ll_value = do.call( paste0(self$distribution, ".expert_ll_exact"), list(y = y, params = private$params) )
        return(ll_value)
      },
      
      
      #' Calculate the non-exact log likelihood of the Expert Distribution
      #'
      #' @param tl (`numeric`)\cr A vector of length N: lower bounds of truncation.
      #' @param yl (`numeric`)\cr A vector of length N: lower bounds of censoring.
      #' @param yu (`numeric`)\cr A vector of length N: upper bounds of censoring.
      #' @param tu (`numeric`)\cr A vector of length N: upper bounds of truncation.
      ll_not_exact = function(tl, yl, yu, tu) {
        result_set = do.call( paste0(self$distribution, ".expert_ll_not_exact"), 
                              list(tl = tl, yl = yl, yu = yu, tu = tu, params = private$params) )
        return(result_set)
      },
      
      #' @description 
      #' Perform the EM optimization with non-exact observations
      #' 
      #' @param expert_old Old expert function
      #' @param tl (`numeric`)\cr A vector of length N: lower bounds of truncation.
      #' @param yl (`numeric`)\cr A vector of length N: lower bounds of censoring.
      #' @param yu (`numeric`)\cr A vector of length N: upper bounds of censoring.
      #' @param tu (`numeric`)\cr A vector of length N: upper bounds of truncation.
      #' @param exposure A vector of length N: exposures for observations
      #' @param z_e_obs Calculated from E-step
      #' @param z_e_lat Calculated from E-step
      #' @param k_e Calculated from E-step
      #' @param penalty T/F: whether penalty is imposed
      #' @param pen_params A vector of penalty parameters
      #' 
      EM_notexact = function(expert_old, tl, yl, yu, tu, exposure,
                             z_e_obs, z_e_lat, k_e,
                             penalty, pen_params) {
        result = do.call( paste0(self$distribution, ".EM_notexact"), 
                          list(expert_old = expert_old, 
                               tl = tl, yl = yl, yu = yu, tu = tu,
                               exposure = exposure, 
                               z_e_obs = z_e_obs, z_e_lat = z_e_lat, k_e = k_e,
                               penalty = penalty, pen_params = pen_params) )
        return(result)
      },
      
      #' @description 
      #' Perform the EM optimization with exact observations
      #' 
      #' @param expert_old Old expert function
      #' @param ye A vector of length N: exact observations
      #' @param exposure A vector of length N: exposures for observations
      #' @param z_e_obs Calculated from E-step
      #' @param penalty T/F: whether penalty is imposed
      #' @param pen_params A vector of penalty parameters
      #' 
      EM_exact = function(expert_old, ye, exposure, z_e_obs, penalty, pen_params) {
        result = do.call( paste0(self$distribution, ".EM_exact"), 
                          list(expert_old = expert_old, ye = ye, exposure = exposure, z_e_obs = z_e_obs, penalty = penalty, pen_params = pen_params) )
        return(result)
      }
    ),
                          
    private = list(
      params = NULL,
      penalty_params = NULL,
      penalty = 0
    )
)
