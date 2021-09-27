#' @title Expert Library Classes
#' @name ExpertMatrix
#'
#' @description 
#' This is the library of all the Expert Functions. In another word, it stored all the metadata of the Expert Functions. 
#' This object is a singleton, which means it should be only created when the package is loaded. 
#' It is also hidden in the global environment. Please do not clear the hidden objects
#' And it you cleared the hidden objects. Please reload the package to fix that.
#' 
#' @export
ExpertLibrary = R6Class("ExpertMatrix", cloneable = TRUE,
  public = list(
    #' @field continuous_experts (`NULL`)\cr
    #' The names of continuous ExpertFunction
    continuous_experts = c(),
    
    #' @field discrete_experts (`NULL`)\cr
    #' The names of discrete ExpertFunction
    discrete_experts = c(),
    
    #' @field distribution_names (`NULL`)\cr
    #' The names of all expert functions
    distribution_names = c(),
    
    #' Register a new Expert Function
    #'
    #' @param expert_name (`character`)\cr
    #' The name of the expert function
    #' @param continuous (`logical`)\cr
    #' Is this new expert function continuous?
    register = function(expert_name, continuous) {
      self$distribution_names = c(self$distribution_names, expert_name)
      if(continuous) { self$continuous_experts = c(self$continuous_experts, expert_name) }
      else{ self$discrete_experts = c(self$discrete_experts, expert_name) }
    },
    
    #' Is given Expert Function continuous?
    #'
    #' @param expert_name (`character`)\cr
    #' The name of the expert function
    is_continous = function(expert_name) {
      return(expert_name %in% self$continuous_experts)
    },
    
    #' Is given Expert Function discrete?
    #'
    #' @param expert_name (`character`)\cr
    #' The name of the expert function
    is_discrete = function(expert_name) {
      return(expert_name %in% self$discrete_experts)
    },
    
    #' Is given Expert Function got a correct name?
    #'
    #' @param expert_name (`character`)\cr
    #' The name of the expert function
    is_valid_expert = function(expert_name) {
      if(expert_name %in% self$distribution_names) {
        return(expert_name)
      }else{
        stop(paste0(expert_name, " was not found in expertlib, check expert name or register in zzz.R"))
      }
    }
  
  ),
  private = list()
)