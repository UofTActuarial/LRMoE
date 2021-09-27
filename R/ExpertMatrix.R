#' @title Expert Matrix Classes
#' @name ExpertMatrix
#'
#' @description 
#' This is the matrix class of Expert Functions. Users are expected to perform all vector operations through ExpertMatrix to reduce index errors.
#' 
#' @export
#' @examples
#' params = matrix(list( list(meanlog = 1, sdlog = 2),  list(meanlog = 1, sdlog  = 2),
#'                       list(meanlog = 3, sdlog = 1),  list(meanlog = 5, sdlog = 1)), nrow = 2)
#' penalty_params = matrix( list( c(2,1,2), c(1,1,1), c(3,1,2), c(1,3,4)), nrow = 2)
#' expert_names = matrix( c("lognormal", "lognormal", "lognormal", "lognormal"), nrow = 2)
#' # ExpertMatrixObj = ExpertMatrix$new(expert_matrix = expert_names, expert_params_matrix = params)
#' # ExpertMatrixObj$set_penalty_params(penalty_params)
#' # ExpertMatrixObj$set_params(params)
ExpertMatrix = R6Class("ExpertMatrix", cloneable = TRUE,
  public = list(
    #' @field nrow  (`numeric`)\cr
    #' The number of rows of ExpertMatrix
    nrow = NULL,
    
    #' @field ncol  (`numeric`)\cr
    #' The number of columns of ExpertMatrix
    ncol = NULL,
    
    #' @field expert_matrix  (`matrix`)\cr
    #' The matrix object to store all the ExpertFunction R6 objects
    expert_matrix = NULL,
    
    #' @field penalty_params  (`matrix`)\cr
    #' The matrix object to store all the penalty parameters for the Expert Functions
    penalty_params = NULL,
    
    #' @field expert_params  (`matrix`)\cr
    #' The matrix object to store all the parameters for the Expert Functions
    expert_params = NULL,
    
    #' Initialze a Expert Matrix Class
    #'
    #' @param expert_matrix (`matrix`)\cr
    #' The name of Expert Functions, stored as a matrix of characters.
    #' @param expert_params_matrix  (`matrix`)\cr
    #' The parameters of Expert Functions, stored as a matrix of lists.
    initialize = function(expert_matrix, expert_params_matrix) {
      #Get the columns and rows of the expert_matrix
      self$ncol = ncol(expert_matrix)
      self$nrow = nrow(expert_matrix)
      
      # Validate the expert matrix.
      checkMatrix(expert_matrix, mode = "character")
      checkMatrix(expert_params_matrix, mode = "list")
      
      # Store all the related parameters and ExpertFunction obj to a matrix
      self$expert_params = expert_params_matrix
      self$expert_matrix = vector(mode = "list", length = self$ncol * self$nrow)
      
      for(index in c(1: (self$ncol * self$nrow) )){
        row_index = ceiling(index/self$ncol)
        col_index = ifelse(index%%self$ncol==0, self$ncol, index%%self$ncol)
        self$expert_matrix[[index]] = ExpertFunction$new( distribution = expert_matrix[row_index, col_index],
                                                          params = expert_params_matrix[row_index, col_index][[1]]
                                                          )
      }
    },
    
    #' select the expert function object inside expert matrix
    #'
    #' @param row_index 
    #' @param col_index 
    select = function(row_index = 0, col_index = 0) {
      # Validate the row index and col index
      assertNumber(row_index, lower = 0, upper = self$nrow)
      assertNumber(col_index, lower = 0, upper = self$ncol)
      mat = self$expert_matrix
      
      if(row_index && col_index) { # Select one entry
        return(self$expert_matrix[[ (row_index-1) * self$ncol + col_index ]])
      }else if(row_index) { # Select one row
        return(mat[ c( ((1-1) * self$ncol + 1) : (1 * self$ncol)) ])
      }else if(col_index) { # Select one column
        return(mat[ seq(col_index, (self$nrow-1) * self$ncol + col_index, self$ncol) ])
      }
    },
    
    #' Exposurize all the expert functions inside the expert matrix
    #'
    #' @param exposure The exposurized value
    exposurize = function(exposure) {
      result = self$clone()
      for(index in c(1: (self$ncol * self$nrow) )){
        result$expert_matrix[[index]] = self$expert_matrix[[index]]$exposurize(exposure)
      }
      return(result)
    },
    
    #' @description
    #' Return the mean matrix of the ExpertMatrix
    get_mean = function() {
      mean_matrix = matrix(0, nrow = self$nrow, ncol = self$ncol)
      for(i in c(1:self$nrow)){
        for(j in c(1:self$ncol)){
          mean_matrix[i,j] = self$select(i,j)$get_mean()
        }
      }
      return(mean_matrix)
    },
    
    #' @description
    #' Return the variance matrix of the ExpertMatrix
    get_variance = function() {
      var_matrix = matrix(0, nrow = self$nrow, ncol = self$ncol)
      for(i in c(1:self$nrow)){
        for(j in c(1:self$ncol)){
          var_matrix[i,j] = self$select(i,j)$get_variance()
        }
      }
      return(var_matrix)
    },
    
    #' set the penalty parameters for all the expert functions
    #'
    #' @param expert_penalty_params_matrix (`matrix`)
    #' To be more specific, it is a matrix of list like this
    #' 
    #' penalty_params = matrix( list( c(2,1,2), c(1,1,1), c(3,1,2), c(1,3,4)), nrow = 2)
    set_penalty_params = function(expert_penalty_params_matrix) {
      for(index in c(1: (self$ncol * self$nrow) )){
        row_index = ceiling(index/self$ncol)
        col_index = ifelse(index%%self$ncol==0, self$ncol, index%%self$ncol)
        self$expert_matrix[[index]]$set_penalty_params(expert_penalty_params_matrix[row_index, col_index][[1]]) 
        self$expert_matrix[[index]]$initialize_penalty()
      }
    },
    
    #' set the parameters for all the expert functions
    #'
    #' @param expert_params_matrix (`matrix`)
    #' To be more specific, it is a matrix of list of named list like this
    #' 
    #' params = matrix(list( list(meanlog = 1, sdlog = 2),  list(meanlog = 1, sdlog  = 2),
    #'                       list(meanlog = 3, sdlog = 1),  list(meanlog = 5, sdlog = 1)), nrow = 2)
    set_params = function(expert_params_matrix) {
      for(index in c(1: (self$ncol * self$nrow) )){
        row_index = ceiling(index/self$ncol)
        col_index = ifelse(index%%self$ncol==0, self$ncol, index%%self$ncol)
        self$expert_matrix[[index]]$set_params(expert_params_matrix[row_index, col_index][[1]]) 
      }
    },
    
    #' Get the sum of penalty value of all the expert functions. Users are required to set the penalty parameters before call this function.
    #' @return value
    #' The total penalty value of Expert Matrix
    get_penalty_value = function() {
      value = 0
      for(expert in self$expert_matrix) {
        value = value + expert$get_penalty()
      }
      return(value)
    },
    
    #' Count the total number of parameters (for the Expert Functions) that exist in the Expert Matrix.
    #' @return result (`integer`)
    #' Total number of parameters
    count_params = function() {
      result = 0
      for(expert in self$expert_matrix){
        result = result + length(expert$get_params())
      }
      return(result)
    }
  )
)
