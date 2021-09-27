##############################################################################
# The S3 Methods for better ExpertFunction and ExpertMatrix summary and print
##############################################################################

# Print the information of the Expert Function
print.ExpertFunction <- function(obj) {
  print(paste("Expert Function Distribution: ", obj$distribution))
  print(paste("Is Continuous?: ", obj$continuous))
  print(paste("Expert Parameters: ", print_nl(obj$get_params()) ))
  print(paste("Expert Penalty Parameters: ", print_vc(obj$get_penalty_params()) ))
  print(paste("Expert Penalty Value: ", obj$get_penalty()))
  print(paste("Expert Mean: ", obj$get_mean(), " Expert Variance: ", obj$get_variance()))
  print(paste("Avaliable methods for current Expert: ", paste(names(obj), collapse=";  " )))
}

# Map the summary(expert_function) to print(expert_function)
summary.ExpertFunction <- function(obj) {
  print(obj)
}

# Print the shorter information of the Expert Function
shorter_summary <- function(obj) {
  print(paste("Expert Function Distribution: ", obj$distribution))
  print(paste("Expert Parameters: ", print_nl(obj$get_params()) ))
  print(paste("Expert Penalty Parameters", print_vc(obj$get_penalty_params()) ))
  print(paste("Expert Penalty Value: ", obj$get_penalty()))
}

# Print the information of the ExpertMatrix
print.ExpertMatrix <- function(obj) {
  limit = 0
  max_print = 60
  for(row in c(1:obj$nrow)){
    for(col in c(1:obj$ncol)){
      print(paste("Expert Summary at position (row=", row, ", col=", col, ")"))
      shorter_summary(obj$select(row, col))
      print("----------------------------------------------------------------")
      limit = limit + 1
      if(limit >= max_print) {
        warning("Reach Maximum Print Limit 60, change limit at R/model_summary.R")
        break
      }
    }
    if(limit >= max_print) { break }
  }
}

# Map summary(expert_matrix) to print(expert_ma)
summary.ExpertMatrix <- function(obj) {
  print(ExpertMatrix)
}