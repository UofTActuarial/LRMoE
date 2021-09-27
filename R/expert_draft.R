# Distribution Expert Function
# 17 Functions need to be implemented
######################################################################################
# Initialize the parameters of the distribution
######################################################################################
distribution.params_init <- function(y){
  # Estimate all the parameters that needed for further calculation
  list() # Tha named list contained all the parameters of distribution
}

distribution.exposurize <- function(params, exposure){
  # Calculate the exposures
  params # The parameters after exposurize
}

distribution.set_params <- function(params){
  #Check the parameters are valid for distribution
  params # The paramsters after validation check
}

######################################################################################
# Calculate the log likelihood and initialize the penalty function
######################################################################################
distribution.expert_ll_exact <- function(y, params){
  # If all the observations are exact, calculate its corresponding likelihood
  NULL # The loglikelihood of distribution based on y
}

distribution.expert_ll_not_exact <- function(tl, tu, yl, yu, params){
  # If some observations are not exact, calculate its corresponding likelihood
  list(expert_ll = NULL, expert_tn = NULL, expert_tn_bar = NULL) # The required variables
}

distribution.penalty <- function(params, penalty_params) {
  # Return the penalty applied on the parameters.
  # Keep in mind to set the default penalty parameters if length(penalty_params) == 0
  if(!length(penalty_params)) { penalty_params = c() }
  NULL # The penal value
}

######################################################################################
# ddistribution, pdistribution, qdistribution and rdistribution implementations.
######################################################################################
distribution.simulation <- function(params, n) {
  # simulated n points based on the distribution params
  NULL # The simulations
}

distribution.mean <- function(params) {
  # Calculate the mean based on the params
  NULL # Mean
}

distribution.variance <- function(params) {
  # Calculate the variance based on the params
  NULL # Variance
}

distribution.logpdf <- function(params, x) {
  # Return the log pdf based on the input x
  NULL # log pdf
}

distribution.pdf <- function(params, x) {
  # Return the pdf based on the input x
  NULL # pdf
}

distribution.logcdf <- function(params, q) {
  # return the log cdf based on the input x
  NULL # log cdf
}

distribution.cdf <- function(params, q) {
  # return the cdf based on the input x
  NULL # cdf
}

distribution.quantile <- function(params, p) {
  # return the percentage points based on the value of p
  NULL # quantile
}

# ATTENTION: ONLY need to implement for continuous function
distribution.lev <- function(params, u) {
 # Limited Expected Value 
}

######################################################################################
# E Step, M Step and EM Optimization steps.
######################################################################################
distribution._EStep <- function() {
  # Perform the E step
  NULL
}

distribution._MStep <- function() {
  # Perform the M step
  NULL
}

distribution.compute_EM <- function() {
  # Perform the EM optimization
  NULL
}
######################################################################################
# Register the distribution at zzz.R to the ExpertLibrary Object (Examples included)
######################################################################################