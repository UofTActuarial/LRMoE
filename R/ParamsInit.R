# Use kmeans to cluster normalized covariates
cluster_covariates <- function(X, n_cluster) {
  means = apply(X, 2, FUN=mean)
  sds = apply(X, 2, FUN=sd)
  X_std = X
  for(col in 1:ncol(X)){
    X_std[,col] = (X[,col] - means[col]) / sds[col]
  }
  X_std[is.na(X_std)] = 0 # Intercept becomes NaN after standardization
  tmp = kmeans(X_std, n_cluster)
  return(list(groups = tmp$cluster, prop = tmp$size/sum(tmp$size)))
}

# Initialize parameters according to the type of observation
params_init_switch <- function(Y_j, type_j) {
  default_expert_continuous = c("lognormal",
                                "gamma",
                                "inversegaussian",
                                "weibull",
                                "burr",
                                "zilognormal",
                                "zigamma",
                                "ziinversegaussian",
                                "ziweibull",
                                "ziburr")
  default_expert_discrete = c("poisson",
                              "negativebinomial",
                              "binomial",
                              "gammacount",
                              "zipoisson",
                              "zinegativebinomial",
                              "zibinomial",
                              "zigammacount")
  if(type_j=="continuous"){
    result = list()
    for(i in 1:length(default_expert_continuous)){
      result[[i]] = list(distribution = default_expert_continuous[i],
                         params = do.call(paste0(default_expert_continuous[i], ".params_init"),
                                          list(y = Y_j)))
    }
    return(result)
  }else if(type_j=="discrete"){
    result = list()
    for(i in 1:length(default_expert_discrete)){
      result[[i]] = list(distribution = default_expert_discrete[i],
                         params = do.call(paste0(default_expert_discrete[i], ".params_init"),
                                          list(y = Y_j)))
    }
    return(result)
  }else{
    stop("Invalid specification of distribution types.")
  }
}


# tmp_continuous = LRMoE:::params_init_switch(Y[,6], "continuous")
# tmp_discrete = LRMoE:::params_init_switch(Y[,2], "discrete")

cmm_transform_inexact_Y <- function(Y) {
  d = ncol(Y)/4
  result = matrix(0, nrow = nrow(Y), ncol = d)
  for(i in 1:d){
    result[,i] = pmin(Y[,4*(i-1)+4], 0.5*(Y[,4*(i-1)+2]+Y[,4*(i-1)+3]))
  }
  return(result)
}

# tmpY = cmm_transform_inexact_Y(Y)

cmm_init_exact <- function(Y, X, n_comp, type) {
  n_dim = ncol(Y)
  n_cov = ncol(X)
  cluster_result = cluster_covariates(X, n_comp)
  label = cluster_result$groups

  # Initialize alpha: a constant term according to prop, last class is reference
  alpha_init = matrix(0, nrow = n_comp, ncol = n_cov)
  alpha_init[,1] = log(cluster_result$prop) - log(cluster_result$prop[n_comp])

  # Summary Statistics
  zero_y = matrix(0, nrow = n_dim, ncol = n_comp)
  mean_y_pos = matrix(0, nrow = n_dim, ncol = n_comp)
  var_y_pos = matrix(0, nrow = n_dim, ncol = n_comp)
  skewness_y_pos = matrix(0, nrow = n_dim, ncol = n_comp)
  kurtosis_y_pos = matrix(0, nrow = n_dim, ncol = n_comp)

  # Params Init
  params_init = list()

  # ll of Init Experts
  ll_init = list()

  # Best init based on ll
  ll_best = list()

  for(d in 1:n_dim){
    params_init[[d]] = list()
    ll_init[[d]] = list()
    ll_best[[d]] = list()
    for(j in 1:n_comp){
      Y_comp = Y[label==j,d]
      zero_y[d,j] = sum(Y_comp==0)/length(Y_comp)
      mean_y_pos[d,j] = mean(Y_comp[Y_comp>0])
      var_y_pos[d,j] = var(Y_comp[Y_comp>0])
      skewness_y_pos[d,j] = skewness(Y_comp[Y_comp>0])
      kurtosis_y_pos[d,j] = kurtosis(Y_comp[Y_comp>0])

      params_init[[d]][[j]] = params_init_switch(Y_comp, type[d])

      ll_init[[d]][[j]] = list()
      ll_init[[d]][[j]][[1]] = rep(-Inf, length(params_init[[d]][[j]]))
      for(k in 1:length(params_init[[d]][[j]])){
        tmp_expert = ExpertFunction$new(params_init[[d]][[j]][[k]]$distribution,
                                        params_init[[d]][[j]][[k]]$params)
        ll_init[[d]][[j]][[1]][k] = sum(tmp_expert$ll_exact(Y_comp))
      }

      max_idx = which.max(ll_init[[d]][[j]][[1]])
      ll_best[[d]][[j]] = params_init[[d]][[j]][[max_idx]]
    }
  }



  return(list(alpha_init = alpha_init,
              params_init = params_init,
              ll_init = ll_init,
              ll_best = ll_best,
              zero_y = zero_y,
              mean_y_pos = mean_y_pos,
              var_y_pos = var_y_pos,
              skewness_y_pos = skewness_y_pos,
              kurtosis_y_pos = kurtosis_y_pos))
}


#' Initialize an LRMoE model based on Clustered Method of Moments
#' @param Y A N by d (\code{exact_Y=T}) or N by 4d (\code{exact_Y=F}) matrix of numerics,
#'          where N is sample size and d is the dimension of each obsevation.
#'          If the size is N by 4d, Each block of four columns should be organized as \code{(tl, yl, yu, tu)}, representing the
#'          truncation lower bound, censoring lower bound, censoring upper bound and truncation upper bound.
#' @param X X A N*P matrix of numerics, where P is the number of covariates.
#'          The first column of \code{X} should be 1, which is the intercept.
#' @param n_comp Numeric, representing how many latent classes/groups to use.
#' @param type A vector of strings of either "continuous" or "exact", representing
#'             whether each dimension of \code{Y} is continuous or discrete.
#' @param exact_Y TRUE/FALSE: whether \code{Y} is observed exactly, or with censoring and/or truncation.
#'
#' @return A list where \code{zero_y}, \code{mean_y_pos}, \code{var_y_pos},
#'         \code{skewness_y_pos} and \code{kurtosis_y_pos} represent the summary statistics of \code{Y}
#'         by dimension and by component.
#'         \code{alpha_init} and \code{experts_init} represents the parameter initializations.
#'         \code{ll_init} contains the loglikelihood of each experts fitted to \code{Y} by
#'         dimension and by component.
#'         \code{ll_best} suggests an initialization of expert functions based on the best loglikelihood.
#' @export
cmm_init <- function(Y, X, n_comp, type, exact_Y = FALSE) {
  if(exact_Y){
    Y_transform = Y
  }else{
    Y_transform = cmm_transform_inexact_Y(Y)
  }
  tmp = cmm_init_exact(Y = Y_transform, X = X, n_comp = n_comp, type = type)

  return(list(zero_y = tmp$zero_y,
              mean_y_pos = tmp$mean_y_pos,
              skewness_y_pos = tmp$skewness_y_pos,
              kurtosis_y_pos = tmp$kurtosis_y_pos,
              alpha_init = tmp$alpha_init,
              experts_init = tmp$params_init,
              ll_init = tmp$ll_init,
              ll_best = tmp$ll_best))
}



# tmp_init = LRMoE::cmm_init(Y = tmpY, X = X,
#                             n_comp = 3, type = c("discrete", "continuous"),
#                             exact_Y = TRUE)
