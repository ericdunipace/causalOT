#' Standardized absolute mean difference calculations
#' 
#' This function will calculate the difference in means between treatment groups standardized by the pooled standard-deviation of the respective covariates.
#'
#' @param x Either a matrix, an object of class [causalOT::dataHolder-class], or an object of class DataSim
#' @param z A integer vector denoting the treatments of each observations. Can be null if `x` is a DataSim object or already of class [causalOT::dataHolder-class].
#' @param weights An object of class [causalWeights][causalOT::causalWeights-class].
#' @param ... Not used at this time.
#'
#' @return A vector of mean balances
#' @export
#'
#' @examples
#' n <- 100
#' p <- 6
#' x <- matrix(stats::rnorm(n * p), n, p)
#' z <- stats::rbinom(n, 1, 0.5)
#' weights <- calc_weight(x = x, z = z, estimand = "ATT", method = "Logistic")
#' mb <- mean_balance(x = x, z = z, weights = weights)
#' print(mb)
mean_balance <- function(x = NULL, z = NULL, weights = NULL, ...) {
  
  is_cw <- inherits(weights, "causalWeights")
  if( is_cw ) {
    w <- NULL
  } else {
    w <- weights
  }
  
  if(is.null(x) && is_cw) {
    dH <- weights@data
  } else {
    dH <- dataHolder(x = x, z = z, weights = w)
  }
  
  x1 <- get_x1(dH)
  x0 <- get_x0(dH)
  
  w  <- get_w(dH)
  z  <- get_z(dH)
  
  if ( is_cw ) {
    w1 <- weights@w1
    w0 <- weights@w0
  } else {
    w1 <- renormalize(w[z==1])
    w0 <- renormalize(w[z==0])
  } 
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  n  <- n0 + n1
  
  sigma_x1 <- matrixStats::colWeightedVars(x1, renormalize(w[z==1]) * n1)
  sigma_x0 <- matrixStats::colWeightedVars(x0, renormalize(w[z==0]) * n0)
  
  if (n0 > 1 && n1 > 1) {
    pool_sd <- sqrt(sigma_x1 * n1/n + sigma_x0 * n0/n)
  } else {
    pool_sd <- sqrt(matrixStats::colWeightedVars(get_x(dH), w * n))
  } 
  
  mean_1 <- matrixStats::colWeightedMeans(x1, w1)
  mean_0 <- matrixStats::colWeightedMeans(x0, w0)
  
  if (any(pool_sd == 0)) {
    pool_sd[pool_sd == 0] <- Inf
  }
  
  return(abs(mean_1 - mean_0) / pool_sd)
}


mean_balance_diagnostic <- function(object) {
  
  stopifnot(inherits(object, "causalWeights"))
  
  if(object@estimand %in% c("ATT", "ATC")) {
    pre <- mean_balance(x = get_x(object@data),
                        z = get_z(object@data),
                        weights = get_w(object@data))
    post <- mean_balance(weights = object)
  } else if (object@estimand == "ATT") {
    x0 <- get_x0(object@data)
    x1 <- get_x1(object@data)
    z  <- get_z(object@data)
    w  <- get_w(object@data)
    w1 <- renormalize(w[z==1])
    w0 <- renormalize(w[z==0])
    n0 <- get_n0(object@data)
    n1 <- get_n1(object@data)
    n  <- n0 + n1
    
    sds <- sqrt(matrixStats::colWeightedVars(x0, w = w0 * n0) * n0/n + 
      n1/n * matrixStats::colWeightedVars(x1, w = w1 * n1))
    
    full.means <- crossprod(x, w)/sds
    
    pre <- crossprod(x0, w0)/sds - crossprod(x1, w1)/sds
    post<- crossprod(x0, object@w0)/sds - crossprod(x1, object@w1)/sds
    
  } else if (object@estimand == "ATE") {
    x <- get_x(object@data)
    x0 <- get_x0(object@data)
    x1 <- get_x1(object@data)
    z  <- get_z(object@data)
    w  <- get_w(object@data)
    w1 <- renormalize(w[z==1])
    w0 <- renormalize(w[z==0])
    n  <- get_n(object@data)
    
    sds <- matrixStats::colWeightedSds(x, w = w * n)
    
    full.means <- crossprod(x, w)/sds
    
    pre <- list()
    pre$controls  <- abs(full.means - crossprod(x0, w0)/sds)
    pre$treated   <- abs(full.means - crossprod(x1, w1)/sds)
    
    post <- list()
    post$controls <- abs(full.means - crossprod(x0, object@w0)/sds)
    post$treated  <- abs(full.means - crossprod(x1, object@w1)/sds)
    
    
  }
  
  
  res <- list(pre = pre, post = post)
  
  return(res)
  
}


