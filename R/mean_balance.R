#' Standardized absolute mean difference calculations
#' 
#' This function will calculate the difference in means between treatment groups standardized by the pooled standard-deviation of the respective covariates.
#'
#' @param data Either a data.frame, matrix, or object of class DataSim
#' @param weights An object of class causalWeights or a list with slots w0 and w1.
#' @param ... Additional arguments passed to the function to know which covariates are to be balanced ("balance.covariates"), and treatment indicator ("treatment.indicator"). These can be column names or numbers. These arguments are not needed if using the causalOT DataSim class.
#'
#' @return A vector of mean balances
#' @export
#'
#' @examples
#' 
#' n <- 100
#' p <- 6
#' x0 <- matrix(rnorm(2*n * p), 2*n, p)
#' x1 <- matrix(rnorm(n * p), n, p)
#' weights <- list(w0 = rep(1/(2*n), 2 * n), w1 = rep(1/n, n))
#' data <- cbind(rbind(x0,x1), z = c(rep(0,2*n), rep(1, n)))
#' colnames(data) <- c(paste0("x", 1:p), "z")
#' mb <- mean_bal(data, weights, balance.covariates = paste0("x", 1:p), 
#'                treatment.indicator = "z")
#' print(mb)
mean_bal <- function(data, weights = NULL, ...) {
  
  
  
  xs <- extract_x(data, ...)
  x1 <- as.matrix(xs$x1)
  x0 <- as.matrix(xs$x0)
  
  if (is.null(weights)) {
    sw <- get_sample_weight(NULL, z = c(rep(1,nrow(x1)), rep(0, nrow(x0))) )
    weights <- list(w0 = sw$a, w1 = sw$b)
  }
  
  if (!inherits(weights, "causalWeights")) {
    if (!any(names(weights) == "w0") & !any(names(weights) == "w1") ) stop("weights must be of class 'causalWeights' or be a list with named slots 'w0' and 'w1'.")
  }
  
  sigma_x1 <- matrixStats::colVars(x1)
  sigma_x0 <- matrixStats::colVars(x0)
  
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  n  <- n0 + n1
  
  if (n0 > 1 && n1 > 1) {
    pool_sd <- sqrt(sigma_x1 * n1/n + sigma_x0 * n0/n)
  } else {
    pool_sd <- sqrt(matrixStats::colVars(rbind(x1, x0)))
  } 
  
  mean_1 <- matrixStats::colWeightedMeans(x1, weights$w1)
  mean_0 <- matrixStats::colWeightedMeans(x0, weights$w0)
  
  if (any(pool_sd == 0)) {
    pool_sd[pool_sd == 0] <- Inf
  }
  
  return(abs(mean_1 - mean_0) / pool_sd)
}


