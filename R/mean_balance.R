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
  
  sigma_x1 <- colVar(x1)
  sigma_x2 <- colVar(x0)
  
  pool_sd <- sqrt(sigma_x1 * 0.5 + sigma_x2 * 0.5)
  
  mean_1 <- matrixStats::colWeightedMeans(x1, weights$w1)
  mean_0 <- matrixStats::colWeightedMeans(x0, weights$w0)
  
  if (any(pool_sd == 0)) {
    pool_sd[pool_sd == 0] <- Inf
  }
  
  return(abs(mean_1 - mean_0) / pool_sd)
}


