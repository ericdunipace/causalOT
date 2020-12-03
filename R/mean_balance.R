mean_bal <- function(data, weights, ...) {
  
  if(!inherits(weights, "causalWeights")){
    if(!all(names(weights) %in% c("w0","w1"))) stop("weights must be of class 'causalWeights' or be a list with named slots 'w0' and 'w1'.")
  }
  
  xs <- extract_x(data, ...)
  x1 <- as.matrix(xs$x1)
  x0 <- as.matrix(xs$x0)
  
  sigma_x1 <- colVar(x1)
  sigma_x2 <- colVar(x0)
  
  pool_sd <- sqrt(sigma_x1*0.5 + sigma_x2 *0.5)
  
  mean_1 <- c(crossprod(x1, weights$w1))
  mean_0 <- c(crossprod(x0, weights$w0))
  
  return(abs(mean_1 - mean_0) /pool_sd)
}


