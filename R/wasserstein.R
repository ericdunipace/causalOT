wasserstein_p.default <- function(a, b, p = 1, tplan = NULL, cost = NULL, ...) {
  
  if (is.numeric(a)) {
    if (!isTRUE(all.equal(sum(a) ,1))) a <- renormalize(a)
  }
  if (is.numeric(b)) {
    if (!isTRUE(all.equal(sum(b) ,1))) b <- renormalize(b)
  }
  if (is.null(cost)) {
    stop("cost matrix must be specified if only masses are given")
  }
  # cost_a <- list(...)$cost_a
  # cost_b <- list(...)$cost_b
  # if (is.null(tplan)) {
  #   if ( isFALSE(list(...)$neg.weights) ) {
  #     nzero_row <- a > 0
  #     nzero_col <- b > 0
  #   } else {
  #     nzero_row <- a != 0
  #     nzero_col <- b != 0
  #   }
  #   a <- a[nzero_row]
  #   b <- b[nzero_col]
  #   cost <- cost[nzero_row, nzero_col, drop = FALSE]
  #   
  #   if (!is.null(cost_a) && !is.null(cost_b)) {
  #     cost_a <- cost_a[nzero_row, nzero_row, drop = FALSE]
  #     cost_b <- cost_b[nzero_col, nzero_col, drop = FALSE]
  #   }
  # }
  n_a <- length(a)
  n_b <- length(b)  
  if (n_a == 1) {
    return(c((sum(c(cost^p) * c(b)))^(1/p)))
  } else if ( n_b == 1) {
    return(c((sum(c(cost^p) * c(a)))^(1/p)))
  } else {
    # return(transport::wasserstein(a = a, b = b, p = p, tplan = tplan, costm = cost, prob = TRUE,...))
    return(approxOT::wasserstein(a = a, b = b, p = p, ground_p = p, tplan = tplan, cost = cost, 
                                 ...))
  }
}

wasserstein_p.default <- function(a, b, p = 1, tplan = NULL, cost = NULL, ...) {
  
  if (is.numeric(a)) {
    if (!isTRUE(all.equal(sum(a) ,1))) a <- renormalize(a)
  }
  if (is.numeric(b)) {
    if (!isTRUE(all.equal(sum(b) ,1))) b <- renormalize(b)
  }
  if (is.null(cost)) {
    stop("cost matrix must be specified if only masses are given")
  }
  # cost_a <- list(...)$cost_a
  # cost_b <- list(...)$cost_b
  # if (is.null(tplan)) {
  #   if ( isFALSE(list(...)$neg.weights) ) {
  #     nzero_row <- a > 0
  #     nzero_col <- b > 0
  #   } else {
  #     nzero_row <- a != 0
  #     nzero_col <- b != 0
  #   }
  #   a <- a[nzero_row]
  #   b <- b[nzero_col]
  #   cost <- cost[nzero_row, nzero_col, drop = FALSE]
  #   
  #   if (!is.null(cost_a) && !is.null(cost_b)) {
  #     cost_a <- cost_a[nzero_row, nzero_row, drop = FALSE]
  #     cost_b <- cost_b[nzero_col, nzero_col, drop = FALSE]
  #   }
  # }
  n_a <- length(a)
  n_b <- length(b)  
  if (n_a == 1) {
    return(c((sum(c(cost^p) * c(b)))^(1/p)))
  } else if ( n_b == 1) {
    return(c((sum(c(cost^p) * c(a)))^(1/p)))
  } else {
    
    
      return(approxOT::wasserstein(a = a, b = b, p = p, ground_p = p, tplan = tplan, cost = cost,
                                   ...))
    
  }
}

wasserstein_p.causalWeights <- function(a, b = NULL, p = 1, tplan = NULL, cost = NULL,...) {
  mass_a <- as.numeric(a$w0)
  mass_b <- as.numeric(a$w1)
  if((a$estimand == "feasible") & !is.null(a$gamma)){
    idx <- which(a$gamma != 0, arr.ind = TRUE)
    tplan <- data.frame(from = idx[,1], to = idx[,2], mass = a$gamma[idx])
  }
  return(wasserstein_p.default(a = mass_a, b = mass_b, p = p, tplan = tplan, cost = cost, ...))
}

wasserstein_p.matrix <- function(a, b, p = 1, tplan = NULL, cost = NULL, dist = "Lp",...) {
  if(is.null(cost)) {
    cost.calc <- switch(dist, "Lp" = cost_calc_lp,
                        "mahalanobis" = cost_mahalanobis,
                        "sdLp" = cost_calc_sdlp)
    cost <- cost.calc(a, b, p, direction = "rowwise")
  }
  if(!is.null(tplan)) {
    mass_a <- as.numeric(tapply(tplan$mass, factor(tplan$from), sum))
    mass_b <- as.numeric(tapply(tplan$mass, factor(tplan$to), sum))
  } else {
    mass_a <- rep(1/nrow(a), nrow(a))
    mass_b <- rep(1/nrow(b), nrow(b))
  }
  return(wasserstein_p.default(a = mass_a, b = mass_b, p = p, tplan = tplan, cost = cost, ...))
}

setGeneric("wasserstein_p", function(a, b, ...) UseMethod("wasserstein_p"))  
setMethod("wasserstein_p", signature(a = "vector", b= "vector"), wasserstein_p.default)
setMethod("wasserstein_p", signature(a = "causalWeights"), wasserstein_p.causalWeights)
setMethod("wasserstein_p", signature(a = "causalWeights", b= "NULL"), wasserstein_p.causalWeights)
setMethod("wasserstein_p", signature(a = "causalWeights", b= "numeric"), wasserstein_p.causalWeights)
setMethod("wasserstein_p", signature(a = "causalWeights", b= "matrix"), wasserstein_p.causalWeights)
setMethod("wasserstein_p", signature(a = "matrix", b= "matrix"), wasserstein_p.matrix)
