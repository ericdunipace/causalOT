weighted.var <- function(x, weights) {
  weights <- weights/sum(weights)
  v2 <- sum(weights^2)
  mean <- c(crossprod(x, weights))
  wt.ss <- c(crossprod((x-mean)^2, weights))
  var <- wt.ss/(1-v2)
  return(var)
}
colVar <- function(x, ...) {
  dots <- list(...)
  weights <- dots$weights
  if(is.null(weights) & length(dots)>0) weights <- dots[[1]]
  if(is.null(weights)){
    return(apply(x,2,var))
  } else {
    return(apply(x, 2, weighted.var, weights = weights))
  }
}
rowVar <- function(x, ...) {
  dots <- list(...)
  weights <- dots$weights
  if(is.null(weights) & length(dots)>0) weights <- dots[[1]]
  if(is.null(weights)){
    return(apply(x,1,var))
  } else {
    return(apply(x, 1, weighted.var, weights = weights))
  }
}
colSD <- function(x, ...) {
  dots <- list(...)
  weights <- dots$weights
  if(is.null(weights) & length(dots)>0) weights <- dots[[1]]
  if(is.null(weights)){
    return(apply(x,2,sd))
  } else {
    return(sqrt(apply(x, 2, weighted.var, weights = weights)))
  }
}
rowSD <- function(x, ...) {
  dots <- list(...)
  weights <- dots$weights
  if(is.null(weights) & length(dots)>0) weights <- dots[[1]]
  if(is.null(weights)){
    return(apply(x,1,sd))
  } else {
    return(sqrt(apply(x, 1, weighted.var, weights = weights)))
  }
}