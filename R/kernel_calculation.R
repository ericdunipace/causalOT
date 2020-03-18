kernel_calculation <- function(X, z, d = 2, theta = NULL, gamma = NULL, metric = c("Lp","mahalanobis")) {
  
  dir <- match.arg(direction)
  
  if(!is.matrix(X)) X <- as.matrix(X)
  if(!is.matrix(z)) z <- as.matrix(z)
  
  if(nrow(X) != nrow(z)) stop("Observations of X and z must be equal")
  
  met <- match.arg(metric)
  calc_covariance <- isTRUE(met == "mahalanobis")
  
  theta <- kernel_param_check(theta)
  gamma <- kernel_param_check(gamma)
  d <- kernel_power_check(d)
  
  return( kernel_calc_(X_ = X, z_ = z, d = d, 
                       theta_ = theta, gamma_ = gamma,
                       calc_covariance = calc_covariance) )
}

kerenel_param_check <- function(param) {
  
  if(is.null(param)) {
    param <- as.double(c(1.0,1.0))
  }
  if(length(param) != 2) {
    if(length(param) > 2) {
      param <- param[1:2]
    } else if(length(param) == 1) {
      param <- c(param[1], param[2])
    } else if (length(param) == 0) {
      warning("Theta or gamma is of length 0, filling in defaults")
      param <- as.double(c(1.0,1.0))
    }
    warning("Theta or gamma not of length 2. Doubling first element")
  }
  if(!is.double(param)) param <- as.double(param)
  
  return(param)
}

kerenel_power_check <- function(d) {
  
  if(is.missing(d)) d <- 1.0
  if(is.null(d)) d <- 1.0
  
  if (!is.double(d) ) d <- as.double(d)
  return(d)
}