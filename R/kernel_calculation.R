kernel_calculation <- function(X, z, p = 1.0, 
                               theta = NULL, gamma = NULL, sigma_2 = NULL, 
                               metric = c("mahalanobis","Lp"),
                               kernel = c("RBF", "polynomial", "linear"),
                               is.dose = FALSE, 
                               estimand = c("ATE","ATT","ATC")) {
  
  if(!is.matrix(X)) X <- as.matrix(X)
  if(!is.matrix(z)) z <- as.matrix(z)
  
  if(!is.logical(is.dose)) is.dose <- isTRUE(is.dose)
  
  estimand <- match.arg(estimand)
  
  if(nrow(X) != nrow(z)) stop("Observations of X and z must be equal")
  
  met <- match.arg(metric)
  kern <- match.arg(kernel)
  calc_covariance <- isTRUE(met == "mahalanobis")
  
  theta <- kernel_param_check(theta)
  gamma <- kernel_param_check(gamma)
  sigma_2 <- kernel_sigma_check(sigma_2, nrow(X), z, is.dose)
  p <- kernel_power_check(p)
  
  if(is.dose) {
    return(kernel_calc_dose_(X_ = X, z_ = z, p = p, 
                 theta_ = theta, gamma_ = gamma,
                 kernel_ = kern,
                 calc_covariance = calc_covariance))
  } else {
    return(kernel_calc_(X_ = X, z = z, p = p, 
                      theta_ = theta, gamma_ = gamma,
                      sigma_2_ = sigma_2,
                      kernel_ = kern,
                      calc_covariance = calc_covariance,
                      estimand = estimand))
  }
  
}

kernel_param_check <- function(param) {
  
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

kernel_power_check <- function(p) {
  
  if(missing(p)) p <- 1.0
  if(is.null(p)) p <- 1.0
  if(is.na(p)) p <- 1.0
  
  if (!is.double(p) ) p <- as.double(p)
  return(p)
}

kernel_sigma_check <- function(s, n, z, is.dose) {
  if(!is.null(s)) {
    if(length(s) == 1) s <- rep(s, n)
    if(length(s) == 2) {
      s_temp <- s
      s <- rep(0, n)
      s[z==1] <- s_temp[[2]]
      s[z==0] <- s_temp[[1]]
    }
    return(as.double(s))
  } else if (!is.dose) {
    stop("sigma^2 must be specified if using binary treatment kernel")
  } else {
    return(s)
  }
}

calc_similarity <- function( X, z, metric = c("mahalanobis","Lp"), 
                             kernel = c("RBF","polynomial","linear"),
                             is.dose = FALSE,
                             estimand = c("ATE","ATT","ATC")) {
  met <- match.arg(metric)
  estimand <- match.arg(estimand)
  kernel <- match.arg(kernel)
  
  if(!is.matrix(X)) X <- as.matrix(X)
  
  if(!is.logical(is.dose)) is.dose <- isTRUE(is.dose)
  
  calc_covariance <- isTRUE(met == "mahalanobis")
  
  if(is.dose) {
    if(!is.matrix(z)) z <- as.matrix(z)
    
    if(nrow(X) != length(z)) stop("Observations of X and z must be equal")
    if(estimand != "ATE") warning("RKHS dose method can only do ATE. ATE will be returned.")
    return( similarity_calc_dose_(X_ = X, z_ = z,
                       calc_covariance = calc_covariance) )
  } else {
    if(!is.integer(z)) z <- as.integer(z)
    
    if(kernel == "polynomial" | kernel == "linear") {
      return( similarity_calc_(X_ = X, z = z,
                             calc_covariance = calc_covariance,
                             estimand = estimand) )
    } else if (kernel == "RBF") {
      if(calc_covariance) {
        return(cost_mahalanobis(X,X, ground_p = 2, direction = "rowwise")^2)
      } else {
        return(cost_calc_lp(X,X, ground_p = 2, direction = "rowwise")^2)
      }
    }
  }
}



ot_kernel_calculation <- function(X, z, p = 1.0, 
                                  theta = NULL, gamma = NULL, 
                                  # sigma_2 = NULL, 
                                  kernel = c("RBF","polynomial","linear"),
                                  metric = c("mahalanobis","Lp"),
                                  is.dose = FALSE, 
                                  estimand = c("ATE","ATT","ATC")) {
  
  if(!is.matrix(X)) X <- as.matrix(X)
  if(!is.matrix(z)) z <- as.matrix(z)
  
  if(!is.logical(is.dose)) is.dose <- isTRUE(is.dose)
  
  estimand <- match.arg(estimand)
  kernel <- match.arg(kernel)
  
  if(nrow(X) != nrow(z)) stop("Observations of X and z must be equal")
  
  met <- match.arg(metric)
  calc_covariance <- isTRUE(met == "mahalanobis")
  
  theta <- kernel_param_check(theta)
  gamma <- kernel_param_check(gamma)
  # sigma_2 <- kernel_sigma_check(sigma_2, nrow(X), is.dose)
  p <- kernel_power_check(p)
  
  if(is.dose) {
    if(kernel != "polynomial") {
      kernel <- "polynomial"
      warning("Kernel needs to be polynomial for dose method.")
    }
    return(kernel_calc_dose_(X_ = X, z_ = z, p = p, 
                             theta_ = theta, gamma_ = gamma,
                             kernel_ = kernel,
                             calc_covariance = calc_covariance)$cov_kernel)
  } else {
    orders <- order(z)
    X <- X[orders,,drop=FALSE]
    z <- z[orders,,drop=FALSE]
    return(kernel_calc_ot_(X_ = X, z = z, p = p, 
                        theta_ = theta, gamma_ = gamma,
                        # sigma_2 = sigma_2,
                        kernel_ = kernel,
                        calc_covariance = calc_covariance,
                        estimand = estimand))
  }
  
}