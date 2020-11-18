cost_fun <- function(x, y, power = 2, metric = c("mahalanobis","Lp","RKHS"), rkhs.args = NULL, estimand = "ATE", ...) {
  direction <- "rowwise"
  metric <- match.arg(metric)
  
  dist <- switch(metric, 
         "Lp" = causalOT::cost_calc_lp(x,y,ground_p = power, direction = direction),
         "mahalanobis" = causalOT::cost_mahalanobis(x,y, ground_p = power, direction = direction),
         "RKHS" = causalOT::cost_RKHS(X=x, Y=y, rkhs.args = rkhs.args, estimand = estimand, ...))
  
  return(dist)
  
}

cost_calc_lp <- function(X, Y, ground_p = 2, direction = c("rowwise", "colwise")) {
  
  dir <- match.arg(direction)
  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if(!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if(dir == "rowwise") {
    if (ncol(X) != ncol(Y)) {
      stop("ncol X and Y should be equal. They can different numbers of observations")
    }
    X <- t(X)
    Y <- t(Y)
  } else {
    if (nrow(X) != nrow(Y)) {
      stop("Number of covariates of X and Y should be equal.")
    }
  }
  stopifnot(ground_p > 0)
  
  return(causalOT::cost_calculation_(X,Y,as.double(ground_p))) 
}

cost_mahalanobis <- function(X, Y, ground_p = 2, direction = c("rowwise", "colwise")) {
  
  dir <- match.arg(direction)
  
  if(dir == "rowwise") {
    if (ncol(X) != ncol(Y)) {
      stop("ncol X and Y should be equal. They can different numbers of observations")
    }
    X <- t(X)
    Y <- t(Y)
  }
  
  if (!is.double(ground_p) ) ground_p <- as.double(ground_p)
  if (nrow(X) != nrow(Y)) {
    stop("Number of covariates of X and Y should be equal.")
  }
  
  return( cost_mahal_(X, Y, ground_p) )
}

cost_RKHS <- function(X, Y, 
                      kernel = c("RBF", "polynomial", "linear"),
                      rkhs.args = list(p = 1.0, 
                                        theta = c(1,1), 
                                        gamma = c(1,1), 
                                        is.dose = FALSE),
                      estimand = "ATE", ...
                      ){
  if(is.null(rkhs.args)) stop("rkhs args must be specified")
  kernel <- match.arg(kernel)
  # dir <- match.arg(direction)
  # 
  # if(dir == "rowwise") {
  #   if (ncol(X) != ncol(Y)) {
  #     stop("ncol X and Y should be equal. They can different numbers of observations")
  #   }
  #   X <- t(X)
  #   Y <- t(Y)
  # }
  
  if ( ncol(X) != ncol(Y) ) {
    stop("Number of covariates of X and Y should be equal.")
  }
  
  covars <- rbind(X,Y)
  z <- c(rep(0,nrow(X)), rep(1, nrow(Y)))
  
  dist <- ot_kernel_calculation(X = covars, z = z, 
                                p = rkhs.args$p, 
                               theta = rkhs.args$theta, 
                               gamma = rkhs.args$gamma,
                               kernel = kernel,
                               metric = "mahalanobis",
                               is.dose = rkhs.args$is.dose, 
                               estimand = estimand)
  if(estimand != "ATE") {
    dist <- dist[[1]]
  #   dist <- dist[[1]][1:nrow(X), (nrow(X) + 1):nrow(covars)]
  }
  return( dist )
}