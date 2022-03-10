.cost_fun_deprecated <- function(x, y, power = 2, metric = c("mahalanobis","Lp","RKHS"), 
                     rkhs.args = NULL, estimand = "ATE", ...) {
  direction <- "rowwise"
  metric <- match.arg(metric)
  
  dist <- switch(metric, 
         "Lp" = cost_calc_lp(x,y,ground_p = power, direction = direction),
         "mahalanobis" = cost_mahalanobis(x,y, ground_p = power, direction = direction),
         "RKHS" = cost_RKHS(X = x, Y = y, rkhs.args = rkhs.args, estimand = estimand, ...))
  
  return(dist)
  
}

#' Calculate cost matrix for a given estimand
#'
#' @param x An object of class `matrix`
#' @param z A treatment indicator with values in 0 and 1. Should be of class
#' `integer` or `vector`
#' @param power The power used to calculate the the cost matrix: \eqn{\{(x-y)^{power}\}^{(1/{power})}}
#' @param metric One of the values in [dist.metrics][dist.metrics()].
#' @param estimand The estimand desired for the weighting estimator. See details
#' @param ... Arguments passed to the RKHS calculating function including
#' \itemize{
#' \item `kernel`, one of "RBF", "polynomial", "linear"
#' \item `rkhs.args` The arguments used to construct the kernel
#' }
#' `...` can also be used to handle extra arguments passed by mistake so that
#' an error is not thrown.
#' 
#' @details 
#' If the estimand is "ATT" or "ATC", `cost_fun` will calculate 
#' the cost matrix where the rows are the control
#' and the columns are the treated. If "ATE" will calculate to cost matrices
#' with the first having the rows corresponding to the control individual and the
#' second having rows correspond to the treated individuals. For both matrices,
#' the columns will correspond to the full sample. The dimensions of the output will
#'  depend on the estimand. For reference, let \eqn{n_1 = \sum_i z_i}, 
#'  \eqn{n_0 = \sum_i (1-z_i)}, and \eqn{n = n_1 + n_0}.
#' 
#'
#' @return Output depends on the estimand.
#' * For ATT and ATC: a matrix of dimension \deqn{n_0 \times n_1}{n_0 * n_1}.
#' * For ATE: a list of two matrices of dimension \eqn{n_0 \times n}{n_0 * n} and \eqn{n_1 \times n}.
#' See details for more information.
#' 
#' @export
#'
#' @examples
#' n0 <- 100
#' n1 <- 55
#' d <- 5
#' x1 <- matrix(stats::rnorm(n1*d), n1, d)
#' x0 <- matrix(stats::rnorm(n0*d), n0, d)
#' 
#' x <- rbind(x0,x1)
#' z <- c(rep(0,n0), rep(1,n1))
#' power <- 2.0
#' 
#' # ATT
#' estimand <- "ATT"
#' metric <- "Lp"
#' cost_ATT <- cost_fun(x, z, power = power, metric = metric, estimand = estimand)
#' print(dim(cost_ATT))
#' 
#' # ATE
#' estimand <- "ATE"
#' cost_ATT <- cost_fun(x, z, power = power, metric = metric, estimand = estimand)
#' length(cost_ATT)
cost_fun <- function(x, z, power = 2, metric = dist.metrics(), 
                     estimand = "ATE", ...) {
  metric <- match.arg(metric)
  
  dist <- switch(metric, 
                 "Lp" = cost_metric_calc(x,z,ground_p = power, metric = metric, estimand = estimand),
                 "mahalanobis" = cost_metric_calc(x,z,ground_p = power,  metric = metric, estimand = estimand),
                 "sdLp" = cost_metric_calc(x,z,ground_p = power,  metric = metric, estimand = estimand),
                 "RKHS" = cost_RKHS(X = x, z = z, estimand = estimand, ...))
  
  return(dist)
  
}

cost_metric_calc <- function(x,z,ground_p, metric = c("mahalanobis","Lp","sdLp"), estimand) {
  direction <- "rowwise"
  metric <- match.arg(metric)
  
  cost_function <- switch(metric,
                          "Lp" = "cost_calc_lp",
                          "mahalanobis" = "cost_mahalanobis",
                          "sdLp" = "cost_calc_sdlp")
  
  if (estimand == "ATT" | estimand == "ATC") {
    return(match.fun(cost_function)(X = x[z == 0, , drop = FALSE], 
                                    Y = x[z == 1, , drop = FALSE], ground_p = ground_p, direction = direction,
                                    estimand = estimand))
  } else if (estimand == "ATE") {
    return(list(
      match.fun(cost_function)(X = x[z == 0, , drop = FALSE], Y = x, ground_p = ground_p, 
                               direction = direction,
                               estimand = estimand),
      match.fun(cost_function)(X = x[z == 1, , drop = FALSE], Y = x, ground_p = ground_p, 
                               direction = direction,
                               estimand = estimand)
    ))
  }
}

cost_calc_lp <- function(X, Y, ground_p = 2, direction = c("rowwise", "colwise"), ...) {
  
  dir <- match.arg(direction)
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if (dir == "rowwise") {
    if (ncol(X) != ncol(Y)) {
      stop("ncol X and Y should be equal. They can have different numbers of observations")
    }
    X <- t(X)
    Y <- t(Y)
  } else {
    if (nrow(X) != nrow(Y)) {
      stop("Number of covariates of X and Y should be equal.")
    }
  }
  stopifnot(ground_p > 0)
  
  return(cost_calculation_(A_ = X, B_ = Y, p = as.double(ground_p))) 
}

cost_calc_sdlp <- function(X, Y, ground_p = 2, direction = c("rowwise", "colwise"), estimand = "ATE") {
  dir <- match.arg(direction)
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if (dir == "rowwise") {
    if (ncol(X) != ncol(Y)) {
      stop("ncol X and Y should be equal. They can have different numbers of observations")
    }
    X <- t(X)
    Y <- t(Y)
  } else {
    if (nrow(X) != nrow(Y)) {
      stop("Number of covariates of X and Y should be equal.")
    }
  }
  stopifnot(ground_p > 0)
  
  if (estimand == "ATE") {
    scale <- 1/matrixStats::rowSds(Y)
    center <- rowMeans(Y)
  } else  {
    total <- cbind(X,Y)
    scale <- 1/matrixStats::rowSds(total)
    center <- rowMeans(total)
  } 
  
  X <-  scale * (X - center)
  Y <-  scale * (Y - center)
  
  if (any(is.infinite(scale)) ) {
    nonzero.idx <- is.finite(scale)
    X <- X[nonzero.idx, , drop = FALSE]
    Y <- Y[nonzero.idx, , drop = FALSE]
  }
  
  return(cost_calculation_(A_ = X, B_ = Y, p = as.double(ground_p))) 
}

cost_mahalanobis <- function(X, Y, ground_p = 2, direction = c("rowwise", "colwise"),
                             estimand = "ATT") {
  
  dir <- match.arg(direction)
  
  if (dir == "rowwise") {
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
  
  return( cost_mahal_(A_ = X, B_ = Y, p = ground_p, estimand = estimand) )
}

cost_RKHS <- function(X, z, 
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
  
  # if ( ncol(X) != ncol(Y) ) {
  #   stop("Number of covariates of X and Y should be equal.")
  # }
  
  # covars <- rbind(X,Y)
  # z <- c(rep(0,nrow(X)), rep(1, nrow(Y)))
  
  dist <- ot_kernel_calculation(X = X, z = z, 
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


cost_factor <- function(x, z) {
  x0 <- x[z==0]
  x1 <- x[z==1]
  
  return(sapply(x1, function(y1) sapply(x0, function(y0) as.numeric(y0 == y1))))
}
