#' Barycentric Projections
#'
#' @param data Data used to perform projections. Should be a matrix, data.frame, or DataSim class.
#' @param weight an object of class [causalWeights][causalOT::causalWeights-class] output from the [calc_weight()][calc_weight()] function
#' @param ... optional arguments such as "cost", the cost matrix, or the "estimand", "metric", or "p" if not provided by weights. Arguments "balance.covariates" and "treatment.indicator" must be provided if data is of class data.frame or matrix.
#'
#' @return a list containing the barycentric projections with slots "control", "treated", and "observed.treatment"
#' 
#' @keywords internal
barycentric_projection <- function(data, weight, 
                                   ...) {
  
  stopifnot(inherits(weight, "causalWeights"))
  
  est <- weight$estimand
  pow <- weight$args[["power"]]
  met <- weight$args$metric
  rkhs.args <- weight$args$rkhs.args
  gamma <- weight$gamma
  dots <- list(...)
  
  if (is.null(est)) {
    est <- match.arg(dots$estimand, choices = c("ATT","ATC","ATE"))
    warning("Weights didn't have estimand flag. Make sure weights aligned with desired estimand")
  }
  if (is.null(met)) {
    met <- match.arg(dots$metric, choices = c("mahalanobis","sdLp","Lp", "RKHS"))
  }
  if (is.null(pow)) {
    pow <- as.double(dots[["p"]])
  }
  if (is.null(rkhs.args)) {
    rkhs.args <- as.double(dots$rkhs.args)
  }
  if (is.null(pow) || length(pow) == 0) stop("weights must have `p` parameter (power) or must be specified in function arguments")
  if (is.null(met)) stop("weights must have `metric` parameter or must be specified in function arguments")
  if (is.null(est)) stop("weights must have `estimand` parameter or must be specified in function arguments")
  
  # met <- "Lp"
  prep.data <- prep_data(data, ...)
  
  z <- prep.data$z
  y <- prep.data$df$y
  y0 <- y[z == 0]
  y1 <- y[z == 1]
  
  x <- as.matrix(prep.data$df[,-which(colnames(prep.data$df) == "y")])
  x0 <- x[z == 0,]
  x1 <- x[z == 1,]
  
  n <- length(z)
  n1 <- sum(z)
  n0 <- n - n1
  d <- ncol(x)
  
  y_out <- list(control = rep(NA,n),
                treated = rep(NA,n),
                observed.treatment = z)
  y_out$control[z==0] <- y0
  y_out$treated[z==1] <- y1
  
  if(is.null(gamma)) {
    if(is.null(dots$cost)) {
      cost <- causalOT::cost_fun(x, z,
                                 power = pow,
                                 metric = met, rkhs.args = rkhs.args,
                                 estimand = est)
    } else {
      cost <- dots$cost
    }
    
    if(est == "ATE" ) {
      # cost <- list(causalOT::cost_fun(x[z==0,], x,
      #                            power = pow,
      #                            metric = met, rkhs.args = rkhs.args,
      #                            estimand = est),
      #              causalOT::cost_fun(x[z==1,], x,
      #                                         power = pow,
      #                                         metric = met, rkhs.args = rkhs.args,
      #                                         estimand = est
      #              )
      # )
      if (is.list(weight$gamma)) {
        gam0 <- weight$gamma[[1]]
        gam1 <- weight$gamma[[2]]
      } else {
        gam0 <- gam1 <- NULL
      }
      
      gamma <- list(calc_gamma(weights = list(w0 = weight$w0, w1 = rep(1/n,n),
                                              gamma = gam0), 
                               cost = cost[[1]], p = pow, ...),
                    calc_gamma(weights = list(w0 = weight$w1, w1 = rep(1/n,n),
                                              gamma = gam1), 
                               cost = cost[[2]], p = pow, ...))
    } else if(est == "ATT" | est == "ATC") {
      
      gamma <- calc_gamma(weight, cost = cost, p = pow, ...)
    }
    
  }
  # remove any NA values
  y0[weight$w0 == 0] <- 0
  y1[weight$w1 == 1] <- 0
  
  
  if(est == "ATC" | est == "ATT") {
    be_args <- list(gamma = gamma,x0 = x0,
                    x1 = x1, y0 = y0, y1 = y1, estimand = est,
                    metric = met,power = pow, ...)
    be_args <- be_args[!duplicated(names(be_args))]
    be_n  <- lapply(names(be_args), as.name)
    names(be_n) <- names(be_args)
    f.call <- as.call(c(list(as.name("barycenter_estimation")),  be_n))
    
    bproj <- eval(f.call, envir = be_args)
    
    y_out$control[z==1] <- bproj$y0
    y_out$treated[z==0] <- bproj$y1
  } else if (est == "ATE") {
    bproj <- list(y0 = NULL, y1 = NULL)
    be_args <- list(gamma = gamma[[1]], x0 = x0,
                    x1 = x, y0 = y0, y1 = y, estimand = "ATT",
                    metric = met,power = pow, ...)
    be_args <- be_args[!duplicated(names(be_args))]
    be_n  <- lapply(names(be_args), as.name)
    names(be_n) <- names(be_args)
    f.call <- as.call(c(list(as.name("barycenter_estimation")),  be_n))
    
    bproj$y0 <- eval(f.call, envir = be_args)$y0
    
    be_args$gamma <- gamma[[2]]
    be_args$x0 <- x1
    be_args$y0 <- y1
    
    bproj$y1  <- eval(f.call, envir = be_args)$y0
    
    y_out$control <- bproj$y0
    y_out$treated <- bproj$y1
  }
  
  
  
  return(y_out)
  
}


# function to actually estimate barycenters after checks above
barycenter_estimation <- function(gamma,x0,x1,y0,y1,
                                  estimand = c("ATT","ATC"),
                                  metric = c("mahalanobis","sdLp","Lp", "RKHS"),
                                  power = 2, ...) {
  
  stopifnot(power >= 1)
  estimand <- match.arg(estimand)
  metric <- match.arg(metric)
  if (metric == "RKHS")  metric <- "Lp"
  
  # check if all outcomes are the same
  
  if (estimand == "ATT" | estimand == "ATE") {
    if ( length(unique(y0)) == 1)  {
      y_out <- list(y0 = rep(NA_real_, length(y1)),
                    y1 = rep(NA_real_, length(y0)))
      y_out$y0 <- y0[1]
      return(y_out)
    }
  }
  if (estimand == "ATC" | estimand == "ATE") {
    if (length(unique(y1)) == 1) {
      y_out <- list(y0 = rep(NA_real_, length(y1)),
                    y1 = rep(NA_real_, length(y0)))
      y_out$y1 <- y1[1]
      return(y_out)
    }
  }
  
  if (power == 2) {
    return(barycenter_pow2(gamma,x0,x1,y0,y1,estimand,metric))
  } else if (power == 1) {
    return(barycenter_pow1(gamma,x0,x1,y0,y1,estimand,metric))
  }
  n0 <- length(y0)
  n1 <- length(y1)
  
  data <- list()
  data$p <- power
  
  if (metric == "Lp") {
    y0t <- as.matrix(y0)
    y1t <- as.matrix(y1)
  } else if (metric == "mahalanobis") {
    d   <- ncol(x0)
    d_0 <- cbind(x0,y0)
    d_1 <- cbind(x1,y1)
    
    cov <- 0.5*(cov(d_1) + cov(d_0))
    U    <- sqrt_mat(cov)
    U_inv <- solve(U)
    
    y0t <- d_0 %*% U_inv

    y1t <- d_1 %*% U_inv
  } else if (metric == "sdLp") {
    s0 <- sd(y0)
    s1 <- sd(y1)
    y0t <- as.matrix(y0)/s0
    y1t <- as.matrix(y1)/s1
  }
  
  if (estimand == "ATT") {
    data$N <- n0
    data$M <- n1
    data$y <- y0t
    data$D <- ncol(y0t)
    data$gamma <- apply(gamma,2,renormalize)
  } else if (estimand == "ATC") {
    data$M <- n0
    data$N <- n1
    data$y <- y1t
    data$D <- ncol(y1t)
    data$gamma <- apply(gamma,1,renormalize)
  }
  
  dots <- list(...)
  iter <- dots$iter
  precision <- dots$tol
  method <- dots$method
  
  if(is.null(iter)) iter <- 2000L
  if(is.null(precision)) precision <- as.double(1e-16)
  if(is.null(method)) method <- "optim"
  method <- match.arg(method, c("optim","rstan"))
  
  if(method == "optim"){
    data$y <- t(data$y)
    lp_loss <- function(x,y,p,w,n,d) {
      z<- matrix(x,d,n)
      return(sum(w * sapply(1:n, function(i) colSums(abs(z[,i] - y)^p))))
    }
    lp_grad <- function(x,y,p,w,n,d) {
      z <- matrix(x, d,n)
      return(sapply(1:n, function(i) p*(sign(z[,i] - y)*abs(z[,i] - y)^(p-1))%*% w[,i]))
    }
    res <- optim(par = c(x=rep(0,data$M* data$D)), fn = lp_loss, gr = lp_grad,
                 y = data$y, p = power, w = data$gamma, n = data$M, d = data$D,
                 control = list(maxit = iter, reltol = precision), method = "BFGS")
    iter.used <- res$counts[1] 
    convergence <- res$convergence
    output.message <- res$message
    
    zmat <- matrix(res$par, nrow = data$M, ncol = data$D, byrow = TRUE)
  } else if (method == "rstan") {
    arguments <- list(
      data = data,
      iter = iter,
      check_data = FALSE,
      ...,
      as_vector = FALSE
    )
    formals.stan <- c("iter", "save_iterations", 
                      "refresh", "init_alpha", "tol_obj", "tol_grad", "tol_param", 
                      "tol_rel_obj", "tol_rel_grad", "history_size",
                      "object", "data", "seed", "init", "check_data",
                      "sample_file", "algorithm", "verbose", "hessian", 
                      "as_vector", "draws", "constrained", "importance_resampling")
    arguments <- arguments[!duplicated(names(arguments))]
    arguments <- arguments[names(arguments) %in% formals.stan]
    if(is.null(arguments$object)) {
      if(metric == "Lp") {
        arguments$object <- stanbuilder("barycenter_projection")
      } else if (metric == "mahalanobis") {
        arguments$object <- stanbuilder("barycenter_projection")
      }
    }
    
    argn <- lapply(names(arguments), as.name)
    names(argn) <- names(arguments)
    f.call <- as.call(c(list(call("::", as.name("rstan"), 
                                  as.name("optimizing"))), argn))
    res <- eval(f.call, envir = arguments)
    zmat <- res$par$z
    convergence <- res$return_code
    output.message <- NULL
    iter.used <- 0
  }
  
  if(convergence != 0) warning("Convergence code ", convergence,"in method ",method,".Algorithm not converged for Lp minimization in barycenter projection. ", 
                               output.message)
  
  if(iter.used >= iter ) warning("Maximum number of iterations hit")
  y_out <- list(y0 = rep(NA_real_, n1),
                y1 = rep(NA_real_, n0))
  
  if (metric == "mahalanobis") {
    
    if (estimand == "ATT") {
      y_out$y0 <- c((zmat %*% U)[,d+1])
    } else if (estimand == "ATC") {
      y_out$y1 <- c((zmat %*% U)[,d+1])
    }
  } else if (metric == "Lp") {
    if (estimand == "ATT") {
      y_out$y0 <- c(zmat)
    } else if (estimand == "ATC") {
      y_out$y1 <- c(zmat)
    }
  } else if (metric == "sdLp") {
    if (estimand == "ATT") {
      y_out$y0 <- c(zmat) * s0
    } else if (estimand == "ATC") {
      y_out$y1 <- c(zmat) * s1
    }
  }
  return(y_out)
}

# barycenters if power = 2
barycenter_pow2 <- function(gamma,x0,x1,y0,y1,estimand,metric) {
  y_out <- list(y0 = rep(NA, length(y1)),
               y1 = rep(NA, length(y0)))
  if (metric != "Lp") warning("With p = 2 the mahalanobis, sdLp, and Lp metrics give the same answer.")
  # metric <- "Lp"
  if(metric == "Lp" | metric == "sdLp") {
    if(estimand == "ATT") {
      y_out$y0 <-  c(crossprod(gamma,y0) * 1/colSums(gamma))
    } else if (estimand == "ATC") {
      y_out$y1 <- c((gamma %*% y1) * 1/rowSums(gamma))
    }
  } else if (metric == "mahalanobis") {
    n0 <- length(y0)
    n1 <- length(y1)
    d <- ncol(x0)
    stopifnot(ncol(x1)==ncol(x0))
    d_0 <- cbind(x0,y0)
    d_1 <- cbind(x1,y1)
    cov <- 0.5*(cov(d_1) + cov(d_0))
    # U    <- chol(cov)
    # U_inv <- solve(U)
    eigen.vals <- eigen(cov, symmetric = TRUE)
    U    <- tcrossprod(eigen.vals$vectors %*% diag(sqrt(abs(eigen.vals$values))), eigen.vals$vectors)
    U_inv <- tcrossprod(eigen.vals$vectors %*% diag(1/sqrt(abs(eigen.vals$values))), eigen.vals$vectors)
    # grand_mean <- colMeans(rbind(d_0,d_1))
    # mean_0 <- colMeans(rbind(d_0))
    # mean_1 <- colMeans(rbind(d_1))
    
    if(estimand == "ATT") {
      # dt_0 <- scale(d_0, scale = FALSE, center = mean_0) %*% U_inv
      dt_0 <- d_0 %*% U_inv
      dt_0_counterfactual <- crossprod(gamma,dt_0) * 1/colSums(gamma)
      # d_0_counterfactual <- dt_0_counterfactual %*% U + matrix(mean_0, n1,d+1)
      d_0_counterfactual <- dt_0_counterfactual %*% U
      y_0_counterfactual <- d_0_counterfactual[,d+1]
      y_out$y0 <- c(y_0_counterfactual)
    } else if (estimand == "ATC") {
      # mean_1 <- colMeans(d_1)
      # dt_1 <- scale(d_1, scale = FALSE, center = mean_1) %*% U_inv
      dt_1 <- d_1 %*% U_inv
      dt_1_counterfactual <- (gamma %*% dt_1) * 1/rowSums(gamma)
      # d_1_counterfactual <- dt_1_counterfactual %*% U + matrix(mean_1, n0,d+1)
      d_1_counterfactual <- dt_1_counterfactual %*% U
      y_1_counterfactual <- d_1_counterfactual[,d+1]
      y_out$y1 <- c(y_1_counterfactual)
    }
  }
  return(y_out)
}

# barycenters if power = 1
barycenter_pow1 <- function(gamma,x0,x1,y0,y1,estimand,metric) {
  n0 <- length(y0)
  n1 <- length(y1)
  y_out <- list(y0 = rep(NA, n1),
               y1 = rep(NA, n0))
  
  if(metric == "Lp") {
    if(estimand == "ATT") {
      y_out$y0 <- c(apply(gamma,2,function(w) matrixStats::weightedMedian(x=y0, w=w)))
    } else if (estimand == "ATC") {
      y_out$y1 <- c(apply(gamma,1,function(w) matrixStats::weightedMedian(x=y1, w=w)))
    }
  } else if (metric == "sdLp") {
    
    if (estimand == "ATT") {
      s0 <- sd(y0)
      y_out$y0 <- c(apply(gamma,2,function(w) matrixStats::weightedMedian(x = y0/s0, w = w))) * s0
    } else if (estimand == "ATC") {
      s1 <- sd(y1)
      y_out$y1 <- c(apply(gamma,1,function(w) matrixStats::weightedMedian(x = y1/s1, w = w))) * s1
    }
  } else if (metric == "mahalanobis") {
    d <- ncol(x0)
    d_0 <- cbind(x0,y0)
    d_1 <- cbind(x1,y1)
    cov <- 0.5*(cov(d_1) + cov(d_0))
    eigen.vals <- eigen(cov, symmetric = TRUE)
    U    <- tcrossprod(eigen.vals$vectors %*% diag(sqrt(abs(eigen.vals$values))), eigen.vals$vectors)
    U_inv <- tcrossprod(eigen.vals$vectors %*% diag(1/sqrt(abs(eigen.vals$values))), eigen.vals$vectors)
    
    if(estimand == "ATT") {
      # mean_0 <- colMeans(d_0)
      # dt_0 <- scale(d_0, scale = FALSE, center = FALSE) %*% U_inv
      dt_0 <- d_0 %*% U_inv
      dt_0_counterfactual <- apply(gamma,2,function(w) matrixStats::colWeightedMedians(x=dt_0, w=w))
      d_0_counterfactual <- crossprod(dt_0_counterfactual, U)# + matrix(mean_0, n1,d+1)
      y_0_counterfactual <- d_0_counterfactual[,d+1]
      y_out$y0 <- c(y_0_counterfactual)
    } else if (estimand == "ATC") {
      # mean_1 <- colMeans(d_1)
      # dt_1 <- scale(d_1, scale = FALSE, center = FALSE) %*% U_inv
      dt_1 <- d_1 %*% U_inv
      dt_1_counterfactual <- apply(gamma,1,function(w) matrixStats::colWeightedMedians(x=dt_1, w=w))
      d_1_counterfactual <- crossprod(dt_1_counterfactual, U)# + matrix(mean_1, n0,d+1)
      y_1_counterfactual <- d_1_counterfactual[,d+1]
      y_out$y1 <- c(y_1_counterfactual)
    }
  }
  return(y_out)
}

