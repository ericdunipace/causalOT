barycentric_projection <- function(data, weight, 
                                   ...) {
  
  stopifnot(inherits(weight, "causalWeights"))
  
  est <- weight$estimand
  pow <- weight$args$p
  met <- weight$args$metric # seems that Lp performs fine in every case. Which makes sense since it's a linear transform
  rkhs.args <- weight$args$rkhs.args
  gamma <- weight$gamma
  dots <- list(...)
  
  if(is.null(est)) {
    est <- match.arg(dots$estimand)
    warning("Weights didn't have estimand flag. Make sure weights aligned with desired estimand")
  }
  if(is.null(met)) {
    met <- match.arg(dots$metric)
  }
  if(is.null(pow)) {
    pow <- as.double(dots$p)
  }
  if(is.null(rkhs.args)) {
    rkhs.args <- as.double(dots$rkhs.args)
  }
  if(is.null(pow)) stop("weights must have `p` parameter (power) or must be specified in function arguments")
  if(is.null(met)) stop("weights must have `metric` parameter (power) or must be specified in function arguments")
  met <- "Lp"
  prep.data <- prep_data(data,...)
  
  z <- prep.data$z
  y <- prep.data$df$y
  y0 <- y[z==0]
  y1 <- y[z==1]
  
  x <- as.matrix(prep.data$df[,-which(colnames(prep.data$df)=="y")])
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  
  n <- length(z)
  n1 <- sum(z)
  n0 <- n - n1
  d <- ncol(x)
  
  y_out <- list(control = rep(NA,n),
                treated = rep(NA,n))
  y_out$control[z==0] <- y0
  y_out$treated[z==1] <- y1
  
  if(is.null(gamma)) {
    cost <- causalOT::cost_fun(x[z==0,], x[z==1,],
                               power = pow,
                               metric = met, rkhs.args = rkhs.args,
                               estimand = est
                               )
    if(estimand == "ATE" ) {
      gamma <- list(calc_gamma(weight, cost = cost[[1]], p = pow),
                    calc_gamma(weight, cost = cost[[2]], p = pow))
    } else if(estimand == "ATT" | estimand == "ATC") {
      gamma <- calc_gamma(weight, cost = cost, p = pow)
    }
    
  }
  if(estimand == "ATC" | estimand == "ATT") {
    bproj <- barycenter_estimation(gamma,x0,x1,y0,y1,est,met,pow,...)
    y_out$control[z==1] <- bproj$y0
    y_out$treated[z==0] <- bproj$y1
  } else if (estimand == "ATE") {
    bproj <- list(y0 = barycenter_estimation(gamma[[1]],x0,x,y0,y,"ATT",met,pow,...)$y0,
                  y1 = barycenter_estimation(gamma[[1]],x1,x,y1,y,"ATT",met,pow,...)$y0)
    y_out$control <- bproj$y0
    y_out$treated <- bproj$y1
  }
  
  
  
  return(y_out)
  
}


barycenter_estimation <- function(gamma,x0,x1,y0,y1,
                                  estimand = c("ATT","ATC"),
                                  metric = c("Lp", "mahalanobis"),
                                  power = 2, ...) {
  
  stopifnot(power >=1)
  estimand <- match.arg(estimand)
  metric <- match.arg(metric)
  # metric <- "Lp"
  
  if(power == 2) {
    return(barycenter_pow2(gamma,x0,x1,y0,y1,estimand,metric))
  } else if (power == 1) {
    return(barycenter_pow1(gamma,x0,x1,y0,y1,estimand,metric))
  }
  n0 <- length(y0)
  n1 <- length(y1)
  
  data <- list()
  data$p <- power
  
  if(metric == "Lp") {
    y0t <- y0
    y1t <- y1
  } else if (metric == "mahalanobis") {
    d_0 <- cbind(x0,y0)
    d_1 <- cbind(x1,y1)
    
    cov <- 0.5*(cov(d_1) + cov(d_0))
    U    <- chol(cov)
    U_inv <- solve(U)
    
    y0t <- d_0 %*% U_inv

    y1t <- d_1 %*% U_inv
  }
  
  if(estimand == "ATT") {
    data$N <- n0
    data$M <- n1
    data$y <- y0t
    data$gamma <- gamma
  } else if (estimand == "ATC"){
    data$N <- n0
    data$M <- n1
    data$y <- y1t
    data$gamma <- t(gamma)
  }
  arguments <- list(
    data = data,
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
      arguments$object <- stanmodels$barycenter_projection_
    } else if (metric == "mahalanobis") {
      arguments$object <- stanmodels$barycenter_projection_mahalanobis_
    }
  }
  
  argn <- lapply(names(arguments), as.name)
  names(argn) <- names(arguments)
  f.call <- as.call(c(list(call("::", as.name("rstan"), 
                                as.name("optimizing"))), argn))
  res <- eval(f.call, envir = arguments)
  
  y_out <- list(y0 = rep(NA_real_, n1),
                y1 = rep(NA_real_, n0))
  
  if (metric == "mahalanobis") {
    if(estimand == "ATT") {
      y_out$y0 <- c((res$par$z %*% U)[,d+1])
    } else if (estimand == "ATC") {
      y_out$y1 <- c((res$par$z %*% U)[,d+1])
    }
  } else if (metric == "Lp") {
    if(estimand == "ATT") {
      y_out$y0 <- c(res$par$z)
    } else if (estimand == "ATC") {
      y_out$y1 <- c(res$par$z)
    }
  }
  return(y_out)
}

barycenter_pow2 <- function(gamma,x0,x1,y0,y1,estimand,metric) {
  y_out <- list(y0 = rep(NA, length(y1)),
               y1 = rep(NA, length(y0)))
  if(metric != "Lp") warning("With p = 2 the mahalanobis and Lp metrics give the same answer.")
  metric <- "Lp"
  if(metric == "Lp") {
    if(estimand == "ATT") {
      y_out$y0 <-  c(crossprod(gamma,y0) * 1/colSums(gamma))
    } else if (estimand == "ATC") {
      y_out$y1 <- c((gamma %*% y1) * 1/rowSums(gamma))
    }
  } else if (metric == "mahalanobis") {
    d <- ncol(x0)
    stopifnot(ncol(x1)==ncol(x0))
    d_0 <- cbind(x0,y0)
    d_1 <- cbind(x1,y1)
    cov <- 0.5*(cov(d_1) + cov(d_0))
    U    <- chol(cov)
    U_inv <- solve(U)
    grand_mean <- colMeans(rbind(d_0,d_1))
    mean_0 <- colMeans(rbind(d_0))
    mean_1 <- colMeans(rbind(d_1))
    
    if(estimand == "ATT") {
      dt_0 <- scale(d_0, scale = FALSE, center = mean_0) %*% U_inv
      dt_0_counterfactual <- crossprod(gamma,dt_0) * 1/colSums(gamma)
      d_0_counterfactual <- dt_0_counterfactual %*% U + matrix(mean_0, n1,d+1)
      y_0_counterfactual <- d_0_counterfactual[,d+1]
      y_out$y0 <- c(y_0_counterfactual)
    } else if (estimand == "ATC") {
      # mean_1 <- colMeans(d_1)
      dt_1 <- scale(d_1, scale = FALSE, center = mean_1) %*% U_inv
      dt_1_counterfactual <- (gamma %*% dt_1) * 1/rowSums(gamma)
      d_1_counterfactual <- dt_1_counterfactual %*% U + matrix(mean_1, n0,d+1)
      y_1_counterfactual <- d_1_counterfactual[,d+1]
      y_out$y1 <- c(y_1_counterfactual)
    }
  }
  return(y_out)
}

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
  } else if (metric == "mahalanobis") {
    d <- ncol(x0)
    d_0 <- cbind(x0,y0)
    d_1 <- cbind(x1,y1)
    cov <- 0.5*(cov(d_1) + cov(d_0))
    U    <- chol(cov)
    U_inv <- solve(U)
    
    if(estimand == "ATT") {
      # mean_0 <- colMeans(d_0)
      dt_0 <- scale(d_0, scale = FALSE, center = FALSE) %*% U_inv
      dt_0_counterfactual <- apply(gamma,2,function(w) matrixStats::colWeightedMedians(x=dt_0, w=w))
      d_0_counterfactual <- crossprod(dt_0_counterfactual, U)# + matrix(mean_0, n1,d+1)
      y_0_counterfactual <- d_0_counterfactual[,d+1]
      y_out$y0 <- c(y_0_counterfactual)
    } else if (estimand == "ATC") {
      mean_1 <- colMeans(d_1)
      dt_1 <- scale(d_1, scale = FALSE, center = FALSE) %*% U_inv
      dt_1_counterfactual <- apply(gamma,1,function(w) matrixStats::colWeightedMedians(x=dt_1, w=w))
      d_1_counterfactual <- crossprod(dt_1_counterfactual, U)# + matrix(mean_1, n0,d+1)
      y_1_counterfactual <- d_1_counterfactual[,d+1]
      y_out$y1 <- c(y_1_counterfactual)
    }
  }
  return(y_out)
}

