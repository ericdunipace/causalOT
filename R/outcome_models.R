# calculates predictions from gaussian process
gp_pred <- function(formula = NULL, data, weights=NULL,
                    param, estimand = c("ATE","ATT","ATC","cATE"), ...) {
  
  test_pos_def_inv <- function(x,y) {
    e <- eigen(x, symmetric = TRUE)
    if (any(e$values <= 0)) {
      min.e <- min(e$values)
      e$values <- e$values - min.e
    }
    p <- length(e$values)
    return(e$vectors %*% diag(1 / e$values, nrow = p, ncol = p) %*% 
             crossprod(e$vectors, y))
  }
  # w0 <- weights$w0
  # w1 <- weights$w1
  estimand <- match.arg(estimand)
  prep.data <- prep_data(data,...)
  
  z <- prep.data$z
  y <- prep.data$df$y
  x <- as.matrix(prep.data$df[,-which(colnames(prep.data$df) == "y")])
  
  n <- length(z)
  n1 <- sum(z)
  n0 <- n - n1
  
  if(is.null(param)) {
    param <- RKHS_param_opt(x,y,z,...)
  }
  
  
  if(estimand == "cATE"){
    att <- gp_pred(formula, data, weights, param, estimand = "ATT",...)
    atc <- gp_pred(formula, data, weights, param, estimand = "ATC",...)
    return(weighted.mean(c(att,atc), w = c(n1,n0)))
  }
  
  if(param$is.standardized) {
    m_y0   <- mean(y[z == 0])
    m_y1   <- mean(y[z == 1])
    sd_y0   <- sd(y[z == 0])
    sd_y1   <- sd(y[z == 1])
    y[z == 0] <- c(scale(y[z == 0]))
    y[z == 1] <- c(scale(y[z == 1]))
  } else {
    m_y0   <- 0
    m_y1   <- 0
    sd_y0   <- 1
    sd_y1   <- 1
  }
  
  # theta <- kernel_param_check(param$theta)
  # gamma <- kernel_param_check(param$gamma)
  # sigma2 <- as.double(param$sigma_2)
  # p <- as.double(param$p)
  # kernel <-  param$kernel
  # is.dose <- param$is.dose
  # 
  # calc_covariance <- isTRUE(param$metric == "mahalanobis")
  
  
  # Rcpp::List kernel_calc_pred_(const Rcpp::NumericMatrix & X_,  //confounders
  #                              const Rcpp::NumericMatrix & X_test_, //test points
  #                              const Rcpp::IntegerVector & z,  //tx, a vector but easier if matrix
  #                              const double p,
  #                              const Rcpp::NumericVector  & theta_,
  #                              const Rcpp::NumericVector & gamma_,
  #                              const Rcpp::NumericVector & sigma_2_,
  #                              const std::string & kernel_,
  #                              const bool calc_covariance,
  #                              const std::string & estimand)
  # Kernel_full <- kernel_calc_pred_(X_=as.matrix(x), 
  #                                  X_test_ = as.matrix(x),
  #                                  z = as.integer(z), p = p, 
  #                                  theta_ = theta, 
  #                                  gamma_ = gamma, 
  #                                  sigma_2_ = sigma2,
  #                                  calc_covariance = calc_covariance, 
  #                                  kernel = as.character(kernel), 
  #                                  estimand = as.character(estimand))
  # pred0 <- if(estimand == "ATE" | estimand == "ATT") {
  #   crossprod(Kernel_full[[1]]$cross, 
  #             solve(Kernel_full[[1]]$cov, y[z == 0]))
  # } else if (estimand == "ATC") {
  #   y[sel]
  # }
  # pred1 <- if(estimand == "ATE" | estimand == "ATC") {
  #   crossprod(Kernel_full[[2]]$cross, 
  #             solve(Kernel_full[[2]]$cov, y[z == 1]))
  # } else if (estimand == "ATT") {
  #   y[sel]
  # }
  if(param$metric == "mahalanobis") {
    cov <- cov(x)
    A <- scale(x, center = TRUE, scale = FALSE) %*% solve(chol(cov))
    
  } else {
    A <- x
  }
  A0 <- A[z == 0,]
  A1 <- A[z == 1,]
  
  if(param$kernel == "polynomial") {
    kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + 
      param$gamma[1] * (1+ param$theta[1] * tcrossprod(A0))^param$p
    kernel_cross0 <- param$gamma[1] * (1+ param$theta[1] * tcrossprod(A0,A))^param$p
    
    kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + 
      param$gamma[2] * (1+ param$theta[2] * tcrossprod(A1))^param$p
    kernel_cross1 <- param$gamma[2] * (1+ param$theta[2] * tcrossprod(A1,A))^param$p
  } else if (param$kernel == "RBF") {
    if(estimand == "ATE" | estimand == "ATT") {
      kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + 
      param$gamma[1] * exp(-0.5 * param$theta[1] * 
                             cost_calc_lp(A0,A0,ground_p = 2, direction = "rowwise" )^2)
    kernel_cross0 <- param$gamma[1] * exp(-0.5 *  param$theta[1] * 
                                            cost_calc_lp(A0,A,ground_p = 2, direction = "rowwise" )^2)
    }
    if(estimand == "ATE" | estimand == "ATC") {
      kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + 
        param$gamma[2] * exp(-0.5 * param$theta[2] * 
                               cost_calc_lp(A1,A1,ground_p = 2, direction = "rowwise" )^2)
      kernel_cross1 <- param$gamma[2] * exp(-0.5 *  param$theta[2] * 
                                              cost_calc_lp(A1,A,ground_p = 2, direction = "rowwise" )^2)
    }
  } else if (param$kernel == "linear") {
    kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + tcrossprod(A0)
    kernel_cross0 <- tcrossprod(A0,A)
    
    kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + tcrossprod(A1)
    kernel_cross1 <- tcrossprod(A1,A)
  }
  
  tau <- if(estimand == "ATE") {
    pred0 <- crossprod(kernel_cross0, test_pos_def_inv(kernel_cov0, y[z == 0])) * sd_y0 + m_y0
    pred1 <- crossprod(kernel_cross1, test_pos_def_inv(kernel_cov1, y[z == 1])) * sd_y1 + m_y1
    
    mean(pred1 - pred0)
  } else if (estimand == "ATT") {
    pred0 <- crossprod(kernel_cross0, test_pos_def_inv(kernel_cov0, y[z == 0])) * sd_y0 + m_y0
    
    mean((y[z == 1] * sd_y1 + m_y1) - pred0[z == 1])
  } else if (estimand == "ATC") {
    pred1 <- crossprod(kernel_cross1, test_pos_def_inv(kernel_cov1, y[z == 1])) * sd_y1 + m_y1
    
    mean(pred1[z == 0] - (y[z == 0] * sd_y0 + m_y0))
  }
  
  return(tau)
}

# gp model predict
predict.gp <- function(object, newdata = NULL, na.action = na.pass, ...) {
  test_pos_def_inv <- function(x,y) {
    e <- eigen(x, symmetric = TRUE)
    if (any(e$values <= 0)) {
      min.e <- min(e$values)
      e$values <- e$values - min.e
    }
    p <- length(e$values)
    return(e$vectors %*% diag(1 / e$values, nrow = p, ncol = p) %*% 
             crossprod(e$vectors, y))
  }
  tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    X <- object$x
  } else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action, 
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  }
  param <- object$param
  n <- nrow(object$x)
  
      kernel_cov <- diag(param$sigma_2,n,n) + 
        param$gamma * exp(-0.5 * param$theta * 
                               cost_calc_lp(object$x, object$x, ground_p = 2, direction = "rowwise" )^2)
      kernel_cross <- param$gamma * exp(-0.5 *  param$theta * 
                                              cost_calc_lp(object$x, X, ground_p = 2, direction = "rowwise" )^2)
  
    pred <- crossprod(kernel_cross, test_pos_def_inv(kernel_cov, object$y))
    if(object$is.standardized) pred <- pred * object$sd.y + object$m.y
  
  return(pred)
}

# get reisdual from gp object
residuals.gp <- function(object, pred = NULL, ...) {
  
  if(!is.null(pred)) {
    y <- object$y * object$sd.y + object$m.y
    
    stopifnot(length(y) == length(pred))
    return(y - pred)
  } else {
    return(object$residuals)
  }
  
}

#gp model fit
gp <- function(formula = NULL, data, weights = NULL, ...) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  
  x <- model.matrix(mt, mf, contrasts)
  
  dist_mat <- cost_calc_lp(x, x, ground_p = 2)^2
  
  y_std <- scale(y)
  arguments <- list(
    data = list(N = as.integer(length(y)),
                y = c(y_std),
                discrep = dist_mat,
                p = 0.0,
                kernel = 2L,
                a = 0.,
                b = 0.),
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
      arguments$object <- stanbuilder("gp_hyper")
  }
  if(is.null(arguments$algorithm)) {
    arguments$algorithm <- "Newton"
  }
  
  argn <- lapply(names(arguments), as.name)
  names(argn) <- names(arguments)
  f.call <- as.call(c(list(call("::", as.name("rstan"), 
                                as.name("optimizing"))), argn))
  
  param <- list()
  tune.fun <- function(x, l, u) {
    abs(0.01 - pgamma(l, exp(x[1]), exp(x[2]))) + abs(0.01 - pgamma(u, exp(x[1]), exp(x[2]), lower.tail = FALSE))
  }
  
  #gp
  l = sqrt(min(arguments$data$discrep[lower.tri(arguments$data$discrep)]))
  u = sqrt(max(arguments$data$discrep))
  tuned <- optim(par = runif(2), fn = tune.fun, l = l, u = u)
  if (tuned$convergence == 0) {
    arguments$data$a <- exp(tuned$par[1])
    arguments$data$b <- exp(tuned$par[2])
  }
  if (arguments$algorithm != "Newton") {
    res <- lapply(1:10, function(i) eval(f.call, envir = arguments))
    res <- res[[which.max(sapply(res, function(r) r$par$marg_lik))]]
  } else {
    res <- eval(f.call, envir = arguments)
  }
  if (!is.null(res$par$theta_half)) {
    param$theta <- c(1/res$par$theta_half^2)
  }
  if (!is.null(res$par$gamma_half)) param$gamma <- c(res$par$gamma_half^2)
  param$sigma_2   <- c(res$par$sigma_half^2)
  
  output <- list()
  class(output) <- "gp"
  
  output$param <- param
  output$residuals <- NA_real_
  output$fitted.values <- NA_real_
  output$df.residual <- length(y) - length(param) + 1
  output$xlevels <- .getXlevels(mt, mf)
  output$call <- cl
  output$terms <- mt
  output$model <- mf
  output$is.standardized <- TRUE
  output$sd.y <- sd(y)
  output$m.y  <- mean(y)
  output$x <- x
  output$y <- y_std
  
  output$fitted.values <- predict.gp(object = output)
  output$residuals <- y - output$fitted.values 
  output$marg_lik <- res$par$marg_lik
  
  return(output)
}

setClass("gp")
setMethod("predict", signature = c(object = "gp"),
          definition = predict.gp)
setMethod("residuals", signature = c(object = "gp"),
          definition = residuals.gp)
  
# use barycentric mapping
mapping <- function(data, z, weights, estimand, f1, f0, sw, ...) {
  n  <- length(z)
  n1 <- sum(z)
  n0 <- n - n1
  runMap <- TRUE # isTRUE(weights$method %in% ot.methods())
  
  if (runMap) {
    # browser()
    bal.cov    <- colnames(data)[colnames(data) != "y"]
    data$z     <- z
    data$f1    <- f1
    data$f0    <- f0
    bproj_y    <- barycentric_projection(data, weights, 
                                      treatment.indicator = "z", 
                                      outcome = "y",
                                      balance.covariates = bal.cov,
                                      estimand = estimand,
                                      ...)
    
    if (estimand == "ATE" | estimand == "ATC") {
      # the model projection for the controls and then predicted value of E(Y|X for controls)
      bproj_f1 <- barycentric_projection(data, weights, 
                                        treatment.indicator = "z", 
                                        outcome = "f1",
                                        balance.covariates = bal.cov,
                                        estimand = "ATC",
                                        ...)
      
      
    } else {
      bproj_f1 <- list(control = rep(NA_real_,n),
                       treated = rep(NA_real_,n))
    }
    if (estimand == "ATE" | estimand == "ATT") {
      bproj_f0 <- barycentric_projection(data, weights, 
                                           treatment.indicator = "z", 
                                           outcome = "f0",
                                           balance.covariates = bal.cov,
                                           estimand = "ATT",
                                           ...)
    } else {
      bproj_f0 <- list(control = rep(NA_real_,n),
                       treated = rep(NA_real_,n))
    }
    # y0         <- y1 <- rep(NA_real_, n)
    y1 <- (bproj_y$treated - ifelse(z == 1, bproj_f0$treated, bproj_f1$treated))
    y0 <- (bproj_y$control - ifelse(z == 1, bproj_f0$control, bproj_f1$control))
    
    keep <- !is.na(y1) & !is.na(y0)
    y1 <- y1[keep]
    y0 <- y0[keep]
    new_weights <- renormalize(sw$total[keep])
    
  } else {
    n          <- length(z)
    y0         <- y1 <- rep(NA_real_, n)
    y0[z == 0]   <- data$y[z == 0] - f1[z == 0]
    y1[z == 1]   <- data$y[z == 1] - f0[z == 1]
    if (estimand == "ATT") {
      y0[z == 1] <- (data$y[z == 0] -  f0[z == 0]) %*% weights$w0
      y0 <- y0[z == 1]
      y1 <- y1[z == 1]
      new_weights <- sw$b
    } else if (estimand == "ATC") {
      y1[z == 0] <- (data$y[z == 1] -  f1[z == 1]) %*% weights$w1
      y0 <- y0[z == 0]
      y1 <- y1[z == 0]
      new_weights <- sw$a
    } else if (estimand == "ATE" | estimand == "feasible") {
      y0[z == 1] <- (data$y[z == 0] -  f0[z == 0]) %*% weights$w0
      y1[z == 0] <- (data$y[z == 1] -  f1[z == 1]) %*% weights$w1
      # y0[z == 0][weights$w0 == 0] <- 
      # y1[z == 0][weights$w0 == 0] <- NA
      # y0[z == 1][weights$w1 == 0] <- 
      # y1[z == 1][weights$w1 == 0] <- NA
      
      # y0 <- y0[!is.na(y0)]
      # y1 <- y1[!is.na(y1)]
      new_weights <- sw$total
    }
    
  }
  # n1 <- length(y1)
  # n0 <- length(y0)
  # new_weights <- list(w0 = rep(1/n0, n0),
  #                     w1 = rep(1/n1, n1))
  # new_weights <- list(w0 = rep(1/n0, n0),
  #                     w1 = rep(1/n1, n1))
  return(list(y0 = y0,
              y1 = y1,
              weights = new_weights))
}

.outcome_calc_deprecated <- function(data, z, weights, formula, model.fun, matched, estimand) {
  
  w0 <- weights$w0
  w1 <- weights$w1
  t_ind <- z == 1
  c_ind <- z == 0
  
  fit_1 <- model.fun(formula$treated, data[t_ind,,drop=FALSE])
  fit_0 <- model.fun(formula$control, data[c_ind,,drop=FALSE])
  f_1   <- predict(fit_1, data)
  f_0   <- predict(fit_0, data)
  
  if (matched) {
    # if(is.null(weights$gamma)) {
    #   stop("Transport matrix must be specified for matched estimator")
    # }
    # rad   <- (2*z - 1)
    # maps <- mapping(data, z, weights, estimand, f1 = f_1, f0 = f_0, ...)
    
    e_1   <- (data$y - f_1)
    e_0   <- (data$y - f_0)
    
    # idx   <- which(weights$gamma !=0, arr.ind = TRUE)
    # 
    # gamma_vec <- c(weights$gamma[idx])
    
    e_1_t <- e_1[t_ind]
    e_1_c <- e_1[c_ind]
    e_0_t <- e_0[t_ind]
    e_0_c <- e_0[c_ind]
    
    tau_t <- 0
    tau_c <- 0
    if(estimand == "ATT" | estimand == "ATE" | estimand == "feasible") {
      tau_t <- c(mean(e_0_t) - sum(e_0_c*w0))
    }
    if(estimand == "ATC" | estimand == "ATE" | estimand == "feasible") {
      tau_c <- c(sum(e_1_t*w1) - mean(e_1_c))
    }
    if(estimand == "ATT"){
      tx_effect <- tau_t
    } else if (estimand == "ATC") {
      tx_effect <- tau_c
    } else if(estimand == "ATE" | estimand == "feasible") {
      if(estimand == "ATE") {
        t_w   <- sum(t_ind) #1/(sum(w1^2))
        c_w   <- sum(c_ind) #1/(sum(w0^2))
      } else if (estimand == "feasible") {
        t_w   <- 1/(sum(w1^2))
        c_w   <- 1/(sum(w0^2))
      }
      
      tx_effect <- ( t_w * tau_t + c_w * tau_c)/(c_w + t_w)
    }
    
  } else {
    mu_1  <- mean(f_1)
    mu_0  <- mean(f_0)
    e_1   <- fit_1$residuals
    e_0   <- fit_0$residuals
    
    y_1   <- mu_1 + e_1 %*% w1
    y_0   <- mu_0 + e_0 %*% w0
    
    tx_effect <- c(y_1 - y_0)
  }
  
  
  return(tx_effect)
}

# model without predictors
IDmodel <- function(formula, data, weights) {
  
  y <- model.response( model.frame(as.formula(formula), data), "numeric")
  m <- weighted.mean(x = y, w = weights)
  
  out <- list(fit = m,
       residuals = y - m)
  
  class(out) <- "IDmodel"
  
  return(out)
}

predict.IDmodel <- function(object, newdata, ...) {
  rep(object$fit, nrow(newdata))
}

setClass("IDmodel", slots = c(fit = "numeric", residuals = "numeric"))
setMethod("predict", signature = c(object = "IDmodel"), definition = predict.IDmodel)

# calculate treatment effects from the outcome
# also return means and predicted values for CI otherwise recalculate all of this stuff...
outcome_calc <- function(data, z, weights, formula, model.fun, matched, estimand,
                         sample_weight, ...) {
  
  w0 <- weights$w0
  w1 <- weights$w1
  t_ind <- z == 1
  c_ind <- z == 0
  
  n0 <- length(c_ind)
  n1 <- length(t_ind)
  
  environment(formula$treated) <- environment(formula$control) <- environment()
  
  f1w <- sample_weight$total
  f0w <- sample_weight$total
  
  nonmap_sw <- sample_weight$total
  if (estimand == "ATT") {
    f0w[t_ind] <- sample_weight$b
    f0w[c_ind] <- w0
    f0w <- renormalize(f0w)
    nonmap_sw[c_ind] <- 0
    nonmap_sw <- renormalize(nonmap_sw)
    
    fit_1 <- IDmodel(formula$treated, data[t_ind,,drop = FALSE],
                     weights = sample_weight$b)
    # fit_0 <- model.fun(formula$control, data[c_ind,,drop = FALSE], 
    #                    weights = sample_weight$a)
    fit_0 <- model.fun(formula$control, data[c_ind,,drop = FALSE], 
                       weights = w0)
    fit_0$weights <- w0
  } else if (estimand == "ATC") {
    # fit_1 <- model.fun(formula$treated, data[t_ind,,drop = FALSE], weights = sample_weight$b)
    fit_1 <- model.fun(formula$treated, data[t_ind,,drop = FALSE], weights = w1)
    fit_0 <- IDmodel(formula$control, data[c_ind,,drop = FALSE], weights = sample_weight$a)
    
    f1w[t_ind] <- w1
    f1w[c_ind] <- sample_weight$a
    f1w <- renormalize(f1w)
    
    nonmap_sw[t_ind] <- 0
    nonmap_sw <- renormalize(nonmap_sw)
    
    fit_1$weights <- w1
    
  } else if (estimand == "ATE") {
    
    # fit_1 <- model.fun(formula$treated, data[t_ind,,drop = FALSE], weights = sample_weight$b)
    # fit_0 <- model.fun(formula$control, data[c_ind,,drop = FALSE], weights = sample_weight$a)
    
    fit_1 <- model.fun(formula$treated, data[t_ind,,drop = FALSE], weights = w1)
    fit_0 <- model.fun(formula$control, data[c_ind,,drop = FALSE], weights = w0)
    
    fit_1$weights <- w1
    fit_0$weights <- w0
  }
  
  # fit_1 <- model.fun(formula$treated, data[t_ind,,drop = FALSE], weights = sample_weight$b)
  # fit_0 <- model.fun(formula$control, data[c_ind,,drop = FALSE], weights = sample_weight$a)
  f_1   <- predict(object = fit_1, newdata = data, weights = f1w)
  f_0   <- predict(object = fit_0, newdata = data, weights = f0w)
  
  if (matched) {

    maps <- mapping(data = data, z = z, weights = weights, 
                    estimand = estimand, f1 = f_1, f0 = f_0,
                    sw = sample_weight,
                    ...)
    # browser()
    # mw        <- switch(estimand,
    #                     "ATT" = sample_weight$b,
    #                     "ATC" = sample_weight$a,
    #                     "ATE" = sample_weight$total)
    #                     
    # tx_effect <- weighted.mean(x = maps$y1 - maps$y0, 
    #                            w = mw[mw != 0])
    y_1  <- weighted.mean(x = maps$y1, 
                          w = maps$weights)
    y_0  <- weighted.mean(x = maps$y0, 
                          w = maps$weights)
    tx_effect <- weighted.mean(x = maps$y1 - maps$y0, 
                               w = maps$weights)
      
  } else {
    mu_1  <- weighted.mean(f_1, nonmap_sw) 
    mu_0  <- weighted.mean(f_0, nonmap_sw)
    e_1   <- residuals(fit_1, 
                       pred = f_1[t_ind])
    e_0   <- residuals(fit_0, 
                       pred = f_0[c_ind])
    
    y_1   <- mu_1 + weighted.mean(x = e_1, w = w1)
    y_0   <- mu_0 + weighted.mean(x = e_0, w = w0)
    
    tx_effect <- c(y_1 - y_0)
  }
  
  
  # return(tx_effect)
  
  
  # get E(Y(z)|X)
  e_y1_x <- NA
  e_y0_x <- NA
  if (estimand == "ATT") {
    # set to NA if not using outcome model...
    if(formula$control != "y~1" && formula$control != "y ~ 1") e_y0_x <- f_0
  } else if (estimand == "ATC") {
    if(formula$treated != "y~1" && formula$treated != "y ~ 1") e_y1_x <- f_1
  } else if (estimand == "ATE") {
    # otherwise return predictions
    if(formula$control != "y~1" && formula$control != "y ~ 1") e_y0_x <- f_0
    if(formula$treated != "y~1" && formula$treated != "y ~ 1") e_y1_x <- f_1
  } else {
    stop("estimand not recognized when saving predicted values")
  }
  
  # returns estimate and quantities needed for asymptotic variance
  return(list(tx_effect = tx_effect,
              outcome.model.fit = NULL, #no need to return it here
              E_Y1 = y_1, #estimates of means
              E_Y0 = y_0,
              E_Y1_X = e_y1_x, #conditional mean estimates if specifying an outcome, otherwise should be NA
              E_Y0_X = e_y0_x)) #should also be a vector of values for each 
}

# calculate treatment effects using a model only
outcome_calc_model <- function(data, z, weights, formula, model.fun, 
                               matched, estimand, 
                               ...) {
  w0 <- weights$w0
  w1 <- weights$w1

  n <- length(z)
  w <- rep(NA,n)
  w[z == 1] <- w1
  w[z == 0] <- w0
  

  environment(formula) <- environment()
  
  fit <- model.fun(formula, 
                   data =  cbind(data, 
                                 z = z), 
                   weights = w, 
                   ...)

  tx_effect <- coef(fit, tx.name = "z", estimand = estimand)["z"]

  # return(tx_effect)
  e_y1_x <- predict(fit, newdata = cbind(data, z = 1))
  e_y0_x <- predict(fit, newdata = cbind(data, z = 0))
  
  y_1 <- mean(e_y1_x)
  y_0 <- mean(e_y0_x)
  
  stopifnot(isTRUE(all.equal(y_1 - y_0,  tx_effect, check.attributes = FALSE)))
  
  return(list(tx_effect = tx_effect,
              outcome.model.fit = fit, #return the model so we can use HW se if possible
              E_Y1 = y_1, #estimates of means
              E_Y0 = y_0,
              E_Y1_X = e_y1_x, #conditional mean estimates from model
              E_Y0_X = e_y0_x))
  
}

# set-up formula for treatment effects
calc_form <- function(formula, doubly.robust, target, split.model) {
  if (!is.null(formula)) {
    if ( isTRUE(split.model) ) {
      if(!(isTRUE(all(names(formula) %in% c("treated","control"))) &&
                  isTRUE(!is.null(names(formula))))) {
        stop(
          "`formula` must be specified as list with slots 'treated' and 'control'"
        )
      }
      formula$treated <- as.formula(formula$treated)
      formula$control <- as.formula(formula$control)
      
    } else {
      formula <- as.formula(formula)
    }
    return(formula)
  }
  form_mod <- formula(y~.)
  form_obs <- formula(y~1)
  
  if (isTRUE(split.model) ) {
    formula$treated <- form_obs
    formula$control <- form_obs
    
    if (doubly.robust) {
      if (target == "ATT") {
        formula$control <- form_mod
      } else if (target == "ATC") {
        formula$treated <- form_mod
      } else if (target == "ATE" | target == "feasible") {
        formula$treated <- form_mod
        formula$control <- form_mod
      }
    }
  } else {
    formula <- formula(y ~ z + .)
  } 
  
  return(formula)
}

# setup model for treatment effects
calc_model <- function(model.fun) {
  if (!is.null(model.fun)) {
    if (is.character(model.fun)) {
      model.fun <- match.fun(model.fun)
    }
    stopifnot(is.function(model.fun))
  } else {
    model.fun <- lm
  }
  return(model.fun)
}

# setup hajek weights
calc_hajek <- function(weights, target, hajek) {
  if (hajek) {
    weights$w0 <- renormalize(weights$w0)
    weights$w1 <- renormalize(weights$w1)
  }
  if (inherits(weights, "causalWeights")) {
    estimand <- switch(weights$estimand,
                       "ATT" = "ATT",
                       "ATC" = "ATC",
                       "ATE" = "ATE",
                       "cATE" = "ATE",
                       "feasible" = "feasible")
    if (estimand != target) stop("Weights do not estimate the target estimand")
  }
  return(weights)
}


#' causalEffect class
#'
#' @slot estimate The estimated treatment effect.
#' @slot data The original data as a `data.frame`.
#' @slot model The function used as the outcome model.
#' @slot formula The formula for the outcome model.
#' @slot weights The weights as an object of class [causalWeights][causalOT::causalWeights-class]
#' @slot estimand A character denoting the estimand targeted by the weights. One of "ATT","ATC", or "ATE". 
#' @slot variance.components Objects for the asymptotic variance calculation designed so expensive models
#' don't have to be re-fit.
#' @slot options A list with the arguments from the [estimate_effect][causalOT::estimate_effect()] function. See details.
#' @slot call The call from the [estimate_effect()][causalOT::estimate_effect] function.
#' 
#' @details The `variance.components` slot is a list with slots
#' \itemize{
#' \item `E_Y1`: The mean if the target population had all been treated.
#' \item `E_Y0`: The mean if the target population had all received control
#' \item `E_Y1_X`: The predicted conditional mean if the target population had all been treated.
#' \item `E_Y0_X`: The predicted conditional mean if the target population had all received control.
#' }
#' Note that for "ATT" and "ATC" estimands, `E_Y1_X` or `E_Y0_X` will be NA, respectively.
#' 
#' Meanwhile, the `options` slot is a list with slots
#' \itemize{
#' \item `hajek`: Were weights normalized to sum to 1 (TRUE/FALSE)
#' \item `doubly.robust`: Was an augmented estimator used? (TRUE/FALSE)
#' \item `matched`: Wass barycentric projection estimator used? (TRUE/FALSE)
#' \item `split.model` Was the outcome model calculated separately in each
#' treatment group? (TRUE/FALSE)
#' \item `balance.covariates`: The covariates selected for balance or in the outcome model in slot `data`
#' \item `treatment.indicator`: The column that is the treatment indicator in slot `data`
#' \item `outcome`: The columns that is the outcome in  slot `data`
#' \item `addl.args`: Any additional arguments passed in the dots (`...`) 
#' of [estimate_effect()][causalOT::estimate_effect].
#' }
#' 
#' @docType class
#' 
#' @rdname causalEffect
#'
#' @export
setClass("causalEffect", slots = c(estimate = "numeric", 
                                   data = "data.frame",
                                   model = "function", 
                                   formula = "formula",
                                   weights = "causalWeights",
                                   estimand = "character",
                                   variance.components = "list",
                                   options = "list",
                                   call = "call"))


#' Estimate treatment effects
#'
#' @param data A `data.frame`, a `list`, or a [DataSim][causalOT::DataSim] object
#' @param formula the outcome model formula
#' @param weights An object of class [causalWeights][causalOT::causalWeights-class]
#' @param hajek Should the weights be normalized to sum to 1 (TRUE/FALSE)
#' @param doubly.robust Should an augmented estimator be used? (TRUE/FALSE)
#' @param matched Should a matched or barycentric project 
#' estimator be used? (TRUE/FALSE)
#' @param estimand Estimand to use. Should agree with estimand in the weights or can
#' be left blank. One of "ATT", "ATC", or "ATE".
#' @param model The outcome model as a character referring to a function or function
#' @param split.model Should the outcome model be calculated separately in each
#' treatment group? (TRUE/FALSE)
#' @param sample_weight The sample weights. Either NULL or an object of class
#' [sampleWeights][causalOT::sampleWeights-class]
#' @param ... Pass additional arguments to the outcome modeling functions like `lm`. Arguments "balance.covariates" and "treatment.indicator" must be provided if data is of class data.frame or matrix.
#'
#' @return an object of class [causalEffect][causalOT::causalEffect-class]
#' @export
#'
#' @examples
#' # set-up data
#' data <- Hainmueller$new()
#' data$gen_data()
#' 
#' # calculate quantities
#' weight <- calc_weight(data, method = "Logistic")
#' tx_eff <- estimate_effect(data = data, weights = weight)
#' 
#' # get estimate
#' print(tx_eff$estimate)
estimate_effect <- function(data, formula = NULL, weights, 
                                  hajek = TRUE, 
                                  doubly.robust = TRUE,
                                  matched = FALSE,
                                  estimand = c("ATT", "ATC", "ATE", "feasible"),
                                  model = NULL, 
                                  split.model = TRUE, 
                                  sample_weight = NULL,
                                  ...) {
  #get args
  # dots <- list(...)
  
  # hajek <- match.arg(hajek)
  dr <- isTRUE(doubly.robust[1])
  matched <- isTRUE(matched[1])
  hajek   <- isTRUE(hajek)
  split.model <- isTRUE(split.model)
  
  if (missing(estimand) || is.null(estimand)) {
    # if (!is.null(target)) estimand <- target
    if (inherits(weights, "causalWeights")) {
      estimand <- switch(weights$estimand,
                       "ATT" = "ATT",
                       "ATC" = "ATC",
                       "ATE" = "ATE",
                       "cATE" = "ATE",
                       "feasible" = "feasible")
    } else {
      stop("Estimand must be specified or weight object of class causalWeights given.")
    }
  # } else if (missing(estimand) && !missing(target)) {
  #   estimand <- target
  #   estimand <- switch(estimand,
  #                      "ATT" = "ATT",
  #                      "ATC" = "ATC",
  #                      "ATE" = "ATE",
  #                      "cATE" = "ATE",
  #                      "feasible" = "feasible")
  #   estimand <- match.arg(estimand)
  # } else if (missing(estimand)) {
  #   if (inherits(weights, "causalWeights")) {
  #     estimand <- switch(weights$estimand,
  #                        "ATT" = "ATT",
  #                        "ATC" = "ATC",
  #                        "ATE" = "ATE",
  #                        "cATE" = "ATE",
  #                        "feasible" = "feasible")
  #   }
  } else {
    estimand <- switch(estimand,
                       "ATT" = "ATT",
                       "ATC" = "ATC",
                       "ATE" = "ATE",
                       "cATE" = "ATE",
                       "feasible" = "feasible")
    estimand <- match.arg(estimand)
  }
  if(inherits(weights, "causalWeights")) {
    if(!(weights$method %in% c("Logistic","None","CBPS"))) {
      hajek <- TRUE
    }
  }
  
  #set up model structures
  #gets model function, or returns "lm" by default
  model.fun <- calc_model(model) 
  #gets specified formula, or returns "~1" (no formula) by default
  formula <- calc_form(formula, dr, estimand, split.model) 
  # changes weights to hajek weights (sum to 1) if desired 
  weights <- calc_hajek(weights, estimand, hajek)
  
  #setup data
  prep.data <- prep_data(data,...) #separates outcome, tx indicator, covariates
  sample_weight <- get_sample_weight(sample_weight, prep.data$z)
  
  #get estimate
  estimation.return <- if (split.model) {
    #estimates separate models on each treatment group
    outcome_calc(prep.data$df, prep.data$z, weights, 
                 formula, model.fun,
                 matched, estimand, sample_weight,
                 ...)
  } else {
    #estimates one model, such as a linear regression, on data jointly
    #i.e., all treatment groups thrown into the same model
    outcome_calc_model(prep.data$df, prep.data$z, weights, 
                       formula, model.fun,
                       matched, estimand, ...)
  }
  
  estimate <- estimation.return$tx_effect
  
  # get model fun if not null
  if (is.null(model)) {
    model.fun <- NULL
  }
  
  #save covariates + tx indicator
  output.df <- cbind(prep.data$df, z = prep.data$z)
  
  #set attributes for other functions such as ci
  attr(output.df, "balance.covariates") <- attributes(prep.data$df)$balance.covariates
  attr(output.df, "outcome") <- attributes(prep.data$df)$outcome
  attr(output.df, "treatment.indicator") <- "z"
  
  #setup output list
  output <- list(estimate = estimate,
                 data = output.df,
                 model = model.fun,
                 formula = formula,
                 weights = weights,
                 estimand = estimand,
                 outcome.model.fit = estimation.return$outcome.model.fit,
                 variance.components = list(E_Y1 = estimation.return$E_Y1,
                   E_Y0 = estimation.return$E_Y0,
                   E_Y1_X = estimation.return$E_Y1_X,
                   E_Y0_X = estimation.return$E_Y0_X),
                 options = list(hajek = hajek,
                                doubly.robust = doubly.robust,
                                matched = matched,
                                split.model = split.model,
                                balance.covariates = attr(output.df, "balance.covariates"),
                                treatment.indicator = attr(output.df, "treatment.indicator"),
                                outcome = attr(output.df, "outcome"),
                                addl.args = list(...)
                                ),
                 call = match.call())
  
  #set class
  class(output) <- "causalEffect"
  
  #return output list
  return(output)
}

#' Confidence Intervals for Causal Effects
#'
#' @param object An object of class [causalEffect][causalOT::causalEffect-class]
#' @param parm Unused. Included to match forms of other `confint` functions
#' @param level Confidence level. Should be between 0 and 1. Default is 0.95.
#' @param method How to calculate the confidence interval. Choices are "bootstrap" for
#' a bootstrap confidence interval, "asymptotic" for "asymptotic" confidence intervals, and
#' "jacknife" for jacknife confidence intervals. Default is "asymptotic" since it is faster.
#' @param ... Additional arguments if method is "bootstrap". Can include 
#' \itemize{
#' \item `n.boot`. How many bootstrap samples should be used. Default is 1000.
#' \item `boot.method`. One of "n-out-of-n" or "m-out-of-n". Optimal transport methods default to 
#' "m-out-of-n".
#' \item `verbose`. Should a progress bar be printed? (TRUE/FALSE) Defaults to FALSE.
#' }
#'
#' @return A list with slots "CI" giving the confidence bounds and "SD" giving estimates of the standard
#' error of the causal effects. If method is "bootstrap" and `boot.method` is "m-out-of-n", then there will
#' also be a slot named "unadjusted" giving the unadjusted confidence interval and standard error estimate
#' for reference.
#' @method confint causalEffect
#' @export
#'
#' @examples
#' # set-up data
#' set.seed(1234)
#' data <- Hainmueller$new()
#' data$gen_data()
#' 
#' # calculate quantities
#' weight <- calc_weight(data, method = "Logistic")
#' tx_eff <- estimate_effect(data = data, weights = weight)
#' 
#' # get asymptotic C.I.
#' confint(tx_eff, model = "lm", method = "asymptotic",
#'     formula = list(treated = "y ~ .", control = "y ~ ."))
confint.causalEffect <- function(object, parm, level = 0.95, 
                                 method = c("asymptotic", "bootstrap",  "jackknife"), 
                                 ...) {
  method <- match.arg(method)
  if (method == "bootstrap") {
    return(ci_boot_ce(object, parm, level, ...))
  } else if (method == "jackknife") {
    return(ci_jack_ce(object, parm, level, ...))
  } else if (method == "asymptotic") {
    return(ci_asympt(object, parm, level, ...))
  } else {
    stop("confidence interval method not currently implemented")
  }
  
  
}

# default vcovHC in sandwich doesn't give correct results for models with 0 weights for some observations!!!
# meatHC <- function (x, type = c("HC3", "const", "HC", "HC0", "HC1", "HC2",
#                                 "HC4", "HC4m", "HC5"), omega = NULL, ...)
# {
#   if (is.list(x) && !is.null(x$na.action))
#     class(x$na.action) <- "omit"
#   X <- model.matrix(x)
#   if (any(alias <- is.na(coef(x))))
#     X <- X[, !alias, drop = FALSE]
#   attr(X, "assign") <- NULL
#   n <- NROW(X)
#   diaghat <- try(hatvalues(x), silent = TRUE)
#   df <- n - NCOL(X)
#   ef <- sandwich::estfun(x, ...)
#   res <- rowMeans(ef/X, na.rm = TRUE)
#   all0 <- apply(abs(ef) < .Machine$double.eps, 1L, all)
#   res[all0] <- 0
#   if (any(all0) && substr(type, 1L, 1L) == "c") {
#     if (inherits(x, "glm")) {
#       res <- as.vector(residuals(x, "working")) * weights(x,
#                                                           "working")
#       if (!(substr(x$family$family, 1L, 17L) %in% c("poisson",
#                                                     "binomial", "Negative Binomial"))) {
#         res <- res * sum(weights(x, "working"), na.rm = TRUE)/sum(res^2,
#                                                                   na.rm = TRUE)
#       }
#     }
#     else if (inherits(x, "lm")) {
#       res <- as.vector(residuals(x))
#       if (!is.null(weights(x)))
#         res <- res * weights(x)
#     }
#   }
#   pos <- weights(x) > 0
#   if(sum(pos) < n) {
#     diaghat_full <- rep(0, n)
#     diaghat_full[pos] <- diaghat
#     diaghat <- diaghat_full
#   }
#   if (is.null(omega)) {
#     type <- match.arg(type)
#     if (type == "HC")
#       type <- "HC0"
#     switch(type, const = {
#       omega <- function(residuals, diaghat, df) rep(1,
#                                                     length(residuals)) * sum(residuals^2)/df
#     }, HC0 = {
#       omega <- function(residuals, diaghat, df) residuals^2
#     }, HC1 = {
#       omega <- function(residuals, diaghat, df) residuals^2 *
#         length(residuals)/df
#     }, HC2 = {
#       omega <- function(residuals, diaghat, df) residuals^2/(1 -
#                                                                diaghat)
#     }, HC3 = {
#       omega <- function(residuals, diaghat, df) residuals^2/(1 -
#                                                                diaghat)^2
#     }, HC4 = {
#       omega <- function(residuals, diaghat, df) {
#         n <- length(residuals)
#         p <- as.integer(round(sum(diaghat), digits = 0))
#         delta <- pmin(4, n * diaghat/p)
#         residuals^2/(1 - diaghat)^delta
#       }
#     }, HC4m = {
#       omega <- function(residuals, diaghat, df) {
#         gamma <- c(1, 1.5)
#         n <- length(residuals)
#         p <- as.integer(round(sum(diaghat), digits = 0))
#         delta <- pmin(gamma[1], n * diaghat/p) + pmin(gamma[2],
#                                                       n * diaghat/p)
#         residuals^2/(1 - diaghat)^delta
#       }
#     }, HC5 = {
#       omega <- function(residuals, diaghat, df) {
#         k <- 0.7
#         n <- length(residuals)
#         p <- as.integer(round(sum(diaghat), digits = 0))
#         delta <- pmin(n * diaghat/p, pmax(4, n * k *
#                                             max(diaghat)/p))
#         residuals^2/sqrt((1 - diaghat)^delta)
#       }
#     })
#   }
#   if (is.function(omega))
#     omega <- omega(res, diaghat, df)
#   rval <- sqrt(omega) * X
#   rval <- crossprod(rval)/n
#   return(rval)
# }
# 
# 
# vcovHC.lm <- function (x, type = c("HC3", "const", "HC", "HC0", "HC1", "HC2",
#                                    "HC4", "HC4m", "HC5"), omega = NULL, sandwich = TRUE, ...)
# {
#   type <- match.arg(type)
#   rval <- causalOT:::meatHC(x, type = type, omega = omega)
#   if (sandwich)
#     rval <- sandwich::sandwich(x, meat. = rval, ...)
#   return(rval)
# }

ci_asympt <- function(object, parm, level, ...) {
  
  if (object$options$split.model == FALSE && 
      !is.null(object$outcome.model.fit) && 
      "vcov" %in% attributes(methods(class = class(object$outcome.model.fit))  )$info$generic) {
    if(!is.null(object$outcome.model.fit$weights) && any(object$outcome.model.fit$weights == 0)) {
    #   # E_Y1 <- object$variance.component$E_Y1 #expectation of Y(1)
    #   # E_Y0 <- object$variance.component$E_Y0 #expectation of Y(0)
    #   # E_Y1_X <- object$variance.component$E_Y1_X  #estimated values of E(Y(1) | X)
    #   # E_Y0_X <- object$variance.component$E_Y0_X 
    #   # object$formula <-
    #   # object$outcome.model.fit <- NULL
    #   # return(ci_semiparm_eff(object, parm, level = level, ...))
      SDS <- sqrt(diag(sandwich::vcovBS(object$outcome.model.fit,...)))
      SD <- SDS["z"]
    } else {
      vcov_mat <- tryCatch(sandwich::vcovCL(object$outcome.model.fit,...),
                           error = sandwich::vcovBS(object$outcome.model.fit,...),
                           warning = sandwich::vcovBS(object$outcome.model.fit,...))
      SDS <- sqrt(diag(vcov_mat))
      SD <- SDS["z"]
    }
    
    
    # get levels
    quant.lower <- (1 - level) * 0.5
    quant.upper <- 1 - quant.lower
    CI <- c(qnorm(c(quant.lower, quant.upper),
                mean = object$estimate, sd = SD))
    
    return(list(CI = CI,
                SD = SD
                ))
  } else {
    return(ci_semiparm_eff(object, parm, level, ...))
  }
  
}

ci_semiparm_eff <- function(object, parm, level, ...) {
  
  # calculates variance of estimating y(z)
  semipar_var <- function(y, z, yhat, w, e_y, n) {
    resid <- w * y * z - e_y - yhat * (w * z - 1)
    return( mean( resid^2 ) / denom )
  }
  
  # mean(w^2/n^2 *z *(y - yhat)) + mean((w * z * yhat - tau)^2)
  
  # semipar_var_ate <- function(y, z, yhat_1, yhat_0, w, tau, n) {
  #   E_var_y1_given_x <- mean(w * z * (y - yhat_1)^2)
  #   E_var_y0_given_x <- mean(w * (1-z) * (y - yhat_0)^2)
  #   var_tau  <- mean( ((yhat_1 - yhat_0) - tau)^2 )
  #   return((E_var_y1_given_x + E_var_y0_given_x + var_tau)/n)
  # }
  
  semipar_var_ate <- function(y, z, yhat_1, yhat_0, w, tau, n) {
    # E_var_y1_given_x <- mean(w * z * (y - yhat_1)^2)
    # E_var_y0_given_x <- mean(w * (1-z) * (y - yhat_0)^2)
    # var_tau  <- mean( ((yhat_1 - yhat_0) - tau)^2 )
    return(mean((w * z * (y - yhat_1) - w * (1 - z) * (y - yhat_0) + (yhat_1 - yhat_0) - tau)^2) / n )
  }
  
  # semipar_var_ate <- function(y, z, yhat_1, yhat_0, w, tau) {
  #   return(mean((w * z * (y - yhat_1) - w * (1-z) * (y - yhat_0) + (yhat_1 - yhat_0) - tau)^2 ))
  # }
  
  # from jose paper
  # semipar_var <- function(y, z, yhat, w, e_y) {
  #   resid <- w * y * z - e_y - yhat * (w*z - 1)
  #   return(mean(resid^2))
  # }
  
  # semipar_var <- function(y, z, yhat, w, e_y) {
  #   # resid <- w * y * z - e_y - yhat * (w*z - 1)
  #   # return(mean(resid^2))
  #   mean(w * z * (y - yhat)^2) + mean((yhat - e_y)^2)
  # }
  
  model.fun <- object$model
  form <- object$formula
  estimand <- object$estimand
  dots <- list(...)
  # b_names <- base::...names
  
  if(is.null(model.fun)) {
    if("model" %in% names(dots)) {
      
      model.fun <- calc_model(dots$model)
    } else {
      stop("Must specify an argument 'model' as input to this function or in the effect estimation function")
    }
    form <- NULL
  }
  
  if(is.null(form)) {
    if("formula" %in% names(dots)) {
      form <- calc_form(formula = dots$formula, 
                        doubly.robust = object$options$doubly.robust, 
                        target = object$estimand,
                        split.model = object$options$split.model
                        )
      environment(form$treated) <- environment(form$control) <- environment()
    } else {
      stop("Must specify an argument 'formula' as input to this function or in the effect estimation function")
    }
  } else {
    environment(form$treated) <- environment(form$control) <- environment()
  }
  
  # treatment effect
  tau  <- object$estimate #estimated treatment effect
  
  # get variance components
  E_Y1 <- object$variance.component$E_Y1 #expectation of Y(1)
  E_Y0 <- object$variance.component$E_Y0 #expectation of Y(0)
  E_Y1_X <- object$variance.component$E_Y1_X  #estimated values of E(Y(1) | X)
  E_Y0_X <- object$variance.component$E_Y0_X  #estimated values of E(Y(0) | X)
  
  #get observed values of outcome
  data.obj <- object$data
  z       <- data.obj[,object$options$treatment.indicator]
  y       <- data.obj[,object$options$outcome]
  n       <- length(y)
  n1      <- sum(z)
  n0      <- n-n1
  
  weights <- object$weights
  # get conditional mean predictions for each unit
  if ( estimand == "ATE" && ( any(is.na(E_Y1_X)) || any(is.na(E_Y0_X)) ) ) {
    x     <- data.obj[, object$options$balance.covariates, drop = FALSE]
    fit_0 <- model.fun(form$control, as.data.frame(cbind(y, x)[z==0,,drop = FALSE]), weights = weights$w0)
    fit_1 <- model.fun(form$treated, as.data.frame(cbind(y, x)[z==1,,drop = FALSE]), weights = weights$w1)
    
    E_Y0_X <- predict(fit_0, newdata = x)
    E_Y1_X <- predict(fit_1, newdata = x)
    
  } else if (estimand == "ATT" && any(is.na(E_Y0_X)) ) {
    x    <- data.obj[, object$options$balance.covariates, drop = FALSE]
    
    fit_0 <- model.fun(form$control, as.data.frame(cbind(y, x)[z==0,,drop = FALSE]), weights = weights$w0)
    
    E_Y0_X <- predict(fit_0, newdata = x)
  } else if (estimand == "ATC" && any(is.na(E_Y1_X)) ) {
    x    <- data.obj[, object$options$balance.covariates, drop = FALSE]
    
    fit_1 <- model.fun(form$treated, as.data.frame(cbind(y, x)[z==1,,drop = FALSE]), weights = weights$w1)
    
    E_Y1_X <- predict(fit_1, newdata = x)
  } 
  
  # reorder data
  orders  <- order(z)
  w       <- c(weights$w0, weights$w1)
  
  # if (object$options$hajek == TRUE) {
  #   if (estimand == "ATE") {
  #     w <- w * n
  #     denom <- n
  #   } else if (estimand == "ATT") {
  #     w <- w * n1
  #     denom <- n1
  #   } else if (estimand == "ATC") {
  #     w <- w * n0
  #     denom <- n0
  #   }
  # } else {
  #   if (estimand == "ATE") {
  #     denom <- n
  #   } else if (estimand == "ATT") {
  #     denom <- n1
  #   } else if (estimand == "ATC") {
  #     denom <- n0
  #   }
  # }
  if (estimand == "ATE") {
    w <- w * n
    denom <- n
  } else if (estimand == "ATT") {
    w <- w * n1
    denom <- n1
  } else if (estimand == "ATC") {
    w <- w * n0
    denom <- n0
  }

  y       <- y[orders]
  z       <- sort(z)
  
  # get variance estimates
  # semipar_var(y, z, yhat, w, e_y)
  if (estimand == "ATE" ) {
    E_Y0_X <- E_Y0_X[orders]
    E_Y1_X <- E_Y1_X[orders]
    # VAR_Y1 <- semipar_var(y, z = z,
    #                       yhat = E_Y1_X, w = w, e_y = E_Y1,
    #                       n = denom)
    # VAR_Y0 <- semipar_var(y, z = 1 - z, 
    #                       yhat = E_Y0_X, w = w, e_y = E_Y0,
    #                       n = denom)
    VAR    <- semipar_var_ate(y, z = z, yhat_1 = E_Y1_X,
                              yhat_0 = E_Y0_X, w, tau = tau,
                              n = denom)
    # VAR <- (n1/n * VAR_Y1 + n0/n * VAR_Y0)
    # VAR <- (VAR_Y1 + VAR_Y0)
    
  } else if (estimand == "ATT" ) {
    E_Y0_X <- E_Y0_X[orders]
    VAR_Y1 <- var(y[z==1]) # assuming var(Y | Z) = var(Y(Z))
    VAR_Y0 <- semipar_var(y, z = 1 - z, yhat = E_Y0_X, w = w, e_y = E_Y0, n = denom)
    VAR <- (VAR_Y1 + VAR_Y0)/n1
  } else if (estimand == "ATC"  ) {
    E_Y1_X <- E_Y1_X[orders]
    VAR_Y1 <- semipar_var(y, z = z,     yhat = E_Y1_X, w = w, e_y = E_Y1, n = denom)
    VAR_Y0 <- var(y[z==0]) # assuming var(Y | Z) = var(Y(Z))
    VAR <- (VAR_Y1 + VAR_Y0)/n0
  }
  
  
  SD  <- sqrt(VAR)
  
  # get levels
  quant.lower <- (1 - level) * 0.5
  quant.upper <- 1 - quant.lower
  
  # calculate CI
  LWR <- qnorm(quant.lower, mean = tau, sd = SD)
  UPR <- qnorm(quant.upper, mean = tau, sd = SD)
  
  output <- list(
                  CI = c(LWR, UPR),
                  SD = SD
                )
  
  return(output)
  
}

ci_jack_ce <- function(object, parm, level, ...) {
    est_jk <- function(i, dat, weight, n0, n1, dr) {
      idx_i <- 1:n0
      idx_j <- 1:n1
      if (i <= n0) {
        idx_i <- idx_i[-i]
      } else if (i > n0) {
        idx_j <- idx_j[-(i - n0)]
      }
      
      dat.new <- dat[-i,,drop = FALSE]
      
      weight.new <- weight
      estimand <- weight.new$estimand
      method <- weight.new$method
      
      if (method %in% ot.methods()) {
        if (estimand != "ATE") {
          weight.new$gamma <- matrix(renormalize(weight.new$gamma[idx_i, idx_j]),
                                     length(idx_i), length(idx_j))
          weight.new$w0 <- rowSums(weight.new$gamma)
          weight.new$w1 <- colSums(weight.new$gamma)
        } else {
          weight.new$gamma <- list(matrix(renormalize(weight.new$gamma[[1]][idx_i, -i]),
                                          length(idx_i), nrow(dat.new)),
                                   matrix(renormalize(weight.new$gamma[[2]][idx_j, -i]),
                                          length(idx_j), nrow(dat.new)))
          weight.new$w0 <- rowSums(weight.new$gamma[[1]])
          weight.new$w1 <- rowSums(weight.new$gamma[[2]])
        }
      } else {
        weight.new$w0 <- renormalize(weight.new$w0[idx_i])
        weight.new$w1 <- renormalize(weight.new$w1[idx_j])
      }
      
      
      return(estimate_effect(data = dat.new, weights = weight.new,
                             estimand = weight.new$estimand, doubly.robust = dr,
                             outcome = 1, treatment.indicator = 2,
                             matched = FALSE,
                             balance.covariates = 3:ncol(dat.new))$estimate)
    }
    
    data <- object$data
    z    <- data[,data$options$treatment.indicator]
    x    <- data[,data$options$balance.covariates, drop = FALSE]
    y    <- data[,data$options$outcome]
    x0 <- x[z == 0, ]
    x1 <- x[z == 1, ]
    
    y0 <- y[z == 0]
    y1 <- y[z == 1]
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    n <- n0 + n1
    x  <- rbind(x0, x1)
    
    
    dat <- data.frame(y = c(y0,y1), z = c(rep(0, n0), rep(1,n1)), x)
    
    tau <- object$estimate
    
    weight <- object$weights
    
    dr <- object$options$doubly.robust
    
    jk.estimates <- vapply(1:nrow(dat), FUN = est_jk, FUN.VALUE = 1, dat = dat,
                           n0 = n0, n1 = n1,
                           weight = weight, dr = dr)
    var.jk <- (n - 1) / n * sum( (jk.estimates - tau)^2)
    
    upr <- tau + qnorm( 1 - (1 - level)/2) * sqrt(var.jk)
    lwr <- tau + qnorm((1 - level)/2)      * sqrt(var.jk)
    
    out <- list(CI = c(lwr, upr),
                SD =  sqrt(var.jk),
                jk = jk.estimates,
                tau = tau)
    
    return(out)
}

ci_boot_ce <- function(object, parm = NULL, level, n.boot = 1000, 
                       boot.method = NULL,
                       verbose = FALSE, ...) {
  
  boot.fun <- function(idx, object,  n0, n1, ...) {
    weight <- object$weights
    balance.covariates <- attributes(object$data)$balance.covariates
    outcome <- attributes(object$data)$outcome
    treatment.indicator <- attributes(object$data)$treatment.indicator
    # wt.args <- c(list(data = object$data, constraint = weight$args$constraint,
    #                   estimand = object$estimand, method = weight$method,
    #                   transport.matrix = !is.null(weight$gamma),
    #                   grid.search = isTRUE(weight$args$grid.search), 
    #                   formula = weight$args$formula,
    #                   balance.covariates = balance.covariates,
    #                   treatment.indicator = treatment.indicator), 
    #              weight$args,
    #              ...)
    # 
    # wt.args <- wt.args[!duplicated(names(wt.args))]
    # wtargn <- names(wt.args)
    # wf.call <- as.call(c(list(as.name("calc_weight")), wtargn))
    # 
    # weight <- eval(wf.call, envir = wt.args) 
    weight.boot <- ci_idx_to_wt(idx = idx,
                            estimand = object$estimand, 
                            method = weight$method, 
                            weight = weight, 
                            object = object,
                            balance.covariates = balance.covariates, 
                            treatment.indicator = treatment.indicator,
                            outcome = outcome,
                            n0 = n0, n1 = n1, ...)
    if (all(is.na(weight.boot$w0)) || all(is.na(weight.boot$w1))) {
      return(NULL)
    }
    # est.args <- c(list(data = data, formula = object$formula,
    #                    weights = weight,
    #                    model = object$model,
    #                    hajek = object$options$hajek,
    #                    doubly.robust = object$options$doubly.robust,
    #                    matched = object$options$matched,
    #                    estimand = object$estimand,
    #                    split.model = object$split.model,
    #                    balance.covariates = balance.covariates,
    #                    treatment.indicator = treatment.indicator,
    #                    outcome = outcome),
    #               object$options$addl.args,
    #               ...)
    # 
    # est.args <- est.args[!duplicated(names(est.args))]
    # estargn <- names(est.args)
    # estf.call <- as.call(c(list(as.name("estimate_effect")), estargn))
    # 
    # est <- eval(estf.call, envir = est.args)$estimate
    est  <- ci_idx_to_est(idx = idx, 
                              object = object,
                              weight = weight.boot,
                              balance.covariates = balance.covariates, 
                              treatment.indicator = treatment.indicator,
                              outcome = outcome,
                              n0 = n0, n1 = n1,
                              ...)
    # est <- eval(ef$call, envir = ef$args)$estimate
    # clear_subset(data)
    return(est)
  }
  
  if (is.null(boot.method)) {
    boot.method <- switch(object$weights$method,
                          "None" = "n-out-of-n",
                          "SBW" = "n-out-of-n",
                          "Logistic" = "n-out-of-n",
                          "NNM" = "m-out-of-n",
                          "Wasserstein" = "m-out-of-n",
                          "Constrained Wasserstein" = "m-out-of-n",
                          "m-out-of-n"
                          )
  } else {
    boot.method <- match.arg(boot.method, choices = c("m-out-of-n", "n-out-of-n"))
  } 
  
  ns <- get_n(object$data, balance.covariates = object$options$balance.covariates,
              treatment.indicator = object$options$treatment.indicator, ...)
  # N  <- sum(ns)
  
  boot.idx <- cot_boot_samples(n.boot = n.boot,
                               boot.method = boot.method,
                               estimand = object$estimand,
                               method = object$weights$method,
                               matched = object$options$matched,
                               n0 = ns[1],
                               n1 = ns[2]
                               )
  # replicate(n = n.boot, sample.int(N,N, replace = TRUE))
  if (verbose) {
    pbapply::pboptions(type = "timer", style = 3, char = "=")
  } else {
    pbapply::pboptions(type = "none")
  }
  # if (boot.method == "m-out-of-n") {
  #   nSamp <- c(adjust_m_of_n_btstrp(ns[1]),
  #              adjust_m_of_n_btstrp(ns[2]))
  # } else {
  #   nSamp <- ns
  # }
  boots <- unlist(pbapply::pbsapply(boot.idx, boot.fun, data = data, object = object, 
                  n0 = ns[1], n1 = ns[2],
                  ...))
  pbapply::pboptions(type = "none")
  
  ci <- quantile(boots, probs = c((1 - level)/2, 1 - (1 - level)/2))
  # names(ci) <- c("lower", "upper")
  sd_boot <- sd(boots)
  
  output <- list(CI = ci, SD = sd_boot)
  if (boot.method == "m-out-of-n") {
    ci_raw <- ci
    sd_raw <- sd_boot
    
    m <- sum(c(adjust_m_of_n_btstrp(ns[1]),
             adjust_m_of_n_btstrp(ns[2])))
    n <- sum(ns)
    
    theta_hat <- object$estimate
    theta_hat_m <- mean(boots)
    l <- ci[1]
    u <- ci[2]
    ci_names <- names(ci)
    output$CI <- c(theta_hat + sqrt(m / n) * (l - theta_hat_m),
            theta_hat + sqrt(m / n) * (u - theta_hat_m))
    names(output$CI) <- ci_names
    output$SD <- sd_boot * sqrt(m / n)
    output$unadjusted = list(CI = ci_raw,
                             SD = sd_raw)
  }
  
  return(output)
  
}

adjust_m_of_n_btstrp <- function(n, scale = 0.75) {
  return(ceiling(n * 0.75))
}

cot_boot_samples <- function(n.boot, boot.method, estimand, method, matched, n0, n1, ...) {
  n <- n0 + n1
  if (boot.method == "m-out-of-n") {
    n0samp <- adjust_m_of_n_btstrp(n0)
    n1samp <- adjust_m_of_n_btstrp(n1)
    nsamp  <- n0samp + n1samp
  } else {
    n0samp <- n0
    n1samp <- n1
    nsamp  <- n
  }
  matched <- isTRUE(matched)
  # boot.idx <- replicate(n = n.boot, sample.int(n,n,replace = TRUE),
  #                       simplify = FALSE)
  # return(boot.idx)
  # if that doesn't work then: 
  if (!matched && !(method %in% ot.methods() )) {
    boot.idx <- replicate(n = n.boot, sample.int(n = n, size = nsamp,
                                                 replace = TRUE),
                          simplify = FALSE)
  } else {
    if (estimand == "ATE") {
      boot.idx <- replicate(n = n.boot,
                            list(control = c(sample.int(n0, 
                                                        n0samp, replace = TRUE)),
                                 treated = c(sample.int(n1, 
                                                        n1samp, replace = TRUE)),
                                 total   = sample.int(n, 
                                                      nsamp, replace = TRUE)), 
                            simplify = FALSE)
    } else if (estimand == "ATT" | estimand == "ATC") {
      boot.idx <- replicate(n = n.boot, 
        list(sample.int(n0,n0samp, replace = TRUE), 
          sample.int(n1,n1samp, replace = TRUE))
      , simplify = FALSE)
    }
  }
  
  
  return(boot.idx)
}

ci_idx_to_wt <- function(idx, estimand, method, weight, object,
                              balance.covariates, treatment.indicator,
                              outcome,
                              n0, n1, ...) {
  weight$args$grid.search <- FALSE
  # wt.args <- c(list(data = object$data[idx,, drop = FALSE], 
  #                   constraint = weight$args$constraint,
  #                   estimand = object$estimand, method = weight$method,
  #                   transport.matrix = !is.null(weight$gamma),
  #                   grid.search = isTRUE(weight$args$grid.search), 
  #                   formula = weight$args$formula,
  #                   balance.covariates = balance.covariates,
  #                   treatment.indicator = treatment.indicator), 
  #              weight$args,
  #              ...)
  # wt.args <- wt.args[!duplicated(names(wt.args))]
  # wtargn <- lapply(names(wt.args), as.name)
  # names(wtargn) <- names(wt.args)
  # wf.call <- as.call(c(list(as.name("calc_weight")), wtargn))
  # return(eval(wf.call, envir = wt.args))
  #if this doesn't work, then try the following
  if (object$options$matched == FALSE &&
    (method == "SBW" | method == "Logistic" | method == "None") ) {
    wt.args <- c(list(data = object$data[idx,, drop = FALSE], 
                      constraint = weight$args$constraint,
                      estimand = object$estimand, method = weight$method,
                      transport.matrix = !is.null(weight$gamma),
                      grid.search = isTRUE(weight$args$grid.search), 
                      formula = weight$args$formula,
                      outcome = outcome,
                      balance.covariates = balance.covariates,
                      treatment.indicator = treatment.indicator,
                      ...), 
                 weight$args
                 )
    wt.args <- wt.args[!duplicated(names(wt.args))]
    wt.args <- wt.args[!sapply(wt.args, is.null)]
    wtargn <- lapply(names(wt.args), as.name)
    names(wtargn) <- names(wt.args)
    wf.call <- as.call(c(list(as.name("calc_weight")), wtargn))
    return( tryCatch(eval(wf.call, envir = wt.args),
                    error = function(e) {
                      warning(e$message)
                      return(calc_weight_error())
                    }) )
    
  } else {
    if (estimand == "ATE") {
      tab.0 <- tabulate(idx[[1]], n0)
      tab.1 <- tabulate(idx[[2]], n1)
      tab   <- tabulate(idx[[3]], n0 + n1)
      # tab   <- tabulate(unlist(idx[1:2]), 
      #                   n0 + n1)
      z     <- object$data[, treatment.indicator]
      dat.0 <-  rbind(object$data[z == 0, ],
                      object$data)
      dat.1 <-  rbind(object$data[z == 1, ],
                      object$data)
      dat.0[, treatment.indicator] <- c(rep(0, n0), 
                                        rep(1, n1 + n0))
      dat.1[, treatment.indicator] <- c(rep(0, n1), 
                                        rep(1, n1 + n0))
      
      wt.args <- list(
                  c(list(data = dat.0, constraint = weight$args$constraint[[1]],
                        estimand = "ATT", method = weight$method,
                        transport.matrix = !is.null(weight$gamma),
                        grid.search = isTRUE(weight$args$grid.search), 
                        formula = weight$args$formula,
                        outcome = outcome,
                        balance.covariates = balance.covariates,
                        treatment.indicator = treatment.indicator,
                        sample_weight = renormalize(c(tab.0, tab)),
                        ...), 
                   weight$args),
                   c(list(data = dat.1, constraint = weight$args$constraint[[2]],
                          estimand = "ATT", method = weight$method,
                          transport.matrix = !is.null(weight$gamma),
                          grid.search = isTRUE(weight$args$grid.search), 
                          formula = weight$args$formula,
                          outcome = outcome,
                          balance.covariates = balance.covariates,
                          treatment.indicator = treatment.indicator, 
                     sample_weight = renormalize(c(tab.1, tab)),
                     ...),
                     weight$args)
      )
      wt.args <- list(wt.args[[1]][!duplicated(names(wt.args[[1]]))],
                      wt.args[[2]][!duplicated(names(wt.args[[2]]))])
      wtargn  <- lapply(names(wt.args[[1]]), as.name)
      names(wtargn)  <- names(wt.args[[1]])
      
      wf.call <- as.call(c(list(as.name("calc_weight")), wtargn))
      wts <- tryCatch(eval(wf.call, envir = wt.args[[1]]),
                      error = function(e) {
                        warning(e$message)
                        return(calc_weight_error()) }  )
      wts_1 <- tryCatch(eval(wf.call, envir = wt.args[[2]]),
                        error = function(e) {
                          warning(e$message)
                          return(calc_weight_error()) }  )
      wts$w1 <- wts_1$w0
      wts$gamma <- list(wts$gamma, wts_1$gamma)
      wts$estimand <- "ATE"
      return(wts)
    } else {
      tab.0 <- tabulate(idx[[1]], n0)
      tab.1 <- tabulate(idx[[2]], n1)
      sample.wt <- list(a = renormalize(tab.0),
                        b = renormalize(tab.1))
      
      # dat <- object$data[idx,]
      # dat[,treatment.indicator] <- c(rep(0, n0), rep(1,n1))
      wt.args <- c(list(data = object$data, constraint = weight$args$constraint,
                        estimand = object$estimand, method = weight$method,
                        transport.matrix = !is.null(weight$gamma),
                        grid.search = isTRUE(weight$args$grid.search), 
                        formula = weight$args$formula,
                        outcome = outcome,
                        balance.covariates = balance.covariates,
                        treatment.indicator = treatment.indicator, 
                        sample_weight = sample.wt,
                        ...),
                   weight$args)
    }
    wt.args <- wt.args[!duplicated(names(wt.args))]
    wt.args <- wt.args[!sapply(wt.args, is.null)]
    wtargn <- lapply(names(wt.args), as.name)
    names(wtargn) <- names(wt.args)
    wf.call <- as.call(c(list(as.name("calc_weight")), wtargn))
    return(tryCatch(eval(wf.call, envir = wt.args),
                    error = function(e) {
                      warning(e$message)
                      return(calc_weight_error())
                    })
           )
  }
  
}
  
ci_idx_to_est <- function(idx, 
                               object,
                               weight, 
                               balance.covariates,
                               treatment.indicator,
                               outcome,
                               n0, n1, ...) {
  # if (!(weight$method %in% c("SBW", "Logistic")) & object$estimand == "ATE") {
  #   idx <- c(idx[[1]], idx[[2]])
  # }
  # dat <- object$data[idx,, drop = FALSE]
  environment(object$formula) <- environment()
  
  if (!object$options$matched && (weight$method == "SBW" | weight$method == "Logistic" | weight$method == "None") ){
    sample.wt <- rep(1 / (n0 + n1), n0 + n1)
    
    est.args <- c(list(data = object$data[idx,], 
                       formula = object$formula,
                       weights = weight,
                       model = object$model,
                       hajek = object$options$hajek,
                       doubly.robust = object$options$doubly.robust,
                       matched = object$options$matched,
                       estimand = object$estimand,
                       split.model = object$options$split.model,
                       sample_weight = sample.wt,
                       balance.covariates = balance.covariates,
                       treatment.indicator = treatment.indicator,
                       outcome = outcome,
                       ...),
                  object$options$addl.args)
  } else {
    sample.wt <- list(
                      a = renormalize(tabulate(idx[[1]], n0)),
                      b = renormalize(tabulate(idx[[2]], n1))
                      )
    
    est.args <- c(list(data = object$data, 
                       formula = object$formula,
                       weights = weight,
                       model = object$model,
                       hajek = object$options$hajek,
                       doubly.robust = object$options$doubly.robust,
                       matched = object$options$matched,
                       estimand = object$estimand,
                       split.model = object$options$split.model,
                       sample_weight = sample.wt,
                       balance.covariates = balance.covariates,
                       treatment.indicator = treatment.indicator,
                       outcome = outcome,
                       ...),
                  object$options$addl.args)
    
  }
  
  est.args <- est.args[!duplicated(names(est.args))]
  estargn <- lapply(names(est.args), as.name)
  names(estargn) <- names(est.args)
  estf.call <- as.call(c(list(as.name("estimate_effect")), estargn))
  
  return(eval(estf.call, envir = est.args)$estimate)
}

setMethod("confint", c(object = "causalEffect"), confint.causalEffect)

  