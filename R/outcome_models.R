gp_pred <- function(formula = NULL, data, weights=NULL,
                    param, estimand = c("ATE","ATT","ATC","cATE"),...) {
  
  test_pos_def_inv <- function(x,y) {
    e <- eigen(x)
    if(any(e$values <=0)) {
      min.e <- min(e$values)
      e$values <- e$values - min.e
    }
    return(e$vectors %*% diag(1/e$values) %*% crossprod(e$vectors, y))
  }
  # w0 <- weights$w0
  # w1 <- weights$w1
  estimand <- match.arg(estimand)
  prep.data <- prep_data(data,...)
  
  z <- prep.data$z
  y <- prep.data$df$y
  x <- as.matrix(prep.data$df[,-which(colnames(prep.data$df)=="y")])
  
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
    m_y0   <- mean(y[z==0])
    m_y1   <- mean(y[z==1])
    sd_y0   <- sd(y[z==0])
    sd_y1   <- sd(y[z==1])
    y[z==0] <- c(scale(y[z==0]))
    y[z==1] <- c(scale(y[z==1]))
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
  #             solve(Kernel_full[[1]]$cov, y[z==0]))
  # } else if (estimand == "ATC") {
  #   y[sel]
  # }
  # pred1 <- if(estimand == "ATE" | estimand == "ATC") {
  #   crossprod(Kernel_full[[2]]$cross, 
  #             solve(Kernel_full[[2]]$cov, y[z==1]))
  # } else if (estimand == "ATT") {
  #   y[sel]
  # }
  if(param$metric == "mahalanobis") {
    cov <- cov(x)
    A <- scale(x, center = TRUE, scale = FALSE) %*% solve(chol(cov))
    
  } else {
    A <- x
  }
  A0 <- A[z==0,]
  A1 <- A[z==1,]
  
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
  }
  
  tau <- if(estimand == "ATE") {
    pred0 <- crossprod(kernel_cross0, test_pos_def_inv(kernel_cov0, y[z==0])) * sd_y0 + m_y0
    pred1 <- crossprod(kernel_cross1, test_pos_def_inv(kernel_cov1, y[z==1])) * sd_y1 + m_y1
    
    mean(pred1 - pred0)
  } else if (estimand == "ATT") {
    pred0 <- crossprod(kernel_cross0, test_pos_def_inv(kernel_cov0, y[z==0])) * sd_y0 + m_y0
    
    mean(y[z==1] * sd_y1 + m_y1 - pred0[z==1])
  } else if (estimand == "ATC") {
    pred1 <- crossprod(kernel_cross1, test_pos_def_inv(kernel_cov1, y[z==1])) * sd_y1 + m_y1
    
    mean(pred1[z==0]* - y[z==0]*sd_y0 - m_y0)
  }
  
  return(tau)
}
  

outcome_calc <- function(data, z, weights, formula, model.fun, matched, estimand) {
  
  w0 <- weights$w0
  w1 <- weights$w1
  t_ind <- z==1
  c_ind <- z==0
  
  fit_1 <- model.fun(formula$treated, data[t_ind,,drop=FALSE])
  fit_0 <- model.fun(formula$control, data[c_ind,,drop=FALSE])
  f_1   <- predict(fit_1, data) 
  f_0   <- predict(fit_0, data)
  
  if (matched) {
    # if(is.null(weights$gamma)) {
    #   stop("Transport matrix must be specified for matched estimator")
    # }
    # rad   <- (2*z - 1)
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


outcome_calc_model <- function(data, z, weights, formula, model.fun, matched, estimand,...) {
  w0 <- weights$w0
  w1 <- weights$w1

  n <- length(z)
  w <- rep(NA,n)
  w[z==1] <- w1
  w[z==0] <- w0

  fit <- model.fun(formula, data = cbind(data,z=z), weights = w, ...)


  tx_effect <- coef(fit)["z"]

  return(tx_effect)
}

calc_form <- function(formula, doubly.robust, target) {
  if(!is.null(formula)){
    stopifnot(isTRUE(all(names(formula) %in% c("treated","control"))))
    formula$treated <- as.formula(formula$treated)
    formula$control <- as.formula(formula$control)
    return(formula)
  }
  form_mod <- formula(y~.)
  form_obs <- formula(y~1)
  
  formula$treated <- form_obs
  formula$control <- form_obs
  
  if (doubly.robust) {
      if(target == "ATT") {
        formula$control <- form_mod
      } else if(target == "ATC") {
        formula$treated <- form_mod
      } else if (target == "ATE" | target == "feasible") {
        formula$treated <- form_mod
        formula$control <- form_mod
      }
  }
  return(formula)
}

calc_model <- function(model.fun) {
  if(!is.null(model.fun)) {
    model.fun <- model.fun
    stopifnot(is.function(model.fun))
  } else {
    model.fun <- lm
  }
  return(model.fun)
}

calc_hajek <- function(weights, target, hajek) {
  if(hajek) {
    weights$w0 <- renormalize(weights$w0)
    weights$w1 <- renormalize(weights$w1)
  }
  
  if(inherits(weights, "causalWeights")) {
    estimand <- switch(weights$estimand,
                       "ATT" = "ATT",
                       "ATC" = "ATC",
                       "ATE" = "ATE",
                       "cATE" = "ATE",
                       "feasible" = "feasible")
    if(estimand != target) stop("Weights not estimating the target")
  }
  return(weights)
}

estimate_effect <- function(data, formula = NULL, weights, 
                                  hajek = TRUE, 
                                  doubly.robust = TRUE,
                                  matched = FALSE,
                                  target = c("ATT", "ATC", "ATE", "feasible"), 
                                  model = NULL, 
                                  split.model = TRUE, ...) {
  #get args
  dots <- list(...)
  
  # hajek <- match.arg(hajek)
  dr <- isTRUE(doubly.robust[1])
  matched <- isTRUE(matched[1])
  
  split.model <- isTRUE(split.model)
  
  if(is.null(target)){
    if(inherits(weights, "causalWeights")){
      target <- switch(weights$estimand,
                       "ATT" = "ATT",
                       "ATC" = "ATC",
                       "ATE" = "ATE",
                       "cATE" = "ATE",
                       "feasible" = "feasible")
    }
  }
  target <- switch(target,
                   "ATT" = "ATT",
                   "ATC" = "ATC",
                   "ATE" = "ATE",
                   "cATE" = "ATE",
                   "feasible" = "feasible")
  targ <- match.arg(target)
  
  #set up model structures
  model.fun <- calc_model(model)
  formula <- calc_form(formula, dr, target)
  weights <- calc_hajek(weights, target, hajek)
  
  #setup data
  prep.data <- prep_data(data)
  
  #get estimate
  estimate <- if(split.model) {
    outcome_calc(prep.data$df, prep.data$z, weights, formula, model.fun,
                           matched, target)
  } else {
    outcome_calc_model(prep.data$df, prep.data$z, weights, formula, model.fun,
                       matched, target)
  }

  return(estimate)
}
  