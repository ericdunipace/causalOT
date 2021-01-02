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
  } else if (param$kernel == "linear") {
    kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + tcrossprod(A0)
    kernel_cross0 <- tcrossprod(A0,A)
    
    kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + tcrossprod(A1)
    kernel_cross1 <- tcrossprod(A1,A)
  }
  
  tau <- if(estimand == "ATE") {
    pred0 <- crossprod(kernel_cross0, test_pos_def_inv(kernel_cov0, y[z==0])) * sd_y0 + m_y0
    pred1 <- crossprod(kernel_cross1, test_pos_def_inv(kernel_cov1, y[z==1])) * sd_y1 + m_y1
    
    mean(pred1 - pred0)
  } else if (estimand == "ATT") {
    pred0 <- crossprod(kernel_cross0, test_pos_def_inv(kernel_cov0, y[z==0])) * sd_y0 + m_y0
    
    mean((y[z==1] * sd_y1 + m_y1) - pred0[z==1])
  } else if (estimand == "ATC") {
    pred1 <- crossprod(kernel_cross1, test_pos_def_inv(kernel_cov1, y[z==1])) * sd_y1 + m_y1
    
    mean(pred1[z==0] - (y[z==0] * sd_y0 + m_y0))
  }
  
  return(tau)
}
  

mapping <- function(data, z, weights, estimand, f1, f0, ...) {
  n  <- length(z)
  n1 <- sum(z)
  n0 <- n - n1
  runMap <- isTRUE(weights$method %in% ot.methods())
  
  if (runMap) {
    bal.cov <- colnames(data)[colnames(data) != "y"]
    data$z  <- z
    data$f1 <- f1
    data$f0 <- f0
    bproj_y <- barycentric_projection(data, weights, 
                                      treatment.indicator = "z", 
                                      outcome = "y",
                                      balance.covariates = bal.cov,
                                      estimand = estimand,
                                      ...)
    
    if(estimand == "ATE" | estimand == "ATC") {
      # the model projection for the controls and then predicted value of E(Y|X for controls)
      bproj_f1<- barycentric_projection(data, weights, 
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
      bproj_f0<- barycentric_projection(data, weights, 
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
    
    y1 <- y1[!is.na(y1)]
    y0 <- y0[!is.na(y0)]
    
  } else {
    n          <- length(z)
    y0         <- y1 <- rep(NA_real_, n)
    y0[z==0]   <- data$y[z==0] - f1[z==0]
    y1[z==1]   <- data$y[z==1] - f0[z==1]
    if(estimand == "ATE" | estimand == "ATT") {
      y0[z==1] <- (data$y[z==0] -  f0[z==0]) %*% weights$w0
      y0 <- y0[z==1]
      y1 <- y1[z==1]
    }
    if (estimand == "ATE" | estimand == "ATC") {
      y1[z==0] <- (data$y[z==1] -  f1[z==1]) %*% weights$w1
      y0 <- y0[z==0]
      y1 <- y1[z==0]
      
    }
    
    
  }
  n1 <- length(y1)
  n0 <- length(y0)
  new_weights <- list(w0 = rep(1/n0, n0),
                      w1 = rep(1/n1, n1))
  new_weights <- list(w0 = rep(1/n0, n0),
                      w1 = rep(1/n1, n1))
  return(list(y0 = y0,
              y1 = y1,
              weights = new_weights))
}

.outcome_calc_deprecated <- function(data, z, weights, formula, model.fun, matched, estimand) {
  
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

    maps <- mapping(data = data, z = z, weights = weights, estimand = estimand, f1 = f_1, f0 = f_0)
    
    # tau_t <- 0
    # tau_c <- 0
    # if(estimand == "ATT" | estimand == "ATE" | estimand == "feasible") {
    #   tau_t <- c(mean(maps$y1[z==1]) - mean(maps$y0[z==1]))
    # }
    # if(estimand == "ATC" | estimand == "ATE" | estimand == "feasible") {
    #   tau_c <- c(mean(maps$y1[z==0]) - mean(maps$y0[z==0]))
    # }
    # if(estimand == "ATT"){
    #   tx_effect <- tau_t
    # } else if (estimand == "ATC") {
    #   tx_effect <- tau_c
    # } else if(estimand == "ATE" | estimand == "feasible") {
    #   if(estimand == "ATE") {
    #     t_w   <- sum(t_ind) #1/(sum(w1^2))
    #     c_w   <- sum(c_ind) #1/(sum(w0^2))
    #   } else if (estimand == "feasible") {
    #     t_w   <- 1/(sum(w1^2))
    #     c_w   <- 1/(sum(w0^2))
    #   }
    #   
    #   tx_effect <- ( t_w * tau_t + c_w * tau_c)/(c_w + t_w)
    # }
    tx_effect <- mean(maps$y1 - maps$y0)
      
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
  w[z == 1] <- w1
  w[z == 0] <- w0
  

  environment(formula) <- environment()
  
  fit <- model.fun(formula, 
                   data =  cbind(data, 
                                 z = z), 
                   weights = w, 
                   ...)

  tx_effect <- coef(fit)["z"]

  return(tx_effect)
}

calc_form <- function(formula, doubly.robust, target, split.model) {
  if (!is.null(formula)) {
    if ( isTRUE(split.model) ) {
      stopifnot(isTRUE(all(names(formula) %in% c("treated","control"))))
      formula$treated <- as.formula(formula$treated)
      formula$control <- as.formula(formula$control)
      
    } else {
      formula <= as.formula(formula)
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
    if (estimand != target) stop("Weights do not estimating the target estimand")
  }
  return(weights)
}

estimate_effect <- function(data, formula = NULL, weights, 
                                  hajek = TRUE, 
                                  doubly.robust = TRUE,
                                  matched = FALSE,
                                  estimand = c("ATT", "ATC", "ATE", "feasible"),
                                  target = c("ATT", "ATC", "ATE", "feasible"), 
                                  model = NULL, 
                                  split.model = TRUE, ...) {
  #get args
  dots <- list(...)
  
  # hajek <- match.arg(hajek)
  dr <- isTRUE(doubly.robust[1])
  matched <- isTRUE(matched[1])
  
  split.model <- isTRUE(split.model)
  
  if (is.null(estimand)) {
    if (!is.null(target)) estimand <- target
    if (inherits(weights, "causalWeights")) {
      estimand <- switch(weights$estimand,
                       "ATT" = "ATT",
                       "ATC" = "ATC",
                       "ATE" = "ATE",
                       "cATE" = "ATE",
                       "feasible" = "feasible")
    }
  } else {
    estimand <- switch(estimand,
                       "ATT" = "ATT",
                       "ATC" = "ATC",
                       "ATE" = "ATE",
                       "cATE" = "ATE",
                       "feasible" = "feasible")
    estimand <- match.arg(estimand)
  }
  
  # targ <- match.arg(target)
  
  #set up model structures
  model.fun <- calc_model(model)
  formula <- calc_form(formula, dr, estimand, split.model)
  weights <- calc_hajek(weights, estimand, hajek)
  
  #setup data
  prep.data <- prep_data(data,...)
  
  #get estimate
  estimate <- if (split.model) {
    outcome_calc(prep.data$df, prep.data$z, weights, formula, model.fun,
                           matched, estimand)
  } else {
    outcome_calc_model(prep.data$df, prep.data$z, weights, formula, model.fun,
                       matched, estimand)
  }
  output.df <- cbind(prep.data$df, z = prep.data$z)
  attr(output.df, "balance.covariates") <- attributes(prep.data$df)$balance.covariates
  attr(output.df, "outcome") <- attributes(prep.data$df)$outcome
  attr(output.df, "treatment.indicator") <- "z"
  
  output <- list(estimate = estimate,
                 data = output.df,
                 model = model.fun,
                 formula = formula,
                 weights = weights,
                 estimand = estimand,
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
  
  class(output) <- "causalEffect"
  return(output)
}

confint.causalEffect <- function(object, parm, level = 0.95, method = c("bootstrap", "asymptotic"), ...) {
  
  if (method == "bootstrap") {
    return(ci_boot_ce(object, parm, level, ...))
  } else {
    stop("asymptotic confidence intervals not currently implemented")
  }
  
  
}

ci_boot_ce <- function(object, parm, level, n.boot = 1000, balance.covariates = NULL,
                       treatment.indicator = NULL, outcome = NULL, 
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
                            n0 = n0, n1 = n1, ...)
    
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
  
  ns <- get_n(object$data, balance.covariates = object$options$balance.covariates,
              treatment.indicator = object$options$treatment.indicator, ...)
  # N  <- sum(ns)
  
  boot.idx <- cot_boot_samples(n.boot = n.boot,
                               estimand = object$estimand,
                               method = object$weights$method,
                               n0 = ns[1],
                               n1 = ns[2]
                               )
  # replicate(n = n.boot, sample.int(N,N, replace = TRUE))
  if (verbose) {
    pbapply::pboptions(type = "timer", style = 3, char = "=")
  } else {
    pbapply::pboptions(type = "none")
  }
  boots <- pbapply::pbsapply(boot.idx, boot.fun, data = data, object = object, 
                  n0 = ns[1], n1 = ns[2],
                  ...)
  pbapply::pboptions(type = "none")
  
  ci <- quantile(boots, probs = c(level/2, (1 - level/2)))
  # names(ci) <- c("lower", "upper")
  
  return(list(CI = ci, SD = sd(boots)))
  
}

cot_boot_samples <- function(n.boot, estimand, method, n0, n1, ...) {
  n <- n0 + n1
  boot.idx <- replicate(n = n.boot, sample.int(n,n,replace = TRUE),
                        simplify = FALSE)
  return(boot.idx)
  # if that doesn't work then: 
  if (method == "SBW" | method == "Logistic") {
    boot.idx <- replicate(n = n.boot, sample.int(n,n,replace = TRUE),
                          simplify = FALSE)
  } else {
    if (estimand == "ATE") {
      boot.idx <- replicate(n = n.boot,
                            function() {
                              total.idx <- sample.int(n, n, replace = TRUE)
                              return(list(control = c(sample.int(n0,n0, replace = TRUE), total.idx),
                                          treated = c(sample.int(n1,n1, replace = TRUE), total.idx)))
                            }, simplify = FALSE)
    } else if (estimand == "ATT" | estimand == "ATC") {
      boot.idx <- replicate(n = n.boot, function() {
        c(sample.int(n0,n0, replace = TRUE), 
          sample.int(n1,n1, replace = TRUE))
      }, simplify = FALSE)
    }
  }
  
  
  return(boot.idx)
}

ci_idx_to_wt <- function(idx, estimand, method, weight, object,
                              balance.covariates, treatment.indicator,
                              n0, n1, ...) {
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
  if (method == "SBW" | method == "Logistic") {
    wt.args <- c(list(data = object$data[idx,, drop = FALSE], 
                      constraint = weight$args$constraint,
                      estimand = object$estimand, method = weight$method,
                      transport.matrix = !is.null(weight$gamma),
                      grid.search = isTRUE(weight$args$grid.search), 
                      formula = weight$args$formula,
                      balance.covariates = balance.covariates,
                      treatment.indicator = treatment.indicator), 
                 weight$args,
                 ...)
    wt.args <- wt.args[!duplicated(names(wt.args))]
    wtargn <- names(wt.args)
    
  } else {
    if (estimand == "ATE") {
      dat.0 <-  object$data[idx[[1]], ]
      dat.1 <-  object$data[idx[[2]], ]
      dat.0[, treatment.indicator] <- c(rep(0, n0), rep(1, n1 + n0))
      dat.1[, treatment.indicator] <- c(rep(0, n1), rep(1, n1 + n0))
      wt.args <- list(
                  c(list(data = dat.0, constraint = weight$args$constraint,
                        estimand = object$estimand, method = weight$method,
                        transport.matrix = !is.null(weight$gamma),
                        grid.search = isTRUE(weight$args$grid.search), 
                        formula = weight$args$formula,
                        balance.covariates = balance.covariates,
                        treatment.indicator = treatment.indicator), 
                   weight$args,
                   ...),
                   c(list(data = dat.1, constraint = weight$args$constraint,
                          estimand = object$estimand, method = weight$method,
                          transport.matrix = !is.null(weight$gamma),
                          grid.search = isTRUE(weight$args$grid.search), 
                          formula = weight$args$formula,
                          balance.covariates = balance.covariates,
                          treatment.indicator = treatment.indicator), 
                     weight$args,
                     ...)
      )
      wt.args <- list(wt.args[[1]][!duplicated(names(wt.args[[1]]))],
                      wt.args[[2]][!duplicated(names(wt.args[[2]]))])
      wtargn <- names(wt.args[[1]])
    } else {
      dat <- object$data[idx,]
      dat[,treatment.indicator] <- c(rep(0, n0), rep(1,n1))
      wt.args <- c(list(data = dat, constraint = weight$args$constraint,
                        estimand = object$estimand, method = weight$method,
                        transport.matrix = !is.null(weight$gamma),
                        grid.search = isTRUE(weight$args$grid.search), 
                        formula = weight$args$formula,
                        balance.covariates = balance.covariates,
                        treatment.indicator = treatment.indicator), 
                   weight$args,
                   ...)
    }
    wt.args <- wt.args[!duplicated(names(wt.args))]
    wtargn <- names(wt.args)
    
  }
  wf.call <- as.call(c(list(as.name("calc_weight")), wtargn))
  return(eval(wf.call, envir = wt.args))
}
  
ci_idx_to_est <- function(idx, 
                               object,
                               weight, 
                               balance.covariates,
                               treatment.indicator,
                               outcome,
                               n0, n1, ...) {
  if (!(weight$method %in% c("SBW", "Logistic")) & object$estimand == "ATE") {
    idx <- c(idx[[1]][1:n0], idx[[2]][1:n1])
  }
  dat <- object$data[idx,, drop = FALSE]
  environment(object$formula) <- environment()
  est.args <- c(list(data = dat, 
                     formula = object$formula,
                     weights = weight,
                     model = object$model,
                     hajek = object$options$hajek,
                     doubly.robust = object$options$doubly.robust,
                     matched = object$options$matched,
                     estimand = object$estimand,
                     split.model = object$options$split.model,
                     balance.covariates = balance.covariates,
                     treatment.indicator = treatment.indicator,
                     outcome = outcome),
                object$options$addl.args,
                ...)
  
  est.args <- est.args[!duplicated(names(est.args))]
  estargn <- lapply(names(est.args), as.name)
  names(estargn) <- names(est.args)
  estf.call <- as.call(c(list(as.name("estimate_effect")), estargn))
  
  return(eval(estf.call, envir = est.args)$estimate)
}

setMethod("confint", c(object = "causalEffect"), confint.causalEffect)

  