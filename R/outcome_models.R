outcome_model_DRH <- function(.self, weights, ...) {
  dots <- list(...)
  est <- dots$est
  if(is.null(est)) est <- "ATE"
  design <- model.matrix(.self$formula, data=.self$x)
  n <- .self$n
  # design0 <- design[z==0,,drop=FALSE]
  # design1 <- design[z==1,,drop=FALSE]
  # y0 <- y[z==0]
  # y1 <- y[z==1]
  w0 <- weights$z0
  w1 <- weights$z1
  w <- rep(NA, n)
  w[.self$z==0] <- w0
  w[.self$z==1] <- w1
  w <- w/sum(w)
  
  #(w * (1-z))/sum((w * (1-z)))
  
  if(est == "ATE") {
    model0 <- lm(.self$y ~ design, subset = .self$z == 0)
    f0 <- predict(model0, data.frame(design))
    mu_0 <- mean(f0)
    e_0 <- model0$resid #c(y - f0)
    y_0 <- mu_0 + e_0 %*% w0 
    
    model1 <- lm(.self$y ~ design, subset = .self$z == 1)
    f1 <- predict(model1, data.frame(design)) #model1$fitted.values#
    mu_1 <- mean(f1)
    e_1 <- model1$resid #c(y - f1)
    y_1 <- mu_1 + e_1 %*% w1 #(w * z)/sum(w*z)
  } else if (est == "ATT") {
    y_1 <- mean(.self$y[.self$z==1])
  } else if (est == "ATC") {
    model1 <- lm(y ~ design, subset = z == 1)
    f1 <- predict(model1, data.frame(design)) #model1$fitted.values#
    mu_1 <- mean(f1)
    e_1 <- model1$resid #c(y - f1)
    y_1 <- mu_1 + e_1 %*% w1 #(w * z)/sum(w*z)
    y_0 <- mean(.self$y[.self$z==0])
  }
  
  tx_effect <- y_1 - y_0
  
  return(tx_effect)
}

outcome_model_H <- function(data, weights, ...) {
  dots <- list(...)
  est <- dots$est
  if(is.null(est)) est <- "ATE"
  design <- .data$get_x()
  ns <- data$get_n()
  n <- sum(ns)
  n1 <- ns["n1"]
  n0 <- ns["n0"]
  # design0 <- design[z==0,,drop=FALSE]
  # design1 <- design[z==1,,drop=FALSE]
  y <- data$get_y()
  z <- data$get_z()
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights$w0
  w1 <- weights$w1
  if(is.null(w1)) w1 <- rep(1/n1,n1)
  w <- rep(NA, n)
  w[z==0] <- w0
  w[z==1] <- w1
  w  <- renormalize(w)
  w1 <- renormalize(w1)
  w0 <- renormalize(w0)
  
  if(est == "ATE") {
    y_1 <- y1 %*% w1
    y_0 <- y0 %*% w0
  } else if (est == "ATT") {
    y_1 <- mean(y1)
    y_0 <- y0 %*% w0
  } else if (est == "ATC") {
    y_0 <- mean(y0)
    y_1 <- y1 %*% w1
  }
  
  tx_effect <- y_1 - y_0
  
  return(tx_effect)
}

# outcome_model <- function(.self, weights, method = c("Hajek"), doubly.robust =c(TRUE, FALSE),...) {
#   meth <- match.arg(method)
#   dr <- match.arg(doubly.robust)
#   estimator <- ifelse(dr, paste0("DR.",meth), meth)
#   fun <- switch(estimator, Hajek = outcome_model_H,
#          DR.Hajek = outcome_model_DRH)
#   return(fun(.self, weights, ...))
# }

outcome_calc <- function(data, z, weights, formula, model.fun) {
  
  w0 <- weights$w0
  w1 <- weights$w1
  
  
  fit_1 <- model.fun(formula$treated, data, subset = z==1)
  fit_0 <- model.fun(formula$control, data, subset = z==0)
  
  f_1   <- predict(fit_1, data) 
  f_0   <- predict(fit_0, data) 
  mu_1  <- mean(f_1)
  mu_0  <- mean(f_0)
  e_1   <- fit_1$resid 
  e_0   <- fit_0$resid 
  
  y_1   <- mu_1 + e_1 %*% w1
  y_0   <- mu_0 + e_0 %*% w0
  
  tx_effect <- y_1 - y_0
  
  return(c(tx_effect))
}

calc_form <- function(formula, doubly.robust, target) {
  if(!is.null(formula)){
    stopifnot(names(formula) %in% c("treated","control"))
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
    if(weights$estimate != target) stop("Weights not estimating the target")
  }
  return(weights)
}

prep_data.data.frame <- function(data,...) {
  #create data.frame
  dots <- list(...)
  tx_ind <- dots$treatment.indicator
  cov_bal <- dots$balance.covariates
  outcome <- dots$outcome
  tx.var <- if(is.character(tx_ind)) {
    match(tx_ind, colnames(data))
  } else {
    tx_ind
  }
  x.vars <- if(is.character(cov.bal)) {
    match(cov_bal, colnames(data))
  } else {
    cov_bal
  }
  y.var <- if(is.character(y.val)) {
    match(outcome, colnames(data))
  } else {
    outcome
  }
  
  df <- data[c(y.var, x.vars)]
  z <- data[[tx.var]]
  
  return(list(df = df, z = z))
}

prep_data.list <- function(data, ...) {
  #create data.frame
  df <- data.frame(data$y, data$x)
  z <- data$z
  
  return(list(df = df, z = z))
}

prep_data.DataSim <- function(data, ...) {
  #create data.frame
  df <- data.frame(data$get_y(), data$get_x())
  z <- data$get_z()
  
  return(list(df = df, z = z))
}

outcome_model <- function(data, formula = NULL, weights, 
                                  hajek = TRUE
                          , 
                                  doubly.robust = c(TRUE, FALSE), 
                                  target = c("ATT", "ATC", "ATE", "feasible"), 
                                  model = NULL, ...) {
  #get args
  dots <- list(...)
  # hajek <- match.arg(hajek)
  dr <- doubly.robust
  if(is.null(target)){
    if(inherits(weights, "causalWeights")){
      target <- weights$estimate
    }
  }
  targ <- match.arg(target)
  
  #set up model structures
  model.fun <- calc_model(model)
  formula <- calc_form(formula, dr, target)
  weights <- calc_hajek(weights, target, hajek)
  
  #setup data
  prep.data <- prep_data(data)
  
  #get estimate
  estimate <- outcome_calc(prep.data$df, prep.data$z, weights, formula, model.fun)
  
  return(estimate)
}

setGeneric("prep_data", function(data, ...) UseMethod("prep_data"))
setMethod("prep_data", signature(data = "DataSim"), prep_data.DataSim)
setMethod("prep_data", signature(data = "data.frame"), prep_data.data.frame)
setMethod("prep_data", signature(data = "list"), prep_data.list)
  