outcome_calc <- function(data, z, weights, formula, model.fun, matched, estimate) {
  
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
    if(estimate == "ATT" | estimate == "ATE" | estimate == "feasible") {
      tau_t <- c(mean(e_0_t) - sum(e_0_c*w0))
    }
    if(estimate == "ATC" | estimate == "ATE" | estimate == "feasible") {
      tau_c <- c(sum(e_1_t*w1) - mean(e_1_c))
    }
    if(estimate == "ATT"){
      tx_effect <- tau_t
    } else if (estimate == "ATC") {
      tx_effect <- tau_c
    } else if(estimate == "ATE" | estimate == "feasible") {
      if(estimate == "ATE") {
        t_w   <- sum(t_ind) #1/(sum(w1^2))
        c_w   <- sum(c_ind) #1/(sum(w0^2))
      } else if (estimate == "feasible") {
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

outcome_model <- function(data, formula = NULL, weights, 
                                  hajek = TRUE, 
                                  doubly.robust = TRUE,
                                  matched = FALSE,
                                  target = c("ATT", "ATC", "ATE", "feasible"), 
                                  model = NULL, ...) {
  #get args
  dots <- list(...)
  
  # hajek <- match.arg(hajek)
  dr <- isTRUE(doubly.robust[1])
  matched <- isTRUE(matched[1])
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
  estimate <- outcome_calc(prep.data$df, prep.data$z, weights, formula, model.fun,
                           matched, target)
  
  return(estimate)
}
  