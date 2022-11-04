#' Estimate treatment effects
#'
#' @param causalWeights An object of class [causalWeights][causalOT::causalWeights-class]
#' @param x A dataHolder, matrix, data.frame, or object of class DataSim. See [calc_weight][causalOT::calc_weight] for more details how to input the data. If `NULL`, will use the info in the `causalWeights` argument.
#' @param y The outcome vector.
#' @param sampleWeights An optional vector of the sample weights prior to any adjustment. By default, each observation is assumed to have the same sampling weights.
#' @param model.function The modeling function to use, if desired. Must take arguments "formula", "data", and "weights". Other arguments passed via `...`, the dots.
#' @param estimate.separately Should the outcome model be estimated separately in each treatment group? TRUE or FALSE.
#' @param augmented.model Should an augmented, doubly robust estimator be used?
#' @param normalize.weights Should the weights in the `causalWeights` argument be normalized to sum to one prior to effect estimation?
#' @param ... Pass additional arguments to the outcome modeling functions.
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
#' weight <- calc_weight(data, method = "COT", options = list(lambda = 0))
#' tx_eff <- estimate_effect(data = data, weights = weight)
#' 
#' # get estimate
#' print(tx_eff$estimate)
estimate_effect <- function(causalWeights,
                            x = NULL, y = NULL,  
                            model.function,
                            estimate.separately = TRUE,
                            augmented.model = FALSE,
                            normalize.weights = TRUE,
                            ...) {
  # save the call
  mc     <- match.call()
  
  # process the data
  data <- causalWeights@data
  if(!is.null(x)) {
    data <- update_x(data, x)
  }
  if(!is.null(y)) {
    data <- update_y(data, y)
  }
  if(all(is.na(data@y))) stop("Must provide an outcome vector!")
  
  # nickname for causalWeights
  cw <- causalWeights
  
  # make sure logicals are logical
  augmented.model <- isTRUE(augmented.model)
  normalize.weights <- isTRUE(normalize.weights)
  estimate.separately <- isTRUE(estimate.separately)
  
  # normalize weights for the likleihood methods
  if (normalize.weights && cw@method %in% likelihood_methods()) {
    cw@w0 <- renormalize(cw@w0)
    cw@w1 <- renormalize(cw@w1)
  }
  
  # use an outcome model and get relevant quantities
  if (!missing(model.function) && !is.null(model.function)) {
    model.outputs <- estimate_model(data, cw, model.function,
                                    estimate.separately, ...)
    
  } else {
    model.outputs <- NULL
  }
  
  # estimate treatement effect
  res <- causalEffect(data, cw, model.outputs, augmented.model, mc)
  
  return(res)
  
}

#' causalEffect class
#'
#' @slot estimate The estimated treatment effect.
#' @slot weights The weights as an object of class [causalWeights][causalOT::causalWeights-class]
#' @slot augmentedData The data as a `data.frame` with variables `weights`, `y_obs`, `y_0`, `y_1`, `y_hat_0`, `y_hat_1`, `x`, and `z`. See details for more info.
#' @slot fit The fitted model if present. See details.
#' @slot call The call from the [estimate_effect()][causalOT::estimate_effect] function.
#' 
#' @details The variables in slot `augmentedData` are 
#' \itemize{
#' \item `weights`: The [causalWeights][causalOT::causalWeights-class] targeting the causal estimand.
#' \item `y_obs`: The vector of the observed outcomes for each observation
#' \item `y_0`: The outcome under the control condition. Missingness respects the design of the experiment. i.e., \eqn{Y(0) | Z = 1} = `NA`.
#' \item `y_hat_0`: The conditional mean outcome under the control condition. Estimated from a model.
#' \item `y_hat_1`: The conditional mean outcome under the treatment condition. Estimated from a model.
#' \item `x`: The columns denoting the covariates.
#' \item `z`: The treatment indicator.
#' }
#'
#' The slot `fit` is a list with slots `control`, `treated`, and `overall_sample`. Control and treated will be filled if `estimate.separately` is TRUE in estimate_effect. `overall_sample` will be filled if `estimate.separately` is FALSE.
#' 
#' @docType class
#' 
#' @rdname causalEffect
#'
#' @export
setClass("causalEffect", slots = c(estimate = "numeric", 
                                   weights = "causalWeights",
                                   augmentedData = "data.frame",
                                   fit = "list",
                                   call = "call"))


#' causalEffect constructor function
#'
#' @param data an object of class [dataHolder][causalOT::dataHolder-class]
#' @param causalWeights  an object of class [causalWeights][causalOT::causalWeights-class]
#' @param model.outputs Outputs of the [estimate_model()][causalOT::estimate_model] function
#' @param augmented.estimate Is the estimate to be the augmented (doubly robust) estimator? TRUE/FALSE
#' @param call the call used to calculate the treatment effects
#'
#' @return  an object of class [causalEffect][causalOT::causalEffect-class]
causalEffect <- function(data, causalWeights, model.outputs, augmented.estimate, call) {
  
  # get outcome and data length (terms from dataHolder)
  n        <- get_n(data)
  y_0_obs  <- get_y0(data)
  y_1_obs  <- get_y1(data)
  z        <- get_z(data)
  
  # causalWeights
  w0       <- causalWeights@w0
  w1       <- causalWeights@w1
  w        <- rep(NA_real_, n)
  w[ z==1] <- w1
  w[ z==0] <- w0
  w_init   <- get_w(causalWeights@data)
  
  # target estimand
  estimand <- causalWeights@estimand
  
  # if an outcome model is used, get relevant quantities
  if (!is.null(model.outputs)) {
    y_hat_0<- model.outputs$y_hat_0
    y_hat_1<- model.outputs$y_hat_1
    delta  <- y_hat_1 - y_hat_0
    
    tau    <- switch(estimand,
                     "ATT" = sum(delta[z==1] * w1), # these weights aren't adjusted
                     "ATC" = sum(delta[z==0] * w0), # these weights aren't adjusted
                     "ATE" = sum(delta       * w_init) # these weights aren't adjusted
              )
    
    # adjust outcome with the observed outcome
    if (augmented.estimate) {
      tau <- tau + sum(w1 * ( y_1_obs - y_hat_1[z==1]) ) - 
                   sum(w0 * ( y_0_obs - y_hat_0[z==0]) )
    }
    
  } else {
    tau <- sum( w1 * y_1_obs ) - sum( w0 * y_0_obs ) # these weights auto-target the estimand
    if (!is.list(model.outputs)) model.outputs <- list(model.outputs)
    
    model.outputs$fit <- list(NULL)
    
    y_hat_0<- rep(NA_real_, n)
    y_hat_1<- rep(NA_real_, n)
  }
  
  # fill out the "science"
  y_0      <- rep(NA_real_, n)
  y_1      <- rep(NA_real_, n)
  
  y_0[z==0]<- y_0_obs
  y_1[z==1]<- y_1_obs
  
  # create causalEffect object
  res <- new("causalEffect",
      estimate = tau,
      weights = causalWeights,
      augmentedData = data.frame(
        weights = renormalize(w),
        y_obs = data@y,
        y_0 = y_0,
        y_1 = y_1,
        y_hat_0 = y_hat_0,
        y_hat_1 = y_hat_1,
        z = data@z,
        data@x
      ),
      fit = model.outputs$fit,
      call = call
    )
  return(res)
  
}

estimate_model <- function(data, causalWeights, model.function,
                           separate.estimation, ...) {
  
  x       <- get_x(data)
  
  if(is.character(model.function)) {
    mod <- get(model.function)
  } else {
    mod <- model.function
  }
  
  if ( isTRUE(all.equal(model.function, barycentric_projection)) ) {
    separate.estimation <- FALSE
  }
  
  if (separate.estimation) { # fit models separately in each treatment group
    y_1_obs <- get_y1(data)
    y_0_obs <- get_y0(data)
    x_0     <- get_x0(data)
    x_1     <- get_x1(data)
    
    oform   <- formula("y ~ . ")
    environment(oform) <- list2env(list(w0=causalWeights@w0,
                                        w1=causalWeights@w1), 
                                   parent=environment(oform))
    fit_0 <- mod(formula = oform, 
                 data = data.frame(x_0, y = y_0_obs), 
                 weights = w0, ...)
    fit_1 <- mod(formula = oform, 
                 data = data.frame(x_1, y = y_1_obs), 
                 weights = w1, ...)
    
    # save models
    fit   <- list(control = fit_0, treated = fit_1, overall_sample = NULL)
    
    # save predicted mean outcomes
    y_hat_0 <- predict(object = fit_0, newdata = data.frame(x))
    y_hat_1 <- predict(object = fit_1, newdata = data.frame(x))
    
  } else { #fit model jointly on all the data
    
    # get overall data
    y <- get_y(data)
    z <- get_z(data)
    w <- rep( NA_real_, get_n(data) )
    w[ z == 0 ] <- causalWeights@w0
    w[ z == 1 ] <- causalWeights@w1
    
    oform   <- formula("y ~ . ")
    environment(oform) <- list2env(list(w = w), 
                                   parent=environment(oform))
    
    # fit model
    fit_overall <- mod(formula = oform, 
                       data = data.frame(x, z = z, y = y), 
                       weights = w, ...)
    
    # save fit
    fit   <- list(control = NULL, treated = NULL, 
                  overall_sample = fit_overall)
    
    # get predicted mean outcomes
    y_hat_0 <- predict(object = fit_overall, data.frame(x, z = 0L, y = y),
                       source.sample = z)
    y_hat_1 <- predict(object = fit_overall, data.frame(x, z = 1L, y = y),
                       source.sample = z)
    
  }
  # save outputs
  output <- list(
    y_hat_0 = y_hat_0,
    y_hat_1 = y_hat_1,
    fit = fit)
  
  return(output)
}


#' Confidence Intervals for Causal Effects
#'
#' @param object An object of class [causalEffect][causalOT::causalEffect-class]
#' @param parm Unused. Included to match forms of other `confint` functions
#' @param level Confidence level. Should be between 0 and 1. Default is 0.95.
#' @param ... Currently unused.
#'
#' @return A list with slots "CI" giving the confidence bounds and "SD" giving estimates of the standard
#' error of the causal effects.
#' 
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
#' confint(tx_eff, model = "lm",
#'     formula = list(treated = "y ~ .", control = "y ~ ."))
confint.causalEffect <- function(object, parm, level = 0.95, 
                                 ...) {
 
  if (object$options$split.model == FALSE && 
      !is.null(object$outcome.model.fit) && 
      "vcov" %in% attributes(methods(class = class(object$fit))  )$info$generic) {
    if(!is.null(object$outcome.model.fit$weights) && any(object$outcome.model.fit$weights == 0)) {
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
  
  
  semipar_var_ate <- function(y, z, yhat_1, yhat_0, w, tau, n) {
    return(mean((w * z * (y - yhat_1) - w * (1 - z) * (y - yhat_0) + (yhat_1 - yhat_0) - tau)^2) / n )
  }
  
  # pull out relevant terms from augmentedData
  
  
  # run appropriate functions
  
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

