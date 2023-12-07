#' Estimate treatment effects
#'
#' @param causalWeights An object of class [causalWeights][causalOT::causalWeights-class]
#' @param x A dataHolder, matrix, data.frame, or object of class DataSim. See [calc_weight][causalOT::calc_weight] for more details how to input the data. If `NULL`, will use the info in the `causalWeights` argument.
#' @param y The outcome vector.
#' @param model.function The modeling function to use, if desired. Must take arguments "formula", "data", and "weights". Other arguments passed via `...`, the dots.
#' @param estimate.separately Should the outcome model be estimated separately in each treatment group? TRUE or FALSE.
#' @param augment.estimate Should an augmented, doubly robust estimator be used?
#' @param normalize.weights Should the weights in the `causalWeights` argument be normalized to sum to one prior to effect estimation?
#' @param ... Pass additional arguments to the outcome modeling functions.
#'
#' @return an object of class [causalEffect][causalOT::causalEffect-class]
#' @export
#'
#' @examples
#' if ( torch::torch_is_installed() ){
#' # set-up data
#' data <- Hainmueller$new()
#' data$gen_data()
#' 
#' # calculate quantities
#' weight <- calc_weight(data, method = "COT", 
#'                       estimand = "ATT",
#'                       options = list(lambda = 0))
#' tx_eff <- estimate_effect(causalWeights = weight)
#' 
#' # get estimate
#' print(tx_eff@estimate)
#' all.equal(coef(tx_eff), c(estimate = tx_eff@estimate))
#' }
estimate_effect <- function(causalWeights,
                            x = NULL, y = NULL,  
                            model.function,
                            estimate.separately = TRUE,
                            augment.estimate = FALSE,
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
  augmented.model <- isTRUE(augment.estimate)
  normalize.weights <- isTRUE(normalize.weights)
  estimate.separately <- isTRUE(estimate.separately)
  
  # normalize weights for the likleihood methods
  if (normalize.weights && cw@method %in% likelihood_methods()) {
    cw@w0 <- renormalize(cw@w0)
    cw@w1 <- renormalize(cw@w1)
  }
  
  # use an outcome model and get relevant quantities
  if (!missing(model.function) && !is.null(model.function)) {
    model.outputs <- estimate_model(data = data, causalWeights = cw, 
                                    model.function = model.function,
                                    separate.estimation = estimate.separately,
                                    ...)
    
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
#' @slot estimand The estimand of interest
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
                                   estimand = "character",
                                   weights = "causalWeights",
                                   augmentedData = "data.frame",
                                   fit = "list",
                                   call = "call"))


#' causalEffect constructor function
#'
#' @param data an object of class [dataHolder][causalOT::dataHolder-class]
#' @param causalWeights  an object of class [causalWeights][causalOT::causalWeights-class]
#' @param model.outputs Outputs of the [estimate_model()][causalOT::estimate_model] function
#' @param augment.estimate Is the estimate to be the augmented (doubly robust) estimator? TRUE/FALSE
#' @param call the call used to calculate the treatment effects
#'
#' @return  an object of class [causalEffect][causalOT::causalEffect-class]
#' @keywords internal
causalEffect <- function(data, causalWeights, model.outputs, augment.estimate, call) {
  
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
    if (augment.estimate) {
      tau <- tau + sum(w1 * ( y_1_obs - y_hat_1[z==1]) ) - 
                   sum(w0 * ( y_0_obs - y_hat_0[z==0]) )
    }
    
  } else {
    # these weights auto-target the estimand
    mu_0 <- sum( w0 * y_0_obs )
    mu_1 <- sum( w1 * y_1_obs )
    
    # tx effect
    tau <- mu_1 - mu_0 
    if (!is.list(model.outputs)) model.outputs <- list(model.outputs)
    
    model.outputs$fit <- list(NULL)
    
    y_hat_0<- rep(mu_0, n)
    y_hat_1<- rep(mu_1, n)
  }
  
  # fill out the "science"
  y_0      <- rep(NA_real_, n)
  y_1      <- rep(NA_real_, n)
  
  y_0[z==0]<- y_0_obs
  y_1[z==1]<- y_1_obs
  
  # create causalEffect object
  res <- methods::new("causalEffect",
      estimate = tau,
      estimand = estimand,
      weights = causalWeights,
      augmentedData = data.frame(
        weights = w,
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

#' Function to estimate outcome models
#'
#' @param data A [causalOT::dataHolder()] object
#' @param causalWeights A [causalOT::causalWeights-class] object
#' @param model.function The model function passed by the user
#' @param separate.estimation TRUE or FALSE, should models be estimated separately in each group?
#' @param ... Extra agruments passed to the predict functions
#'
#' @return a list with slots `y_hat_0`, `y_hat_1`, and `fit`.
#' @keywords internal
estimate_model <- function(data, causalWeights, model.function,
                           separate.estimation, ...) {
  
  x       <- get_x(data)
  
  if(is.character(model.function)) {
    mod <- get(model.function)
  } else {
    mod <- model.function
  }
  
  if ( identical(model.function, barycentric_projection) ) {
    separate.estimation <- FALSE
  }
  
  if (separate.estimation) { # fit models separately in each treatment group
    y_1_obs <- get_y1(data)
    y_0_obs <- get_y0(data)
    x_0     <- get_x0(data)
    x_1     <- get_x1(data)
    
    oform   <- formula("y ~ . ")
    w0      <- causalWeights@w0
    w1      <- causalWeights@w1
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
    
    # save stats::predicted mean outcomes
    y_hat_0 <- stats::predict(object = fit_0, newdata = data.frame(x), ...)
    y_hat_1 <- stats::predict(object = fit_1, newdata = data.frame(x), ...)
    
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
    
    # get stats::predicted mean outcomes
    y_hat_0 <- stats::predict(object = fit_overall, data.frame(x, z = 0L, y = y),
                       source.sample = z, ...)
    y_hat_1 <- stats::predict(object = fit_overall, data.frame(x, z = 1L, y = y),
                       source.sample = z, ...)
    
  }
  # save outputs
  output <- list(
    y_hat_0 = y_hat_0,
    y_hat_1 = y_hat_1,
    fit = fit)
  
  return(output)
}

#' Extract treatment effect estimate
#'
#' @param object An object of class [causalOT::causalEffect-class]
#' @param ... Not used
#'
#' @return A number corresponding to the estimated treatment effect
#' @export
#' @method coef causalEffect
#'
#' @examples
#' # set-up data
#' set.seed(1234)
#' data <- Hainmueller$new()
#' data$gen_data()
#' 
#' # calculate quantities
#' weight <- calc_weight(data, method = "Logistic", estimand = "ATE")
#' tx_eff <- estimate_effect(causalWeights = weight)
#' 
#' all.equal(coef(tx_eff), c(estimate = tx_eff@estimate))
coef.causalEffect <- function(object, ...) {
  ret <- object@estimate
  names(ret) <- "estimate"
  return(ret)
}

#' Get the variance of a causalEffect
#'
#' @param object An object of class [causalEffect][causalOT::causalEffect-class]
#' @param ... Passed on to the sandwich estimator if there is a model fit that supports one
#'
#' @return The variance of the treatment effect as a matrix
#' @export
#' 
#' @method vcov causalEffect
#'
#' @examples
#' # set-up data
#' set.seed(1234)
#' data <- Hainmueller$new()
#' data$gen_data()
#' 
#' # calculate quantities
#' weight <- calc_weight(data, estimand = "ATT", method = "Logistic")
#' tx_eff <- estimate_effect(causalWeights = weight)
#' 
#' vcov(tx_eff)
vcov.causalEffect <- function(object, ...) {
  
  run.sandwich <- sandwich_avail_check(object)
  
  VAR <- if(run.sandwich) {
    sandwich_vars(object, ...)
  } else {
    semiparm_eff_var(object, ...)
  }
  
  VAR <- as.matrix(VAR)
  rownames(VAR) <- colnames(VAR) <- "estimate"
  
  return(VAR)
}

semiparm_eff_var <- function(object, ...) {
  
  # calculates variance of estimating y(z)
  semipar_var <- function(y, z, yhat, w, denom) {
    # estimate of expected y(z)
    e_y   <- sum(w * y * z - yhat * (w * z - 1))/denom
    
    # centered data
    resid <- w * y * z - e_y - yhat * (w * z - 1)
    
    # weighted variance
    return( sum( resid^2)/(denom - 1) )
  }
  
  semipar_var_ATT_ATC <- function(y, z, yhat_0, yhat_1, w, tau, denom, estimand) {
    delta       <- w * (y_hat_1 - y_hat_0 )
    
    delta_model <- switch(estimand,
                          "ATT" = (delta - tau) * z,
                          "ATC" = (delta - tau) * (1-z))
    
    
    # weights targeted to ATC or ATT
    resid       <-   w * (y - y_hat_1) *      z  - # treated mean estimate
                     w * (y - y_hat_0) * (1 - z) + # control mean estimate
                     delta_model
    
    # variance
    return( sum( resid^2 ) / (denom - 1) )
    
  }
  
  # ATE specific variance calculation
  semipar_var_ate <- function(y, z, yhat_1, yhat_0, w, w_init, tau, denom) {
    resid <- w * z * (y - yhat_1) - # treated model residual
             w * (1 - z) * (y - yhat_0) + # control model residual
             (yhat_1 - yhat_0) - # model stats::prediction for all data
             tau # estimate of treatment effect
    
    # weighted variance
    return( sum( w_init * resid^2 )/(denom - 1)  )
  }
  
  # pull out relevant terms from augmentedData
  data    <- object@augmentedData
  z       <- data$z
  w       <- data$weights #adjusted weights
  y       <- data$y_obs
  y_hat_0 <- data$y_hat_0
  y_hat_1 <- data$y_hat_1
  n       <- length(y)
  n1      <- sum(z)
  n0      <- n - n1
  
  # treatment effect
  tau  <- object@estimate #estimated treatment effect
  
  # get causalWeights object
  cw   <- object@weights
  
  # adjust stats::predicted outcomes if not using a model
  if (all(is.na(y_hat_0))) y_hat_0 <- rep(0.0, n)
  if (all(is.na(y_hat_1))) y_hat_1 <- rep(0.0, n)
  
  # get proper denominator
  estimand <- object@estimand
  
  denom <- switch(estimand,
                  "ATE" = n,
                  "ATT" = n1,
                  "ATC" = n0)
  
  # rescale weight vector to frequency weights
  w      <- w * denom
  
  # get sample variance estimates
  if (estimand == "ATE" ) {
    
    # rescale weight vectors to frequency weights
    w_init <- cw@data@weights
    w_init <- w_init * denom
    
    VAR_sample    <- semipar_var_ate(y, z = z, yhat_1 = y_hat_1,
                              yhat_0 = y_hat_0, w = w, 
                              w_init = w_init,
                              tau = tau,
                              denom = denom)
  } else {
    
    # total var assuming no correlation between obs
    VAR_sample <- semipar_var_ATT_ATC(y, z = z, yhat_1 = y_hat_1,
                                      yhat_0 = y_hat_0, w = w, 
                                      tau = tau,
                                      denom = denom,
                                      estimand = estimand)
  }
  
  # variance of mean estimate
  VAR <- VAR_sample/denom
  
  return(VAR)
  
}

sandwich_avail_check <- function(object) {
  fit <- object@fit$overall_sample
  return(!is.null(fit) && "vcov" %in% attributes(methods(class = class(fit))  )$info$generic)
}

sandwich_vars <- function(object, parm, level, ...) {
  fit <- object@fit$overall_sample
  vcov_mat <- if(any(fit$weights == 0)) {
    sandwich::vcovBS(fit,...)
  } else {
    tryCatch(sandwich::vcovCL(fit,...),
             error = sandwich::vcovBS(fit,...),
             warning = sandwich::vcovBS(fit,...))
    
    
  }
  VARS <- diag(vcov_mat)
  return( VARS["z"] )
}
