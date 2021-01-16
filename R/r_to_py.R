check_python_modules <- function(method = "auto", conda = "auto") {
  
  if(!reticulate::py_module_available("scipy")) {
    stop("Required python package not installed. Package missing: ","scipy")
    # reticulate::py_install("scipy", method = method, conda = conda)
  }
  
  if(!reticulate::py_module_available("numpy")) {
    stop("Required python package not installed. Package missing: ","numpy")
    # reticulate::py_install("numpy", method = method, conda = conda)
  }
  
  if(!reticulate::py_module_available("torch")) {
    stop("Required python package not installed. Package missing: ","torch")
    # reticulate::py_install("torch", method = method, conda = conda)
  }
  
  if(!reticulate::py_module_available("geomloss")) {
    stop("Required python package not installed. Package missing: ","geomloss")
    # reticulate::py_install("geomloss", method = method, conda = conda)
  }
  # if(!reticulate::py_module_available("utils")) {
  #   stop("Required python package not installed. Package missing: ","utils")
  #   # reticulate::py_install("utils", method = method, conda = conda)
  # }
  if(!reticulate::py_module_available("logging")) {
    stop("Required python package not installed. Package missing: ","logging")
    # reticulate::py_install("logging", method = method, conda = conda)
  }
}

ot_imputer <- function(formula, data, subset,
                       weights = NULL,
                       method = "OT", ...) {
  
  # TODO: add in ability to just return mean of target group in only provide intercept!!!!!
  method <- match.arg(method, c("OT","RR"))
  
  formula  <- as.formula(formula)
  Terms    <- terms(formula, data = data)
  varnames <- all.vars(Terms)
  
  
  if (attributes(Terms)$intercept == 1 && all(varnames == varnames[1])) {
    # return mean only...
    intercept.only = TRUE
    newform <- formula
  } else {
    newform  <- paste0(varnames[1], "~ 0 +",
                       paste0(varnames[-1], collapse = " + "))
    
    intercept.only = FALSE
  }
  
  out.dat <- model.frame(newform, data)
  
 
  
  if (is.null(weights)) weights <- rep(1/nrow(out.dat), nrow(out.dat))
  
  out <- list(data = out.dat, 
              formula = formula,
              weights = weights,
              method = method,
              intercept.only = intercept.only,
              dots = list(...))
  
  class(out) <- c("ot_imputer")
  
  return(out)
}

predict.ot_imputer <- function(object, newdata, subset,
                       weights,
                       models,
                       niter = 2000L, 
                       epsilon = NULL,
                       learning.rate = 0.01,
                       batch.size = 128L,
                       # python.path,
                       verbose = FALSE) {
  
  # if (is.null(python.path)) python.path <- "/usr/bin/python3"
  
  # reticulate::use_python(python = python.path, required = TRUE)
  
  check_python_modules()
  
  reticulate::source_python(file = ot_fun_path)
  
  if (missing(weights) ) {
    weights <- rep(1/nrow(newdata), nrow(newdata))
  } else if (is.null(weights)) {
    weights <- rep(1/nrow(newdata), nrow(newdata))
  }
  
  # Setup newdata
  Terms <- terms(object$data)
  
  y <- model.response(object$data, "numeric")
  mm_orig <- cbind(y, 
                   model.matrix(Terms, object$data))
  
  outcome.col <- attr(Terms, "response")
  m <- model.frame(delete.response(Terms), newdata)
  mm_pred <- cbind(rep(NA_real_, nrow(m)), 
                   model.matrix(delete.response(Terms), m))
  colnames(mm_pred) <- c("y", colnames(mm_orig)[-1])
  
  weights_orig     <- object$weights
  weights_new      <- weights
  weights_combined <- renormalize(c(weights_orig, 
                  weights_new))
  
  X_dat <- rbind(mm_orig, 
                 mm_pred)
  
  if (object$intercept.only ) {
    X_dat <- as.matrix(X_dat[,1,drop = FALSE])
  }
  
  nz.idx  <- which(weights_combined > 0)
  orig.idx <- 1:nrow(mm_orig)
  new.idx <-  (nrow(mm_orig) + 1):(nrow(mm_orig) + nrow(mm_pred))
  # nz.new.idx <- new.idx[new.idx %in% nz.indx]
  
  #convert to python datatypes
  np <- reticulate::import("numpy", convert = FALSE)
  X_np <- np$array(X_dat[nz.idx, , drop = FALSE])
  X_miss <- torch$DoubleTensor(X_np)
  weights_py <- torch$DoubleTensor(weights_combined[nz.idx])
  
  # setup epsilon
  if(is.null(epsilon)) {
    epsilon <- reticulate::r_to_py(pick_epsilon(X_miss))
  }
  
  report.interval <- reticulate::r_to_py(floor(niter/4))
  niter <- reticulate::r_to_py(as.integer(niter))
  batch.size <- reticulate::r_to_py(as.integer(batch.size))
  learning.rate <- reticulate::r_to_py(as.double(learning.rate))
  
  imputer <- switch(object$method,
                    "OT" = OTimputer(eps = epsilon, batchsize = batch.size, 
                                     lr = learning.rate, niter = niter),
                    "RR" = RRimputer(models = models,
                                     eps = epsilon, batchsize = batch.size, 
                                     lr = learning.rate, niter = niter))
  
  #impute
  if (verbose) {
    message("Beginning python code to do OT imputation...")
  }
  otimp = imputer$fit_transform(X = X_miss, 
                                mass = weights_py,
                                verbose = reticulate::r_to_py(isTRUE(verbose)), 
                                report_interval = report.interval, 
                                X_true = reticulate::r_to_py(NULL)
    )
  if (verbose) {
    message("Imputation complete!")
  }
  # extract imputations
  y_imp <- otimp$data$numpy()[, outcome.col]
  
  y_out <- rep(NA_real_, nrow(X_dat))
  y_out[nz.idx] <- y_imp
  y_out[orig.idx] <- y
  
  return(y_out[new.idx])
  
}

residuals.ot_imputer <- function(object, pred, ...) {
  
  y <- c(model.response(object$data, "numeric"))
  
  stopifnot(length(y) == length(pred))
  return(y - pred)
}

coef.ot_imputer <- function(object, tx.name, estimand, ...) {
  
  # if (!is.null(list(...)treatment.indicator)) {
  #   tx.col <- list(...)treatment.indicator)
  # } else {
  #   stop("need treatment indicator column to proceed")
  # }
  
  # if (is.numeric(tx.col)) {
  #   tx.col <- data.names[tx.col]
  # }
  
  Terms <- terms(object$data)
  balnames <- all.vars(delete.response(Terms))
  z  <- object$data[,tx.name]
  # sw <- get_sample_weight(sample_weight, z)
  tx.pos <- which(balnames == tx.name)
  newTerms <- drop.terms(Terms, dropx = tx.pos, keep.response = TRUE)
  
  obj0 <- obj1 <- object
  
  obj0$data <- model.frame(newTerms, object$data[ z == 0, , drop = FALSE])
  obj0$weights <- renormalize(object$weights[z == 0])
  
  obj1$data <- model.frame(newTerms, object$data[z == 1, , drop = FALSE])
  obj1$weights <- renormalize(object$weights[z == 1]) 
  
  
  
  keep <- if (estimand == "ATE" | estimand == "feasible") {
    1:nrow(object$data)
  } else if (estimand == "ATT") {
    z == 1
  } else if (estimand == "ATC") {
    z == 0
  }
  
  sw <- rep(1.0, nrow(object$data))
  
  if (estimand == "ATT") {
    sw[z == 0] <- 0
  } else if (estimand == "ATC") {
    sw[z == 1] <- 0
  }
  
  sw <- renormalize(sw[keep])
  newdata <- delete.response(object$data[keep, , drop = FALSE])
  # newweights <- rep(1/nrow(newdata), nrow(newdata))
  
  y0 <- predict.ot_imputer(obj0, newdata = newdata, weights = sw)
  y1 <- predict.ot_imputer(obj1, newdata = newdata, weights = sw)
  
  effect <- weighted.mean(x = y1 - y0, w = sw)
  
  return(c(z = effect))
  
}

setClass("ot_imputer", slots = c(data = "matrix"))
setMethod("predict", signature = c(object = "ot_imputer"),
          definition = predict.ot_imputer)
setMethod("residuals", signature = c(object = "ot_imputer"),
          definition = residuals.ot_imputer)
