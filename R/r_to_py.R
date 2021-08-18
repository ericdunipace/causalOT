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
  torch <- reticulate::import("torch", convert = TRUE)
  
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
  reticulate::py_set_seed(sample.int(.Machine$integer.max, size = 1), 
                          disable_hash_randomization = TRUE)
  
  # np$random$seed(sample.int(.Machine$integer.max, size = 1))
  torch$manual_seed(sample.int(.Machine$integer.max, size = 1))
  
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

#' Sinkhorn Loss
#'
#' @param power power of the optimal transport distance.
#' @param blur The finest level of detail that should be handled by the loss function - in order to prevent overfitting on the samples’ locations.
#' @param reach specifies the typical scale associated to the constraint strength
#' @param diameter A rough indication of the maximum distance between points, which is used to tune the espilon-scaling descent and provide a default heuristic for clustering multiscale schemes. If None, a conservative estimate will be computed on-the-fly.
#' @param scaling specifies the ratio between successive values of sigma in the epsilon-scaling descent. This parameter allows you to specify the trade-off between speed (scaling < .4) and accuracy (scaling > .9).
#' @param truncate If backend is "multiscale", specifies the effective support of a Gaussian/Laplacian kernel as a multiple of its standard deviation
#' @param metric Set the metric. One of "Lp","sdLp", or "mahalanobis".
#' @param cost specifies the cost function that should be used instead of the default
#' @param kernel 
#' @param cluster_scale If backend is "multiscale", specifies the coarse scale at which cluster centroids will be computed. If None, a conservative estimate will be computed from diameter and the ambient space’s dimension, making sure that memory overflows won’t take place.
#' @param debias specifies if we should compute the unbiased Sinkhorn divergence instead of the classic, entropy-regularized “SoftAssign” loss.
#' @param potentials When this parameter is set to True, the SamplesLoss layer returns a pair of optimal dual potentials 𝐹 F  and 𝐺 G , sampled on the input measures, instead of differentiable scalar value. These dual vectors (𝐹(𝑥𝑖)) ( F ( x i ) )  and (𝐺(𝑦𝑗)) ( G ( y j ) )  are encoded as Torch tensors, with the same shape as the input weights (𝛼𝑖) ( α i )  and (𝛽𝑗) ( β j )
#' @param verbose if backend is "multiscale", specifies whether information on the clustering and epsilon-scaling descent should be displayed in the standard output.
#' @param backend one of "auto", "tensorized", "online", or "multiscale"
#'
#' @description This function serves as an R wrapper to the Python function SamplesLoss in the 
#' GeomLoss package 
#' <http://www.kernel-operations.io/geomloss/api/pytorch-api.html?highlight=samplesloss#geomloss.SamplesLoss>
#'
#' @return a list with slots "loss", "f", "g". "loss" is the Sinkhorn distance,
#' "f" is the potential corresponding to data `x`, and "g" is the potential
#' corresponding to data `y`.
#' 
#' @details Serves as an `R` interface for the `SamplesLoss` function from the 
#' geomloss library in Python.
#' 
#' @export
#'
#' @examples
#' x <- stats::rnorm(100, 100, 10)
#' a <- rep(1/100, 100)
#' 
#' y <- stats::rnorm(50, 50, 10)
#' b <- rep(1/50, 50)
#' 
#' sink <- sinkhorn_geom(x = x, y = y, a = a, b = b, power = 2,
#'               metric = "Lp")
#' 
#' # sinkhorn distance, de-biased
#' print(sink$loss)
#' 
#' # potentials for first 5 obs in each group
#' print(sink$f[1:5])
#' print(sink$g[1:5])
sinkhorn_geom <- function(x, y, a, b, power = 2, 
                          blur = 0.05, reach = NULL, diameter = NULL,
                          scaling = 0.5, truncate = 5,
                          metric = "sdLp", kernel = NULL,
                          cluster_scale=NULL, 
                          debias=TRUE, 
                          verbose=FALSE, backend='auto', ... ) {
  
  # retrieve python packages
  np <- reticulate::import("numpy", convert = TRUE)
  torch <- reticulate::import("torch", convert = TRUE)
  geomloss <- reticulate::import("geomloss", convert = TRUE)
  # cmake <- reticulate::import("cmake", convert = TRUE)
  
  # get data dimensions
  n <- nrow(x)
  m <- nrow(y)
  d <- ncol(x)
  
  # adjust data if metric is "sdLp" or "mahalanobis"
  if (metric == "mahalanobis") {
    total <- rbind(x,y)
    U <- inv_sqrt_mat(cov(total), symmetric = TRUE)
    
    update <- (total - matrix(colMeans(total), nrow = n+m,
                                       ncol = d, byrow = TRUE)) %*% U
    
    x <- update[1:n,,drop = FALSE]
    y <- update[-(1:n),,drop = FALSE]
    
  } else if (metric == "sdLp") {
    total <- rbind(x,y)
    update <- scale(total)
    x <- update[1:n,,drop = FALSE]
    y <- update[-(1:n),,drop = FALSE]
  }
  
  # set up appropriate cost function for backend
  if (backend == "tensorized" || (n*m <= 5000^2 && backend != "multiscale" && backend != "online")) {
    if (power == 2) {
      cost <- geomloss$utils$squared_distances
    } else if (power == 1) {
      reticulate::source_python(file = lp_python_path)
      cost <- l1_loss
    } else {
      # reticulate::source_python(file = lp_python_path)
      # cost <- lp_loss
      cost <- paste0("Sum(Pow(X-Y,", power,"))")
      backend <- "online"
      power <- 1L
    }
  } else if (n*m > 5000^2 || backend == "multiscale" || backend == "online") {
    if (power == 2) {
      cost <- "SqDist(X,Y)"
    } else if (power == 1) {
      cost <- "Sum(Abs(X - Y))"
    } else {
      cost <- paste0("Sum(Pow(X-Y,", power,"))")
      power <- 1L
    }
    pykeops <- reticulate::import("pykeops", convert = TRUE)
    pykeops$clean_pykeops()
  } else {
    cost <- NULL
  }
  
  
  # sets up python data types
  xt <- torch$DoubleTensor(np$array(x))$contiguous()
  yt <- torch$DoubleTensor(np$array(y))$contiguous()
  at <- torch$DoubleTensor(a)$contiguous()
  bt <- torch$DoubleTensor(b)$contiguous()
  
  
  # sets up python function
  Loss <- geomloss$SamplesLoss("sinkhorn", p = power, blur = blur, reach = reach,
                      diameter = diameter, 
                      scaling = scaling, 
                      cost = cost, kernel = kernel,
                      cluster_scale = cluster_scale,
                      debias = debias,
                      potentials = TRUE,
                      verbose = verbose,
                      backend = backend)
  
  # run python function from geomloss library
  potentials.torch <- Loss(at, xt, bt, yt)
  
  # get potentials
  f  <- c(potentials.torch[[1]]$float()$numpy())
  g  <- c(potentials.torch[[2]]$float()$numpy())
  
  # use potentials to calculate the sinkhorn loss
  loss <- sum(f * a) + sum(g * b)
  
  # setup output matrix
  retVal <- list(f = f,
                 g = g,
                 loss = loss)
  
  return(retVal)
  
}



setClass("ot_imputer", slots = c(data = "matrix"))
setMethod("predict", signature = c(object = "ot_imputer"),
          definition = predict.ot_imputer)
setMethod("residuals", signature = c(object = "ot_imputer"),
          definition = residuals.ot_imputer)
