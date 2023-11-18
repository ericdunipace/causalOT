# this file creates a barycentric projection function to use as an outcome function.


#' Barycentric Projection outcome estimation
#'
#' @param formula A formula object specifying the outcome and covariates. 
#' @param data A data.frame of the data to use in the model.
#' @param weights Either a vector of weights, one for each observations, or an object of class [causalWeights][causalOT::causalWeights-class].
#' @param separate.samples.on The variable in the data denoting the treatment indicator. How to separate samples for the optimal transport calculation
#' @param penalty The penalty parameter to use in the optimal transport calculation. By default it is \eqn{1/\log(n)}{1/log(n)}.
#' @param cost_function A user supplied cost function. If supplied, must take arguments `x1`, `x2`, and `p`.
#' @param p The power to raise the cost function. Default is 2.0. For user supplied cost functions, the cost will not be raised by this power unless the user so specifies.
#' @param debias Should debiased barycentric projections be used? See details.
#' @param cost.online Should an online cost algorithm be used? Default is "auto", which selects an online cost algorithm when the sample size in each group specified by `separate.samples.on`, \eqn{n_0}{n0} and \eqn{n_1}{n1}, is such that \eqn{n_0 \cdot n_1 \geq 5000^2}{n_0 * n_1 >= 5000^2.} Must be one of "auto", "online", or "tensorized". The last of these is the offline option.
#' @param diameter The diameter of the covariate space, if known.
#' @param niter The maximum number of iterations to run the optimal transport problems
#' @param tol The tolerance for convergence of the optimal transport problems
#' @param ... Not used at this time.
#'
#' @return An object of class "bp" which is a list with slots:
#' \itemize{
#' \item `potentials` The dual potentials from calculating the optimal transport distance
#' \item `penalty` The value of the penalty parameter used in calculating the optimal transport distance
#' \item `cost_function` The cost function used to calculate the distances between units.
#' \item `cost_alg` A character vector denoting if an \eqn{L_1}{L1} distance, a squared euclidean distance, or other distance metric was used.
#' \item `p` The power to which the cost matrix was raised if not using a user supplied cost function.
#' \item `debias` Whether barycentric projections should be debiased.
#' \item `tensorized` TRUE/FALSE denoting wether to use offline cost matrices.
#' \item `data` An object of class [causalOT::dataHolder-class] with the data used to calculate the optimal transport distance.
#' \item `y_a` The outcome vector in the first sample.
#' \item `y_b` The outcome vector in the second sample.
#' \item `x_a` The covariate matrix in the first sample.
#' \item `x_b` The covariate matrix in the second sample.
#' \item `a` The empirical measure in the first sample.
#' \item `b` The empirical measure in the second sample.
#' \item `terms` The terms object from the formula.
#' }
#' @export
#' 
#' @details 
#' The barycentric projection uses the dual potentials from the optimal transport distance between the two samples to calculate projections from one sample into another. For example, in the sample of controls, we may wish to know their outcome had they been treated. In general, we then seek to minimize 
#' \deqn{\text{argmin}_{\eta} \sum_{ij} cost(\eta_i, y_j) \pi_{ij} }
#' where \eqn{\pi_{ij}} is the primal solution from the optimal transport problem.
#' 
#' These values can also be de-biased using the solutions from running an optimal transport problem of one sample against itself. Details are listed in Pooladian et al. (2022) <https://arxiv.org/abs/2202.08919>.
#'
#' @examples
#' if(torch::torch_is_installed()) {
#' set.seed(23483)
#' n <- 2^5
#' pp <- 6
#' overlap <- "low"
#' design <- "A"
#' estimate <- "ATT"
#' power <- 2
#' data <- causalOT::Hainmueller$new(n = n, p = pp,
#' design = design, overlap = overlap)
#' 
#' data$gen_data()
#' 
#' weights <- causalOT::calc_weight(x = data,
#'   z = NULL, y = NULL,
#'   estimand = estimate,
#'   method = "NNM")
#'   
#'  df <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
#'   
#'  fit <- causalOT::barycentric_projection(y ~ ., data = df, 
#'     weight = weights,
#'     separate.samples.on = "z",
#'     niter = 2)
#'  inherits(fit, "bp")
#'  }
barycentric_projection <- function(formula, data, weights, 
                                   separate.samples.on = "z", 
                                   penalty = NULL, cost_function = NULL,
                                   p = 2, 
                                   debias = FALSE, cost.online = "auto",
                                   diameter = NULL, niter = 1000L,
                                   tol = 1e-7, ...) {
  
  tx_form  <- paste0(separate.samples.on, "~ 0")
  
  if (!(separate.samples.on %in% colnames(data))) stop("must give some way of distinguishing samples in argument 'separate.samples.on'.")
  if (inherits(weights, "causalWeights")) {
    cw <- weights
    weights <- NA_real_
  } else {
    cw <- NULL
  }
  dH       <- df2dataHolder(treatment.formula = tx_form, outcome.formula = formula,
                data = data, weights = weights)
  mt       <- attr(dH, "terms")
  y_a      <- get_y0(dH)
  y_b      <- get_y1(dH)
  x_a      <- get_x0(dH)
  x_b      <- get_x1(dH)
  w        <- get_w(dH)
  z        <- get_z(dH)
  
  if (inherits(cw, "causalWeights")) {
    a      <- renormalize(cw@w0)
    b      <- renormalize(cw@w1)
    # dH@weights <- renormalize(c(a,b))[ranks(dH@z, ties = "first")] # this isn't great because it copies the whole object
  } else {
    a      <- renormalize(w[ z == 0])
    b      <- renormalize(w[ z == 1])
  }
  
  
  # solve for the potentials
  if ( missing(penalty) || is.null(penalty)) penalty <- 1/log(get_n(dH))
  stopifnot(penalty > 0)
  
  ot       <- OT$new(x = x_a, y = x_b, a = a, b = b,
                 penalty = penalty, cost_function = cost_function, p = p, 
                 debias = TRUE, tensorized = cost.online,
                 diameter=diameter)
  ot$sinkhorn_opt(niter, tol)
  
  potentials    <- ot$potentials
  cost_function <- ot$C_xy$fun
  cost_alg      <- ot$C_xy$algorithm
  p             <- ot$C_xy$p
  penalty       <- ot$penalty
  debias        <- as.logical(debias)
  tensorized    <- ot$tensorized
  
  if(cost_alg == "L1" && !tensorized) warning("With an L1 cost and `online.cost` set to 'online', you need to provide a cost function that generates the appropriate cost matrix. Otherwise predictions will generate an error.")
  
  output <- list(
    potentials    = list(f_ab = as.numeric(potentials$f_xy$to(device = "cpu")),
                         g_ba = as.numeric(potentials$g_yx$to(device = "cpu")),
                         f_aa = as.numeric(potentials$f_xx$to(device = "cpu")),
                         g_bb = as.numeric(potentials$g_yy$to(device = "cpu"))),
    penalty       = as.numeric(penalty),
    cost_function = cost_function,
    cost_alg      = cost_alg,
    p             = p,
    debias        = debias,
    tensorized    = tensorized,
    data          = dH,
    y_a           = y_a,
    y_b           = y_b,
    x_a           = x_a,
    x_b           = x_b,
    a             = a,
    b             = b,
    terms         = mt
  )
  
  class(output) <- "bp"
  
  return(output)
}

#' Predict method for barycentric projection models
#' 
#' @param object An object of class "bp"
#'
#' @param newdata a data.frame containing new observations
#' @param source.sample a vector giving the sample each observations arise from
#' @param cost_function a cost metric between observations
#' @param niter number of iterations to run the barycentric projection for powers > 2.
#' @param tol Tolerance on the optimization problem for projections with powers > 2.
#' @param ... Dots passed to the lbfgs method in the torch package.
#'
#' @export
#' 
#' @examples
#' if(torch::torch_is_installed()) {
#' set.seed(23483)
#' n <- 2^5
#' pp <- 6
#' overlap <- "low"
#' design <- "A"
#' estimate <- "ATT"
#' power <- 2
#' data <- causalOT::Hainmueller$new(n = n, p = pp,
#' design = design, overlap = overlap)
#' 
#' data$gen_data()
#' 
#' weights <- causalOT::calc_weight(x = data,
#'   z = NULL, y = NULL,
#'   estimand = estimate,
#'   method = "NNM")
#'   
#'  df <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
#'   
#'  # undebiased
#'  fit <- causalOT::barycentric_projection(y ~ ., data = df, 
#'     weight = weights,
#'     separate.samples.on = "z", niter = 2)
#'     
#'  #debiased
#'  fit_d <- causalOT::barycentric_projection(y ~ ., data = df, 
#'     weight = weights,
#'     separate.samples.on = "z", debias = TRUE, niter = 2)
#'  
#'  # predictions, without new data
#'  undebiased_predictions <- predict(fit,   source.sample = df$z)
#'  debiased_predictions   <- predict(fit_d, source.sample = df$z)
#'  
#'  isTRUE(all.equal(unname(undebiased_predictions), df$y)) # FALSE
#'  isTRUE(all.equal(unname(debiased_predictions), df$y)) # TRUE
#'  }
predict.bp <- function(object, newdata = NULL, source.sample, cost_function = NULL, niter = 1e3, tol = 1e-7, ...) {
  
  stopifnot(inherits(object, "bp"))
  
  y_a      <- object$y_a
  y_b      <- object$y_b
  x_a      <- object$x_a
  x_b      <- object$x_b
  a        <- object$a
  b        <- object$b
  a_log    <- log_weights(a)
  b_log    <- log_weights(b)
  n        <- length(a)
  m        <- length(b)
  
  # options
  debias   <- object$debias
  tensorized<- object$tensorized
  
  # change newdata to matrix
  if( missing(newdata) || is.null(newdata) ) {
    if (debias) {
      return(get_y(object$data))
    } else {
      x_new    <- get_x(object$data)
      y_new    <- get_y(object$data)
      z_source <- c("a","b")[get_z(object$data) + 1]
      z_target <- z_source #rep("both", nrow(x_new))
    }
    
  } else {
    process.newdat <- setup_new_data_bp(object,
                                        newdata, 
                                        source.sample)
    x_new    <- process.newdat$x
    y_new    <- process.newdat$y
    z_source <- process.newdat$z_s
    z_target <- process.newdat$z_t
    
  }
  
  n_new    <- nrow(x_new)
  if(all(is.na(y_new)) || is.null(y_new)) debias <- FALSE
  
  # variables for bp calc
  p        <- object$p
  cost_fun <- object$cost_function
  cost_alg <- object$cost_alg
  lambda   <- object$penalty
  
  if (!is.null(cost_function)) {
    cost_fun <- cost_function
    if(!is.character(cost_fun)) {
      tensorized <- TRUE
    } else {
      tensorized <- FALSE
    }
  }
  
  if(cost_alg == "L1" && !tensorized) stop("With an L1 cost you must either run the 'barycentric_projection' function without an online cost or you must provide a non-online version to 'predict' via the 'cost_function' argument.")
  
  # get dual potential
  f_ab     <- object$potentials$f_ab
  g_ba     <- object$potentials$g_ba
  f_aa     <- object$potentials$f_aa
  g_bb     <- object$potentials$g_bb
  
  # get bp function
  bp_fun   <- switch(cost_alg,
                   "L1" = bp_pow1,
                   "squared.euclidean" = bp_pow2,
                   bp_general)
  
  a_to_a   <- any((z_target == "a" | debias) & z_source == "a")
  b_to_a   <- any(z_target == "a" & z_source == "b")
  a_to_b   <- any(z_target == "b" & z_source == "a")
  b_to_b   <- any((z_target == "b" | debias) & z_source == "b")
  
  y_hat_a  <- rep(NA_real_, n_new)
  y_hat_b  <- rep(NA_real_, n_new)
  
  dots     <- list(...) # collect dots to prevent name conflicts
  
  # target a measure
  if ( a_to_a ) {
    if (debias) {
      z_targ_a <- ifelse(z_source == "a", "a", z_target)
    } else {
      z_targ_a <- z_target
    }
    y_hat_a <- construct_bp_est(bp_fun = bp_fun, lambda = lambda, 
                                source = "a", target = "a",
                                z_source = z_source, z_target = z_targ_a,
                                cost_function = cost_fun, p = p, 
                                x_new = x_new, x_t = x_a, y_hat_t = y_hat_a, y_t = y_a,
                                f_st = f_aa, l_target_measure = a_log,
                                tensorized = tensorized,
                                niter = niter, tol = tol, dots)
  }
  if (b_to_a) {
    y_hat_a <- construct_bp_est(bp_fun = bp_fun, lambda = lambda,  
                     source = "b", target = "a",
                     z_source = z_source, z_target = z_target,
                     cost_function = cost_fun, p = p, 
                     x_new = x_new, x_t = x_a, y_hat_t = y_hat_a, y_t = y_a,
                     f_st = f_ab, l_target_measure = a_log,
                     tensorized = tensorized,
                     niter = niter, tol = tol, dots)
  }
  
  # target b measure
  if (a_to_b) {
    y_hat_b <- construct_bp_est(bp_fun, lambda, 
                                source = "a", target = "b",
                                z_source, z_target,
                                cost_fun, p, 
                                x_new, x_b, y_hat_b, y_b,
                                g_ba, b_log,
                                tensorized, 
                                niter, tol, dots)
  }
  if (b_to_b) {
    if (debias) {
      z_targ_b <- ifelse(z_source == "b", "b", z_target)
    } else {
      z_targ_b <- z_target
    }
    y_hat_b <- construct_bp_est(bp_fun, lambda, 
                                source = "b", target = "b",
                                z_source, z_targ_b,
                                cost_fun, p, 
                                x_new, x_b, y_hat_b, y_b,
                                g_bb, b_log,
                                tensorized, 
                                niter, tol, dots)
  }
  
  # "debias" estimates using self BP
  if (debias) { # only works if has some outcome info already...
    # debiased potential towards self just returns the same outcome
    y_hat_a_debias <- ifelse(z_source == "a", y_new, y_hat_a + y_new - y_hat_b)
    y_hat_b_debias <- ifelse(z_source == "b", y_new, y_hat_b + y_new - y_hat_a)
    
    y_hat_a <- y_hat_a_debias
    y_hat_b <- y_hat_b_debias
  }
  
  # not make sense if both given...
  y_hat    <- ifelse(z_target == "a", y_hat_a, y_hat_b)
  return(y_hat)
  
}

setup_new_data_bp <- function(object, newdata, from.sample) {
  
  if( missing(from.sample)) stop("Must specify which sample the observations come from with argument 'source.sample'.")
  if (!is.vector(from.sample)) stop("'from.sample' must be a vector.")
  
  mt      <- object$terms
  tx_form <- mt$treatment
  out_form<- mt$outcome
  
  mf <- model.frame(tx_form, newdata)
  tx_target <- model.response(mf, type = "any")
  if(length(unique(tx_target)) > 2L) stop("treatment indicator must have only two levels")
  if ( is.factor(tx_target) ) {
    tx_target <- as.integer(tx_target != levels(tx_target)[1L])
  } else {
    if (!all(tx_target %in% c(0L, 1L))) {
      tx_target <- as.integer(tx_target != tx_target[1L])
    }
  }
  
  z_target   <- c("a", "b")[tx_target + 1]
  
  if( !is.null(from.sample) ) {
    newdata[[attr(mf,"names")[1]]] <- from.sample
  }
  
  dHnew   <- df2dataHolder(treatment.formula = tx_form, 
                           outcome.formula = out_form,
                           data = newdata)
  
  x_new  <- get_x(dHnew)
  y_new  <- get_y(dHnew)
  z_new <- get_z(dHnew)
  
  debias <- object$debias
  if(all(is.na(y_new))) debias <- FALSE
  
  if(all(is.na(dHnew@y)))  debias <- FALSE
  
  z_source <- c("a", "b")[z_new + 1]
  
  return(list(x = x_new, y = y_new,
              z_s = z_source, z_t = z_target))
}

construct_bp_est  <- function(bp_fun,lambda, 
                              source, target,
                              z_source, z_target,
                              cost_function, p, 
                              x_new, x_t, y_hat_t, y_t,
                              f_st,l_target_measure,
                              tensorized,
                              niter, tol, dots) {
  s_to_t_sel <- z_source == source & (z_target == target | z_target == "both")
  x_s_to_t   <- x_new[ s_to_t_sel , , drop = FALSE]
  C_s_to_t   <- cost(x_s_to_t, x_t, p = p, tensorized = tensorized, cost_function = cost_function)
  y_s_to_t   <- bp_fun(nrow(x_s_to_t), lambda, C_s_to_t, y_t,
                       f_st, l_target_measure, tensorized, 
                       niter, tol, dots)
  y_hat_t[ s_to_t_sel ] <- y_s_to_t
  return(y_hat_t)
}

bp_pow2 <- function(n, eps, C_yx, y_target, f, a_log, tensorized, niter = 1e3, tol = 1e-7, dots) {
  
  if(tensorized) {
    G <- if(inherits(f, "torch_tensor") ) {
      f/eps + a_log
    } else {
      torch::torch_tensor(f/eps + a_log,
                          dtype = C_yx$data$dtype,
                          device = C_yx$data$device)
    }
    if(!inherits(y_target, "torch_tensor")) {
      y_targ <-  torch::torch_tensor(y_target, dtype = C_yx$data$dtype, device = C_yx$data$device)
    } else {
      y_targ <- y_target
    }
    
    y_source <- torch::torch_matmul((G - C_yx$data/eps)$log_softmax(2)$exp(), y_targ)$to(device = "cpu")
    
  } else {
    G <- f/eps + a_log
    if(inherits(C_yx$data$x, "torch_tensor")) {
      x <- as.matrix(C_yx$data$x$to(device = "cpu"))
      y <- as.matrix(C_yx$data$y$to(device = "cpu"))
    } else {
      x <- as.matrix(C_yx$data$x)
      y <- as.matrix(C_yx$data$y)
    }
    
    d <- ncol(x)
    if(inherits(y_target, "torch_tensor")) {
      y_targ <- y_target$to("cpu")
    } else {
      y_targ <- y_target
    }
    
    sum_data <- list(X = as.matrix(x),
                     Y = as.matrix(y),
                     Outcome = as.numeric(y_targ),
                     G = as.numeric(G),
                     P = as.numeric(1/eps)
                     
    )
    
    if (utils::packageVersion("rkeops") >= 2.0) {
      # online_red <- rkeops::keops_kernel(
      #   formula = paste0("SumSoftMaxWeight_Reduction(G - P *", C_yx$fun,", Outcome, 1)"),
      #   args = c(
      #     paste0("X = Vi(",d,")"),
      #     paste0("Y = Vj(",d,")"),
      #     paste0("Outcome = Vj(",1,")"),
      #     "G = Vj(1)",
      #     "P = Pm(1)")
      # )
      wt_red <- rkeops::keops_kernel(
        formula = paste0("LogSumExp_Reduction(G - P *", C_yx$fun,", 1)"),
        args = c(
          paste0("X = Vi(",d,")"),
          paste0("Y = Vj(",d,")"),
          "G = Vj(1)",
          "P = Pm(1)")
      )
      
      online_red <- rkeops::keops_kernel(
        formula = paste0("Sum_Reduction(Outcome * Exp(G - P *", C_yx$fun," - Norm), 1)"),
        args = c(
          paste0("X = Vi(",d,")"),
          paste0("Y = Vj(",d,")"),
          paste0("Outcome = Vj(",1,")"),
          "G = Vj(1)",
          "P = Pm(1)",
          "Norm = Vi(1)")
      )
      
      wt_norm <- wt_red(list(X = as.matrix(x),
                             Y = as.matrix(y),
                             G = as.numeric(G),
                             P = as.numeric(1/eps)))
      
      sum_data$Norm <- wt_norm
    } else {
      wt_red <- rkeops::keops_kernel(
        formula = paste0("Max_SumShiftExp_Reduction(G - P *", C_yx$fun,", 0)"),
        args = c(
          paste0("X = Vi(",d,")"),
          paste0("Y = Vj(",d,")"),
          "G = Vj(1)",
          "P = Pm(1)")
      )
      
      online_red <- rkeops::keops_kernel(
        formula = paste0("Sum_Reduction(Outcome * Exp(G - P *", C_yx$fun," - Norm), 0)"),
        args = c(
          paste0("X = Vi(",d,")"),
          paste0("Y = Vj(",d,")"),
          paste0("Outcome = Vj(",1,")"),
          "G = Vj(1)",
          "P = Pm(1)",
          "Norm = Vi(1)")
      )
      
      wt_norm <- wt_red(list(X = as.matrix(x),
                             Y = as.matrix(y),
                             G = as.numeric(G),
                             P = as.numeric(1/eps)))
      
      sum_data$Norm <- log(wt_norm[,2]) + wt_norm[,1]
    }
    y_source <- online_red(sum_data)
  }
  
  return(as_numeric(y_source))
}

bp_pow1 <- function(n, eps, C_yx, y_target, f, a_log, tensorized, niter = 1e3, tol = 1e-7, dots) {
  if(!tensorized) stop("L1 norm must have a tensorized cost function")
  
  G <- if(inherits(f, "torch_tensor")) {
    f/eps + a_log
  } else {
    torch::torch_tensor(f/eps + a_log, device = C_yx$data$device,
                        dtype = C_yx$data$dtype)
  }
  wts <- (G - C_yx$data/eps)$log_softmax(2)$exp()$to(device = "cpu")
  if(inherits(y_target, "torch_tensor")) {
    y_targ <- as.numeric(y_target$to("cpu"))
  } else {
    y_targ <- as.numeric(y_target)
  }
  
  y_source <- apply(as.matrix(wts),1,function(w) matrixStats::weightedMedian(x=y_targ, w=c(w)))
  
  return(y_source)
}

bp_general <- function(n, eps, C_yx, y_target, f, a_log, tensorized, niter = 1e3, tol = 1e-7, dots) {
  
  fun <- tensorized_switch_generator(tensorized)
  
  if(tensorized) {
    device = C_yx$data$device
    dtype <- C_yx$data$dtype
  } else {
    if(inherits(C_yx$data$x, "torch_tensor")) {
      device = C_yx$data$x$device
      dtype <- C_yx$data$x$dtype
    } else {
      device <- cuda_device_check(NULL)
      dtype  <- cuda_dtype_check(NULL, device)
      
    }
    
  }
  
  y_source  <- torch::torch_zeros(c(n), dtype = dtype, 
                                  device = device,
                                  requires_grad = TRUE)
  
  # get and set dots args
  opt_args  <- c(list(params = y_source), dots)
  opt_cons  <- torch::optim_lbfgs
  poss_args <- names(formals(opt_cons))
  m         <- match(poss_args, 
                     names(opt_args), 
                     0L)
  # if(is.null(opt_args$line_search_fn) || opt_args$line_search_fn != "strong_wolfe") warning("line_search_fn = 'strong_wolfe' recommended for the torch LBFGS optimizer. You can pass the `line_search_fn` arg as an additional argument to this function.")
  opt       <- do.call(what = opt_cons, 
                      opt_args[m])
  torch_lbfgs_check(opt)
  y_s_old <- y_source$detach()$clone()
  
  if(!tensorized) {
    closure <- function() {
      opt$zero_grad()
      loss <- fun(y_source, y_target, eps, C_yx, f, a_log)
      loss$backward()
      return(loss)
    }
    
    for(i in 1:niter) {
      opt$step(closure)
      if(converged(y_source, y_s_old, tol)) break
      y_s_old <- y_source$detach()$clone()
    }
    
  } else {
    closure <- function() {
      opt$zero_grad()
      C_st <- C_yx$fun(y_source$view(c(n,1)), y_target$view(c(m,1)), C_yx$p)
      loss <- (C_st * wt_mat)$sum()
      loss$backward()
      return(loss)
    }
    if(!inherits(y_target, "torch_tensor")) y_target <- torch::torch_tensor(as.numeric(y_target), dtype = dtype, device = device)
    m <- length(y_target)
    wt_mat <- fun(eps, C_yx, f, a_log)
    for (i in 1:niter) {
      loss <- opt$step(closure)
      if(converged(y_source, y_s_old, tol)) break
      y_s_old <- y_source$detach()$clone()
    }
  }
  torch_cubic_reassign()
  return(as.numeric(y_source$to(device = "cpu")))
}

# functions if tensorized or not
tensorized_switch_generator <- function(tensorized) {
  switch(as.numeric(tensorized) + 1,
         torch::autograd_function(
           forward = function(ctx, y_source, y_target, eps, C_yx, f, a_log) {
             x_form <- C_yx$fun
             y_form <- sub("X", "S", x_form)
             y_form <- sub("Y", "T", y_form)
             
             G <- as_numeric(f/eps + a_log)
             x <- as_matrix(C_yx$data$x)
              
             y <- as_matrix(C_yx$data$y)
             d <- ncol(x)
             ys<- as_numeric(y_source)
             yt<- as_numeric(y_target)
             
             sum_data <- list(X = as.matrix(x),
                              Y = as.matrix(y),
                              S = as.numeric(ys),
                              T = as.numeric(yt),
                              G = as.numeric(G),
                              P = as.numeric(1.0/eps)
             )
             
             if (utils::packageVersion("rkeops") >= 2.0) {
               # online_red <- rkeops::keops_kernel(
               #   formula = paste0("SumSoftMaxWeight_Reduction(  G - P *", x_form,",", y_form,", 1)"),
               #   args = c(
               #     paste0("X = Vi(",d,")"),
               #     paste0("Y = Vj(",d,")"),
               #     paste0("S = Vi(",1,")"),
               #     paste0("T = Vj(",1,")"),
               #     "G = Vj(1)",
               #     "P = Pm(1)")
               # )
               wt_red <- rkeops::keops_kernel(
                 formula = paste0("LogSumExp_Reduction(G - P *", x_form,", 1)"),
                 args = c(
                   paste0("X = Vi(",d,")"),
                   paste0("Y = Vj(",d,")"),
                   "G = Vj(1)",
                   "P = Pm(1)")
               )
               
               online_red <- rkeops::keops_kernel(
                 formula = paste0("Sum_Reduction(", y_form,"*  Exp(G - P *", x_form,"- Norm), 1)"),
                 args = c(
                   paste0("X = Vi(",d,")"),
                   paste0("Y = Vj(",d,")"),
                   paste0("S = Vi(",1,")"),
                   paste0("T = Vj(",1,")"),
                   "G = Vj(1)",
                   "P = Pm(1)",
                   "Norm = Vi(1)")
               )
               
               wt_norm <- wt_red(list(X = as.matrix(x),
                                      Y = as.matrix(y),
                                      G = as.numeric(G),
                                      P = as.numeric(1/eps)))
               
               
               sum_data$Norm <- wt_norm
             } else {
               wt_red <- rkeops::keops_kernel(
                 formula = paste0("Max_SumShiftExp_Reduction(G - P *", x_form,", 0)"),
                 args = c(
                   paste0("X = Vi(",d,")"),
                   paste0("Y = Vj(",d,")"),
                   "G = Vj(1)",
                   "P = Pm(1)")
               )
               
               online_red <- rkeops::keops_kernel(
                 formula = paste0("Sum_Reduction(", y_form,"*  Exp(G - P *", x_form,"- Norm), 0)"),
                 args = c(
                   paste0("X = Vi(",d,")"),
                   paste0("Y = Vj(",d,")"),
                   paste0("S = Vi(",1,")"),
                   paste0("T = Vj(",1,")"),
                   "G = Vj(1)",
                   "P = Pm(1)",
                   "Norm = Vi(1)")
               )
               
               wt_norm <- wt_red(list(X = as.matrix(x),
                                      Y = as.matrix(y),
                                      G = as.numeric(G),
                                      P = as.numeric(1/eps)))
               
               
               sum_data$Norm <- log(wt_norm[,2]) + wt_norm[,1]
             }
             
             reds <- online_red(sum_data)
             # print(y_source$device)
             # print(y_source$dtype)
             ctx$save_for_backward(data = sum_data, 
                                   kernel_op = online_red,
                                   dtype = y_source$dtype,
                                   device = y_source$device)
             loss <- sum(reds)
             return(loss)
           },
           backward = function(ctx, grad_output) {
             # browser()
             s <- ctx$saved_variables
             s_grad <- rkeops::keops_grad(s$kernel_op, var="S")
             # browser()
             eta <- as.matrix(rep(as_numeric(grad_output), length(s$data$Norm)))
             grad_data <- if(utils::packageVersion("rkeops") >= 2.0) {
               c(s$data, list(varReNPP = eta))
             } else {
               c(s$data, list(eta = eta))
             }
             grad_data <- unname(grad_data)
             cpu_grad <- as_numeric(s_grad(grad_data))
             grad <- list(y_source = torch::torch_tensor(cpu_grad,
                                                         device = s$device,
                                                         dtype = s$dtype))
             # print(grad$y_source$device)
             # print(grad$y_source$dtype)
             return(grad)
             
           }),
         function(eps, C_yx, f, a_log) {
           G <- f/eps + a_log
           return( (G - C_yx$data/eps)$log_softmax(2)$exp() )
         }
  )
} 

