converged <- function(new, old, tol) {
  diff = abs(old - new)
  error = diff / abs(old + .Machine$double.eps)
  conv_check = as.logical(sum(error) < tol)
  
  return (conv_check)
}

# make vector sum to 1
renormalize <- function(x) {
  if (all(is.na(x))) return(x)
  
  if (isTRUE(any(x < 0))) {
    warning("Negative weights found! Normalizing to sum to 1 with less accurate function. Make sure negative weights make sense for your problem")
    return(x/sum(x, na.rm = TRUE))
  }
  if (isTRUE(all(x == 0)) ) return(rep(0, length(x)))
  l_x <- log_weights(x)
  return(exp(l_x - logSumExp(l_x)))
}

osqp_R6_solve <- function(model, delta, delta_idx, w = NULL, normalize = TRUE) {
  
  if (!missing(delta) && !is.null(delta)) {
    l <- model$GetData(element = "l")
    u <- model$GetData(element = "u")
    u[delta_idx] <- delta
    l[delta_idx] <- -delta
    model$Update(l = l, u = u)
  }
  if(!is.null(w) && length(w)>0)  model$WarmStart(x = w)
  
  res <- model$Solve()
  w <- res$x
  if(normalize) {
    w[w<0] <- 0
    w <- renormalize(w)
  }
  
  if (res$info$status_val != 1 && res$info$status_val != 2 ) {
    # browser()
    warning("Algorithm did not converge!!! OSQP solver message: ", res$info$status)
  }
  if (res$info$status_val == -3 || res$info$status_val == -4) {
    stop("Problem infeasible")
  }
  return(w)
}

lbfgs3c_R6_solve <- function(init, options, 
                             bounds,
                             objective,
                             gradient,
                             ...) {
  
  fit <- lbfgsb3c::lbfgsb3c(par = init,
                            fn = objective,
                            gr = gradient,
                            lower = bounds[,1],
                            upper = bounds[,2],
                            control = options,
                            ...
  )
  
  if(is.null(fit$convergence) || fit$convergence != 0) warning(fit$message)
  
  return(fit$par)
  
}

lbfgs3c_control <- function(...) {
  control <- list(...)
  control.names <- c("trace","factr",
                     "pgtol", "abstol",
                     "reltol", "lmm",
                     "maxit", "info")
  if(is.null(control$maxit) && !is.null(control$niter)) control$maxit <- control$niter
  control <- control[names(control) %in% control.names]
  
  if (length(control) == 0 || is.null(control) || !is.list(control)) {
    control <- list(trace = 0L,
                    factr = 1e7,
                    pgtol = 0,
                    abstol = 0,
                    reltol = 0,
                    lmm = 5,
                    maxit = 1000L,
                    info = FALSE
    )
  }
  
  if (is.null(control$trace))         control$trace <- 0
  if (is.null(control$factr))         control$factr <- 1e7
  if (is.null(control$pgtol))         control$pgtol <- 0
  if (is.null(control$abstol))        control$abstol <- 0
  if (is.null(control$reltol))        control$reltol <- 0
  if (is.null(control$lmm))           control$lmm <- 5
  if (is.null(control$maxit))         control$maxit <- 1000L
  if (is.null(control$info))          control$info <- FALSE
  
  control$trace <- as.numeric(control$trace)
  control$factr <- as.numeric(control$factr)
  control$pgtol <- as.numeric(control$pgtol)
  control$abstol <- as.numeric(control$abstol)
  control$reltol <- as.numeric(control$reltol)
  control$lmm   <-   as.integer(control$lmm)
  control$maxit <- as.integer(control$maxit)
  control$info   <- isTRUE(control$info)
  
  return(control)
}

arg_not_used <- function(arg) {
  return(missing(arg) || all(is.na(arg)) || is.null(arg))
}

vec_to_row_constraints <- function(rows, cols) {
  # #linear algebra formulation, slower
  # ones_cols <- Matrix::Matrix(data = 1, nrow = 1, ncol = cols)
  # diag_rows <- Matrix::Diagonal(rows, x = 1)
  # return(Matrix::kronecker(ones_cols, diag_rows))
  
  col_idx <- seq(1,rows*cols,rows)
  
  return(Matrix::sparseMatrix(i = rep(1:rows, each = cols),
                              j = c(sapply(0:(rows - 1), function(i) i + col_idx)),
                              x = rep.int(1,rows * cols),
                              dims = c(rows, rows * cols), repr = "C"))
}

vec_to_col_constraints <- function(rows, cols) {
  # ones_rows <- Matrix::Matrix(data = 1, nrow = 1, ncol = rows)
  # diag_cols <- Matrix::Diagonal(cols, x = 1)
  # return(Matrix::kronecker(ones_rows, diag_cols))
  
  
  return( Matrix::sparseMatrix(i = rep(1:cols, each = rows),
                               j = c(sapply(0:(cols - 1), function(i) i * rows + 1:rows)),
                               x = rep(1,rows * cols),
                               dims = c(cols, rows * cols), repr = "C")
          
  )
}

check_weights = function(a, x) {
  if(missing(a) || is.null(a) || all(is.na(a)) ) {
   a <- rep(1.0/nrow(x), nrow(x))
  } else {
   a <- renormalize(a)
  }
  return(a)
}

mirror_softmax <- torch::autograd_function( # for mirror descent
  forward = function(ctx, param) {
    return(param$log_softmax(1)$exp())
  },
  backward = function(ctx, grad_output) {
    # browser()
    # grad_output[1L] <- 0.0 # set first gradient to 0 so is identified
    return(list(param = grad_output))
  }
)

.cubic_interpolate <- function (x1, f1, g1, x2, f2, g2, bounds = NULL) 
{
  if (!is.null(bounds)) {
    xmin_bound <- bounds[1]
    xmax_bound <- bounds[2]
  }
  else if (x1 <= x2) {
    xmin_bound <- x1
    xmax_bound <- x2
  }
  else {
    xmin_bound <- x2
    xmax_bound <- x1
  }
  d1 <- g1$item() + g2$item() - 3 * (f1 - f2)/(x1 - x2)
  d2_square <- d1^2 - g1 * g2
  if (d2_square$item() >= 0) {
    d2 <- sqrt(d2_square)
    min_pos <- if (x1 < x2) {
      x2 - (x2 - x1) * ((g2 + d2 - d1)/(g2 - g1 + 2 * 
                                          d2))
    } else if (x1 == x2) {
      x2
    } else {
      x1 - (x1 - x2) * ((g1 + d2 - d1)/(g1 - g2 + 2 * 
                                          d2))
    }
    as.numeric(min(max(min_pos, xmin_bound), xmax_bound))
  }
  else {
    as.numeric((xmin_bound + xmax_bound)/2)
  }
}

cot_torch_bincount <- function (self, weights = list(), minlength = 0L) 
{
  args <- mget(x = c("self", "weights", "minlength"))
  args$self <- torch_sub(args$self, 1L)
  expected_types <- list(self = "Tensor", weights = "Tensor", 
                         minlength = "int64_t")
  nd_args <- "self"
  return_types <- list(list("Tensor"))
  call_c_function(fun_name = "bincount", args = args, expected_types = expected_types, 
                  nd_args = nd_args, return_types = return_types, fun_type = "namespace")
}

torch_check <- function() {
  testthat::skip_if_not_installed("torch")
  if(!torch::torch_is_installed()) {
    testthat::skip("Torch is not installed")
  }
}

R6_bootStrap <- function() {
  n <- self$n
  m <- self$m
  a <- self$a
  b <- self$b
  
  a_tilde <- rmultinom(1, n, prob = a)/n
  b_tilde <- rmultinom(1, m, prob = b)/m
  
  return(list(a = a_tilde, b = b_tilde))
}

R6_boot <- function(w_list) {
  masses <- private$bootStrap()
  a <- masses$a
  b_tilde <- masses$b
  n <- self$n
  m <- self$m
  
  means <- rep(NA_real_, length(w_list))
  for(i in seq_along(w_list)) {
    w <- w_list[i]
    stopifnot(length(w) == n)
    w_tilde <- renormalize(w *rmultinom(1, n, prob = a))
    means[i] <- self$eval(w_tilde, b_tilde)
  }
  
  return(means)
}