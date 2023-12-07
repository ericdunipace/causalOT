converged <- function(new, old, tol) {
  diff = abs(old - new)
  error = diff / abs(old + .Machine$double.eps)
  rel_conv <- as.logical(as_numeric(sum(error)) < tol)
  abs_conv <- as.logical(as_numeric(sum(diff)) < tol * tol)
  conv_check = rel_conv || abs_conv

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
  
  not_conv <- is.null(fit$convergence) || isFALSE(fit$convergence == 0) || is.na(fit$convergence)
  if(not_conv) warning(fit$message)
  
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
  d2_square <- d1^2 - g1$item() * g2$item()
  # browser()
  if (!is.nan(d2_square) && is.finite(d2_square) && d2_square >= 0) {
    d2 <- sqrt(d2_square)
    min_pos <- if(x1 == x2) {
      x2
    } else if (x1 < x2 ) {
      x2 - (x2 - x1) * ((g2$item() + d2 - d1)/(g2$item() - g1$item() + 2 * d2))
    } else {
      x1 - (x1 - x2) * ((g1$item() + d2 - d1) / (g1$item() - g2$item() + 2 * d2))
    }
    return(as.numeric(min(max(min_pos, xmin_bound), xmax_bound)))
  } else {
    return(as.numeric((xmin_bound + xmax_bound)/2))
  }
}

# cot_torch_bincount <- function (self, weights = list(), minlength = 0L) 
# {
#   args <- mget(x = c("self", "weights", "minlength"))
#   args$self <- torch_sub(args$self, 1L)
#   expected_types <- list(self = "Tensor", weights = "Tensor", 
#                          minlength = "int64_t")
#   nd_args <- "self"
#   return_types <- list(list("Tensor"))
#   call_c_function(fun_name = "bincount", args = args, expected_types = expected_types, 
#                   nd_args = nd_args, return_types = return_types, fun_type = "namespace")
# }

torch_check <- function() {
  testthat::skip_if_not_installed("torch")
  if(!torch::torch_is_installed()) {
    testthat::skip("Torch is not installed")
  }
}

rkeops_check <- function() {
  testthat::skip_if_not_installed("rkeops")
  
  if (utils::packageVersion("rkeops") >= 2.0 && rlang::is_installed("reticulate")) {
  } else if (utils::packageVersion("rkeops") < 2.0 ){
    cmake <- tryCatch(rkeops::check_cmake(system("which cmake", intern = TRUE)),
                      error = function(e) {0L}
    )
    testthat::skip("error in cmake for rkeops")
    
    # from rkeops help pages
    formula = "Sum_Reduction(Exp(-s * SqNorm2(x - y)) * b, 0)"
    
    # input arguments
    args = c("x = Vi(3)",      # vector indexed by i (of dim 3)
             "y = Vj(3)",      # vector indexed by j (of dim 3)
             "b = Vj(6)",      # vector indexed by j (of dim 6)
             "s = Pm(1)")      # parameter (scalar)
    
    # compilation of the corresponding operator
    op <- tryCatch(rkeops::keops_kernel(formula, args),
             error = function(e) {FALSE})
    
    # if an error during compilation, skip
    if(is.logical(op) && isFALSE(op)) {
      testthat::skip("error in compilation for rkeops")
    }
  }
  
}

check_weights_torch <- function(a, x, device) {
  if(missing(a) || is.null(a) || all(is.na(a))) {
    a <- rep(1.0/nrow(x), nrow(x))
  }
  stopifnot("x must be a torch_tensor" = inherits(x, "torch_tensor"))
  return(torch::torch_tensor(a, dtype = x$dtype, device = device)$contiguous())
}

check_weights_torch_torch <- function(a, x, device) {
  if(all(as.logical(a$isnan()$to(device = "cpu")))) {
    a <- rep(1.0/nrow(x), nrow(x))
  }
  stopifnot("x must be a torch_tensor" = inherits(x, "torch_tensor"))
  return(torch::torch_tensor(a, dtype = x$dtype, device = device)$contiguous())
}

check_weights_numeric <- function(a, x, ...) {
  if(missing(a) || is.null(a) || all(is.na(a)) ) {
    a <- rep(1.0/nrow(x), nrow(x))
  } 
  return(a)
}

setOldClass(c("torch_tensor","R7"))
setGeneric("check_weights", function(a, x, ...) standardGeneric("check_weights"))

setMethod("check_weights", 
          signature(a = "ANY", x = "torch_tensor"),
          check_weights_torch)

setMethod("check_weights", 
          signature(a = "torch_tensor", x = "torch_tensor"),
          check_weights_torch_torch)

setMethod("check_weights", 
          signature(a = "ANY", x = "matrix"),
          check_weights_numeric)

#based on scipy implementation
scalar_search_armijo <- function(phi, phi0, derphi0, x, dx, c1=1e-4, alpha0=1, amin=0) {
  #   phi is eval fun scalar
  #   phi0 is eval at starting values, scalar
  #   derphi0 is a scalar starting sum of gradients times (proposal - original)
  #   x is original value, vector
  #   dx are gradients, vector
  # 
  #   Minimize over alpha, the function ``phi(alpha)``.
  #   Uses the interpolation algorithm (Armijo backtracking) as suggested by
  #   Wright and Nocedal in 'Numerical Optimization', 1999, pp. 56-57
  #   alpha > 0 is assumed to be a descent direction.
  #   Returns
  #   -------
  #   alpha
  #   phi1
  
  phi_a0 = phi(x, dx, alpha0)
  if (as.logical((phi_a0 <= phi0 + c1 * alpha0 * derphi0)$to(device = "cpu"))) {
    return(list(alpha = alpha0, phi1 = phi_a0))
  }
  
  # Otherwise, compute the minimizer of a quadratic interpolant:
  alpha1 = -(derphi0) * alpha0 ^ 2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
  phi_a1 = phi(x, dx, alpha1)
  if (as.logical( (phi_a1 <= (phi0 + c1 * alpha1 * derphi0) )$to(device = "cpu"))  && as.logical((alpha1 >= 0)$to(device = "cpu"))  ) {
    return(list(alpha = alpha1, phi1 = phi_a1))
  }
  if (as.logical((alpha1 < 0)$to(device = "cpu")) ) alpha1 <- alpha0 - 0.01  #avoids the negative step size
  # Otherwise, loop with cubic interpolation until we find an alpha which
  # satisfies the first Wolfe condition (since we are backtracking, we will
  # assume that the value of alpha is not too small and satisfies the second
  # condition.
  a <- b <- alpha2 <- phi_a2 <- NULL
  while (as.logical((alpha1 > amin)$to(device = "cpu")) ) {      # we are assuming alpha>0 is a descent direction
    factor = alpha0^2 * alpha1^2 * (alpha1 - alpha0)
    if (as.logical((factor == 0)$to(device = "cpu")) ) break
    a = alpha0^2 * (phi_a1 - phi0 - derphi0 * alpha1) - 
      alpha1^2 * (phi_a0 - phi0 - derphi0 * alpha0)
    a = a / factor
    b = -alpha0^3 * (phi_a1 - phi0 - derphi0 * alpha1) + 
      alpha1^3 * (phi_a0 - phi0 - derphi0 * alpha0)
    b = b / factor
    
    alpha2 = (-b + sqrt(abs(b^2 - 3 * a * derphi0))) / (3.0 * a)
    phi_a2 = phi(x,  dx, alpha2)
    if ( as.logical((phi_a2 <= phi0 + c1 * alpha2 * derphi0)$to(device = "cpu")) && as.logical((alpha2 >= 0)$to(device = "cpu")) ) {
      return(list(alpha = alpha2, phi1 = phi_a2))
    }
    
    if (as.logical( ((alpha1 - alpha2) > alpha1 / 2.0)$to(device = "cpu")) || as.logical(((1 - alpha2/alpha1) < 0.96)$to(device = "cpu") ) ) {
      alpha2 = alpha1 / 2.0
      alpha0 = alpha1
      alpha1 = alpha2
      phi_a0 = phi_a1
      phi_a1 = phi_a2
    }
  }
  
  # ## brute force
  # phi_a3 <- NULL
  # for (alpha3 in seq(alpha0, amin, by = -0.01)) {
  #   phi_a3 = phi(x + dx * alpha3)
  #   if (phi_a3 <= (phi0 + c1 * alpha3 * derphi0) ) {
  #     return(list(alpha = alpha3, phi1 = phi_a3))
  #   }
  # }
  
  # Failed to find a suitable step length
  return(list(alpha = NULL, phi1 = phi_a1))
}

original_cubic <- get(".cubic_interpolate", envir = asNamespace("torch"))

torch_lbfgs_check <- function(opt){
  if (inherits(opt, "optim_lbfgs")) {
    cb <- .cubic_interpolate
    tmpfun <- get(".cubic_interpolate", envir = asNamespace("torch"))
    environment(cb) <- environment(tmpfun)
    attributes(cb) <- attributes(tmpfun)
    utils::assignInNamespace(".cubic_interpolate", cb, "torch" )
    
    ln_srch <- opt$defaults$line_search_fn
    no_ls <- (is.null(ln_srch) || is.na(ln_srch) || ln_srch != "strong_wolfe")
    if(no_ls ) {
      warning(" Torch's LBFGS doesn't work well without 'strong_wolfe' line search on this problem. Specify it with line_search_fn = 'strong_wolfe' in the appropriate options argument.")
    }
  }
}

torch_cubic_reassign <- function() {
  utils::assignInNamespace(".cubic_interpolate", original_cubic, "torch" )
}

lr_reduce = function(sched, loss) {
  if (is.null(sched)) return(TRUE)
  
  check <- TRUE
  
  if (inherits(sched, "lr_reduce_on_plateau")) {
    old_mode <- sched$threshold_mode
    if (as.logical(loss == sched$best) ) sched$threshold_mode <- "abs"
    
    init_lr <- sched$optimizer$defaults$lr
    lr <- tryCatch(
      sched$get_lr(),
      error = function(e) sapply(sched$optimizer$state_dict()$param_groups, function(p) p[["lr"]])
    )
    min_lr <- sched$min_lrs[[1]]
    
    # if ( sched$num_bad_epochs == sched$patience && init_lr != lr) {
    if(abs(lr - min_lr)/lr < 1e-3) {
      check <- TRUE 
    } else {
      check <- FALSE
    }
    sched$step(loss)
    sched$threshold_mode <- old_mode
  } else {
    sched$step(loss)
  }
    
  
  return(check)
}

#' @export
is.na.torch_tensor <- function(x) {
  is.nan.torch_tensor(x)
}

#' @export
is.nan.torch_tensor <- function(x) {
  as.logical(x$isnan()$to(device = "cpu"))
}

is_torch_tensor <- function(x) {
  inherits(x, "torch_tensor")
}

as_numeric <- function(x) {
  return(as.numeric(switch(is_torch_tensor(x) + 1L,
                    x,
                    x$to(device = "cpu"))))
}


as_matrix <- function(x) {
  return(as.matrix(switch(is_torch_tensor(x) + 1L,
                           x,
                           x$to(device = "cpu"))))
}

as_logical <- function(x) {
  return(as.logical(switch(is_torch_tensor(x) + 1L,
                          x,
                          x$to(device = "cpu"))))
}

get_device <- function(...) {
  args <- list(...)
  nargs <- ...length()
  device <- vector("list", nargs) |> setNames(...names())
  for(i in 1:nargs) {
   if (inherits(args[[i]], "torch_tensor"))  {
     device[[i]] <- args[[i]]$device
   } else {
     device[[i]] <- torch::torch_device("cpu")
   }
  }
  return(device)
}

get_dtype <- function(...) {
  args <- list(...)
  nargs <- ...length()
  dtype <- vector("list", nargs) |> setNames(...names())
  for(i in 1:nargs) {
    if (inherits(args[[i]], "torch_tensor"))  {
      dtype[[i]] <- args[[i]]$dtype
    } else {
      dtype[[i]] <- torch::torch_double()
    }
  }
  return(dtype)
}

cuda_device_check <- function(device) {
  if (is.null(device)) {
    cuda_opt <- torch::cuda_is_available() && torch::cuda_device_count() >= 1
    if (cuda_opt) {
      device <-  torch::torch_device("cuda")
    } else {
      device <-  torch::torch_device("cpu")
    }
  } 
  stopifnot("device argument must be NULL or an object of class 'torch_device'" = torch::is_torch_device(device))
  return(device)
}

cuda_dtype_check <- function(dtype, device = NULL) {
  #dtype
  stopifnot("device not set" = !is.null(device))
  if ( is.null(dtype) ) {
    if (grepl("cuda", capture.output(print(device)) ) ) {
      dtype <- torch::torch_float()
    } else {
      dtype <- torch::torch_double()
    }
  }
  stopifnot("Argument 'dtype' must be of class 'torch_dtype'. Please see '?torch_dtype' for more info." = torch::is_torch_dtype(dtype))
  
  return(dtype)
}

jacobian_torch <- function(vector_function, parameters) {
  stopifnot("parameters must be a list"= is.list(parameters))
  # stopifnot("losses must be a list" = is.list(losses))
  
  n_param <- length(parameters)
  n_fun   <- length(vector_function)
  nc      <- ncol(vector_function)
  stopifnot("input must be a vector not a matrix" = is.na(nc) || nc==1)
  
  deriv_list <- vector("list", n_fun)

  for (i in 1:n_fun) {
    deriv_list[[i]] <- torch::autograd_grad(vector_function[i],
                         inputs = parameters,
                         retain_graph = TRUE,
                         create_graph = TRUE,
                         allow_unused = TRUE)
  }
  # deriv_list <- lapply(1:n_fun, function(i) torch::autograd_grad(vector_function[i], 
  #                                                                           inputs = parameters,
  #                                                                           retain_graph = TRUE,
  #                                                                           create_graph = TRUE,
  #                                                                           allow_unused = TRUE))
  
  first_derivatives <- lapply(X = deriv_list, 
                              FUN = .none_deriv_to_zeros,
                              parameters = parameters
                              )
  
  return(
    torch::torch_vstack(first_derivatives)$transpose(-1,1)
  )
}

.none_deriv_to_zeros <- function(derivs, parameters) {
  
  deriv_out <- vector("list", length(parameters))
  
  for (p in seq_along(parameters) ) {
    deriv_out[[p]] <- if(derivs[[p]]$numel() == 0) {
      torch::torch_zeros_like(parameters[[p]])
    } else {
      derivs[[p]]
    }
  }
  
  return( torch::torch_hstack(deriv_out) )
}

hessian_torch <- function(loss, parameters) {
  stopifnot("parameters must be a list"= is.list(parameters))
  
  n_param <- length(parameters)
  
  param_lengths <- sapply(parameters, length)
  total_n       <- sum(param_lengths)
  param_idx     <- .get_parameter_indices(param_lengths)
  hessian_out  <- torch::torch_zeros(c(total_n, total_n),
                                     dtype = parameters[[1]]$dtype,
                                     device = parameters[[1]]$device)
  
  first_deriv <- torch::autograd_grad(loss, 
                               parameters, 
                               retain_graph = TRUE, 
                               create_graph = TRUE,
                               allow_unused = TRUE)
  
  indexes_def <- .check_undefined_grad(first_deriv)
  hessian_list <- vector("list", length(indexes_def))
  cur_deriv <- NULL
  cur_idx   <- NULL
  
  for (p in seq_along(indexes_def)) {
    cur_idx <- indexes_def[[p]]
    cur_deriv <- first_deriv[[cur_idx]]
    hessian_list[[p]] <- .second_deriv(cur_deriv, parameters)
  }
  
  used_idx <- unlist(param_idx[indexes_def])
  target_hessian <- torch::torch_vstack(hessian_list)$detach()
  hessian_out[used_idx][, used_idx] <- target_hessian
  
  return(
    hessian_out
  )
  
}

.get_parameter_indices <- function(p_lengths) {
  
  n          <- sum(p_lengths)
  idx        <- 1:n
  cur_idx    <- 1L
  cur_length <- NULL
  out_idx    <- vector("list", length(p_lengths))
  
  for (p in seq_along(p_lengths) ) {
    cur_length   <- p_lengths[[p]]
    out_idx[[p]] <- cur_idx:(cur_idx +  cur_length - 1L)
    cur_idx      <- cur_idx + cur_length 
  }
  
  return( out_idx )
  
}

.check_undefined_grad <- function(derivatives) {
  
  # is_none <- sapply(derivatives, function(p) p$numel() == 0)
  # none_idx <- which(is_none)
  
  not_none <- sapply(derivatives, function(p) p$numel() != 0)
  
  return(which(not_none))
  
}

.second_deriv <- function(first_param, parameters) {
  
  length_first <- length(first_param)
  selection_vector <- torch::torch_zeros(length_first,
                                         dtype = torch::torch_long(),
                                         device = first_param$device)
  ag_output <- NULL
  ag_vectors <- vector("list", length_first)
  
  for (i in 1:length_first) {
    selection_vector$copy_(0L)
    selection_vector[i] <- 1L

    ag_output <- torch::autograd_grad(outputs = first_param, 
                               inputs = parameters, 
                               grad_outputs = selection_vector,
                               retain_graph = TRUE, 
                               create_graph = FALSE,
                               allow_unused = TRUE)
    ag_vectors[[i]] <- .none_deriv_to_zeros(ag_output, parameters)$detach()
    
  }
  
  return(torch::torch_vstack(ag_vectors))
}

.idx_2d_to_1d <- function(idx, nrow, ncol) {
  rep_idx <- rep(idx,length(idx))
  out_idx <-  rep_idx + nrow * (rep(idx, each = length(idx)) - 1L)
  stopifnot(max(out_idx) <= nrow*ncol)
  return(out_idx)
}
  
# R6_bootStrap <- function() {
#   n <- self$n
#   m <- self$m
#   a <- self$a
#   b <- self$b
#   
#   a_tilde <- rmultinom(1, n, prob = a)/n
#   b_tilde <- rmultinom(1, m, prob = b)/m
#   
#   return(list(a = a_tilde, b = b_tilde))
# }
# 
# R6_boot <- function(w_list) {
#   masses <- private$bootStrap()
#   a <- masses$a
#   b_tilde <- masses$b
#   n <- self$n
#   m <- self$m
#   
#   means <- rep(NA_real_, length(w_list))
#   for(i in seq_along(w_list)) {
#     w <- w_list[i]
#     stopifnot(length(w) == n)
#     w_tilde <- renormalize(w *rmultinom(1, n, prob = a))
#     means[i] <- self$eval(w_tilde, b_tilde)
#   }
#   
#   return(means)
# }
