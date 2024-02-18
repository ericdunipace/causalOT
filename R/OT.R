
# setOldClass("torch_tensor")
# setOldClass(c("OT","R6"))
OT <- R6::R6Class("OT",
  public = list(
    C_xy = "cost",
    C_yx = "cost",
    C_xx = "cost",
    C_yy = "cost",
    n = "integer",
    m = "integer",
    penalty = "torch_tensor",
    softmin = "function",
    debias = "logical",
    device = "torch_device",
    diameter = "torch_tensor",
    dtype = "torch_dtype",
    tensorized = "logical",
    initialize = function(x, y, a = NULL, b = NULL, penalty, 
                          cost_function = NULL, p = 2, debias = TRUE, tensorized = "auto",
                          diameter=NULL, device = NULL, dtype = NULL) {
      # browser()
      if(missing(penalty) || is.null(penalty) || is.na(penalty) || penalty < 0) {
        stop("Must specify a penalty > 0!")
      } else {
        penalty <- as.double(penalty)
      }
      
      # check if should run online version in keops
      tensorized <- private$is_tensorized(tensorized, x, y)
      
      # should setup debiased potentials
      debias <- isTRUE(debias)
      
      # device check
      if(is.null(device) || !torch::is_torch_device(device)) {
        use_cuda <- torch::cuda_is_available() && torch::cuda_device_count()>=1
        if (use_cuda) {
          self$device <-  torch::torch_device("cuda")
        } else {
          self$device <-  torch::torch_device("cpu")
        }
      } else {
        self$device <- device
        attempt_cuda <- grepl("cuda",  capture.output(print(self$device)))
        use_cuda <- attempt_cuda && torch::cuda_device_count()>=1
        if(attempt_cuda && !use_cuda) {
          warning("CUDA not available even though you tried. Switching to CPU") 
          self$device <- torch::torch_device("cpu")
        }
      }
      
      # dtype
      if (!is.null(dtype) && !torch::is_torch_dtype(dtype)) {
        warning("Provided dtype is not or class torch_dtype. Using automatic selection.")
        dtype <- NULL
      }
      if (is.null(dtype)) {
        self$dtype <- switch(as.integer(use_cuda) + 1L,
                             torch::torch_double(),
                             torch::torch_float32()
        )
      } else {
        self$dtype <- dtype
      }
      
      # if(self$dtype == torch::torch_double()) {
      #   if(packageVersion("rkeops") >= pkg_vers_number("2.0")) {
      #     rkeops::rkeops_use_float64()
      #   } else {
      #     rkeops::compile4float64()
      #   }
      # }
      
      # setup data
      if ( ! inherits(x, "torch_tensor")) {
        x <- torch::torch_tensor(x, dtype = self$dtype)$contiguous()
      }
      if ( ! inherits(y, "torch_tensor")) {
        y <- torch::torch_tensor(y, dtype = self$dtype)$contiguous()
      }
      d <- ncol(x)
      
      # setup masses
      a <- check_weights(a, x, self$device)
      b <- check_weights(b, y, self$device)
      a_log <- log_weights(a$detach())$contiguous()$to(device = self$device)
      b_log <- log_weights(b$detach())$contiguous()$to(device = self$device)
      
      if(!tensorized) self$device <- torch::torch_device("cpu") # avoids copying matrices more than needed
      
      # setup costs
      C_xy <- cost(x, y$detach(), p = p, tensorized = tensorized, cost_function = cost_function)
      C_yx <- cost(y, x$detach(), p = p, tensorized = tensorized, cost_function = cost_function)
      C_xy$to_device <- self$device
      C_yx$to_device <- self$device
      if(debias) {
        C_xx <- cost(x, x$detach(), p = p, tensorized = tensorized, cost_function = cost_function)
        C_yy <- cost(y, y$detach(), p = p, tensorized = tensorized, cost_function = cost_function)
        C_xx$to_device <- self$device
        C_yy$to_device <- self$device
      } else {
        C_xx <- C_yy <- NULL
      }
      
      
      if(tensorized) {
        softmin <- softmin_tensorized
      } else {
        if ( rkeops_installed() ) {
          if(capture.output(self$dtype) == "torch_Double") {
            if(utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
              rkeops::rkeops_use_float64()
            } else {
              rkeops::compile4float64()
            }
          } else {
            if(utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
              rkeops::rkeops_use_float32()
            } else {
              rkeops::compile4float32()
            }
            
          }
          if (utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
            reduction <- rkeops::keops_kernel(
              formula = paste0("LogSumExp_Reduction( G - P *", C_xy$fun, ", 1)"),
              args = c(
                paste0("X = Vi(",d,")"),
                paste0("Y = Vj(",d,")"),
                "G = Vj(1)",
                "P = Pm(1)")
            )
          } else {
            reduction <- rkeops::keops_kernel(
              formula = paste0("Max_SumShiftExp_Reduction( G - P *", C_xy$fun, ", 0)"),
              args = c(
                paste0("X = Vi(",d,")"),
                paste0("Y = Vj(",d,")"),
                "G = Vj(1)",
                "P = Pm(1)")
            )
          }
          
          C_xy$reduction <- reduction
          C_yx$reduction <- reduction
          if (debias) {
            C_xx$reduction <- reduction
            C_yy$reduction <- reduction
          }
          softmin <- softmin_online
        }
      }
      
      # setup diameter
      if(missing(diameter) || is.null(diameter) || is.na(diameter) || diameter < 0 || is.infinite(diameter)) {
        diameter <- private$diam_check(x$detach(), y$detach(), 
                                       p = p, tensorized = tensorized, 
                                       cost = C_xy$fun,
                                       C_yx)
      }
      
      self$C_xy = C_xy
      self$C_yx = C_yx
      self$C_xx = C_xx
      self$C_yy = C_yy
      private$a_ = a
      private$b_ = b
      self$n = length(a)
      self$m  = length(b)
      private$a_log = a_log
      private$b_log = b_log
      self$penalty = as.double(penalty)
      self$softmin = softmin
      self$diameter = diameter
      self$debias = debias
      self$tensorized = tensorized
      # private$pot = list(f_xy = torch::torch_zeros(self$n),
      #                 g_yx = torch::torch_zeros(self$m),
      #                 f_xx = torch::torch_zeros(self$n),
      #                 g_yy = torch::torch_zeros(self$m))
      private$pot = vector("list", 0)
      
      if (self$tensorized) {
        
        torch_tens <- torch::torch_zeros_like(self$b,
                                              dtype = self$dtype)
        pen <- torch::torch_tensor(self$penalty, dtype = self$dtype)
        softmin_jit <- torch::jit_compile(
          "def softmin(eps: float, C_xy: Tensor, y_potential: Tensor, b_log: Tensor):
            return -eps * (y_potential/eps + b_log - C_xy/eps).logsumexp(1)
          "
        )
          
        # assignInNamespace("softmin_jit", softmin_jit, "causalOT")
        self$softmin <- function(eps, C_xy, y_potential, b_log) {
          softmin_jit$softmin(eps, C_xy$data, y_potential, b_log)
        }
        
      }
      
      return(invisible(self))
      
    }
  ),
  active = list(
    potentials = function(value) {
      # browser()
      if(missing(value)) return(private$pot)
      if(all(is.character(names(value)))) {
        if(!all(names(value) %in% c("f_xy","g_yx", "f_xx","g_yy"))) {
          stop("Names of potential list, if given, must be in f_xy, g_yx, f_xx, g_yy.")
        }
      } 
      stopifnot(length(value) >= 2 )
      if(is.null(names(value))) {
        warning("No names given but assuming that the first two potentials are f_xy, g_yx")
        names(value)[1:2] <- c("f_xy","g_yx")
      }
      f_xy <- value$f_xy
      g_yx <- value$g_yx
      
      if(length(f_xy) != self$n) stop("f_xy must have length of a")
      if(length(g_yx) != self$m) stop("g_yx must have length of b")
      pot <- list(f_xy=f_xy, g_yx=g_yx)
      if(self$debias) {
        stopifnot(length(value) == 4)
        if(any(is.na(names(value)))) {
          warning("No names given but assuming that the last two potentials are f_xx, g_yy")
          names(value)[3:4] <- c("f_xx","g_yy")
        }
        f_xx <- value$f_xx
        g_yy <- value$g_yy
        if(length(f_xx) != self$n) stop("f_xx must have length of a")
        if(length(g_yy) != self$m) stop("g_yy must have length of b")
        pot <- c(pot, list(f_xx=f_xx, g_yy=g_yy))
      }
      private$pot <- pot
      return(invisible(self))
    },
    a = function(value) {
      if(missing(value)) return(private$a_)
      
      if(length(value) != self$n) stop("Assignment measure to mass a must have same length as original problem")
      
      if(inherits(value, "torch_tensor")) {
        private$a_ <- value$to(device = private$a_$device)
      } else {
        private$a_ <- torch::torch_tensor(value, dtype = self$dtype, device = private$a_$device)
      }
      
      private$a_log <- log_weights(private$a_$detach())$to(device = private$a_$device)
    },
    b = function(value) {
      if(missing(value)) return(private$b_)
      
      if(length(value) != self$m) stop("Assignment measure to mass a must have same length as original problem")
      
      if(inherits(value, "torch_tensor")) {
        private$b_ <- value$to(device = private$b_$device)
      } else {
        private$b_ <- torch::torch_tensor(value, dtype = private$b_$device)
      }
      private$b_log <- log_weights(private$b_$detach())$to(device = private$b_$device)
    }
  ),
  private = list(
    a_ = "torch_tensor",
    b_ = "torch_tensor",
    a_log = "torch_tensor",
    b_log = "torch_tensor",
    pot = "list",
    # softmin_jit = "script_function",
    is_tensorized = function(tensorized, x, y) {
      # browser()
      if(tensorized == "auto") {
        n <- nrow(x)
        m <- nrow(y)
        cutoff <- log(n) + log(m) >= 17.03439
        if ( rkeops_installed() && (is.na(cutoff) || isTRUE(cutoff)) ) {
          tensorized <- FALSE
        } else {
          tensorized <- TRUE
        }
      } else if (tensorized == "tensorized") {
        tensorized <- TRUE
      } else if (tensorized == "online") {
        tensorized <- FALSE
        if (!rkeops_installed()) {
          warning("Package 'rkeops' must be installed to use online cost calculations.")
          tensorized <- TRUE
        }
      } else {
        stop("tensorized must be one of auto, tensorized, or online")
      }
      return(tensorized)
    },
    diam_check = function(x, y, p, tensorized, cost_function, C_yx) {
      # browser()
      diam_cost <- cost(torch::torch_vstack(list(x$max(1)[[1]]$view(c(1,-1)), x$min(1)[[1]]$view(c(1,-1)))), 
                        torch::torch_vstack(list(y$max(1)[[1]]$view(c(1,-1)), y$min(1)[[1]]$view(c(1,-1)))), 
                        p = p, tensorized = tensorized, cost_function = cost_function)
      if (tensorized) {
        diameter <- max(diam_cost$data)$item() 
      } else {
        exp_sums <- C_yx$reduction( list( as.matrix(diam_cost$data$x$to(device = "cpu")),  as.matrix(diam_cost$data$y$to(device = "cpu")), rep(0.,2), -1.) )
        diameter <- max(exp_sums[,1])
      }
      return(diameter)
    }
  )
)


# log_weights function
# setGeneric("log_weights", function(a) standardGeneric("log_weights"))
# setMethod("log_weights", signature(a = "torch_tensor"),
# function(a) {
#   min_val <- as.double(-1e5)
#   torch::with_no_grad({
#     a_log <- torch::torch_log(a)
#     a_log[a_log < min_val] <- min_val
#   })
#   return(a_log)
# }
# ) 
# setMethod("log_weights", signature(a = "numeric"),
#           function(a) {
#             min_val <- as.double(-1e5)
#             a_log <- log(a)
#             a_log[a_log < min_val] <- min_val
#             return(a_log)
#           }
# ) 

log_weights <- function(a) UseMethod("log_weights")
log_weights.torch_tensor <- function(a) {
  min_val <- as.double(-1e5)
  torch::with_no_grad({
    a_log <- torch::torch_log(a)
    a_log[a_log < min_val] <- min_val
  })
  return(a_log)
}

log_weights.numeric <- function(a) {
  min_val <- as.double(-1e5)
  a_log <- log(a)
  a_log[a_log < min_val] <- min_val
  return(a_log)
}

softmin_tensorized <- function(eps, C_xy, y_potential, b_log) {
  return ( -eps * ( torch::torch_add(y_potential / eps, b_log) - C_xy$data / eps )$logsumexp(2) )
  # softmin_jit(eps, C_xy$data, y_potential, b_log)
}
  
# softmin_jit <- torch::jit_trace(
#     function(eps, C_xy, y_potential, b_log) {
#       return ( -eps * ( torch::torch_add(y_potential / eps, b_log) - C_xy / eps )$logsumexp(2) )
#     },
#     torch::torch_tensor(10., dtype = self$dtype),
#     torch::torch_tensor(matrix(0, 2,3), dtype = self$dtype),
#     torch::torch_tensor(rep(0,3), dtype = self$dtype),
#     torch::torch_tensor(rep(0,3), dtype = self$dtype)
#     )

softmin_online <- function(eps, C_xy, y_potential, b_log) {
  x <- C_xy$data$x
  y <- C_xy$data$y
  out <- softmin_keops(eps, x, y, 
                       y_potential$detach(), b_log$detach(), 
                       C_xy$reduction)
  
  return(out)
}

softmin_keops <- torch::autograd_function(
  forward = function(ctx, eps, x, y, y_potential, b_log, reduction) {
    xmat <- as_matrix(x)
    ymat <- as_matrix(y)
    G <- as_numeric(b_log + y_potential / eps)
    one_over_eps <- as_numeric(1.0 / eps)
    sums <- reduction( input = list(X = xmat,  
                                        Y = ymat,
                                        G = G,
                                        P = one_over_eps) 
    )
    if (rkeops_installed() && utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
      out <- torch::torch_tensor(-eps * c(sums), 
                                 dtype = b_log$dtype, device = b_log$device)
    } else {
      out <- torch::torch_tensor(-eps * (log(sums[,2]) + sums[,1]), 
                                 dtype = b_log$dtype, device = b_log$device)
    }
    
    
    ctx$save_for_backward(x = xmat,
                          y = ymat,
                          G = G,
                          one_over_eps = one_over_eps,
                          forward_op = reduction,
                          dtype = b_log$dtype,
                          device = b_log$device)
    
    return(out)
  },
  backward = function(ctx, grad_output) {
    grads <- list(x = NULL)
    saved_var <- ctx$saved_variables
    if (ctx$needs_input_grad$x) {
      cost_grad <- rkeops::keops_grad(op = saved_var$forward_op,
                                      var = "X")
      grads$x <- grad_output * 
        torch::torch_tensor(cost_grad(list(X = saved_var$x,
                                           Y = saved_var$y,
                                           G = saved_var$G,
                                           P = saved_var$one_over_eps)
        ), 
        dtype = saved_var$dtype,
        device = saved_var$device)
    }
    return(grads)
  }
)

epsilon_select <- function(diameter, eps, tot_iter, cur_iter){
  
  if (eps$item() > diameter$item()) {
    return (eps)
  }
  
  if (tot_iter  == cur_iter){
    return( eps)
  }

  eps_new = (diameter - eps) * exp(-0.9*(cur_iter - 1)) + eps
  
  return (eps_new)
}

epsilon_trajectory <- function(diameter, eps, tot_iter){
  
  if (eps > diameter) {
    return (rep(eps, tot_iter))
  }
  
  eps_new = (diameter - eps) * exp(-0.9*(1:tot_iter - 1)) + eps
  
  return (eps_new)
}

OT$set("public", 
"sinkhorn_opt",
function(niter = 1e3, tol = 1e-7) {
  
  fg_list <- private$sinkhorn_loop(niter, tol)
  
  if ( self$debias ) {
    f_xx <- private$sinkhorn_self_loop("x", niter, tol)
    g_yy <- private$sinkhorn_self_loop("y", niter, tol)
    self$potentials <- list(f_xy = fg_list$f_xy,
                            g_yx = fg_list$g_yx,
                            f_xx = f_xx,
                            g_yy = g_yy)
  } else {
    self$potentials <- list(f_xy = fg_list$f_xy,
                            g_yx = fg_list$g_yx)
  }
  return(invisible(self))
})

OT$set("public", 
       "sinkhorn_cot",
function(niter = 1e3, tol = 1e-7) {
  # torch::autograd_set_grad_mode(enabled = FALSE)
  torch::with_no_grad({
  fg_list <- private$sinkhorn_loop(niter, tol)
  
  if ( self$debias ) {
    g_yy <- private$pot$g_yy
    if(is.null(g_yy)) g_yy <- torch::torch_zeros_like(fg_list$g_yx, dtype = self$dtype)
    f_xx <- private$sinkhorn_self_loop("x", niter, tol)
    self$potentials <- list(
         f_xy = fg_list$f_xy,
         g_yx = fg_list$g_yx,
         f_xx = f_xx,
         g_yy = g_yy
    )
  } else {
    self$potentials <- list(f_xy = fg_list$f_xy,
                            g_yx = fg_list$g_yx)
  }
  })
  # torch::autograd_set_grad_mode(enabled = TRUE)
  return(invisible(self))
}
)

OT$set("private", 
"sinkhorn_self_loop",
function(which.margin = "x", niter, tol) {
  
  if (!self$debias) stop("Self potentials can only be run if debias == TRUE")
  # torch::autograd_set_grad_mode(enabled = FALSE)
  
  which.margin <- match.arg(which.margin, c("x","y"))
  if(which.margin == "x") {
    f_xx = private$pot$f_xx
    C_xx = self$C_xx
    a = private$a_$detach()
    a_log = private$a_log$detach()
  } else if (which.margin == "y") {
    f_xx = private$pot$g_yy
    C_xx = self$C_yy
    a = private$b_$detach()
    a_log = private$b_log$detach()
  } else {
    stop("Wrong margin given. Must be one of x or y")
  }
  eps      <- torch::jit_scalar(self$penalty)
  diameter <- torch::jit_scalar(self$diameter)
  softmin  <- self$softmin
  n        <- length(a)
  
  print_period <- 1L
  
  eps_log_switch = diameter/round(log(.Machine$double.xmax))
  
  torch::with_no_grad({
  if(is.null(f_xx)) {
    f_xx <- torch::torch_zeros_like(a_log, dtype = self$dtype)
    if (as.logical(eps > diameter)) {
      f_xx <- softmin(eps, C_xx, f_xx, a_log)
    } else {
      f_xx <- softmin(diameter, C_xx, f_xx, a_log)
    }
  } 
  
  loss <- loss_1 <- loss_2 <- a$dot(f_xx)$item() * 2.0
  ft_1 = f_xx$detach()$clone()
  # f_01 <- f_02 <- f_xx
  # f_1 <- f_2 <- f_xx
  
  epsilons <- epsilon_trajectory(diameter, eps, niter)
  eps_cur <- NULL
  
  for (eps_cur in epsilons) {
    # eps_cur = epsilons[i]
    # eps_cur = epsilon_select(diameter, eps, niter, i)
    # ft_1 = softmin(eps_cur, C_xx, f_2,  a_log)
    # ft_2 = softmin(eps_cur, C_xx, ft_1, a_log)
    
    if (eps_cur > eps_log_switch) {
      # Anderson acceleration
      ft_1 = softmin(torch::jit_scalar(eps_cur), C_xx, f_xx, a_log)
      f_xx$add_(ft_1)$mul_(0.5) #OT(a,a)
      #f_xx = ft_1 + f_xx; f_xx = f_xx * 0.5;  # OT(a,a)
    } else {
      f_xx = softmin(torch::jit_scalar(eps_cur), C_xx, f_xx, a_log)
    }
    # loss = sum((f_1 + f_2) * a)
    loss = a$dot(f_xx)$item() * 2.0
    # if ((i %% print_period) == 0) {
      if (abs(loss - loss_1) < tol) break
      if (abs(loss - loss_2) < tol) break
    # }
    
    loss_2 = loss_1
    loss_1 = loss
  }
  })
  # ft_1 = softmin(eps, C_xx, f_2, a_log)
  # ft_2 = softmin(eps, C_xx, f_1, a_log)
  # f_xx = 0.5 * ft_1 + 0.5 * ft_2
  # torch::autograd_set_grad_mode(enabled = TRUE) # get last step for grad if needed
  f_xx = softmin(eps, C_xx, f_xx$detach(), a_log)
  return(f_xx)
})

OT$set("private", 
       "sinkhorn_loop",
function(niter, tol) {
  # ttorch::autograd_set_grad_mode(enabled = FALSE)
  eps      <- torch::jit_scalar(self$penalty)
  diameter <- torch::jit_scalar(self$diameter)
  softmin  <- self$softmin
  a        <- private$a_$detach()
  b        <- private$b_$detach()
  a_log    <- private$a_log$detach()
  b_log    <- private$b_log$detach()
  C_xy     <- self$C_xy
  C_yx     <- self$C_yx
  n        <- self$n
  m        <- self$m
  f_xy     <- private$pot$f_xy
  g_yx     <- private$pot$g_yx
  
  print_period <- 1L
  
  eps_log_switch = diameter/round(log(.Machine$double.xmax))
  missing_pot <- is.null(f_xy) || is.null(g_yx)
  nan_f <- nan_g <- FALSE
  if(!missing_pot) {
    nan_f <- as.logical(f_xy$isnan()$any()$to(device = "cpu"))
    nan_g <- as.logical(g_yx$isnan()$any()$to(device = "cpu"))
  }
  
  torch::with_no_grad({
  if (missing_pot || nan_f || nan_g) {
    if (as.logical(eps > diameter)) {
      f_xy <- softmin(eps, 
                      C_xy, 
                      torch::torch_zeros_like(b_log, dtype = self$dtype), 
                      b_log)
      g_yx <- softmin(eps, 
                      C_yx, 
                      torch::torch_zeros_like(a_log, dtype = self$dtype), 
                      a_log)
    } else {
      f_xy <- softmin(diameter, 
                      C_xy, 
                      torch::torch_zeros_like(b_log,dtype = self$dtype), 
                      b_log)
      g_yx <- softmin(diameter, 
                      C_yx, 
                      torch::torch_zeros_like(a_log, dtype = self$dtype), 
                      a_log)
    }
  } 
  
  loss <- loss_0 <- (a$dot(f_xy) + b$dot(g_yx))$item()
  norm <- NULL
  ft = f_xy$detach()$clone()
  gt = g_yx$detach()$clone()
  epsilons <- epsilon_trajectory(diameter, eps, niter)
  eps_cur <- NULL
  
  for (e in epsilons) {
    eps_cur = torch::jit_scalar(e)
    # eps_cur = epsilons[i]
    # eps_cur = epsilon_select(diameter, eps, niter, i)
    
    if (eps_cur > eps_log_switch) {
      gt = softmin(eps_cur, C_yx, f_xy, a_log)
      ft = softmin(eps_cur, C_xy, g_yx, b_log)
      # Anderson acceleration
      #g_yx = g_yx + gt;  g_yx = g_yx * 0.5;  # OT(b,a)
      #f_xy = f_xy + ft;  f_xy = f_xy * 0.5; # OT(a,b)
      g_yx$add_(gt)$mul_(0.5);  # OT(b,a)
      f_xy$add_(ft)$mul_(0.5); # OT(a,b)
    } else {
      g_yx = softmin(eps_cur, C_yx, f_xy, a_log) # OT(b,a)
      f_xy = softmin(eps_cur, C_xy, g_yx, b_log) # OT(a,b)
    }
    
    norm <- b$dot(g_yx)
    g_yx$sub_(norm)
    f_xy$add_(norm)
    loss = ( a$dot(f_xy) )$item() #+ b$dot(g_yx))$item()
    # cat(paste0(loss$item(),", "))
    # if ( (i %% print_period) == 0 ) {
      if (abs(loss - loss_0) < tol) break
    # }
    
    loss_0 = loss
    
  }
  
  gt <- g_yx$detach()$clone()
  ft <- f_xy$detach()$clone()
  })
  # torch::autograd_set_grad_mode(enabled = TRUE)
  
  return(list(f_xy = softmin(eps, C_xy, gt, b_log), 
              g_yx = softmin(eps, C_yx, ft, a_log)))
})

OT$set("public", 
       "primal",
 function() {
   if(!self$tensorized) {
     stop("Not implemented for online calculation")
   }
   
   eps      <- self$penalty
   a_log    <- private$a_log
   b_log    <- private$b_log
   C_xy     <- self$C_xy
   n        <- self$n
   m        <- self$m
   pot      <- self$potentials
   if(is.null(pot) || length(pot) == 0) {
     stop("Must run sinkhorn optimization first")
   }
   f        <- pot$f_xy
   g        <- pot$g_yx
   
   pi <- ((f$view(c(n,1)) + g$view(c(1,m)) - C_xy$data )/eps + a_log$view(c(n,1)) + b_log$view(c(1,m)))$exp()
   
   if(self$debias) {
     C_yy     <- self$C_yy
     C_xx     <- self$C_xx
     p <- pot$f_xx
     q <- pot$g_yy
     pi_xx <- (( p$view(c(n,1)) + p$view(c(1,n)) - C_xx$data )/eps + a_log$view(c(n,1)) + a_log$view(c(1,n)) )$exp()
     pi_yy <- (( q$view(c(m,1)) + q$view(c(1,m)) - C_yy$data )/eps + b_log$view(c(m,1)) + b_log$view(c(1,m)) )$exp()
   }  else {
     pi_xx <- pi_yy <- NULL
   }
   
   output <- list(xy = pi,
                  xx = pi_xx,
                  yy = pi_yy)
   return(output)
 })

round_pi <- function(raw_pi, a, b) {
  n <- length(a)
  m <- length(b)
  
  x <- a/rowSums(raw_pi)
  x[x > 1] <- 1
  
  y <- b/colSums(raw_pi)
  y[y > 1] <- 1
  
  X <- diag(x)
  Y <- diag(y)
  
  pi_prime <-  matrix(x, n, m) *  raw_pi
  pi_2_prime <- pi_prime * matrix(y, n, m, byrow = TRUE)
  err_row <- a - rowSums(pi_2_prime)
  err_col <- b - colSums(pi_2_prime)
  
  return(pi_2_prime + matrix(err_row, n, m) * matrix(err_col,n,m,byrow=TRUE)/ sum(abs(err_row)))
  
}

sinkhorn_dist <- function(OT) {
  if(!inherits(OT, "OT")) stop("Must be an OT object")
  pot <- OT$potentials
  a   <- OT$a
  b   <- OT$b
  if(is.null(pot) || length(pot) == 0) {
    stop("Must run sinkhorn optimization first")
  }
  loss <- a$dot(pot$f_xy) + b$dot(pot$g_yx)
  
  # if(loss < 0) {
  #   raw_pi <- ot$primal()
  #   if(self$C_xy)
  #   loss <- sum(round_pi(raw_pi$xy, rep(1,self$n), rep(1,self$m)) *
  #     a$view(c(self$n,1)) * b$view(c(1,self$m)) * selfC_xy$data)
  # }
  
  if (OT$debias) {
    loss <- loss - a$dot(pot$f_xx) - b$dot(pot$g_yy)
  }
  return(loss)
}

sinkhorn_loss <- function(OT) { 
  if(!inherits(OT, "OT")) stop("Must be an OT object")
  linear_terms <- sinkhorn_dist(OT)
  
  eps <- OT$penalty
  C_xy <- OT$C_xy
  pot  <- OT$potentials
  f_xy <- pot$f_xy
  g_yx <- pot$g_yx
  
  if (OT$tensorized) {
    n <- nrow(C_xy)
    m <- ncol(C_xy)
    a_log <- log_weights(OT$a)
    b_log <- log_weights(OT$b)
    
    K_xy <- (g_yx$view(c(1,m)) + f_xy$view(c(n,1)) - C_xy$data +
               a_log$view(c(n,1)) + b_log$view(c(1,m)))/eps
    exponential_terms <-  -eps * K_xy$view(-1)$logsumexp(1)
    if (OT$debias) {
      C_yy <- OT$C_yy$data
      C_xx <- OT$C_xx$data
      
      f_xx <- pot$f_xx
      g_yy <- pot$g_yy
      
      K_xx <- (f_xx$view(c(1,n)) + f_xx - C_xx + a_log$view(c(n,1)) + a_log$view(c(1,n)))/eps
      K_yy <- (g_yy$view(c(1,m)) + g_yy - C_yy + b_log$view(c(m,1)) + b_log$view(c(1,m)))/eps
      
      exponential_terms <- exponential_terms + 
        0.5 * eps * K_xx$view(-1)$logsumexp(1) +
        0.5 * eps * K_yy$view(-1)$logsumexp(1)
    }
  } else {
    x <- C_xy$data$x
    y <- C_xy$data$y
    a_log <- OT$a_log
    b_log <- OT$b_log
    
    n <- nrow(C_xy)
    m <- ncol(C_xy)
    
    exp_sums <- C_xy$reduction( list(C_xy$data$x,  C_xy$data$y, 
                          as.numeric(b_log), as.numeric(g_yx), 
                          1.0 / eps) )
    l_x <- torch::torch_tensor(log(exp_sums[,2]) + exp_sums[,1])
    exponential_terms <-  -eps * (l_x + a_log + f_xy/eps)$logsumexp(1)$exp()
    
    if (OT$debias) {
      f_xx <- pot$f_xx
      g_yy <- pot$g_yy
      
      exp_sums_xx <- C_xy$reduction( list(C_xy$data$x,  C_xy$data$x, 
                            as.numeric(a_log), as.numeric(f_xx), 
                            1.0 / eps))
      l_xx <- torch::torch_tensor(log(exp_sums_xx[,2]) + exp_sums_xx[,1])
      
      exp_sums_yy <- C_xy$reduction( list (C_xy$data$y,  C_xy$data$y, 
                               as.numeric(b_log), as.numeric(g_yy), 
                               1.0 / eps) )
      l_yy <- torch::torch_tensor(log(exp_sums_yy[,2]) + exp_sums_yy[,1])
      exponential_terms <- exponential_terms + 
        (l_xx + a_log + f_xx/eps)$logsumexp(1)$exp() * 0.5 * eps +
        (l_yy + b_log + g_yy/eps)$logsumexp(1)$exp() * 0.5 * eps
    }
  }
  
  loss <- linear_terms + exponential_terms
  return(loss)
}

energy_dist <- function(OT) {
  if(!inherits(OT, "OT")) stop("Must be an OT object")
  if (!OT$debias) {
    stop("Must have option debias set to TRUE for energy distance")
  }
  if (OT$tensorized) {
    a <- OT$a
    b <- OT$b
    loss_cross <- if (OT$C_yx$data$requires_grad) {
      b$dot(OT$C_yx$data$matmul(a))
    } else {
      a$dot(OT$C_xy$data$matmul(b))
    }
    loss <- loss_cross - 
      a$dot(OT$C_xx$data$matmul(a)) * 0.5 -
      b$dot(OT$C_yy$data$matmul(b)) * 0.5
  } else {
    loss <- energy_dist_online(OT$C_xy$data$x,
                               OT$C_xy$data$y,
                               OT$a,
                               OT$b,
                               OT$C_xy$fun)
    
  }
  
  return(loss)
}

energy_dist_online <- torch::autograd_function(
  forward = function(ctx, x, y, a, b, formula) {
    xmat <- as_matrix(x)
    ymat <- as_matrix(y)
    a_vec <- as_numeric(a)
    b_vec <- as_numeric(b)
    
    d <- ncol(xmat)
    
    device <- get_device(x = x, y = y, a = a, b = b)
    dtype <- get_dtype(x = x, y = y, a = a, b = b)
    if (rkeops_installed() && utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
      sumred <- rkeops::keops_kernel(
        formula = paste0("Sum_Reduction( B* ", formula, ", 1)"),
        args = c(
          paste0("X = Vi(",d,")"),
          paste0("Y = Vj(",d,")"),
          "B = Vj(1)")
      )
    } else {
      sumred <- rkeops::keops_kernel(
        formula = paste0("Sum_Reduction( B* ", formula, ", 0)"),
        args = c(
          paste0("X = Vi(",d,")"),
          paste0("Y = Vj(",d,")"),
          "B = Vj(1)")
      )
    }
    
    a_cross_deriv <- sumred(list(xmat,ymat, b_vec))
    # b_cross_deriv <- sumred(list(y,x, a))
    a_self_deriv <- sumred(list(xmat,xmat, a_vec))
    b_self_deriv <- sumred(list(ymat,ymat, b_vec))
    loss <- sum(a_vec * a_cross_deriv) - 
      0.5 * sum(a_vec * a_self_deriv) -
      0.5 * sum(b_vec * b_self_deriv)
    
    ctx$save_for_backward(a_deriv = c(a_cross_deriv - a_self_deriv),
                          b_deriv = c(b_self_deriv),
                          a = a_vec,
                          x = xmat,
                          b = b_vec,
                          y = ymat,
                          forward_op = sumred,
                          dtype = dtype,
                          device = device)
    return(torch::torch_tensor(loss, dtype = dtype$a,
                               device = device$a))
  },
  backward = function(ctx, grad_output) {
    grads <- list(x = NULL,
                  y = NULL,
                  a = NULL,
                  b = NULL)
    sv    <- ctx$saved_variables
    go    <- as_numeric(grad_output)
    if (ctx$needs_input_grad$a) {
      grads$a <- torch::torch_tensor(go * sv$a_deriv,
      device = sv$device$a,
      dtype = sv$dtype$a)
    }
    
    if (ctx$needs_input_grad$b) {
      grads$b <- torch::torch_tensor(go * 
                          c(sv$sumred(list(sv$y,sv$x, sv$a)) - sv$b_deriv),
                                     device = sv$device$b,
                                     dtype = sv$dtype$b)
    }
    
    if (ctx$needs_input_grad$x) {
      
      cost_grad_xy <- rkeops::keops_grad(op = sv$forward_op,
                                           var = "X")
      
      
      grad_data <- list(X = sv$x,
                        Y = sv$y,
                        B = sv$b,
                        eta = matrix(sv$a))
      grad_data2 <- list(X = sv$x,
                         Y = sv$x,
                         B = sv$a,
                         eta = matrix(sv$a))
      if(rkeops_installed() && utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
        grad_data <- unname(grad_data)
        grad_data2 <- unname(grad_data2)
      } 
      grads$x <- cost_grad_xy(grad_data)
      
      cost_grad_xx <- rkeops::keops_grad(op = sv$forward_op,
                                         var = "X")
      
      
      grads$x <- torch::torch_tensor(go * c(grads$x - cost_grad_xx(grad_data2)),
      device = sv$device$x,
      dtype = sv$dtype$x)
    }
    if (ctx$needs_input_grad$y) {
      
      cost_grad <- rkeops::keops_grad(op = sv$forward_op,
                                      var = "X")
      
      grad_data <- list(X = sv$y,
                        Y = sv$x,
                        B = sv$a,
                        eta = matrix(sv$b))
      
      grad_data2 <- list(X = sv$y,
                        Y = sv$y,
                        B = sv$b,
                        eta = matrix(sv$b))
      
      if(rkeops_installed() && utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
        grad_data <- unname(grad_data)
        grad_data2 <- unname(grad_data2)
      } 
      
      grads$y <-  cost_grad(grad_data)
      
      cost_grad_yy <- rkeops::keops_grad(op = sv$forward_op,
                                      var = "X")
      grads$y <- torch::torch_tensor(go *c(grads$y - cost_grad_yy(grad_data2)),
      device = sv$device$y,
      dtype = sv$dtype$y)
    }
    return(grads)
  }
)


inf_sinkhorn_dist <- function(OT) {
  if(!inherits(OT, "OT")) stop("Must be an OT object")
  if(OT$debias) return(energy_dist(OT))
  if (OT$tensorized) {
    a <- OT$a
    b <- OT$b
    loss <- if (OT$C_yx$data$requires_grad) {
      b$dot(OT$C_yx$data$matmul(a))
    } else {
      a$dot(OT$C_xy$data$matmul(b))
    }
  } else {
    loss <- inf_sinkhorn_online(x = OT$C_xy$data$x,
                                y = OT$C_xy$data$y,
                                a = OT$a,
                                b = OT$b,
                                formula = OT$C_xy$fun)
    
  }
  
  return(loss)
}

inf_sinkhorn_online <- torch::autograd_function(
  forward = function(ctx, x, y, a, b, formula) {
    
    device <- get_device(x = x, y = y, a = a, b = b)
    dtype <- get_dtype(x = x, y = y, a = a, b = b)
    
    x_mat <- as_matrix(x)
    y_mat <- as_matrix(y)
    a_vec <- as_numeric(a)
    b_vec <- as_numeric(b)
    
    d <- ncol(x_mat)
    
    if(utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
      sumred <- rkeops::keops_kernel(
        formula = paste0("Sum_Reduction( B* ", formula, ", 1)"),
        args = c(
          paste0("X = Vi(",d,")"),
          paste0("Y = Vj(",d,")"),
          "B = Vj(1)")
      )
    } else {
      sumred <- rkeops::keops_kernel(
        formula = paste0("Sum_Reduction( B* ", formula, ", 0)"),
        args = c(
          paste0("X = Vi(",d,")"),
          paste0("Y = Vj(",d,")"),
          "B = Vj(1)")
      )
    }
    
    a_deriv <- sumred(list(x_mat,y_mat, b_vec))
    # b_deriv <- sumred(list(y,x, a))
    loss <- sum(a_vec * a_deriv)
    
    ctx$save_for_backward(a_deriv = c(a_deriv),
                          # b_deriv = c(b_deriv),
                          a = a_vec,
                          b = b_vec,
                          x = x_mat,
                          y = y_mat,
                          forward_op = sumred,
                          device = device,
                          dtype = dtype)
    return(torch::torch_tensor(loss,
                               dtype = dtype$a,
                               device = device$a))
  },
  backward = function(ctx, grad_output) {
    grads <- list(x = NULL,
                  y = NULL,
                  a = NULL,
                  b = NULL)
    saved_var <- ctx$saved_variables
    go <- as_numeric(grad_output)
    if (ctx$needs_input_grad$a) {
      grads$a <- 
        torch::torch_tensor(go *saved_var$a_deriv, 
                             dtype = saved_var$dtype$a,
                             device = saved_var$device$a)
    }
    if (ctx$needs_input_grad$b) {
      grads$b <- 
        torch::torch_tensor(go * saved_var$forwar_op(list(saved_var$y,
                                                     saved_var$x, 
                                                     saved_var$a)), 
                            dtype = saved_var$dtype$b,
                            device = saved_var$device$b)
    }
    if (ctx$needs_input_grad$x) {
      cost_grad <- rkeops::keops_grad(op = saved_var$forward_op,
                                      var = "X")
      grad_data <- list(X = saved_var$x,
                        Y = saved_var$y,
                        B = saved_var$b,
                        eta = matrix(saved_var$a))
      if(utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
        grad_data <- unname(grad_data)
      } 
      grads$x <- 
        torch::torch_tensor(go * cost_grad(grad_data), 
                  dtype = saved_var$dtype$x,
                  device = saved_var$device$x)
    }
    if (ctx$needs_input_grad$y) {
      cost_grad <- rkeops::keops_grad(op = saved_var$forward_op,
                                      var = "X")
      grad_data <- list(X = saved_var$y,
                        Y = saved_var$x,
                        B = saved_var$a,
                        eta =  matrix(saved_var$b))
      
      if(rkeops_installed() && utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
        grad_data <- unname(grad_data)
      } 
      
      grads$x <- 
        torch::torch_tensor(go * cost_grad(grad_data),
            dtype = saved_var$dtype$y,
            device = saved_var$device$y)
    }
    return(grads)
  }
)

# semi_dual <- function(OT,pot,debias = FALSE) {
#   
#   a <- OT$a
#   a_log <- log(a)
#   eps <- OT$penalty
#   softmin <- OT$softmin
#   
#   if(debias) {
#     #softmin has neg sign already
#     loss <- pot$dot(a) + a$dot(softmin(eps, OT$C_xx, pot, a_log))
#   } else {
#     #softmin has neg sign already
#     b <- OT$b
#     loss <- pot$dot(a) + b$dot(softmin(eps, OT$C_yx, f, a_log))
#   }
#   return(loss)
# }

transportationMatrix <- function(x = NULL, z = NULL, weights = NULL, 
                                 lambda = NULL, p = 2, 
                                 cost = NULL,
                                 debias = FALSE,
                                 diameter=NULL,
                                 niter = 1000, tol = 1e-7) {
  
  #sets up the attribute names
  tm_attr_names <- c("dual", "lambda", "debias",
                     "estimand", "p"
  ) 
  # sets up a blank matrix to plug in if needed
  blank_mat <- list(xy = matrix(0, 0, 0),
                    xx = matrix(0, 0, 0),
                    yy = matrix(0, 0, 0))
  attributes(blank_mat) <- list(dual = numeric(0),
                                lambda = NA_real_,
                                debias = logical(0),
                                estimand = NULL,
                                p = NA_real_)
  if(is.null(weights)) { #escape hatch
    
    tmat <- list(w0 = blank_mat,
                 w1 = blank_mat)
    class(tmat) <- c("transportationMatix", class(tmat))
    return(tmat)
  }
  
  stopifnot(inherits(weights, "causalWeights"))
  
  cw <- weights
  dh <- dataHolder(x = x, z = z)
  n <- get_n(dh)
  
  if(is.null(lambda) || is.na(lambda)) {
    lambda <- 1/log(n)
  }
  stopifnot(lambda > 0.0)
  
  if( cw@estimand == "ATE" ) {
    x0 <- get_x0(dh)
    x1 <- get_x1(dh)
    y <- get_x(dh)
    b <- dh@weights
    
    ot0 <- OT$new(x = x0, y = y, 
                  a = cw@w0, b = b,
                  penalty = lambda, 
                  cost = cost, p = p, debias = debias, 
                  tensorized = "tensorized",
                  diameter = diameter)
    ot1 <- OT$new(x = x1, y = y, 
                  a = cw@w1, b = b,
                  penalty = lambda, 
                  cost = cost, p = p, debias = debias, 
                  tensorized = "tensorized",
                  diameter = diameter)
    
    ot0$sinkhorn_opt(niter, tol)
    ot1$sinkhorn_opt(niter, tol)
    
    pot0 <- ot0$potentials
    f0 <- pot0$f_xy
    g0 <- pot0$g_yx
    
    pot1 <- ot1$potentials
    f1 <- pot1$f_xy
    g1 <- pot1$g_yx
    
    pi_0 <- ot0$primal()
    pi_1 <- ot1$primal()
    
    tmat <- list(w0 = pi_0,
                 w1 = pi_1)
    
    attributes(tmat$w0$xy)[tm_attr_names] <- 
      list(
        dual = list(f = as.numeric(f0), g = as.numeric(g0)),
        lambda = lambda,
        debias = debias,
        estimand = "ATE",
        p = p
      )
    
    attributes(tmat$w1$xy)[tm_attr_names] <- 
      list(
        dual = list(f = as.numeric(f1), g = as.numeric(g1)),
        lambda = lambda,
        debias = debias,
        estimand = "ATE",
        p = p
      )
    
    if (debias) {
      p0 <- pot0$f_xx
      q0 <- pot0$g_yy
      
      p1 <- pot1$f_xx
      q1 <- pot1$g_yy
      
      attributes(tmat$w0$xx)[tm_attr_names] <- 
        list(
          dual = list(f = as.numeric(p0), g = as.numeric(p0)),
          lambda = lambda,
          debias = debias,
          estimand = "ATE",
          p = p
        )
      
      attributes(tmat$w1$xx)[tm_attr_names] <- 
        list(
          dual = list(f = as.numeric(p1), g = as.numeric(p1)),
          lambda = lambda,
          debias = debias,
          estimand = "ATE",
          p = p
        )
      
      attributes(tmat$w0$yy)[tm_attr_names] <- 
        list(
          dual = list(f = as.numeric(q0), g = as.numeric(q0)),
          lambda = lambda,
          debias = debias,
          estimand = "ATE",
          p = p
        )
      
      attributes(tmat$w1$yy)[tm_attr_names] <- 
        list(
          dual = list(f = as.numeric(q1), g = as.numeric(q1)),
          lambda = lambda,
          debias = debias,
          estimand = "ATE",
          p = p
        )
    } 
    
  } else if (cw@estimand != "ATE") {
    if (cw@estimand == "ATT") {
      x <- get_x0(dh)
      y <- get_x1(dh)
      a <- cw@w0
      b <- cw@w1
    } else if (cw@estimand == "ATC") {
      x <- get_x1(dh)
      y <- get_x0(dh)
      a <- cw@w1
      b <- cw@w0
    } else {
      stop("Estimand not found!")
    }
    
    ot <- OT$new(x = x, y = y, 
                 a = a, b = b,
                 penalty = lambda, 
                 cost = cost, p = p, debias = debias, 
                 tensorized = "tensorized",
                 diameter = diameter)
    ot$sinkhorn_opt(niter, tol)
    
    pot <- ot$potentials
    f <- pot$f_xy
    g <- pot$g_yx
    
    mat <- ot$primal()
    attributes(mat$xy)[tm_attr_names] <- 
      list(
        dual = list(f = as.numeric(f), g = as.numeric(g)),
        lambda = lambda,
        debias = debias,
        estimand = cw@estimand,
        p = p
      )
    
    tmat <- switch(cw@estimand,
                   "ATT" = list(w0 = mat,
                                w1 = blank_mat),
                   "ATC" = list(w0 = blank_mat,
                                w1 = mat))
    
    if (debias) {
      attributes(mat$xx)[tm_attr_names] <- 
        list(
          dual = list(f = as.numeric(pot$f_xx), g = as.numeric(pot$f_xx)),
          lambda = lambda,
          debias = debias,
          estimand = cw@estimand,
          p = p
        )
      
      attributes(mat$yy)[tm_attr_names] <- 
        list(
          dual = list(f = as.numeric(pot$g_yy), g = as.numeric(pot$g_yy)),
          lambda = lambda,
          debias = debias,
          estimand = cw@estimand,
          p = p
        )
    }
  }
  class(tmat) <- c("transportationMatrix", class(tmat))
  return(tmat)
}


# function to handle special cases
loss_select <- function(ot, niter, tol) {
  lambda <- ot$penalty
  if (is.finite(lambda)) {
    ot$sinkhorn_opt(niter, tol)
    return(sinkhorn_dist(ot))
  } else if ( is.infinite(lambda) ) {
    return(energy_dist(ot))
  }
}

#' Optimal Transport Distance
#'
#' @param x1 Either an object of class [causalOT::causalWeights-class] or a matrix of the covariates in the first sample
#' @param x2 `NULL` or a matrix of the covariates in the second sample.
#' @param a Empirical measure of the first sample. If NULL, assumes each observation gets equal mass. Ignored for objects of class causalWeights.
#' @param b Empirical measure of the second sample. If NULL, assumes each observation gets equal mass. Ignored for objects of class causalWeights.
#' @param penalty The penalty of the optimal transport distance to use. If missing or NULL, the function will try to guess a suitable value depending if debias is TRUE or FALSE.
#' @param p \eqn{L_p} distance metric power
#' @param cost Supply your own cost function. Should take arguments `x1`, `x2`, and `p`.
#' @param debias TRUE or FALSE. Should the debiased optimal transport distances be used.
#' @param online.cost How to calculate the distance matrix. One of "auto", "tensorized", or "online".
#' @param diameter The diameter of the metric space, if known. Default is NULL.
#' @param niter The maximum number of iterations for the Sinkhorn updates
#' @param tol The tolerance for convergence
#'
#' @return For objects of class matrix, numeric value giving the optimal transport distance. For objects of class causalWeights, results are returned as a list for before ('pre') and after adjustment ('post').
#' @export
#' @rdname ot_distance
#'
#' @examples
#' if ( torch::torch_is_installed()) {
#' x <- matrix(stats::rnorm(10*5), 10, 5)
#' z <- stats::rbinom(10, 1, 0.5)
#' weights <- calc_weight(x = x, z = z, method = "Logistic", estimand = "ATT")
#' ot1 <- ot_distance(x1 = weights, penalty = 100, 
#' p = 2, debias = TRUE, online.cost = "auto", 
#' diameter = NULL)
#' ot2<- ot_distance(x1 = x[z==0, ], x2 = x[z == 1,], 
#' a= weights@w0/sum(weights@w0), b = weights@w1,
#'  penalty = 100, p = 2, debias = TRUE, online.cost = "auto", diameter = NULL)
#'
#'  all.equal(ot1$post, ot2) 
#' }
ot_distance <- function(x1, x2 = NULL, 
         a = NULL, b = NULL,
         penalty, p = 2, 
         cost = NULL, 
         debias = TRUE, online.cost = "auto",
         diameter = NULL,
         niter = 1000, tol = 1e-7) UseMethod("ot_distance")

# setGeneric("ot_distance", function(x1, x2 = NULL, 
#                                    a = NULL, b = NULL,
#                                    penalty, p = 2, 
#                                    cost = NULL, 
#                                    debias = TRUE, online.cost = "auto",
#                                    diameter = NULL,
#                                    niter = 1000, tol = 1e-7) standardGeneric("ot_distance"))


#' @include weightsClass.R
#' @export 
#' @describeIn ot_distance method for causalWeights class
#' @method ot_distance causalWeights
ot_distance.causalWeights <- function(x1, x2 = NULL, a = NULL, b = NULL, penalty, p = 2, 
         cost = NULL, 
         debias = TRUE, online.cost = "auto",
         diameter=NULL,
         niter = 1000, tol = 1e-7) {
  
  
  stopifnot(inherits(x1, "causalWeights"))
  
  
  cw <- x1
  dh <- x1@data
  
  if(missing(penalty) || is.na(penalty) || is.null(penalty)) {
    warning("Penalty parameter not provided. Using estimated cost diameter as a penalty parameter.")
    maxes <- apply(dh@x,2,max)
    mins  <- apply(dh@x,2,min)
    
    penalty <- 1/p * sum((maxes-mins)^p)
    
  }
  
  stopifnot(penalty > 0.0)
  
  if( cw@estimand == "ATE" ) {
    x0 <- get_x0(dh)
    x1 <- get_x1(dh)
    x <- get_x(dh)
    z <- get_z(dh)
    b <- renormalize(get_w(dh))
    a0_init <- renormalize(b[z==0])
    a1_init <- renormalize(b[z==1])
    
    ot0_init <- OT$new(x = x0, y = x, 
                       a = a0_init, b = b,
                       penalty = penalty, 
                       cost = cost, p = p, debias = debias, 
                       tensorized = online.cost,
                       diameter = diameter)
    ot1_init <- OT$new(x = x1, y = x, 
                       a = a1_init, b = b,
                       penalty = penalty, 
                       cost = cost, p = p, debias = debias, 
                       tensorized = online.cost,
                       diameter = diameter)
    
    ot0 <- OT$new(x = x0, y = x, 
                  a = renormalize(cw@w0), b = b,
                  penalty = penalty, 
                  cost = cost, p = p, debias = debias, 
                  tensorized = online.cost,
                  diameter = diameter)
    ot1 <- OT$new(x = x1, y = x, 
                  a = renormalize(cw@w1), b = b,
                  penalty = penalty, 
                  cost = cost, p = p, debias = debias, 
                  tensorized = online.cost,
                  diameter = diameter)
    
    return(list(pre = c(control = as_numeric(loss_select(ot0_init, niter, tol)),
                        treated = as_numeric(loss_select(ot1_init, niter, tol))),
                post =  c(control = as_numeric(loss_select(ot0, niter, tol)),
                          treated = as_numeric(loss_select(ot1, niter, tol)))
    ))
    
    
  } else if (cw@estimand == "ATT" || cw@estimand == "ATC") {
    x0 <- get_x0(dh)
    x1 <- get_x1(dh)
    w <- get_w(dh)
    z <- get_z(dh)
    a_init <- renormalize(w[z==0])
    b_init <- renormalize(w[z==1])
  } else {
    stop("Estimand not found!")
  }
  
  ot_init <- OT$new(x = x0, y = x1, 
                    a = a_init, b = b_init,
                    penalty = penalty, 
                    cost = cost, p = p, debias = debias, 
                    tensorized = online.cost,
                    diameter = diameter)
  
  ot_final <- OT$new(x = x0, y = x1, 
                     a = renormalize(cw@w0), b = renormalize(cw@w1),
                     penalty = penalty, 
                     cost = cost, p = p, debias = debias, 
                     tensorized = online.cost,
                     diameter = diameter)
  
  
  return(list(pre = as_numeric(loss_select(ot_init, niter, tol)),
              post = as_numeric(loss_select(ot_final, niter, tol))))
  
}
# setMethod("ot_distance", signature(x1 = "causalWeights"),
# function(x1, x2 = NULL, a = NULL, b = NULL, penalty, p = 2, 
#          cost = NULL, 
#          debias = TRUE, online.cost = "auto",
#          diameter=NULL,
#          niter = 1000, tol = 1e-7) {
#   
#   
#   stopifnot(inherits(x1, "causalWeights"))
#   
#   
#   cw <- x1
#   dh <- x1@data
#   
#   if(missing(penalty) || is.na(penalty) || is.null(penalty)) {
#     warning("Penalty parameter not provided. Using estimated cost diameter as a penalty parameter.")
#     maxes <- apply(dh@x,2,max)
#     mins  <- apply(dh@x,2,min)
#     
#     penalty <- 1/p * sum((maxes-mins)^p)
#     
#   }
#   
#   stopifnot(penalty > 0.0)
#   
#   if( cw@estimand == "ATE" ) {
#     x0 <- get_x0(dh)
#     x1 <- get_x1(dh)
#     x <- get_x(dh)
#     z <- get_z(dh)
#     b <- renormalize(get_w(dh))
#     a0_init <- renormalize(b[z==0])
#     a1_init <- renormalize(b[z==1])
#     
#     ot0_init <- OT$new(x = x0, y = x, 
#                   a = a0_init, b = b,
#                   penalty = penalty, 
#                   cost = cost, p = p, debias = debias, 
#                   tensorized = online.cost,
#                   diameter = diameter)
#     ot1_init <- OT$new(x = x1, y = x, 
#                   a = a1_init, b = b,
#                   penalty = penalty, 
#                   cost = cost, p = p, debias = debias, 
#                   tensorized = online.cost,
#                   diameter = diameter)
#     
#     ot0 <- OT$new(x = x0, y = x, 
#                   a = renormalize(cw@w0), b = b,
#                   penalty = penalty, 
#                   cost = cost, p = p, debias = debias, 
#                   tensorized = online.cost,
#                   diameter = diameter)
#     ot1 <- OT$new(x = x1, y = x, 
#                   a = renormalize(cw@w1), b = b,
#                   penalty = penalty, 
#                   cost = cost, p = p, debias = debias, 
#                   tensorized = online.cost,
#                   diameter = diameter)
#     
#     return(list(pre = c(control = as_numeric(loss_select(ot0_init, niter, tol)),
#                            treated = as_numeric(loss_select(ot1_init, niter, tol))),
#                 post =  c(control = as_numeric(loss_select(ot0, niter, tol)),
#                            treated = as_numeric(loss_select(ot1, niter, tol)))
#                 ))
#     
#    
#   } else if (cw@estimand == "ATT" || cw@estimand == "ATC") {
#     x0 <- get_x0(dh)
#     x1 <- get_x1(dh)
#     w <- get_w(dh)
#     z <- get_z(dh)
#     a_init <- renormalize(w[z==0])
#     b_init <- renormalize(w[z==1])
#   } else {
#     stop("Estimand not found!")
#   }
#   
#   ot_init <- OT$new(x = x0, y = x1, 
#                      a = a_init, b = b_init,
#                      penalty = penalty, 
#                      cost = cost, p = p, debias = debias, 
#                      tensorized = online.cost,
#                      diameter = diameter)
#   
#   ot_final <- OT$new(x = x0, y = x1, 
#                a = renormalize(cw@w0), b = renormalize(cw@w1),
#                penalty = penalty, 
#                cost = cost, p = p, debias = debias, 
#                tensorized = online.cost,
#                diameter = diameter)
#   
#   
#   return(list(pre = as_numeric(loss_select(ot_init, niter, tol)),
#               post = as_numeric(loss_select(ot_final, niter, tol))))
#   
# }
#           
# )

ot_dist_default <- function(x1, x2, a = NULL, b = NULL, penalty, p = 2, 
                            cost = NULL, 
                            debias = TRUE, online.cost = "auto",
                            diameter=NULL,
                            niter = 1000, tol = 1e-7) {
  
  ot <- OT$new(x = x1, y = x2, 
               a = a, b = b,
               penalty = penalty, 
               cost = cost, p = p, debias = debias, 
               tensorized = online.cost,
               diameter = diameter)
  
  return(as_numeric(loss_select(ot, niter, tol)))
}

#' @export
#' @describeIn ot_distance method for matrices
#' @method ot_distance matrix
ot_distance.matrix <- ot_dist_default
# setMethod("ot_distance", signature(x1 = "matrix"), ot_dist_default)

#' @export 
#' @describeIn ot_distance method for arrays
#' @method ot_distance array
ot_distance.array <- ot_dist_default
# setMethod("ot_distance", signature(x1 = "array"), ot_dist_default)

#' @export 
#' @describeIn ot_distance method for torch_tensors
#' @method ot_distance torch_tensor
ot_distance.torch_tensor <- ot_dist_default
# setMethod("ot_distance", signature(x1 = "torch_tensor"), ot_dist_default)

