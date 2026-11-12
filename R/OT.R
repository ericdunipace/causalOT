# setOldClass("torch_tensor")
# setOldClass(c("OT","R6"))
OT <- R6::R6Class(
  "OT",
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
    sinkhorn_converged = "logical",
    initialize = function(
      x,
      y,
      a = NULL,
      b = NULL,
      penalty,
      cost_function = NULL,
      p = 2,
      debias = TRUE,
      tensorized = "auto",
      diameter = NULL,
      device = NULL,
      dtype = NULL
    ) {
      # browser()
      if (
        missing(penalty) || is.null(penalty) || is.na(penalty) || penalty < 0
      ) {
        stop("Must specify a penalty > 0!")
      } else {
        penalty <- as.double(penalty)
      }

      # check if should run online version in keops
      tensorized <- private$is_tensorized(tensorized, x, y)

      # should setup debiased potentials
      debias <- isTRUE(debias)

      # device check
      if (is.null(device) || !torch::is_torch_device(device)) {
        use_cuda <- torch::cuda_is_available() &&
          torch::cuda_device_count() >= 1
        if (use_cuda) {
          self$device <- torch::torch_device("cuda")
        } else {
          self$device <- torch::torch_device("cpu")
        }
      } else {
        self$device <- device
        attempt_cuda <- grepl("cuda", capture.output(print(self$device)))
        use_cuda <- attempt_cuda && torch::cuda_device_count() >= 1
        if (attempt_cuda && !use_cuda) {
          warning("CUDA not available even though you tried. Switching to CPU")
          self$device <- torch::torch_device("cpu")
        }
      }

      # dtype
      if (!is.null(dtype) && !torch::is_torch_dtype(dtype)) {
        warning(
          "Provided dtype is not or class torch_dtype. Using automatic selection."
        )
        dtype <- NULL
      }
      if (is.null(dtype)) {
        use_float <- use_cuda || identical(self$device$type, "mps")
        self$dtype <- switch(
          as.integer(use_float) + 1L,
          torch::torch_double(),
          torch::torch_float32()
        )
      } else {
        self$dtype <- dtype
      }

      # setup data
      if (!inherits(x, "torch_tensor")) {
        x <- torch::torch_tensor(x, dtype = self$dtype)$contiguous()
      }
      if (!inherits(y, "torch_tensor")) {
        y <- torch::torch_tensor(y, dtype = self$dtype)$contiguous()
      }
      d <- ncol(x)

      # setup masses
      a <- check_weights(a, x, self$device)
      b <- check_weights(b, y, self$device)
      a_log <- log_weights(a$detach())$contiguous()$to(device = self$device)
      b_log <- log_weights(b$detach())$contiguous()$to(device = self$device)

      # if (!tensorized) {
      #   self$device <- torch::torch_device("cpu")
      # } # avoids copying matrices more than needed

      # setup costs
      C_xy <- cost(
        x,
        y$detach(),
        p = p,
        tensorized = tensorized,
        cost_function = cost_function
      )
      C_yx <- cost(
        y,
        x$detach(),
        p = p,
        tensorized = tensorized,
        cost_function = cost_function
      )
      C_xy$to_device <- self$device
      C_yx$to_device <- self$device
      if (debias) {
        C_xx <- cost(
          x,
          x$detach(),
          p = p,
          tensorized = tensorized,
          cost_function = cost_function
        )
        C_yy <- cost(
          y,
          y$detach(),
          p = p,
          tensorized = tensorized,
          cost_function = cost_function
        )
        C_xx$to_device <- self$device
        C_yy$to_device <- self$device
      } else {
        C_xx <- C_yy <- NULL
      }

      softmin <- if (tensorized) softmin_tensorized else softmin_online

      # setup diameter
      if (
        missing(diameter) ||
          is.null(diameter) ||
          is.na(diameter) ||
          diameter < 0 ||
          is.infinite(diameter)
      ) {
        diameter <- private$diam_check(
          x$detach(),
          y$detach(),
          p = p,
          tensorized = tensorized,
          cost = C_xy$fun,
          C_yx
        )
      }

      self$C_xy <- C_xy
      self$C_yx <- C_yx
      self$C_xx <- C_xx
      self$C_yy <- C_yy
      private$a_ <- a
      private$b_ <- b
      self$n <- length(a)
      self$m <- length(b)
      private$a_log <- a_log
      private$b_log <- b_log
      self$penalty <- as.double(penalty)
      self$softmin <- softmin
      self$diameter <- diameter
      self$debias <- debias
      self$tensorized <- tensorized
      private$pot <- vector("list", 0)

      return(invisible(self))
    }
  ),
  active = list(
    potentials = function(value) {
      # browser()
      if (missing(value)) {
        return(private$pot)
      }
      if (is.null(value)) {
        private$pot <- list()
        private$pot_eps <- NULL
        private$g_yy_key <- NULL
        return(invisible(self))
      }
      if (all(is.character(names(value)))) {
        if (!all(names(value) %in% c("f_xy", "g_yx", "f_xx", "g_yy"))) {
          stop(
            "Names of potential list, if given, must be in f_xy, g_yx, f_xx, g_yy."
          )
        }
      }
      stopifnot(length(value) >= 2)
      if (is.null(names(value))) {
        warning(
          "No names given but assuming that the first two potentials are f_xy, g_yx"
        )
        names(value)[1:2] <- c("f_xy", "g_yx")
      }
      f_xy <- value$f_xy
      g_yx <- value$g_yx

      if (length(f_xy) != self$n) {
        stop("f_xy must have length of a")
      }
      if (length(g_yx) != self$m) {
        stop("g_yx must have length of b")
      }
      pot <- list(f_xy = f_xy, g_yx = g_yx)
      if (self$debias) {
        stopifnot(length(value) == 4)
        if (any(is.na(names(value)))) {
          warning(
            "No names given but assuming that the last two potentials are f_xx, g_yy"
          )
          names(value)[3:4] <- c("f_xx", "g_yy")
        }
        f_xx <- value$f_xx
        g_yy <- value$g_yy
        if (length(f_xx) != self$n) {
          stop("f_xx must have length of a")
        }
        if (length(g_yy) != self$m) {
          stop("g_yy must have length of b")
        }
        pot <- c(pot, list(f_xx = f_xx, g_yy = g_yy))
      }
      private$pot <- pot
      private$pot_eps <- as.numeric(self$penalty)
      private$g_yy_key <- NULL
      return(invisible(self))
    },
    a = function(value) {
      if (missing(value)) {
        return(private$a_)
      }

      if (length(value) != self$n) {
        stop(
          "Assignment measure to mass a must have same length as original problem"
        )
      }

      if (inherits(value, "torch_tensor")) {
        private$a_ <- value$to(device = private$a_$device)
      } else {
        private$a_ <- torch::torch_tensor(
          value,
          dtype = self$dtype,
          device = private$a_$device
        )
      }

      private$a_log <- log_weights(private$a_$detach())$to(
        device = private$a_$device
      )
    },
    b = function(value) {
      if (missing(value)) {
        return(private$b_)
      }

      if (length(value) != self$m) {
        stop(
          "Assignment measure to mass a must have same length as original problem"
        )
      }

      if (inherits(value, "torch_tensor")) {
        private$b_ <- value$to(device = private$b_$device)
      } else {
        private$b_ <- torch::torch_tensor(
          value,
          dtype = self$dtype,
          device = private$b_$device
        )
      }
      private$b_log <- log_weights(private$b_$detach())$to(
        device = private$b_$device
      )
    }
  ),
  private = list(
    a_ = "torch_tensor",
    b_ = "torch_tensor",
    a_log = "torch_tensor",
    b_log = "torch_tensor",
    pot = "list",
    pot_eps = NULL, # epsilon at which the stored potentials were computed
    g_yy_key = NULL, # inputs of the stored g_yy, see g_yy_current
    is_tensorized = function(tensorized, x, y) {
      if (tensorized == "auto") {
        n <- nrow(x)
        m <- nrow(y)
        cutoff <- log(n) + log(m) >= 17.03439
        if ((is.na(cutoff) || isTRUE(cutoff))) {
          tensorized <- FALSE
        } else {
          tensorized <- TRUE
        }
      } else if (tensorized == "tensorized") {
        tensorized <- TRUE
      } else if (tensorized == "online") {
        tensorized <- FALSE
      } else {
        stop("tensorized must be one of auto, tensorized, or online")
      }
      return(tensorized)
    },
    diam_check = function(x, y, p, tensorized, cost_function, C_yx) {
      # browser()
      diam_cost <- cost(
        torch::torch_vstack(list(
          x$max(1)[[1]]$view(c(1, -1)),
          x$min(1)[[1]]$view(c(1, -1))
        )),
        torch::torch_vstack(list(
          y$max(1)[[1]]$view(c(1, -1)),
          y$min(1)[[1]]$view(c(1, -1))
        )),
        p = p,
        tensorized = TRUE,
        cost_function = cost_function
      )
      return(max(diam_cost$data)$item())
    },
    xy_converged = "logical",
    xx_converged = list(x = NULL, y = NULL)
  )
)

log_weights <- function(a) {
  if (inherits(a, "torch_tensor")) {
    return(log_weights_torch_tensor(a))
  }
  if (is.numeric(a)) {
    return(log_weights_numeric(a))
  }
  stop(
    "no applicable method for 'log_weights' applied to an object of class ",
    paste0('"', class(a), '"', collapse = ", ")
  )
}

log_weights_torch_tensor <- function(a) {
  min_val <- as.double(-1e5)
  torch::with_no_grad({
    a_log <- a$log()$clamp_min(min_val)
  })
  return(a_log)
}

log_weights_numeric <- function(a) {
  min_val <- as.double(-1e5)
  a_log <- log(a)
  a_log[a_log < min_val] <- min_val
  return(a_log)
}

softmin_tensorized <- function(eps, C_xy, y_potential, b_log) {
  return(
    -eps *
      (torch::torch_add(y_potential / eps, b_log) - C_xy$data / eps)$logsumexp(
        2
      )
  )
}

online_lse_user <- torch::autograd_function(
  forward = function(ctx, x, y, G, P, C) {
    out <- torch::with_no_grad(online_lse_rows(C, x, y, G, P))
    ctx$save_for_backward(x = x, y = y, G = G, P = P, C = C)
    out
  },
  backward = function(ctx, grad_output) {
    s <- ctx$saved_variables
    gx <- torch::torch_zeros_like(s$x)
    gy <- torch::torch_zeros_like(s$y)
    Gv <- s$G$view(c(1L, -1L))
    n <- nrow(s$x)
    rows <- online_rows(s$x, nrow(s$y))
    for (lo in seq.int(1L, n, by = rows)) {
      i <- lo:min(n, lo + rows - 1L)
      xb <- s$x[i, , drop = FALSE]$detach()$requires_grad_(TRUE)
      yb <- s$y$detach()$requires_grad_(TRUE)
      torch::with_enable_grad({
        lse <- (Gv - s$C$fun(xb, yb, s$C$p) * s$P)$logsumexp(2)
        gr <- torch::autograd_grad(
          lse,
          list(xb, yb),
          grad_outputs = grad_output[i]
        )
      })
      gx[i, ] <- gr[[1]]
      gy$add_(gr[[2]])
    }
    list(x = gx, y = gy)
  }
)

softmin_online <- function(eps, C_xy, y_potential, b_log) {
  x <- C_xy$data$x
  y <- C_xy$data$y
  G <- (y_potential / eps + b_log)$detach()
  if (x$requires_grad || y$requires_grad) {
    return(online_lse_user(x, y, G, 1 / eps, C_xy) * (-eps))
  }
  if (identical(C_xy$algorithm, "user")) {
    return(online_lse_rows(C_xy, x, y, G, 1 / eps) * (-eps))
  }
  sinkhorn_compiled()$online_softmin(
    x,
    y,
    G,
    torch::jit_scalar(1 / eps),
    torch::jit_scalar(as.numeric(C_xy$p)),
    torch::jit_scalar(online_rows(x, nrow(y)))
  )
}


# epsilon scaling
epsilon_trajectory <- function(start, eps, tot_iter) {
  path <- if (start > eps) {
    pmax(start * exp(-0.9 * (0:ceiling(log(start / eps) / 0.9))), eps)
  } else {
    eps
  }
  c(path, rep(eps, max(0L, tot_iter - length(path))))
}

# if warm, will restart from the maximum of desired epsilon or epsilon
# used to fit the potentials
sinkhorn_start <- function(warm, pot_eps, diameter, eps) {
  if (!warm) {
    return(max(diameter, eps))
  }
  if (is.null(pot_eps) || !is.finite(pot_eps)) {
    return(max(diameter, eps))
  }
  max(pot_eps, eps)
}

sinkhorn_usable <- function(...) {
  all(vapply(
    list(...),
    function(p) !is.null(p) && !as.logical(p$isnan()$any()$to(device = "cpu")),
    logical(1)
  ))
}

# if cost needs gradient; weights and potentials enter detached
cost_requires_grad <- function(C) {
  d <- C$data
  if (inherits(d, "torch_tensor")) {
    return(d$requires_grad)
  }
  any(vapply(
    d,
    function(t) inherits(t, "torch_tensor") && t$requires_grad,
    logical(1)
  ))
}

sinkhorn_ratio <- TRUE

# Length of Anderson memory
sinkhorn_anderson_m_xy <- 10L
sinkhorn_anderson_m_xx <- 3L


OT$set("public", "sinkhorn_opt", function(niter = 1e3, tol = 1e-7) {
  if (is.finite(self$penalty)) {
    fg_list <- private$sinkhorn_loop(niter, tol)

    if (self$debias) {
      f_xx <- private$sinkhorn_self_loop("x", niter, tol)
      g_yy <- if (private$g_yy_current()) {
        private$pot$g_yy
      } else {
        private$sinkhorn_self_loop("y", niter, tol)
      }
      self$potentials <- list(
        f_xy = fg_list$f_xy,
        g_yx = fg_list$g_yx,
        f_xx = f_xx,
        g_yy = g_yy
      )
      private$g_yy_key <- private$g_yy_inputs()
    } else {
      self$potentials <- list(f_xy = fg_list$f_xy, g_yx = fg_list$g_yx)
    }
  }
  self$sinkhorn_converged <- isTRUE(private$xy_converged) &&
    (!self$debias ||
      (isTRUE(private$xx_converged$x) && isTRUE(private$xx_converged$y)))

  return(invisible(self))
})

# get the inputs needed for g_yy
OT$set("private", "g_yy_inputs", function() {
  if (cost_requires_grad(self$C_yy) || private$b_$requires_grad) {
    return(NULL)
  }
  list(
    b = private$b_$detach()$clone(),
    C = self$C_yy$data,
    eps = as.numeric(self$penalty)
  )
})
# identify if g_yy needs updating
OT$set("private", "g_yy_current", function() {
  key <- private$g_yy_key
  if (is.null(key) || !sinkhorn_usable(private$pot$g_yy)) {
    return(FALSE)
  }
  if (cost_requires_grad(self$C_yy) || private$b_$requires_grad) {
    return(FALSE)
  }
  identical(key$C, self$C_yy$data) &&
    identical(key$eps, as.numeric(self$penalty)) &&
    torch::torch_equal(key$b, private$b_$detach())
})
# only update non-g_yy potentials
OT$set("public", "sinkhorn_cot", function(niter = 1e3, tol = 1e-7) {
  torch::with_no_grad({
    fg_list <- private$sinkhorn_loop(niter, tol)

    if (self$debias) {
      g_yy <- private$pot$g_yy
      if (is.null(g_yy)) {
        g_yy <- torch::torch_zeros_like(fg_list$g_yx, dtype = self$dtype)
      }
      f_xx <- private$sinkhorn_self_loop("x", niter, tol)
      self$potentials <- list(
        f_xy = fg_list$f_xy,
        g_yx = fg_list$g_yx,
        f_xx = f_xx,
        g_yy = g_yy
      )
    } else {
      self$potentials <- list(f_xy = fg_list$f_xy, g_yx = fg_list$g_yx)
    }
  })
  # torch::autograd_set_grad_mode(enabled = TRUE)
  return(invisible(self))
})

OT$set("public", "sinkhorn_ccp_eval", function(niter = 1e3, tol = 1e-7) {
  torch::with_no_grad({
    fg_list <- private$sinkhorn_loop(niter, tol)
    if (!self$debias) {
      stop("CCP Sinkhorn evaluation requires debias = TRUE.")
    }
    g_yy <- private$pot$g_yy
    if (is.null(g_yy)) {
      g_yy <- torch::torch_zeros_like(fg_list$g_yx, dtype = self$dtype)
    }
    f_xx <- private$sinkhorn_self_loop("x", niter, tol)
    self$potentials <- list(
      f_xy = fg_list$f_xy,
      g_yx = fg_list$g_yx,
      f_xx = f_xx,
      g_yy = g_yy
    )
    loss <- self$a$dot(fg_list$f_xy) +
      self$b$dot(fg_list$g_yx) -
      self$a$dot(f_xx) -
      self$b$dot(g_yy)
  })

  list(
    loss = loss,
    f_xy = fg_list$f_xy,
    f_xx = f_xx,
    direction = f_xx - fg_list$f_xy
  )
})


OT$set("public", "sinkhorn_smc", function(niter = 1e3, tol = 1e-7) {
  if (is.finite(self$penalty)) {
    fg_list <- private$sinkhorn_loop(niter, tol)

    if (self$debias) {
      f_xx <- private$pot$f_xx
      if (is.null(f_xx)) {
        f_xx <- private$sinkhorn_self_loop("x", niter, tol)
      }
      g_yy <- private$sinkhorn_self_loop("y", niter, tol)
      self$potentials <- list(
        f_xy = fg_list$f_xy,
        g_yx = fg_list$g_yx,
        f_xx = f_xx,
        g_yy = g_yy
      )
    } else {
      self$potentials <- list(f_xy = fg_list$f_xy, g_yx = fg_list$g_yx)
    }
  }
  return(invisible(self))
})

OT$set(
  "private",
  "sinkhorn_self_loop",
  function(which.margin = "x", niter, tol) {
    if (!self$debias) {
      stop("Self potentials can only be run if debias == TRUE")
    }
    which.margin <- match.arg(which.margin, c("x", "y"))
    if (which.margin == "x") {
      f_xx <- private$pot$f_xx
      C_xx <- self$C_xx
      a_log <- private$a_log$detach()
      a <- a_log$exp()
    } else {
      f_xx <- private$pot$g_yy
      C_xx <- self$C_yy
      a_log <- private$b_log$detach()
      a <- a_log$exp()
    }
    eps <- as.numeric(self$penalty)
    softmin <- self$softmin
    warm <- sinkhorn_usable(f_xx)
    epsilons <- epsilon_trajectory(
      sinkhorn_start(warm, private$pot_eps, as.numeric(self$diameter), eps),
      eps,
      niter
    )

    torch::with_no_grad({
      f_xx <- if (warm) {
        f_xx$detach()$clone()
      } else {
        torch::torch_zeros_like(a_log, dtype = self$dtype)
      }
      converged <- FALSE
      if (self$tensorized) {
        out <- sinkhorn_compiled()$sinkhorn_xx(
          C_xx$data,
          a_log,
          a,
          f_xx,
          torch::torch_tensor(epsilons, dtype = torch::torch_double()),
          torch::jit_scalar(eps),
          torch::jit_scalar(tol),
          torch::jit_scalar(sinkhorn_anderson_m_xx),
          torch::jit_scalar(2),
          torch::jit_scalar(0L),
          torch::jit_scalar(FALSE),
          torch::jit_scalar(sinkhorn_ratio)
        )
        f_xx <- out$f
        iterations <- out$iterations$item()
        converged <- as.logical(out$converged$item())
      } else if (!identical(C_xx$algorithm, "user")) {
        x <- C_xx$data$x$detach()
        out <- sinkhorn_compiled()$sinkhorn_xx(
          x,
          a_log,
          a,
          f_xx,
          torch::torch_tensor(epsilons, dtype = torch::torch_double()),
          torch::jit_scalar(eps),
          torch::jit_scalar(tol),
          torch::jit_scalar(sinkhorn_anderson_m_xx),
          torch::jit_scalar(as.numeric(C_xx$p)),
          torch::jit_scalar(online_rows(x, nrow(x))),
          torch::jit_scalar(TRUE),
          torch::jit_scalar(sinkhorn_ratio)
        )
        f_xx <- out$f
        iterations <- out$iterations$item()
        converged <- as.logical(out$converged$item())
      } else {
        iterations <- 0L
        for (e in epsilons) {
          f_old <- f_xx
          f_xx <- (f_xx + softmin(torch::jit_scalar(e), C_xx, f_xx, a_log)) *
            0.5
          change <- (f_xx - f_old)$abs()$max()$item()
          iterations[[1L]] <- iterations[[1L]] + 1L
          if (e == eps && change < tol * eps) {
            converged <- TRUE
            break
          }
        }
      }
      private$xx_converged[[which.margin]] <- converged
    })
    # if (iterations == niter) {
    #   warning("Max OT iterations hit")
    # }
    # record grads if needed
    if (!cost_requires_grad(C_xx)) {
      return(f_xx)
    }
    softmin(
      eps,
      C_xx,
      f_xx$detach(),
      a_log
    )
  }
)

OT$set("private", "sinkhorn_loop", function(niter, tol) {
  eps <- as.numeric(self$penalty)
  softmin <- self$softmin
  a <- private$a_$detach()
  b <- private$b_$detach()
  a_log <- private$a_log$detach()
  b_log <- private$b_log$detach()
  C_xy <- self$C_xy
  C_yx <- self$C_yx
  f_xy <- private$pot$f_xy
  g_yx <- private$pot$g_yx
  warm <- sinkhorn_usable(f_xy, g_yx)
  epsilons <- epsilon_trajectory(
    sinkhorn_start(warm, private$pot_eps, as.numeric(self$diameter), eps),
    eps,
    niter
  )

  torch::with_no_grad({
    if (warm) {
      f_xy <- f_xy$detach()$clone()
      g_yx <- g_yx$detach()$clone()
    } else {
      f_xy <- torch::torch_zeros_like(a_log, dtype = self$dtype)
      g_yx <- torch::torch_zeros_like(b_log, dtype = self$dtype)
    }
    converged <- FALSE
    if (self$tensorized) {
      out <- sinkhorn_compiled()$sinkhorn_xy(
        C_xy$data,
        C_yx$data,
        a_log,
        b_log,
        a,
        b,
        f_xy,
        g_yx,
        torch::torch_tensor(epsilons, dtype = torch::torch_double()),
        torch::jit_scalar(eps),
        torch::jit_scalar(tol),
        torch::jit_scalar(sinkhorn_anderson_m_xy),
        torch::jit_scalar(2),
        torch::jit_scalar(0L),
        torch::jit_scalar(0L),
        torch::jit_scalar(FALSE),
        torch::jit_scalar(sinkhorn_ratio)
      )
      f_xy <- out$f
      g_yx <- out$g
      converged <- as.logical(out$converged$item())
    } else if (!identical(C_xy$algorithm, "user")) {
      x <- C_xy$data$x$detach()
      y <- C_xy$data$y$detach()
      out <- sinkhorn_compiled()$sinkhorn_xy(
        x,
        y,
        a_log,
        b_log,
        a,
        b,
        f_xy,
        g_yx,
        torch::torch_tensor(epsilons, dtype = torch::torch_double()),
        torch::jit_scalar(eps),
        torch::jit_scalar(tol),
        torch::jit_scalar(sinkhorn_anderson_m_xy),
        torch::jit_scalar(as.numeric(C_xy$p)),
        torch::jit_scalar(online_rows(x, nrow(y))),
        torch::jit_scalar(online_rows(y, nrow(x))),
        torch::jit_scalar(TRUE),
        torch::jit_scalar(sinkhorn_ratio)
      )
      f_xy <- out$f
      g_yx <- out$g
      converged <- as.logical(out$converged$item())
    } else {
      for (e in epsilons) {
        g_yx <- softmin(torch::jit_scalar(e), C_yx, f_xy, a_log) # OT(b,a)
        f_new <- softmin(torch::jit_scalar(e), C_xy, g_yx, b_log) # OT(a,b)
        # center the potentials
        norm <- (b * g_yx)$sum()
        g_yx$sub_(norm)
        f_new$add_(norm)
        change <- (f_new - f_xy)$abs()$max()$item()
        f_xy <- f_new
        if (e == eps && change < tol * eps) {
          converged <- TRUE
          break
        }
      }
    }
    private$xy_converged <- converged
    gt <- g_yx$detach()$clone()
    ft <- f_xy$detach()$clone()
  })

  # A final step records gradients through the cost; without them it would
  # only recompute what the loop has.
  if (!cost_requires_grad(C_xy) && !cost_requires_grad(C_yx)) {
    return(list(f_xy = ft, g_yx = gt))
  }
  eps <- torch::jit_scalar(eps)
  return(list(
    f_xy = softmin(eps, C_xy, gt, b_log),
    g_yx = softmin(eps, C_yx, ft, a_log)
  ))
})

OT$set("public", "primal", function(problem = c("xy", "xx", "yy")) {
  if (!self$tensorized) {
    stop("Not implemented for online calculation")
  }

  problem <- match.arg(problem)
  if (!self$debias && problem %in% c("xx", "yy")) {
    stop("Must have debiased OT problems for these quantities")
  }

  eps <- self$penalty
  a_log <- private$a_log
  b_log <- private$b_log
  n <- self$n
  m <- self$m
  pot <- self$potentials
  if (is.null(pot) || length(pot) == 0) {
    stop("Must run sinkhorn optimization first")
  }
  pi <- pi_xx <- pi_yy <- NULL
  if (problem == "xy") {
    f <- pot$f_xy
    g <- pot$g_yx

    pi <- ((f$view(c(n, 1)) + g$view(c(1, m)) - self$C_xy$data) /
      eps +
      a_log$view(c(n, 1)) +
      b_log$view(c(1, m)))$exp()
  } else if (problem == "xx") {
    C_xx <- self$C_xx
    p <- pot$f_xx
    pi_xx <- ((p$view(c(n, 1)) + p$view(c(1, n)) - C_xx$data) /
      eps +
      a_log$view(c(n, 1)) +
      a_log$view(c(1, n)))$exp()
    pi_xx <- (pi_xx + pi_xx$t()) / 2.0
  } else if (problem == "yy") {
    C_yy <- self$C_yy
    q <- pot$g_yy
    pi_yy <- ((q$view(c(m, 1)) + q$view(c(1, m)) - C_yy$data) /
      eps +
      b_log$view(c(m, 1)) +
      b_log$view(c(1, m)))$exp()
    pi_yy <- (pi_yy + pi_yy$t()) / 2.0
  } else {
    stop("Problem not found")
  }

  output <- list(xy = pi, xx = pi_xx, yy = pi_yy)
  return(output)
})

OT$set(
  "private",
  "kernel",
  function(problem = c("xy", "xx", "yy"), idx1 = NULL, idx2 = NULL) {
    problem <- match.arg(problem)
    if (!self$debias && problem %in% c("xx", "yy")) {
      stop("Must have debiased OT problems for these quantities")
    }

    pot <- self$potentials
    if (is.null(pot) || length(pot) == 0) {
      stop("Must run sinkhorn optimization first")
    }

    parts <- switch(
      problem,
      xy = list(f = pot$f_xy, g = pot$g_yx, C = self$C_xy),
      xx = list(f = pot$f_xx, g = pot$f_xx, C = self$C_xx),
      yy = list(f = pot$g_yy, g = pot$g_yy, C = self$C_yy)
    )

    if (is.null(idx1)) {
      idx1 <- seq_along(parts$f)
    }
    if (is.null(idx2)) {
      idx2 <- seq_along(parts$g)
    }
    if (!all(idx1 %in% seq_along(parts$f))) {
      stop("idx1 is not valid")
    }
    if (!all(idx2 %in% seq_along(parts$g))) {
      stop("idx2 is not valid")
    }

    f <- parts$f[idx1]$view(c(-1L, 1L))
    g <- parts$g[idx2]$view(c(1L, -1L))
    (f + g)$sub_(cost_block(parts$C, idx1, idx2))$div_(as.numeric(
      self$penalty
    ))$exp_()
  }
)

OT$set("private", "kernel_crossprod", function(idx) {
  b <- as_cpu_double_tensor(self$b)
  R <- torch::torch_zeros(
    c(length(idx), length(idx)),
    dtype = torch::torch_double()
  )
  for (lo in seq.int(1L, self$m, by = 512L)) {
    j <- lo:min(self$m, lo + 511L)
    K <- as_cpu_double_tensor(private$kernel("xy", idx, j))
    R <- R + (K * b[j]$view(c(1L, -1L)))$matmul(K$t())
  }
  R
})


OT$set("public", "reset_potentials", function() {
  private$pot[c("f_xy", "g_yx", "f_xx")] <- NULL
  invisible(self)
})

OT$set(
  "public",
  "hessian",
  function(problem = c("xy", "xx", "yy"), reduced = TRUE) {
    # whether to return a full hessian or one only with the nonzero entries
    reduced <- as.logical(isTRUE(reduced))

    problem <- match.arg(problem)
    debias <- self$debias
    if (problem %in% c("xx", "yy") && !debias) {
      stop("xx and yy hessians can only be calculated for debiased OT problems")
    }

    index_a <- which(as.logical(self$a$detach() > 0))
    index_b <- which(as.logical(self$b$detach() > 0))

    n <- length(index_a)
    m <- length(index_b)

    l_b <- private$b_log$detach()[index_b, drop = FALSE]
    eps <- self$penalty

    if (self$tensorized) {
      # if(!debias) {
      # measure <- torch::torch_cat(list(self$a,self$b[index_b]))
      # P_neg   <- -( (f/lambda + l_a)$view(c(n,1)) +
      #                 (g[index_b] / lambda + l_b)$view(c(1,m-1)) -
      #                 C[,index_b] / lambda
      # )$exp()

      # E_d_phi <- torch::torch_diag(-measure)

      # E_d_phi[(n+1):(n + m-1), 1:n]  <- P_neg$transpose(-1,1)
      # E_d_phi[1:n, (n+1):(n + m-1)]  <- P_neg
      # rm(P_neg)
      primal <- self$primal(problem = problem)

      if (problem == "xy") {
        index_subset <- if (m > 1) {
          1:(m - 1)
        } else {
          1
        }
        index_b_1 <- index_b[index_subset]

        P <- primal$xy$to(device = "cpu")$to(dtype = torch::torch_double())[
          index_a,
          index_b,
          drop = FALSE
        ]
        r <- P$sum(2) # row mass (approx a)
        s <- P$sum(1) # col mass (approx b)
        l_b <- l_b$to(device = "cpu")$to(
          dtype = torch::torch_double()
        )
        c <- (l_b[index_subset] - l_b[m])$exp()

        H_ff <- r$diag() * -(1 / eps)

        if (m == 1) {
          if (!reduced) {
            H <- torch::torch_zeros(
              c(self$n + self$m, self$n + self$m),
              dtype = torch::torch_double()
            )
            H[index_a][, index_a] <- H_ff
          } else {
            H <- H_ff
          }
        } else {
          H_fg <- (P[, index_subset] -
            P[, m]$view(c(-1, 1)) * c$view(c(1, -1))) *
            -(1 / eps)
          H_gg <-
            (s[index_subset]$diag() +
              s[m] * c$view(c(-1, 1)) * c$view(c(1, -1))) *
            -(1 / eps)

          H_temp <- torch::torch_vstack(list(
            torch::torch_hstack(list(H_ff, H_fg)),
            torch::torch_hstack(list(H_fg$t(), H_gg))
          ))

          if (!reduced) {
            idx <- c(index_a, self$n + index_b[index_subset])
            H <- torch::torch_zeros(
              c(self$n + self$m, self$n + self$m),
              dtype = torch::torch_double()
            )
            H[idx][, idx] <- H_temp
          } else {
            H <- H_temp
          }
        }
      } else if (problem == "xx") {
        P_xx <- primal$xx$to(device = "cpu")$to(dtype = torch::torch_double())[
          index_a,
          index_a,
          drop = FALSE
        ]
        r_xx <- P_xx$sum(2)
        H_xx <- (r_xx$diag() + P_xx) * -(1 / eps)

        if (!reduced) {
          H <- torch::torch_zeros(
            c(self$n, self$n),
            dtype = torch::torch_double()
          )
          H[index_a][, index_a] <- H_xx
        } else {
          H <- H_xx
        }
      } else if (problem == "yy") {
        P_yy <- primal$yy$to(device = "cpu")$to(dtype = torch::torch_double())[
          index_b,
          index_b,
          drop = FALSE
        ]
        r_yy <- P_yy$sum(2)
        H_yy <- (r_yy$diag() + P_yy) * -(1 / eps)
        if (!reduced) {
          H <- torch::torch_zeros(
            c(self$m, self$m),
            dtype = torch::torch_double()
          )
          H[index_b][, index_b] <- H_yy
        } else {
          H <- H_yy
        }
      }

      return(H)

      # if (debias)
    }
  }
)

OT$set(
  "private",
  "var_phi", #variance of gradient
  function(problem = c("xy", "xx", "yy")) {
    problem <- match.arg(problem)
    if (!self$tensorized) {
      stop("Online cost version of the variance not yet implemented")
    }
    # get basic variables
    n <- self$n
    m <- self$m

    a <- self$a
    b <- self$b
    # check that measures are positive and remove parts that aren't
    nz_a <- which(as.logical(a > 0))
    nz_b <- which(as.logical(b > 0))
    nn <- length(nz_a)
    mm <- length(nz_b)
    a <- a[nz_a, drop = FALSE]$to(device = "cpu")$to(
      dtype = torch::torch_double()
    )
    b <- b[nz_b, drop = FALSE]$to(device = "cpu")$to(
      dtype = torch::torch_double()
    )
    P <- self$primal(problem)[[problem]]$to(device = "cpu")$to(
      dtype = torch::torch_double()
    )
    if (problem == "xy") {
      P <- P[nz_a, nz_b, drop = FALSE]
      phi_f <- a$outer(b)$sub_(P) # nn x mm, driven by y

      if (mm > 1) {
        phi_g <- -P[, 1:(mm - 1)]$t() +
          (b[1:(mm - 1)] / b[mm])$outer(P[, mm]) # (mm-1) x nn, driven by x
      } else {
        # one point in b: g is fixed by the gauge, so no g parameters
        phi_g <- torch::torch_zeros(
          c(0, nn),
          dtype = P$dtype,
          device = P$device
        )
      }
      out <- list(x = phi_g, y = phi_f)
    } else if (problem == "xx") {
      phi_xx <- a$outer(a)$sub_(P[nz_a, nz_a, drop = FALSE]) # nn x nn
      out <- list(x = phi_xx, y = NULL)
    } else {
      phi_yy <- b$outer(b)$sub_(P[nz_b, nz_b, drop = FALSE]) # mm x mm
      out <- list(x = NULL, y = phi_yy)
    }

    lapply(out, function(s) {
      if (is.null(s)) s else s$to(device = "cpu", dtype = torch::torch_double())
    })
  }
)

OT$set("public", "vcov", function() {
  if (!self$tensorized) {
    stop("Online cost version of the variance not yet implemented")
  }
  torch::with_no_grad({
    # get basic variables
    n <- self$n
    m <- self$m
    nz_a <- which(as.logical(self$a > 0))
    nz_b <- which(as.logical(self$b > 0))
    nn <- length(nz_a)
    mm <- length(nz_b)
    nz_b_1 <- if (mm > 1) nz_b[1:(mm - 1)] else integer(0)
    zeros <- function(r, c) {
      torch::torch_zeros(c(r, c), dtype = torch::torch_double())
    }

    s_xy <- private$var_phi("xy")

    if (!self$debias) {
      # f_xy, g_yx
      U_x <- torch::torch_vstack(list(zeros(nn, nn), s_xy$x))
      U_y <- torch::torch_vstack(list(s_xy$y, zeros(mm - 1, mm)))
      H <- self$hessian("xy", reduced = TRUE)
      idxes <- c(nz_a, n + nz_b)
      p <- n + m
    } else {
      # f_xy, g_yx, f_xx, g_yy
      s_xx <- private$var_phi("xx")
      s_yy <- private$var_phi("yy")
      U_x <- torch::torch_vstack(list(
        zeros(nn, nn),
        s_xy$x,
        s_xx$x,
        zeros(mm, nn)
      ))
      U_y <- torch::torch_vstack(list(
        s_xy$y,
        zeros(mm - 1, mm),
        zeros(nn, mm),
        s_yy$y
      ))
      H <- torch::torch_block_diag(list(
        self$hessian("xy", reduced = TRUE),
        self$hessian("xx", reduced = TRUE),
        self$hessian("yy", reduced = TRUE)
      ))
      idxes <- c(nz_a, n + nz_b, n + m + nz_a, 2 * n + m + nz_b)
      p <- 2 * (n + m)
    }
    meat <- U_x$mm(U_x$t()) + U_y$mm(U_y$t())

    L <- torch::linalg_cholesky(-H)
    half <- torch::torch_cholesky_solve(meat, L)
    V <- torch::torch_cholesky_solve(half$t(), L)

    # variance of the very last potential fixed by gauge is complicated
    b_nz <- as.numeric(self$b$detach()$cpu())[nz_b]
    gauge_row <- -b_nz[-mm] / b_nz[mm] # how g_m loads on g_1..g_{m-1}
    G <- torch::torch_vstack(list(
      eye(mm - 1),
      torch::torch_tensor(
        matrix(gauge_row, nrow = 1),
        dtype = torch::torch_double()
      )
    ))
    n_self <- H$size(1) - (nn + mm - 1) # f_xx and g_yy; 0 if not debiased
    J <- torch::torch_block_diag(list(eye(nn), G, eye(n_self)))
    V <- J$mm(V)$mm(J$t())

    VAR <- zeros(p, p)
    VAR[idxes][, idxes] <- V
  })

  VAR$to(device = self$device, dtype = self$dtype)
})

OT$set("public", "hessian_weights", function(index = seq_len(self$n)) {
  a <- as_cpu_double_tensor(self$a)
  av <- a[index]$view(c(1L, -1L))
  I <- eye(length(index))

  R <- private$kernel_crossprod(index) # blocked cross prod
  H <- torch::linalg_solve(I - R * av + av, R)
  if (self$debias) {
    Kxx <- as_cpu_double_tensor(private$kernel("xx", index, index))
    H <- H + torch::linalg_solve(I + Kxx * av, Kxx)
  }
  H <- H * as.numeric(self$penalty)
  (H + H$t()) / 2
})

round_pi <- function(raw_pi, a, b) {
  n <- length(a)
  m <- length(b)

  x <- a / rowSums(raw_pi)
  x[x > 1] <- 1

  y <- b / colSums(raw_pi)
  y[y > 1] <- 1

  X <- diag(x)
  Y <- diag(y)

  pi_prime <- matrix(x, n, m) * raw_pi
  pi_2_prime <- pi_prime * matrix(y, n, m, byrow = TRUE)
  err_row <- a - rowSums(pi_2_prime)
  err_col <- b - colSums(pi_2_prime)

  return(
    pi_2_prime +
      matrix(err_row, n, m) *
        matrix(err_col, n, m, byrow = TRUE) /
        sum(abs(err_row))
  )
}

sinkhorn_dist <- function(OT) {
  if (!inherits(OT, "OT")) {
    stop("Must be an OT object")
  }
  pot <- OT$potentials
  a <- OT$a
  b <- OT$b
  if (is.null(pot) || length(pot) == 0) {
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
  if (!inherits(OT, "OT")) {
    stop("Must be an OT object")
  }
  linear_terms <- sinkhorn_dist(OT)

  eps <- OT$penalty
  C_xy <- OT$C_xy
  pot <- OT$potentials
  f_xy <- pot$f_xy
  g_yx <- pot$g_yx

  if (OT$tensorized) {
    n <- OT$n
    m <- OT$m
    a_log <- log_weights(OT$a)
    b_log <- log_weights(OT$b)

    K_xy <- (g_yx$view(c(1, m)) +
      f_xy$view(c(n, 1)) -
      C_xy$data +
      a_log$view(c(n, 1)) +
      b_log$view(c(1, m))) /
      eps
    exponential_terms <- -eps * K_xy$view(-1)$logsumexp(1)
    if (OT$debias) {
      C_yy <- OT$C_yy$data
      C_xx <- OT$C_xx$data

      f_xx <- pot$f_xx
      g_yy <- pot$g_yy

      K_xx <- (f_xx$view(c(1, n)) +
        f_xx$view(c(n, 1)) -
        C_xx +
        a_log$view(c(n, 1)) +
        a_log$view(c(1, n))) /
        eps
      K_yy <- (g_yy$view(c(1, m)) +
        g_yy$view(c(m, 1)) -
        C_yy +
        b_log$view(c(m, 1)) +
        b_log$view(c(1, m))) /
        eps

      exponential_terms <- exponential_terms +
        0.5 * eps * K_xx$view(-1)$logsumexp(1) +
        0.5 * eps * K_yy$view(-1)$logsumexp(1)
    }
  } else {
    a_log <- log_weights(OT$a)
    b_log <- log_weights(OT$b)
    # softmin(eps, C_xy, y_potential, b_log)
    s_x <- OT$softmin(eps, C_xy, g_yx, b_log)
    log_mass <- ((f_xy - s_x) / eps + a_log)$logsumexp(1)
    exponential_terms <- -eps * log_mass

    if (OT$debias) {
      s_xx <- OT$softmin(eps, OT$C_xx, pot$f_xx, a_log)
      s_yy <- OT$softmin(eps, OT$C_yy, pot$g_yy, b_log)
      exponential_terms <- exponential_terms +
        0.5 * eps * ((pot$f_xx - s_xx) / eps + a_log)$logsumexp(1) +
        0.5 * eps * ((pot$g_yy - s_yy) / eps + b_log)$logsumexp(1)
    }
  }

  loss <- linear_terms + exponential_terms
  return(loss)
}

# each target point j spreads its mass b_j over the sources by a softmin of
# its cost column; the summed source mass is the semi-relaxed starting weight
cot_init_semirelaxed <- function(C, b, lambda) {
  if (inherits(C, "costTensor")) {
    C <- C$data
  }
  if (!inherits(C, "torch_tensor")) {
    C <- torch::torch_tensor(C, dtype = torch::torch_double())
  }
  if (!inherits(b, "torch_tensor")) {
    b <- torch::torch_tensor(b, dtype = C$dtype, device = C$device)
  } else {
    b <- b$to(dtype = C$dtype, device = C$device)
  }
  stopifnot(
    "lambda must be a finite positive scalar." = length(lambda) == 1L &&
      is.finite(lambda) &&
      lambda > 0
  )
  stopifnot(
    "b must have one entry per target point." = length(b) == dim(C)[[2L]]
  )

  log_b <- log_weights(b$detach())
  eta <- -C / lambda
  log_col_norm <- eta$logsumexp(1)$view(c(1L, -1L))
  log_w_terms <- eta - log_col_norm + log_b$view(c(1L, -1L))
  log_w_terms$logsumexp(2)$exp()
}

# skip the b'g term of the self divergence if y is constant
cot_ccp_y_self_term <- function(ot, niter = 1000L, tol = 1e-10) {
  ot_private <- ot$.__enclos_env__$private
  # need to copy g_yy otherwise it's edited by reference
  g_yy_stored <- ot_private$pot$g_yy
  if (!is.null(g_yy_stored)) {
    ot_private$pot$g_yy <- g_yy_stored$detach()$clone()
    on.exit(ot_private$pot$g_yy <- g_yy_stored, add = TRUE)
  }
  torch::with_no_grad({
    g_yy <- ot_private$sinkhorn_self_loop("y", niter, tol)
    ot$b$detach()$dot(g_yy)$item()
  })
}

# scale for stopping tests of weight solvers
cot_conv_scale <- function(ot) {
  diameter <- suppressWarnings(as.numeric(ot$diameter))
  if (length(diameter) == 1L && is.finite(diameter) && diameter > 0) {
    return(diameter)
  }
  max(as.numeric(ot$penalty), .Machine$double.eps)
}

# check final fit convergence
cot_converged <- function(gap, D, atol, tol) {
  values <- list(gap, D, atol, tol)
  if (any(lengths(values) != 1L) || any(!is.finite(unlist(values)))) {
    return(FALSE)
  }
  gap <= atol + tol * max(0, D)
}

cot_fit_valid <- function(weights, loss) {
  if (inherits(weights, "torch_tensor")) {
    weights <- as.numeric(weights$detach()$to(device = "cpu"))
  }
  length(weights) > 0L &&
    all(is.finite(weights)) &&
    all(weights >= -cot_polish_ctl$round_factor * .Machine$double.eps) &&
    sum(weights) > 0 &&
    length(loss) == 1L &&
    is.finite(as.numeric(loss))
}

cot_ccp_kkt <- function(ot, support.tol = 1e-10) {
  if (!inherits(ot, "OT")) {
    stop("Must provide an OT object.")
  }
  if (!ot$debias) {
    stop("CCP KKT diagnostics require debias = TRUE.")
  }
  pot <- ot$potentials
  if (is.null(pot) || length(pot) == 0) {
    stop("Must run Sinkhorn optimization first.")
  }
  w <- ot$a$detach()
  kkt <- (pot$f_xy - pot$f_xx)$detach()
  active <- as.logical((w > support.tol)$to(device = "cpu"))
  kappa <- if (any(active)) {
    (w[active]$dot(kkt[active]) / w[active]$sum())$item()
  } else {
    NA_real_
  }
  kkt_dev <- if (any(active)) {
    (kkt[active] - kappa)$abs()$max()$item()
  } else {
    NA_real_
  }
  fw_gap <- (w$dot(kkt) - kkt$min())$item()
  list(kkt_constant = kappa, kkt_max_deviation = kkt_dev, fw_gap = fw_gap)
}

cot_ccp_diagnostics <- function(ot, support.tol = 1e-10) {
  if (!inherits(ot, "OT")) {
    stop("Must provide an OT object.")
  }
  if (!ot$tensorized) {
    stop("CCP diagnostics require tensorized costs.")
  }
  if (!ot$debias) {
    stop("CCP diagnostics require debias = TRUE.")
  }
  pot <- ot$potentials
  if (is.null(pot) || length(pot) == 0) {
    stop("Must run Sinkhorn optimization first.")
  }
  kkt <- cot_ccp_kkt(ot, support.tol = support.tol)
  plan <- ot$primal()$xy$detach()
  objective <- (ot$a$dot(pot$f_xy) +
    ot$b$dot(pot$g_yx) -
    ot$a$dot(pot$f_xx))$detach()$item() -
    cot_ccp_y_self_term(ot)
  list(
    objective = objective,
    kkt_constant = kkt$kkt_constant,
    kkt_max_deviation = kkt$kkt_max_deviation,
    fw_gap = kkt$fw_gap,
    relative_fw_gap = kkt$fw_gap /
      max(abs(objective), cot_conv_scale(ot)),
    row_marginal_max_error = (plan$sum(2) - ot$a)$abs()$max()$item(),
    column_marginal_max_error = (plan$sum(1) - ot$b)$abs()$max()$item()
  )
}

cost_bilinear <- function(C, u, v) {
  if (inherits(C, "costOnline")) {
    return(online_bilinear(C$data$x, C$data$y, u, v, C))
  }
  u$dot(C$data$matmul(v))
}

energy_dist <- function(OT) {
  if (!inherits(OT, "OT")) {
    stop("Must be an OT object")
  }
  if (!OT$debias) {
    stop("Must have option debias set to TRUE for energy distance")
  }
  a <- OT$a
  b <- OT$b
  loss_cross <- if (cost_requires_grad(OT$C_yx)) {
    cost_bilinear(OT$C_yx, b, a)
  } else {
    cost_bilinear(OT$C_xy, a, b)
  }

  loss <- loss_cross -
    cost_bilinear(OT$C_xx, a, a) * 0.5 -
    cost_bilinear(OT$C_yy, b, b) * 0.5

  return(loss)
}


inf_sinkhorn_dist <- function(OT) {
  if (!inherits(OT, "OT")) {
    stop("Must be an OT object")
  }
  if (OT$debias) {
    return(energy_dist(OT))
  }
  loss <- if (cost_requires_grad(OT$C_yx)) {
    cost_bilinear(OT$C_yx, OT$b, OT$a)
  } else {
    cost_bilinear(OT$C_xy, OT$a, OT$b)
  }

  return(loss)
}

online_rows <- function(x, m) {
  mb <- getOption("causalOT.online_block_mb", 128)
  max(1L, as.integer(floor(mb * 2^20 / (max(1, m) * x$element_size()))))
}


online_by_rows <- function(C, x, y, f) {
  n <- nrow(x)
  rows <- online_rows(x, nrow(y))
  out <- lapply(seq.int(1L, n, by = rows), function(lo) {
    i <- lo:min(n, lo + rows - 1L)
    f(C$fun(x[i, , drop = FALSE], y, C$p), i)
  })
  torch::torch_cat(out, 1)
}

online_lse_rows <- function(C, x, y, G, P) {
  Gv <- G$view(c(1L, -1L))
  online_by_rows(C, x, y, function(Cb, i) (Cb * (-P))$add_(Gv)$logsumexp(2))
}

online_cost_mv <- function(C, x, y, v) {
  online_by_rows(C, x, y, function(Cb, i) Cb$mv(v))
}

online_cost_tmv <- function(C, x, y, u) {
  out <- torch::torch_zeros(nrow(y), dtype = x$dtype, device = x$device)
  rows <- online_rows(x, nrow(y))
  for (lo in seq.int(1L, nrow(x), by = rows)) {
    i <- lo:min(nrow(x), lo + rows - 1L)
    out$add_(C$fun(x[i, , drop = FALSE], y, C$p)$t()$mv(u[i]))
  }
  out
}


online_bilinear <- torch::autograd_function(
  forward = function(ctx, x, y, u, v, C) {
    Cv <- online_cost_mv(C, x, y, v)
    ctx$save_for_backward(x = x, y = y, u = u, v = v, Cv = Cv, C = C)
    u$dot(Cv)
  },
  backward = function(ctx, grad_output) {
    s <- ctx$saved_variables
    nig <- ctx$needs_input_grad
    out <- list(x = NULL, y = NULL, u = NULL, v = NULL)
    if (isTRUE(nig$u)) {
      out$u <- grad_output * s$Cv
    }
    if (isTRUE(nig$v)) {
      out$v <- grad_output * online_cost_tmv(s$C, s$x, s$y, s$u)
    }
    if (isTRUE(nig$x) || isTRUE(nig$y)) {
      g <- online_bilinear_grad(
        s$x$detach(),
        s$y$detach(),
        s$u$detach(),
        s$v$detach(),
        s$C,
        isTRUE(nig$x),
        isTRUE(nig$y)
      )
      if (isTRUE(nig$x)) {
        out$x <- grad_output * g$x
      }
      if (isTRUE(nig$y)) out$y <- grad_output * g$y
    }
    out
  }
)

online_bilinear_grad <- function(x, y, u, v, C, need_x, need_y) {
  if (!identical(C$algorithm, "user") && C$p == 2) {
    # grad_x C_ij = x_i - y_j for C = |x - y|^2 / 2.
    return(list(
      x = if (need_x) {
        u$view(c(-1L, 1L)) * (x * v$sum() - v$view(c(1L, -1L))$matmul(y))
      },
      y = if (need_y) {
        v$view(c(-1L, 1L)) * (y * u$sum() - u$view(c(1L, -1L))$matmul(x))
      }
    ))
  }
  n <- nrow(x)
  gx <- if (need_x) torch::torch_empty_like(x)
  gy <- if (need_y) torch::torch_zeros_like(y)
  rows <- online_rows(x, nrow(y))
  for (lo in seq.int(1L, n, by = rows)) {
    i <- lo:min(n, lo + rows - 1L)
    torch::with_enable_grad({
      xb <- x[i, , drop = FALSE]$clone()$requires_grad_(need_x)
      yb <- y$clone()$requires_grad_(need_y)
      sb <- u[i, drop = FALSE]$dot(C$fun(xb, yb, C$p)$mv(v))
      inputs <- c(if (need_x) list(xb), if (need_y) list(yb))
      gr <- torch::autograd_grad(list(sb), inputs)
    })
    if (need_x) {
      gx[i, ] <- gr[[1]]
    }
    if (need_y) {
      gy <- gy + gr[[length(gr)]]
    }
    online_gc(x)
  }
  list(x = gx, y = gy)
}

online_gc <- function(x) {
  if (!identical(x$device$type, "cpu")) {
    invisible(gc(full = FALSE))
  }
}

transportationMatrix <- function(
  x = NULL,
  z = NULL,
  weights = NULL,
  lambda = NULL,
  p = 2,
  cost = NULL,
  debias = FALSE,
  diameter = NULL,
  niter = 1000,
  tol = 1e-7
) {
  #sets up the attribute names
  tm_attr_names <- c("dual", "lambda", "debias", "estimand", "p")
  # sets up a blank matrix to plug in if needed
  blank_mat <- list(
    xy = matrix(0, 0, 0),
    xx = matrix(0, 0, 0),
    yy = matrix(0, 0, 0)
  )
  attributes(blank_mat) <- list(
    dual = numeric(0),
    lambda = NA_real_,
    debias = logical(0),
    estimand = NULL,
    p = NA_real_
  )
  if (is.null(weights)) {
    #escape hatch

    tmat <- list(w0 = blank_mat, w1 = blank_mat)
    class(tmat) <- c("transportationMatix", class(tmat))
    return(tmat)
  }

  stopifnot(inherits(weights, "causalWeights"))

  cw <- weights
  dh <- dataHolder(x = x, z = z)
  n <- get_n(dh)

  if (is.null(lambda) || is.na(lambda)) {
    lambda <- 1 / log(n)
  }
  stopifnot(lambda > 0.0)

  if (cw@estimand == "ATE") {
    x0 <- get_x0(dh)
    x1 <- get_x1(dh)
    y <- get_x(dh)
    b <- dh@weights

    ot0 <- OT$new(
      x = x0,
      y = y,
      a = cw@w0,
      b = b,
      penalty = lambda,
      cost = cost,
      p = p,
      debias = debias,
      tensorized = "tensorized",
      diameter = diameter
    )
    ot1 <- OT$new(
      x = x1,
      y = y,
      a = cw@w1,
      b = b,
      penalty = lambda,
      cost = cost,
      p = p,
      debias = debias,
      tensorized = "tensorized",
      diameter = diameter
    )

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

    tmat <- list(w0 = pi_0, w1 = pi_1)

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

    ot <- OT$new(
      x = x,
      y = y,
      a = a,
      b = b,
      penalty = lambda,
      cost = cost,
      p = p,
      debias = debias,
      tensorized = "tensorized",
      diameter = diameter
    )
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

    tmat <- switch(
      cw@estimand,
      "ATT" = list(w0 = mat, w1 = blank_mat),
      "ATC" = list(w0 = blank_mat, w1 = mat)
    )

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
  } else if (is.infinite(lambda)) {
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
ot_distance <- function(
  x1,
  x2 = NULL,
  a = NULL,
  b = NULL,
  penalty,
  p = 2,
  cost = NULL,
  debias = TRUE,
  online.cost = "auto",
  diameter = NULL,
  niter = 1000,
  tol = 1e-7
) {
  UseMethod("ot_distance")
}

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
ot_distance.causalWeights <- function(
  x1,
  x2 = NULL,
  a = NULL,
  b = NULL,
  penalty,
  p = 2,
  cost = NULL,
  debias = TRUE,
  online.cost = "auto",
  diameter = NULL,
  niter = 1000,
  tol = 1e-7
) {
  stopifnot(inherits(x1, "causalWeights"))

  cw <- x1
  dh <- x1@data

  if (missing(penalty) || is.na(penalty) || is.null(penalty)) {
    warning(
      "Penalty parameter not provided. Using estimated cost diameter as a penalty parameter."
    )
    maxes <- apply(dh@x, 2, max)
    mins <- apply(dh@x, 2, min)

    penalty <- 1 / p * sum((maxes - mins)^p)
  }

  stopifnot(penalty > 0.0)

  if (cw@estimand == "ATE") {
    x0 <- get_x0(dh)
    x1 <- get_x1(dh)
    x <- get_x(dh)
    z <- get_z(dh)
    b <- renormalize(get_w(dh))
    a0_init <- renormalize(b[z == 0])
    a1_init <- renormalize(b[z == 1])

    ot0_init <- OT$new(
      x = x0,
      y = x,
      a = a0_init,
      b = b,
      penalty = penalty,
      cost = cost,
      p = p,
      debias = debias,
      tensorized = online.cost,
      diameter = diameter
    )
    ot1_init <- OT$new(
      x = x1,
      y = x,
      a = a1_init,
      b = b,
      penalty = penalty,
      cost = cost,
      p = p,
      debias = debias,
      tensorized = online.cost,
      diameter = diameter
    )

    ot0 <- OT$new(
      x = x0,
      y = x,
      a = renormalize(cw@w0),
      b = b,
      penalty = penalty,
      cost = cost,
      p = p,
      debias = debias,
      tensorized = online.cost,
      diameter = diameter
    )
    ot1 <- OT$new(
      x = x1,
      y = x,
      a = renormalize(cw@w1),
      b = b,
      penalty = penalty,
      cost = cost,
      p = p,
      debias = debias,
      tensorized = online.cost,
      diameter = diameter
    )

    return(list(
      pre = c(
        control = as_numeric(loss_select(ot0_init, niter, tol)),
        treated = as_numeric(loss_select(ot1_init, niter, tol))
      ),
      post = c(
        control = as_numeric(loss_select(ot0, niter, tol)),
        treated = as_numeric(loss_select(ot1, niter, tol))
      )
    ))
  } else if (cw@estimand == "ATT" || cw@estimand == "ATC") {
    x0 <- get_x0(dh)
    x1 <- get_x1(dh)
    w <- get_w(dh)
    z <- get_z(dh)
    a_init <- renormalize(w[z == 0])
    b_init <- renormalize(w[z == 1])
  } else {
    stop("Estimand not found!")
  }

  ot_init <- OT$new(
    x = x0,
    y = x1,
    a = a_init,
    b = b_init,
    penalty = penalty,
    cost = cost,
    p = p,
    debias = debias,
    tensorized = online.cost,
    diameter = diameter
  )

  ot_final <- OT$new(
    x = x0,
    y = x1,
    a = renormalize(cw@w0),
    b = renormalize(cw@w1),
    penalty = penalty,
    cost = cost,
    p = p,
    debias = debias,
    tensorized = online.cost,
    diameter = diameter
  )

  return(list(
    pre = as_numeric(loss_select(ot_init, niter, tol)),
    post = as_numeric(loss_select(ot_final, niter, tol))
  ))
}

ot_dist_default <- function(
  x1,
  x2,
  a = NULL,
  b = NULL,
  penalty,
  p = 2,
  cost = NULL,
  debias = TRUE,
  online.cost = "auto",
  diameter = NULL,
  niter = 1000,
  tol = 1e-7
) {
  ot <- OT$new(
    x = x1,
    y = x2,
    a = a,
    b = b,
    penalty = penalty,
    cost = cost,
    p = p,
    debias = debias,
    tensorized = online.cost,
    diameter = diameter
  )

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
