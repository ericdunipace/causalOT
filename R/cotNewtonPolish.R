# Float64 refinement of solutions

# Internal constants
cot_polish_ctl <- list(
  maxit = 20L,
  tol = 1e-8,
  sinkhorn_tol = 1e-13,
  sinkhorn_niter = 10000L,
  max_free = 8000L,
  armijo = 1e-4,
  alpha_min = 1e-8,
  zero_tol = 1e-12,
  lbfgs_m = 10L,
  lbfgs_maxit = 200L,
  round_factor = 64, # multiple of eps treated as rounding noise
  sinkhorn_factor = 16 # multiple of the sinkhorn tol treated as potential noise
)

cot_polish_twin <- function(ot, x, y) {
  b <- renormalize(as_cpu_double_tensor(ot$b))
  OT$new(
    x = as_cpu_double_tensor(x),
    y = as_cpu_double_tensor(y),
    a = renormalize(as_cpu_double_tensor(ot$a)),
    b = b,
    penalty = ot$penalty,
    # the default cost must stay NULL or OT treats it as user-supplied
    cost_function = if (identical(ot$C_xy$algorithm, "user")) ot$C_xy$fun,
    p = ot$C_xy$p,
    debias = TRUE,
    tensorized = if (ot$tensorized) "tensorized" else "online",
    diameter = ot$diameter,
    device = torch::torch_device("cpu"),
    dtype = torch::torch_float64()
  )
}


cot_polish_eval <- function(twin, w, ctl = cot_polish_ctl) {
  lambda <- as.numeric(twin$penalty)
  twin$a <- as_cpu_double_tensor(w)
  # sinkhorn tolerances are relative to epsilon
  twin$sinkhorn_opt(
    niter = ctl$sinkhorn_niter,
    tol = cot_polish_sinkhorn_tol(lambda, ctl) / lambda
  )
  if (!isTRUE(twin$sinkhorn_converged)) {
    stop("Polish Sinkhorn did not converge")
  }
  pot <- twin$potentials
  grad <- pot$f_xy - pot$f_xx
  list(
    loss = (twin$a$dot(grad) + twin$b$dot(pot$g_yx))$item(),
    y_self = twin$b$dot(pot$g_yy)$item(),
    grad = as.numeric(grad)
  )
}

cot_polish_sinkhorn_tol <- function(lambda, ctl) {
  max(ctl$sinkhorn_tol, ctl$round_factor * .Machine$double.eps * lambda)
}

# Active-set solve of the simplex QP min d'Hd/2 + g'd, sum(d) = 0, a + d >= 0.
# The binding set is retained across Levenberg shifts and outer steps.
cot_polish_exact <- function(H, g, a, ctl, state) {
  n <- length(a)
  binding <- if (length(state$active) == n) {
    state$active
  } else {
    # The free set changed, so start from the last binding set restricted to it.
    (a <= ctl$zero_tol & g > sum(a * g) / sum(a)) |
      state$binding[state$indices]
  }
  if (all(binding)) {
    binding[which.min(g)] <- FALSE
  }
  kkt_tol <- ctl$round_factor * .Machine$double.eps * max(1, max(abs(g)))
  # The Householder reflection Q = I - beta v v' maps e_1 to the all-ones
  # direction, so its last s - 1 columns Z span sum(d) = 0. Then Z'AZ is QAQ
  # without its first row and column, and applying Q directly is O(s^2).
  reduced_step <- function(F, rhs) {
    s <- length(F)
    v <- rep(1 / sqrt(s), s)
    v[1L] <- v[1L] + 1
    beta <- 2 / sum(v^2)
    v <- as_cpu_double_tensor(v)
    reflect <- function(u) u - v * (beta * v$dot(u)$item())

    A <- H[F, F, drop = FALSE]
    Av <- A$matmul(v)
    w <- Av * beta - v * (0.5 * beta^2 * v$dot(Av)$item())
    L <- torch::linalg_cholesky(
      (A - v$outer(w) - w$outer(v))[2:s, 2:s, drop = FALSE]
    )

    z <- torch::torch_cholesky_solve(
      -reflect(as_cpu_double_tensor(rhs))[2:s]$view(c(-1L, 1L)),
      L
    )
    zero <- torch::torch_zeros(1L, dtype = torch::torch_double())
    as.numeric(reflect(torch::torch_cat(list(zero, z$view(-1L)))))
  }
  for (attempt in seq_len(20L)) {
    state$active <- binding
    F <- which(!binding)
    d <- -a
    d[F] <- sum(a[binding]) / length(F)
    if (length(F) > 1L) {
      d[F] <- d[F] +
        reduced_step(
          F,
          g[F] +
            as.numeric(H[F, , drop = FALSE]$matmul(as_cpu_double_tensor(d)))
        )
    }
    negative <- F[a[F] + d[F] < -ctl$zero_tol]
    if (length(negative)) {
      binding[negative] <- TRUE
      next
    }
    q <- g + as.numeric(H$matmul(as_cpu_double_tensor(d)))
    release <- which(binding & q < mean(q[F]) - kkt_tol)
    if (!length(release)) {
      state$binding[state$indices] <- binding
      return(d)
    }
    binding[release] <- FALSE
  }
  stop("Polish QP needs a larger Levenberg shift")
}

# Free and binding coordinates. Bertsekas' epsilon shrinks with the
# projected-gradient residual; the inverse cost scale makes the gradient step
# dimensionless on the simplex.
cot_polish_sets <- function(a, g, scale, ctl = cot_polish_ctl) {
  S <- which(a > ctl$zero_tol)
  nu <- sum(a[S] * g[S]) / sum(a[S])
  v <- a - (g - nu) / max(scale, .Machine$double.eps)
  u <- sort(v, decreasing = TRUE)
  theta <- (cumsum(u) - 1) / seq_along(u)
  projected <- pmax(v - theta[max(which(u > theta))], 0)
  epsilon <- min(1e-2 / length(a), sqrt(sum((a - projected)^2)))
  # Avoid declaring a small optimal weight binding because of potential noise.
  binding <- which(a <= epsilon & g > nu + ctl$sinkhorn_factor * ctl$sinkhorn_tol)
  zeros <- which(a <= ctl$zero_tol & g < nu)
  list(
    free = sort(unique(c(setdiff(which(a > 0), binding), zeros))),
    binding = binding
  )
}

cot_polish_lbfgs <- function(g, memory, mean_weight) {
  q <- g - mean(g)
  alpha <- numeric(length(memory))
  for (i in rev(seq_along(memory))) {
    pair <- memory[[i]]
    alpha[i] <- sum(pair$s * q) / pair$sy
    q <- q - alpha[i] * pair$y
  }
  scale <- if (length(memory)) {
    p <- memory[[length(memory)]]
    p$sy / sum(p$y^2)
  } else {
    magnitude <- max(abs(q))
    if (magnitude > 0) mean_weight / magnitude else 1
  }
  r <- scale * q
  for (i in seq_along(memory)) {
    pair <- memory[[i]]
    beta <- sum(pair$y * r) / pair$sy
    r <- r + pair$s * (alpha[i] - beta)
  }
  -r + mean(r)
}

cot_polish_newton <- function(twin, a, ev, base, F, state, ctl) {
  H <- twin$hessian_weights(F)

  # The Hessian gauge is irrelevant on sum(d) = 0. Remove it before ridge
  # scaling and active-set multiplier tests, especially at large lambda.
  rows <- H$mean(dim = 2, keepdim = TRUE)
  H <- H - rows - rows$t() + H$mean()
  model_g <- ev$grad[F] + as.numeric(H$matmul(as_cpu_double_tensor(base[F])))
  centered_g <- model_g - mean(model_g)
  hscale <- max(abs(as.numeric(H$diag())), .Machine$double.eps)
  eye <- torch::torch_eye(length(F), dtype = torch::torch_float64())
  shifts <- c(0, hscale * 10^(-10:0))
  # The last shift that worked is a good guess, so only try one step below it.
  first <- if (is.null(state$shift)) 1L else max(1L, state$shift - 1L)
  for (i in first:length(shifts)) {
    dc <- tryCatch(
      cot_polish_exact(
        H + eye * shifts[[i]],
        centered_g,
        a[F] + base[F],
        ctl,
        state
      ),
      error = function(e) NULL
    )
    if (!is.null(dc) && all(is.finite(dc))) {
      d <- base
      d[F] <- d[F] + dc
      if (sum(ev$grad * d) < 0) {
        state$shift <- i
        return(d)
      }
    }
  }
  base
}


# Monotone Newton polish. Uses exact Newton steps on the free set, or L-BFGS
# when the free set exceeds max_free. Returns float64 weights; if the polish
# cannot improve an unconverged start, the start is returned unchanged.
cot_newton_polish <- function(twin, ctl = cot_polish_ctl) {
  a <- as.numeric(twin$a)
  if (any(!is.finite(a)) || any(a < 0) || sum(a) <= 0) {
    stop("Polish requires probability weights")
  }
  a <- a / sum(a)
  lambda <- as.numeric(twin$penalty)
  scale <- cot_conv_scale(twin)

  # Loss resolution, and the gap resolution: a potential error of size delta
  # moves the computed gap by up to 2 * (delta_xy + delta_xx), and each delta is
  # about the Sinkhorn tolerance over (1 - contraction rate). 16 covers rates up
  # to 3/4.
  floor <- ctl$round_factor * .Machine$double.eps * scale
  atol <- floor + ctl$sinkhorn_factor * cot_polish_sinkhorn_tol(lambda, ctl)
  fw_gap <- function(w, g) max(0, sum(w * g) - min(g))
  evaluate <- function(w) cot_polish_eval(twin, w, ctl)

  initial <- a
  start <- ev <- evaluate(a)
  y_self <- start$y_self
  done <- function(gap, loss) gap <= atol + ctl$tol * max(0, loss - y_self)

  gap_start <- fw_gap(a, ev$grad)
  status <- "maxit"
  accepted <- iter <- 0L
  memory <- list()
  prev_a <- prev_g <- NULL
  last_dir <- ""
  qp_state <- new.env(parent = emptyenv())
  qp_state$binding <- logical(length(a))
  # L-BFGS needs many more iterations than Newton
  maxit <- if (sum(a > 0) > ctl$max_free) ctl$lbfgs_maxit else ctl$maxit
  for (k in seq_len(maxit)) {
    gap <- fw_gap(a, ev$grad)
    if (done(gap, ev$loss)) {
      break
    }

    iter <- k
    sets <- cot_polish_sets(a, ev$grad, scale, ctl)
    F <- sets$free
    if (!identical(qp_state$indices, F)) {
      qp_state$active <- NULL
      qp_state$indices <- F
    }
    dir <- if (length(F) > ctl$max_free) "lbfgs" else "newton"
    if (dir != last_dir) {
      memory <- list()
      prev_a <- prev_g <- NULL
    }
    last_dir <- dir
    base <- numeric(length(a))
    base[sets$binding] <- -a[sets$binding]
    base[F] <- -sum(base) / length(F)
    if (dir == "lbfgs") {
      # Pairs are stored at full length and restricted to the current free set,
      # because the free set changes on most steps.
      if (!is.null(prev_a)) {
        memory <- utils::tail(
          c(memory, list(list(s = a - prev_a, y = ev$grad - prev_g))),
          ctl$lbfgs_m
        )
      }
      pairs <- list()
      for (p in memory) {
        s <- p$s[F] - mean(p$s[F])
        y <- p$y[F] - mean(p$y[F])
        sy <- sum(s * y)
        if (sy > 1e-12 * sqrt(sum(s^2) * sum(y^2))) {
          pairs <- c(pairs, list(list(s = s, y = y, sy = sy)))
        }
      }
      d <- base
      d[F] <- d[F] + cot_polish_lbfgs(ev$grad[F], pairs, mean(a))
      # The first step is scaled to probability-mass units.
      if (!length(pairs) && max(abs(d)) > mean(a)) {
        d <- d * mean(a) / max(abs(d))
      }
      prev_a <- a
      prev_g <- ev$grad
    } else {
      d <- cot_polish_newton(twin, a, ev, base, F, qp_state, ctl)
    }
    if (!all(is.finite(d)) || sum(ev$grad * d) >= 0) {
      status <- "no_descent"
      break
    }

    alpha <- 1
    moved <- FALSE
    while (alpha >= ctl$alpha_min) {
      next_ev <- NULL
      slope <- NA_real_
      cand <- pmax(a + alpha * d, 0)
      cand[cand < ctl$zero_tol] <- 0
      if (sum(cand) > 0) {
        cand <- cand / sum(cand)
        slope <- sum((ev$grad - sum(a * ev$grad)) * (cand - a))
        if (slope < 0) {
          next_ev <- evaluate(cand)
          new_gap <- fw_gap(cand, next_ev$grad)
          armijo_ok <- next_ev$loss <= ev$loss + ctl$armijo * slope &&
            (next_ev$loss < ev$loss - floor || new_gap < gap)
          # flat loss to rounding but gap halved, take step anyways
          flat_ok <- next_ev$loss <= ev$loss + floor && new_gap < 0.5 * gap
          if (is.finite(next_ev$loss) && (armijo_ok || flat_ok)) {
            moved <- TRUE
            break
          }
        }
      }
      shrink <- 0.5
      if (dir == "lbfgs" && !length(pairs) && !is.null(next_ev) && slope < 0) {
        # A rejected probe also supplies a secant curvature estimate, which
        # avoids many halvings on a sharply curved face.
        curvature <- sum((next_ev$grad - ev$grad) * (cand - a))
        if (is.finite(curvature) && curvature > 0) {
          shrink <- max(0.05, min(0.5, -slope / curvature))
        }
      }
      alpha <- alpha * shrink
    }
    if (!moved) {
      status <- "line_search_failed"
      break
    }
    improvement <- ev$loss - next_ev$loss
    a <- cand
    ev <- next_ev

    accepted <- accepted + 1L
    # Stalled only when neither the loss nor the gap made real progress.
    if (improvement <= floor && !(fw_gap(a, ev$grad) < 0.5 * gap)) {
      status <- "stalled"
      break
    }
  }

  gap <- fw_gap(a, ev$grad)
  if (done(gap, ev$loss)) {
    status <- if (accepted > 0L) "refined" else "already_optimal"
  } else if (start$loss - ev$loss <= floor) {
    a <- initial
    ev <- start
    gap <- gap_start
    status <- "no_improvement"
  }
  list(
    weights = a,
    loss = ev$loss,
    divergence = ev$loss - y_self,
    gap = gap,
    changed = !identical(a, initial),
    status = status,
    iterations = iter,
    atol = atol
  )
}
