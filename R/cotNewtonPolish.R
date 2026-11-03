# Float64 CPU refinement of a single, unconstrained primal COT problem.
# All tensors are detached; neither the evaluator nor the optimizer mutates OT.
cot_polish_control <- function(...) {
  defaults <- list(max_free = 2000L, maxit = 20L, tol = 1e-8,
    sinkhorn_tol = 1e-13, sinkhorn_niter = 10000L, armijo = 1e-4,
    alpha_min = 1e-8, zero_tol = 1e-12, release_margin = 0,
    release_max = NULL, lbfgs_m = 10L, osqp = list(), block_size = 512L,
    conv_scale = NULL, active_tol = NULL, roundoff = 64,
    resolution_tol = 1e-4, weight_eps = .Machine$double.eps,
    warm_start = TRUE, marginal_tol = 1e-13)
  args <- list(...)
  if (length(args) && (is.null(names(args)) || anyNA(names(args)) || any(!nzchar(names(args))) || anyDuplicated(names(args))))
    stop("Polish controls must have unique, nonempty names")
  unknown <- setdiff(names(args), names(defaults))
  if (length(unknown)) stop("Unknown polish control: ", paste(unknown, collapse = ", "))
  defaults[names(args)] <- args
  for (nm in c("max_free", "maxit", "sinkhorn_niter", "lbfgs_m", "block_size")) {
    x <- defaults[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x))
      stop("polish.control$", nm, " must be a positive integer")
  }
  for (nm in c("sinkhorn_tol", "armijo", "alpha_min", "zero_tol", "release_margin", "tol", "conv_scale", "release_max", "active_tol", "roundoff", "resolution_tol", "weight_eps")) {
    x <- defaults[[nm]]
    if (!is.null(x) && (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0))
      stop("polish.control$", nm, " must be a finite nonnegative scalar")
  }
  if (defaults$armijo <= 0 || defaults$armijo >= 1 || defaults$alpha_min <= 0 || defaults$alpha_min > 1 || defaults$sinkhorn_tol <= 0 || defaults$zero_tol >= 1)
    stop("Invalid polish line-search or Sinkhorn controls")
  if (!is.null(defaults$release_max) && defaults$release_max != floor(defaults$release_max)) stop("release_max must be an integer")
  if (is.null(defaults$tol)) defaults$tol <- 1e-8
  if (defaults$weight_eps < .Machine$double.eps) stop("weight_eps must be at least double precision epsilon")
  if (defaults$roundoff < 1) stop("roundoff must be at least one")
  if (!is.logical(defaults$warm_start) || length(defaults$warm_start) != 1L || is.na(defaults$warm_start)) stop("warm_start must be TRUE or FALSE")
  if (!is.numeric(defaults$marginal_tol) || length(defaults$marginal_tol) != 1L || !is.finite(defaults$marginal_tol) || defaults$marginal_tol <= 0) stop("marginal_tol must be positive")
  if (!is.list(defaults$osqp)) stop("polish.control$osqp must be a list")
  defaults
}

cot_polish_tensor <- function(x) {
  if (inherits(x, "torch_tensor")) return(x$detach()$to(device = "cpu")$to(dtype = torch::torch_float64()))
  torch::torch_tensor(x, dtype = torch::torch_float64(), device = "cpu")
}

# Cost accessors avoid a full CPU copy of a GPU matrix and materialization of
# online costs. Online custom KeOps expressions are deliberately unsupported.
cot_polish_cost <- function(x, block_size = 512L) {
  if (inherits(x, "cot_polish_cost")) return(x)
  online <- NULL
  if (inherits(x, "costOnline")) {
    if (identical(x$algorithm, "user")) stop("Custom online cost expressions are not supported by polish")
    xx <- cot_polish_tensor(x$data$x); yy <- cot_polish_tensor(x$data$y); p <- x$p
    online <- list(x = xx, y = yy, p = as.numeric(p))
    nr <- nrow(xx); nc <- nrow(yy)
    get <- function(i, j) torch::torch_cdist(xx[i, , drop = FALSE], yy[j, , drop = FALSE], p = p)^p / p
  } else {
    if (inherits(x, "costTensor")) x <- x$data
    nr <- nrow(x); nc <- ncol(x)
    if (is.null(nr) || is.null(nc)) stop("Polish costs must be matrices or cost objects")
    get <- function(i, j) cot_polish_tensor(x[i, j, drop = FALSE])
  }
  dense <- if (!inherits(x, "costOnline") && max(nr, nc) <= block_size) get(seq_len(nr), seq_len(nc)) else NULL
  structure(list(nrow = nr, ncol = nc, get = get, block_size = block_size, dense = dense, online = online), class = "cot_polish_cost")
}

cot_polish_softmin <- function(cost, h, eps, transpose = FALSE) {
  nr <- if (transpose) cost$ncol else cost$nrow
  nc <- if (transpose) cost$nrow else cost$ncol
  out <- torch::torch_empty(nr, dtype = torch::torch_float64())
  for (lo in seq.int(1L, nr, by = cost$block_size)) {
    ii <- lo:min(nr, lo + cost$block_size - 1L)
    # Chunk both axes: memory stays bounded even for a very wide target.
    acc <- torch::torch_full(length(ii), -Inf, dtype = torch::torch_float64())
    for (lj in seq.int(1L, nc, by = cost$block_size)) {
      jj <- lj:min(nc, lj + cost$block_size - 1L)
      C <- if (transpose) cost$get(jj, ii)$t() else cost$get(ii, jj)
      v <- ((h[jj]$view(c(1L, -1L)) - C) / eps)$logsumexp(2)
      acc <- torch::torch_logaddexp(acc, v)
    }
    out[ii] <- -eps * acc
  }
  out
}

cot_polish_diameter <- function(cost) {
  if (!is.null(cost$diameter)) return(cost$diameter)
  ans <- 0
  for (lo in seq.int(1L, cost$nrow, by = cost$block_size)) {
    for (lj in seq.int(1L, cost$ncol, by = cost$block_size)) {
      C <- cost$get(lo:min(cost$nrow, lo + cost$block_size - 1L), lj:min(cost$ncol, lj + cost$block_size - 1L))
      if (!all(is.finite(as.numeric(C)))) stop("Polish costs must be finite")
      ans <- max(ans, C$max()$item())
    }
  }
  ans
}

# Compile the small dense evaluator once. Keeping the inner Sinkhorn loop in
# torch avoids an R/C++ round trip for every update; large/online costs retain
# the bounded-memory block evaluator below.
cot_polish_dense <- local({
  compiled <- NULL
  function() {
    if (is.null(compiled)) compiled <<- torch::jit_compile("
def polish_online_softmin(x: Tensor, y: Tensor, h: Tensor, eps: float, p: float, block: int):
  out = torch.empty(x.size(0), dtype=x.dtype, device=x.device)
  for i in range(0, x.size(0), block):
    xi = x[i:i+block]
    acc = torch.full([xi.size(0)], -float('inf'), dtype=x.dtype, device=x.device)
    for j in range(0, y.size(0), block):
      yj = y[j:j+block]
      if p == 2.0:
        C = torch.clamp((xi * xi).sum(1).view(-1, 1) / 2.0 + (yj * yj).sum(1).view(1, -1) / 2.0 - xi.matmul(yj.t()), min=0.0)
      else:
        C = torch.cdist(xi, yj, p=p).pow(p) / p
      v = ((h[j:j+block].view(1, -1) - C) / eps).logsumexp(1)
      acc = torch.logaddexp(acc, v)
    out[i:i+block] = -eps * acc
  return out

def polish_online(a: Tensor, b: Tensor, x: Tensor, y: Tensor, p: float, path: Tensor, tol: float, maxit: int, block: int):
  la = a.log()
  lb = b.log()
  f = torch.zeros_like(a)
  g = torch.zeros_like(b)
  s = torch.zeros_like(a)
  done_xy = False
  done_xx = False
  iterations = 0
  for k in range(path.numel()):
    eps = float(path[k])
    final = k == path.numel() - 1
    done_xy = False
    done_xx = False
    for it in range(maxit if final else 5):
      iterations += 1
      if not done_xy:
        fn = polish_online_softmin(x, y, g + eps * lb, eps, p, block)
        fn = fn - a.dot(fn)
        gn = polish_online_softmin(y, x, fn + eps * la, eps, p, block)
        done_xy = final and bool((fn - f).abs().max() <= tol)
        f = fn
        g = gn
      if not done_xx:
        sn = (s + polish_online_softmin(x, x, s + eps * la, eps, p, block)) / 2.0
        done_xx = final and bool((sn - s).abs().max() <= tol)
        s = sn
      if done_xy and done_xx:
        break
  assert done_xy and done_xx, 'Polish Sinkhorn evaluator did not converge'
  return {'f_xy': f, 'g_yx': g, 'f_xx': s, 'iterations': torch.tensor(iterations)}

def polish_self_dense(b: Tensor, C: Tensor, path: Tensor, tol: float, maxit: int):
  lb = b.log()
  f = torch.zeros_like(b)
  done = False
  for k in range(path.numel()):
    eps = float(path[k])
    final = k == path.numel() - 1
    for it in range(maxit if final else 5):
      fn = (f - eps * ((f + eps * lb).view(1, -1) / eps - C / eps).logsumexp(1)) / 2.0
      done = final and bool((fn - f).abs().max() <= tol)
      f = fn
      if done:
        break
  assert done, 'Polish target self evaluator did not converge'
  return f

def polish_self_online(b: Tensor, x: Tensor, p: float, path: Tensor, tol: float, maxit: int, block: int):
  lb = b.log()
  f = torch.zeros_like(b)
  done = False
  for k in range(path.numel()):
    eps = float(path[k])
    final = k == path.numel() - 1
    for it in range(maxit if final else 5):
      fn = (f + polish_online_softmin(x, x, f + eps * lb, eps, p, block)) / 2.0
      done = final and bool((fn - f).abs().max() <= tol)
      f = fn
      if done:
        break
  assert done, 'Polish target self evaluator did not converge'
  return f


def polish_warm_dense(a: Tensor, b: Tensor, C: Tensor, D: Tensor, eps: float, f: Tensor, g: Tensor, s: Tensor, tol: float, residual_tol: float, maxit: int, accelerate: bool):
  la = a.log()
  lb = b.log()
  done_xy = False
  done_xx = False
  iterations = 0
  for it in range(maxit):
    iterations += 1
    if not done_xy:
      fn = -eps * ((g + eps * lb).view(1, -1) / eps - C / eps).logsumexp(1)
      fn = fn - a.dot(fn)
      gn = -eps * ((fn + eps * la).view(1, -1) / eps - C.t() / eps).logsumexp(1)
      row_f = -eps * ((gn + eps * lb).view(1, -1) / eps - C / eps).logsumexp(1)
      # Column marginals equal b by construction of gn.
      error = (a * ((fn - row_f) / eps).exp() - a).abs().sum()
      done_xy = bool(error <= tol)
      f = fn
      g = gn
    if not done_xx:
      target = -eps * ((s + eps * la).view(1, -1) / eps - D / eps).logsumexp(1)
      done_xx = bool((target - s).abs().max() <= residual_tol)
      # Remove the scalar oscillation exactly. On smooth, high-epsilon
      # problems the remaining self map contracts much faster without damping.
      s = target + a.dot(s - target) / 2.0 if accelerate else (s + target) / 2.0
    if done_xy and done_xx:
      break
  assert done_xy and done_xx, 'Warm polish evaluator did not converge'
  return {'f_xy': f, 'g_yx': g, 'f_xx': s, 'iterations': torch.tensor(iterations)}

def polish_warm_online(a: Tensor, b: Tensor, x: Tensor, y: Tensor, p: float, block: int, eps: float, f: Tensor, g: Tensor, s: Tensor, tol: float, residual_tol: float, maxit: int, accelerate: bool):
  la = a.log()
  lb = b.log()
  done_xy = False
  done_xx = False
  iterations = 0
  for it in range(maxit):
    iterations += 1
    if not done_xy:
      fn = polish_online_softmin(x, y, g + eps * lb, eps, p, block)
      fn = fn - a.dot(fn)
      gn = polish_online_softmin(y, x, fn + eps * la, eps, p, block)
      row_f = polish_online_softmin(x, y, gn + eps * lb, eps, p, block)
      # Column marginals equal b by construction of gn.
      error = (a * ((fn - row_f) / eps).exp() - a).abs().sum()
      done_xy = bool(error <= tol)
      f = fn
      g = gn
    if not done_xx:
      target = polish_online_softmin(x, x, s + eps * la, eps, p, block)
      done_xx = bool((target - s).abs().max() <= residual_tol)
      # Remove the scalar oscillation exactly. On smooth, high-epsilon
      # problems the remaining self map contracts much faster without damping.
      s = target + a.dot(s - target) / 2.0 if accelerate else (s + target) / 2.0
    if done_xy and done_xx:
      break
  assert done_xy and done_xx, 'Warm polish evaluator did not converge'
  return {'f_xy': f, 'g_yx': g, 'f_xx': s, 'iterations': torch.tensor(iterations)}

def polish_dense(a: Tensor, b: Tensor, C: Tensor, D: Tensor, path: Tensor, tol: float, maxit: int):
  la = a.log()
  lb = b.log()
  f = torch.zeros_like(a)
  g = torch.zeros_like(b)
  s = torch.zeros_like(a)
  done_xy = False
  done_xx = False
  iterations = 0
  for k in range(path.numel()):
    eps = float(path[k])
    final = k == path.numel() - 1
    done_xy = False
    done_xx = False
    steps = maxit if final else 5
    for it in range(steps):
      iterations += 1
      if not done_xy:
        fn = -eps * ((g + eps * lb).view(1, -1) / eps - C / eps).logsumexp(1)
        fn = fn - a.dot(fn)
        gn = -eps * ((fn + eps * la).view(1, -1) / eps - C.t() / eps).logsumexp(1)
        done_xy = final and bool((fn - f).abs().max() <= tol)
        f = fn
        g = gn
      if not done_xx:
        sn = (s - eps * ((s + eps * la).view(1, -1) / eps - D / eps).logsumexp(1)) / 2.0
        done_xx = final and bool((sn - s).abs().max() <= tol)
        s = sn
      if done_xy and done_xx:
        break
  assert done_xy and done_xx, 'Polish Sinkhorn evaluator did not converge'
  return {'f_xy': f, 'g_yx': g, 'f_xx': s, 'iterations': torch.tensor(iterations)}
")
    compiled
  }
})

# The target-only constant is solved once, in the same float64 arithmetic as
# the weight objective. It cannot be inferred from arbitrary cross costs.
cot_polish_self <- function(b, C_yy, lambda, control) {
  b <- cot_polish_tensor(b)
  yy <- cot_polish_cost(C_yy, control$block_size)
  diameter <- cot_polish_diameter(yy)
  control$sinkhorn_tol <- max(control$sinkhorn_tol, 8 * .Machine$double.eps * lambda)
  path <- if (lambda >= diameter) lambda else unique(pmax(c(diameter * exp(-0.9 * (0:ceiling(log(diameter / lambda) / 0.9))), lambda), lambda))
  if (!is.null(yy$dense)) {
    f <- cot_polish_dense()$polish_self_dense(b, yy$dense, cot_polish_tensor(path),
      torch::jit_scalar(control$sinkhorn_tol), torch::jit_scalar(as.integer(control$sinkhorn_niter)))
  } else if (!is.null(yy$online) && identical(as.matrix(yy$online$x), as.matrix(yy$online$y))) {
    f <- cot_polish_dense()$polish_self_online(b, yy$online$x, torch::jit_scalar(yy$online$p),
      cot_polish_tensor(path), torch::jit_scalar(control$sinkhorn_tol),
      torch::jit_scalar(as.integer(control$sinkhorn_niter)), torch::jit_scalar(as.integer(control$block_size)))
  } else {
    f <- torch::torch_zeros_like(b); lb <- b$log(); done <- FALSE
    for (eps in path) {
      final <- eps == lambda
      for (it in seq_len(if (final) control$sinkhorn_niter else 5L)) {
        fn <- (f + cot_polish_softmin(yy, f + eps * lb, eps)) / 2
        done <- final && (fn - f)$abs()$max()$item() <= control$sinkhorn_tol
        f <- fn
        if (done) break
      }
    }
    if (!done) stop("Polish target self evaluator did not converge")
  }
  b$dot(f)$item()
}

cot_polish_eval <- function(a, b, C_xy, C_xx, lambda, control = cot_polish_control(), diameter = NULL) {
  a <- cot_polish_tensor(a); b <- cot_polish_tensor(b)
  a <- a / a$sum(); b <- b / b$sum()
  xy <- cot_polish_cost(C_xy, control$block_size); xx <- cot_polish_cost(C_xx, control$block_size)
  if (is.null(diameter)) diameter <- max(cot_polish_diameter(xy), cot_polish_diameter(xx))
  control$sinkhorn_tol <- max(control$sinkhorn_tol, 8 * .Machine$double.eps * lambda)
  path <- if (lambda >= diameter) lambda else unique(pmax(c(diameter * exp(-0.9 * (0:ceiling(log(diameter / lambda) / 0.9))), lambda), lambda))
  if (!is.null(xy$online) && !is.null(xx$online) &&
      identical(xy$online$p, xx$online$p) &&
      identical(as.matrix(xy$online$x), as.matrix(xx$online$x)) &&
      identical(as.matrix(xx$online$x), as.matrix(xx$online$y))) {
    pot <- cot_polish_dense()$polish_online(a, b, xy$online$x, xy$online$y,
      torch::jit_scalar(xy$online$p), cot_polish_tensor(path),
      torch::jit_scalar(control$sinkhorn_tol), torch::jit_scalar(as.integer(control$sinkhorn_niter)),
      torch::jit_scalar(as.integer(control$block_size)))
    grad <- pot$f_xy - pot$f_xx
    return(list(loss = (a$dot(grad) + b$dot(pot$g_yx))$item(), grad = as.numeric(grad), pot = pot, iterations = as.integer(pot$iterations)))
  }
  if (!is.null(xy$dense) && !is.null(xx$dense)) {
    pot <- cot_polish_dense()$polish_dense(a, b, xy$dense, xx$dense,
      cot_polish_tensor(path), torch::jit_scalar(control$sinkhorn_tol),
      torch::jit_scalar(as.integer(control$sinkhorn_niter)))
    grad <- pot$f_xy - pot$f_xx
    return(list(loss = (a$dot(grad) + b$dot(pot$g_yx))$item(),
      grad = as.numeric(grad), pot = pot, iterations = as.integer(pot$iterations)))
  }
  la <- a$log(); lb <- b$log()
  f <- torch::torch_zeros_like(a); g <- torch::torch_zeros_like(b); s <- torch::torch_zeros_like(a)
  for (eps in path) {
    final <- eps == lambda
    cross_done <- self_done <- FALSE
    for (it in seq_len(if (final) control$sinkhorn_niter else 5L)) {
      if (!cross_done) {
        fn <- cot_polish_softmin(xy, g + eps * lb, eps)
        fn <- fn - a$dot(fn)
        gn <- cot_polish_softmin(xy, fn + eps * la, eps, TRUE)
        err <- (fn - f)$abs()$max()$item()
        f <- fn; g <- gn
        cross_done <- final && err <= control$sinkhorn_tol
      }
      if (!self_done) {
        sn <- (s + cot_polish_softmin(xx, s + eps * la, eps)) / 2
        err <- (sn - s)$abs()$max()$item()
        s <- sn
        self_done <- final && err <= control$sinkhorn_tol
      }
      if (cross_done && self_done) break
    }
  }
  if (!cross_done || !self_done) stop("Polish Sinkhorn evaluator did not converge")
  list(loss = (a$dot(f - s) + b$dot(g))$item(), grad = as.numeric(f - s),
       pot = list(f_xy = f, g_yx = g, f_xx = s))
}

# Exact candidate block. S includes every strictly positive mass so the block
# triangular identity remains exact even below the support selection tolerance.
cot_polish_hessian <- function(a, b, C_xy, C_xx, lambda, pot, candidates = seq_along(a), block_size = 512L, support = which(a > 0)) {
  xy <- cot_polish_cost(C_xy, block_size); xx <- cot_polish_cost(C_xx, block_size)
  S <- support; C <- unique(c(S, candidates)); ns <- length(S); n <- length(C)
  R <- torch::torch_zeros(c(n, n), dtype = torch::torch_float64())
  for (lo in seq.int(1L, xy$ncol, by = block_size)) {
    j <- lo:min(xy$ncol, lo + block_size - 1L)
    K <- ((pot$f_xy[C]$view(c(-1L, 1L)) + pot$g_yx[j]$view(c(1L, -1L)) - xy$get(C, j)) / lambda)$exp()
    R <- R + (K * cot_polish_tensor(b[j])$view(c(1L, -1L)))$matmul(K$t())
  }
  Kxx <- ((pot$f_xx[C]$view(c(-1L, 1L)) + pot$f_xx[C]$view(c(1L, -1L)) - xx$get(C, C)) / lambda)$exp()
  si <- seq_len(ns); av <- cot_polish_tensor(a[S])$view(c(1L, -1L))
  eye <- torch::torch_eye(ns, dtype = torch::torch_float64())
  XS <- torch::linalg_solve(eye - R[si, si, drop = FALSE] * av + av, R[si, , drop = FALSE])
  YS <- torch::linalg_solve(eye + Kxx[si, si, drop = FALSE] * av, Kxx[si, , drop = FALSE])
  X <- XS; Y <- YS
  if (n > ns) {
    zi <- (ns + 1L):n
    X <- torch::torch_cat(list(XS, R[zi, , drop = FALSE] - (-R[zi, si, drop = FALSE] * av + av)$matmul(XS)), dim = 1L)
    Y <- torch::torch_cat(list(YS, Kxx[zi, , drop = FALSE] - (Kxx[zi, si, drop = FALSE] * av)$matmul(YS)), dim = 1L)
  }
  H <- lambda * (X + Y); H <- (H + H$t()) / 2
  idx <- match(candidates, C)
  H[idx, idx, drop = FALSE]
}

cot_polish_basis <- function(n) {
  # Householder reflection maps e_1 to the normalized all-ones vector.
  if (n <= 1L) return(matrix(numeric(), n, 0L))
  v <- rep(1 / sqrt(n), n); v[1L] <- v[1L] + 1
  Q <- diag(n) - 2 * tcrossprod(v) / sum(v^2)
  Q[, -1L, drop = FALSE]
}

# Feasible QP on the epsilon-free model. Retain its binding set across
# Levenberg shifts and outer steps instead of restarting each factorization.
cot_polish_exact <- function(H, g, a, zero_tol, state = new.env(parent = emptyenv())) {
  n <- length(a)
  binding <- if (length(state$active) == n) state$active else a <= zero_tol & g > sum(a * g) / sum(a)
  if (all(binding)) binding[which.min(g)] <- FALSE
  kkt_tol <- 128 * .Machine$double.eps * max(1, max(abs(g)))
  seen <- character()
  for (attempt in seq_len(20L)) {
    key <- paste(which(binding), collapse = ",")
    if (key %in% seen) break
    seen <- c(seen, key)
    state$active <- binding
    F <- which(!binding)
    if (!length(F)) {binding[which.min(g)] <- FALSE; next}
    trial <- -a; trial[F] <- sum(a[binding]) / length(F)
    if (length(F) > 1L) {
      Z <- cot_polish_tensor(cot_polish_basis(length(F)))
      L <- torch::linalg_cholesky(Z$t()$matmul(H[F, F, drop = FALSE])$matmul(Z))
      rhs <- -Z$t()$matmul(cot_polish_tensor(g[F] + as.numeric(H[F, , drop = FALSE]$matmul(cot_polish_tensor(trial))))$view(c(-1L, 1L)))
      trial[F] <- trial[F] + as.numeric(Z$matmul(torch::torch_cholesky_solve(rhs, L)))
      state$n_factor <- if (is.null(state$n_factor)) 1L else state$n_factor + 1L
    }
    negative <- F[a[F] + trial[F] < -1e-12]
    if (length(negative)) {binding[negative] <- TRUE; next}
    q <- g + as.numeric(H$matmul(cot_polish_tensor(trial)))
    release <- which(binding & q < mean(q[F]) - kkt_tol)
    if (!length(release)) {state$active <- binding; return(trial)}
    binding[release] <- FALSE
  }
  if (all(binding)) binding[which.min(g)] <- FALSE
  F <- which(!binding)
  d <- -a; d[F] <- sum(a[binding]) / length(F)
  for (it in seq_len(4L * n + 4L)) {
    state$active <- binding
    F <- which(!binding)
    q <- g + as.numeric(H$matmul(cot_polish_tensor(d)))
    nu <- mean(q[F])
    if (max(abs(q[F] - nu)) <= kkt_tol) {
      release <- which(binding & q < nu - kkt_tol)
      if (!length(release)) {state$active <- binding; return(d)}
      binding[release[which.min(q[release])]] <- FALSE
      next
    }
    p <- numeric(n)
    if (length(F) > 1L) {
      Z <- cot_polish_tensor(cot_polish_basis(length(F)))
      L <- torch::linalg_cholesky(Z$t()$matmul(H[F, F, drop = FALSE])$matmul(Z))
      rhs <- -Z$t()$matmul(cot_polish_tensor(q[F])$view(c(-1L, 1L)))
      p[F] <- as.numeric(Z$matmul(torch::torch_cholesky_solve(rhs, L)))
      state$n_factor <- if (is.null(state$n_factor)) 1L else state$n_factor + 1L
    }
    if (max(abs(p)) <= .Machine$double.eps) stop("Polish QP needs a larger Levenberg shift")
    hits <- which(p < 0 & !binding)
    ratios <- (a[hits] + d[hits]) / -p[hits]
    alpha <- if (length(hits)) min(1, max(0, min(ratios))) else 1
    d <- d + alpha * p
    if (alpha < 1) {
      hit <- hits[which.min(ratios)]
      d[hit] <- -a[hit]; binding[hit] <- TRUE
    }
  }
  state$active <- binding
  stop("Polish QP needs a larger Levenberg shift")
}

cot_polish_sets <- function(a, g, gap, control, scale = 1) {
  S <- which(a > control$zero_tol)
  nu <- sum(a[S] * g[S]) / sum(a[S])
  epsilon_max <- if (is.null(control$active_tol)) 1e-2 / length(a) else control$active_tol
  # Bertsekas' epsilon shrinks with the projected-gradient residual. The
  # inverse cost scale gives a dimensionless gradient step on the simplex.
  v <- a - (g - nu) / max(scale, .Machine$double.eps)
  u <- sort(v, decreasing = TRUE)
  theta <- (cumsum(u) - 1) / seq_along(u)
  rho <- max(which(u > theta))
  projected <- pmax(v - theta[rho], 0)
  epsilon <- min(epsilon_max, sqrt(sum((a - projected)^2)))
  # Avoid declaring a small optimal weight binding because of potential noise.
  margin <- 8 * control$sinkhorn_tol
  binding <- which(a <= epsilon & g > nu + margin)
  zeros <- which(a <= control$zero_tol & g < nu - control$release_margin * gap)
  zeros <- zeros[order(g[zeros])]
  if (!is.null(control$release_max)) zeros <- utils::head(zeros, control$release_max)
  free <- sort(unique(c(setdiff(which(a > 0), binding), zeros)))
  # The iterate and its potentials retain ALL current mass. Bound coordinates
  # are held out of the Hessian; they may move only in the projected step.
  list(free = free, binding = binding, weights = a, nu = nu)
}

cot_polish_threshold <- function(loss, y_self, scale, control) {
  control$tol * max(0, loss - y_self)
}

cot_polish_resolution_threshold <- function(loss, y_self, scale, control) {
  max(control$roundoff * .Machine$double.eps * scale,
      control$resolution_tol * max(0, loss - y_self, control$roundoff * control$weight_eps * scale))
}

# Shared certificate for phase 1 and refinement. A resolution certificate needs
# evidence of a stalled evaluator, except at the double arithmetic floor.
cot_certify <- function(gap, divergence, scale, weight_eps, tol, resolution = FALSE, roundoff = 64, resolution_tol = 1e-4) {
  if (any(!is.finite(c(gap, divergence, scale, weight_eps, tol))) || gap < 0)
    return("uncertified")
  if (gap <= tol * max(0, divergence)) return("certified")
  floor <- roundoff * .Machine$double.eps * scale
  threshold <- max(floor, resolution_tol * max(0, divergence, roundoff * weight_eps * scale))
  if ((resolution || gap <= floor) && gap <= threshold) return("certified_resolution")
  "uncertified"
}

cot_fit_valid <- function(weights, loss) {
  if (inherits(weights, "torch_tensor")) weights <- as.numeric(weights$detach()$to(device = "cpu"))
  length(weights) > 0L && all(is.finite(weights)) &&
    all(weights >= -64 * .Machine$double.eps) && sum(weights) > 0 &&
    length(loss) == 1L && is.finite(as.numeric(loss))
}

cot_polish_convergence <- function(gap, loss, y_self, scale, control, resolution = FALSE) {
  status <- cot_certify(gap, loss - y_self, scale, control$weight_eps, control$tol, resolution, control$roundoff, control$resolution_tol)
  if (status == "uncertified") NULL else status
}

cot_polish_is_converged <- function(status) {
  status %in% c("certified", "certified_resolution")
}

cot_polish_lbfgs <- function(g, memory, mean_weight = 1 / length(g)) {
  q <- g - mean(g); alpha <- numeric(length(memory))
  for (i in rev(seq_along(memory))) {
    pair <- memory[[i]]; alpha[i] <- sum(pair$s * q) / pair$sy
    q <- q - alpha[i] * pair$y
  }
  scale <- if (length(memory)) {p <- memory[[length(memory)]]; p$sy / sum(p$y^2)} else {
    magnitude <- max(abs(q))
    if (magnitude > 0) mean_weight / magnitude else 1
  }
  r <- scale * q
  for (i in seq_along(memory)) {
    pair <- memory[[i]]; beta <- sum(pair$y * r) / pair$sy
    r <- r + pair$s * (alpha[i] - beta)
  }
  -r + mean(r)
}

cot_newton_polish <- function(a, b, C_xy, C_xx, lambda, potentials = NULL,
                              direction = c("exact", "osqp", "lbfgs"), control = cot_polish_control(), C_yy = NULL, cache = NULL) {
  requested_direction <- direction <- match.arg(direction); control <- do.call(cot_polish_control, control)
  started <- proc.time()[["elapsed"]]
  a <- as.numeric(cot_polish_tensor(a)); b <- as.numeric(cot_polish_tensor(b))
  if (any(!is.finite(a)) || any(a < 0) || sum(a) <= 0 || any(!is.finite(b)) || any(b < 0) || abs(sum(b) - 1) > 1e-6 || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0)
    stop("Polish requires probability weights and a finite positive lambda")
  solver_weights <- a
  a <- a / sum(a); b <- b / sum(b)
  quantize <- function(w) {
    if (control$weight_eps > .Machine$double.eps)
      as.numeric(torch::torch_tensor(w, dtype = torch::torch_float32())) else w
  }
  # Preserve the phase-1 point exactly for the cold baseline comparison. Only
  # proposed steps are quantized; normalizing then re-rounding the baseline
  # itself would silently change the answer against which monotonicity is tested.
  initial <- a; initial_solver_weights <- solver_weights
  xy <- cot_polish_cost(C_xy, control$block_size); xx <- cot_polish_cost(C_xx, control$block_size)
  if (xy$nrow != length(a) || xy$ncol != length(b) || xx$nrow != length(a) || xx$ncol != length(a)) stop("Incompatible polish cost dimensions")
  xy$diameter <- cot_polish_diameter(xy); xx$diameter <- cot_polish_diameter(xx)
  diameter <- max(xy$diameter, xx$diameter)
  scale <- if (is.null(control$conv_scale)) max(diameter, lambda * .Machine$double.eps, .Machine$double.eps) else control$conv_scale
  if (is.null(C_yy)) stop("C_yy is required for divergence-relative polish convergence")
  yy <- cot_polish_cost(C_yy, control$block_size)
  if (yy$nrow != length(b) || yy$ncol != length(b)) stop("Incompatible target self-cost dimensions")
  yy$diameter <- cot_polish_diameter(yy)
  scale <- max(scale, yy$diameter)
  ty <- proc.time()[["elapsed"]]
  cache_key <- as.character(lambda)
  y_self <- if (!is.null(cache) && !is.null(cache[[cache_key]])) cache[[cache_key]] else {
    value <- cot_polish_self(b, yy, lambda, control)
    if (!is.null(cache)) cache[[cache_key]] <- value
    value
  }
  y_time <- proc.time()[["elapsed"]] - ty
  threshold <- function(loss) cot_polish_threshold(loss, y_self, scale, control)
  resolution_threshold <- function(loss) cot_polish_resolution_threshold(loss, y_self, scale, control)
  loss_floor <- control$roundoff * .Machine$double.eps * scale
  classify <- function(gap, loss, resolution = FALSE) cot_polish_convergence(gap, loss, y_self, scale, control, resolution)
  evaluate <- function(w) cot_polish_eval(w, b, xy, xx, lambda, control, diameter)
  ts <- proc.time()[["elapsed"]]
  ev <- evaluate(a); start_ev <- ev
  ev_is_cold <- TRUE; readback_evals <- 0L
  start_time <- proc.time()[["elapsed"]] - ts
  gap0 <- max(0, sum(a * ev$grad) - min(ev$grad))
  accepted <- iter <- 0L; max_shift <- 0; status <- "maxit"
  loss_resolution <- FALSE
  memory <- list(); prevF <- integer(); prev_a <- prev_g <- NULL
  osqp_model <- NULL; osqp_n <- 0L
  qp_state <- new.env(parent = emptyenv())
  timings <- c(start = start_time, y_self = y_time, hessian = 0, factor = 0, line_search = 0, readback = 0)
  probe_timings <- numeric(); probe_iterations <- integer()
  direction_history <- character(); free_history <- integer(); binding_history <- integer()
  use_warm <- control$warm_start
  line_search_evals <- 0L; first_step_norm <- NA_real_
  for (k in seq_len(control$maxit)) {
    gap <- max(0, sum(a * ev$grad) - min(ev$grad))
    certified <- classify(gap, ev$loss)
    if (!is.null(certified)) {
      if (!ev_is_cold) {
        tr <- proc.time()[["elapsed"]]
        ev <- evaluate(a)
        ev_is_cold <- TRUE; readback_evals <- readback_evals + 1L
        timings["readback"] <- timings["readback"] + proc.time()[["elapsed"]] - tr
        use_warm <- FALSE
        gap <- max(0, sum(a * ev$grad) - min(ev$grad))
        certified <- classify(gap, ev$loss)
      }
      if (!is.null(certified)) {status <- certified; break}
    }
    iter <- k
    sets <- cot_polish_sets(a, ev$grad, gap, control, scale)
    C <- F <- sets$free
    if (!identical(qp_state$indices, C)) {
      qp_state$active <- NULL; qp_state$indices <- C
    }
    direction <- if (requested_direction != "lbfgs" && length(C) > control$max_free) "lbfgs" else requested_direction
    if (length(direction_history) && direction != utils::tail(direction_history, 1L)) {
      memory <- list(); prevF <- integer(); prev_a <- prev_g <- NULL
    }
    direction_history <- c(direction_history, direction)
    free_history <- c(free_history, length(C)); binding_history <- c(binding_history, length(sets$binding))
    base <- numeric(length(a))
    base[sets$binding] <- -a[sets$binding]
    base[C] <- -sum(base) / length(C)
    d <- base; shift <- 0
    t0 <- proc.time()[["elapsed"]]
    if (direction == "lbfgs") {
      if (!identical(F, prevF)) memory <- list()
      if (identical(F, prevF) && !is.null(prev_a)) {
        s <- a[F] - prev_a[F]; y <- ev$grad[F] - prev_g[F]
        s <- s - mean(s); y <- y - mean(y); sy <- sum(s * y)
        if (sy > 1e-12 * sqrt(sum(s^2) * sum(y^2))) memory <- utils::tail(c(memory, list(list(s = s, y = y, sy = sy))), control$lbfgs_m)
      }
      d[F] <- d[F] + cot_polish_lbfgs(ev$grad[F], memory, mean(a))
      if (!length(memory) && max(abs(d)) > mean(a)) d <- d * mean(a) / max(abs(d))
      prevF <- F; prev_a <- a; prev_g <- ev$grad
    } else {
      # Truncate only the tiny active columns in the Hessian operator. The
      # gradient, potentials, line-search point and objective keep their mass.
      # No full-source solve or extra cold sparse-model evaluation is needed.
      H <- tryCatch(cot_polish_hessian(a, b, xy, xx, lambda, ev$pot, C,
        control$block_size, support = intersect(which(a > 0), C)), error = function(e) NULL)
      timings["hessian"] <- timings["hessian"] + proc.time()[["elapsed"]] - t0
      t0 <- proc.time()[["elapsed"]]
      if (is.null(H)) {status <- "no_descent"; break}
      # The Hessian gauge is irrelevant on sum(d)=0. Remove it before ridge
      # scaling and active-set multiplier tests, especially at large lambda.
      rows <- H$mean(dim = 2, keepdim = TRUE)
      H <- H - rows - rows$t() + H$mean()
      model_g <- ev$grad[C] + as.numeric(H$matmul(cot_polish_tensor(base[C])))
      centered_g <- model_g - mean(model_g)
      hscale <- max(abs(as.numeric(H$diag())), .Machine$double.eps)
      repeat {
        max_shift <- max(max_shift, shift)
        HH <- H + shift * torch::torch_eye(length(C), dtype = torch::torch_float64())
        dc <- tryCatch({
          if (direction == "exact") cot_polish_exact(HH, centered_g, (a[C] + base[C]), control$zero_tol, qp_state) else {
            P <- Matrix::forceSymmetric(Matrix::Matrix(as.matrix(HH), sparse = TRUE), uplo = "U")
            AA <- rbind(Matrix::Matrix(1, 1, length(C), sparse = TRUE), Matrix::Diagonal(length(C)))
            pars <- utils::modifyList(list(verbose = FALSE, eps_abs = 1e-10, eps_rel = 1e-10, polishing = TRUE, max_iter = 100000L), control$osqp)
            modern <- utils::packageVersion("osqp") >= "1.0.0"
            if (!modern && !is.null(pars$polishing)) {pars$polish <- pars$polishing; pars$polishing <- NULL}
            # Reuse the model when the candidate count and sparse pattern agree.
            pattern <- list(P@p, P@i)
            if (is.null(osqp_model) || osqp_n != length(C) || !identical(pattern, osqp_pattern)) {
              osqp_model <- osqp::osqp(P = P, q = centered_g, A = AA, l = c(0, -(a[C] + base[C])), u = c(0, rep(Inf, length(C))), pars = do.call(osqp::osqpSettings, pars))
              osqp_n <- length(C); osqp_pattern <- pattern
            } else {
              update <- if (modern) osqp_model@Update else osqp_model$Update
              update(Px = P@x, q = centered_g, l = c(0, -(a[C] + base[C])), u = c(0, rep(Inf, length(C))))
            }
            solve <- if (modern) osqp_model@Solve else osqp_model$Solve
            ans <- solve()
            if (!(ans$info$status_val %in% c(1L, 2L))) stop("OSQP did not solve the polish QP")
            as.numeric(ans$x)
          }
        }, error = function(e) NULL)
        trial <- base
        if (!is.null(dc)) trial[C] <- trial[C] + dc
        if (!is.null(dc) && all(is.finite(dc)) && sum(ev$grad * trial) < 0) {d <- trial; break}
        shift <- if (shift == 0) 1e-10 * hscale else shift * 10
        if (shift > hscale * 1.01) break
      }
    }
    timings["factor"] <- timings["factor"] + proc.time()[["elapsed"]] - t0
    if (!all(is.finite(d)) || sum(ev$grad * d) >= 0) {status <- "no_descent"; break}
    if (k == 1L) first_step_norm <- max(abs(d))
    alpha <- 1; moved <- FALSE; t0 <- proc.time()[["elapsed"]]
    while (alpha >= control$alpha_min) {
      next_ev <- NULL; slope <- NA_real_
      candidate <- pmax(a + alpha * d, 0)
      candidate[candidate < control$zero_tol] <- 0
      if (sum(candidate) > 0) {
        candidate <- candidate / sum(candidate)
        candidate_solver_weights <- quantize(candidate)
        candidate <- candidate_solver_weights / sum(candidate_solver_weights)
        slope <- sum(ev$grad * (candidate - a))
        next_ev <- if (slope < 0) {
          line_search_evals <- line_search_evals + 1L
          tp <- proc.time()[["elapsed"]]
          probe_is_cold <- !use_warm
          failed_warm_work <- 0L
          probe <- tryCatch(if (use_warm)
            cot_polish_warm(candidate, b, xy, xx, lambda, ev$pot, control) else evaluate(candidate),
            error = function(e) {
              probe_is_cold <<- TRUE
              if (grepl("Warm polish evaluator did not converge", conditionMessage(e), fixed = TRUE))
                failed_warm_work <<- control$sinkhorn_niter
              evaluate(candidate)
            })
          probe_work <- failed_warm_work + if (is.null(probe$iterations)) 0L else probe$iterations
          # Close to resolution, the cold evaluator decides both Armijo and
          # certification; marginal residual noise must not choose the support.
          if (use_warm && abs(ev$loss - probe$loss) <= 1e-11 * scale) {
            probe <- evaluate(candidate)
            probe_work <- probe_work + if (is.null(probe$iterations)) 0L else probe$iterations
            probe_is_cold <- TRUE
            use_warm <- FALSE
          }
          probe_timings <- c(probe_timings, proc.time()[["elapsed"]] - tp)
          probe_iterations <- c(probe_iterations, as.integer(probe_work))
          probe
        } else NULL
        if (!is.null(next_ev) && is.finite(next_ev$loss) && next_ev$loss <= ev$loss &&
            next_ev$loss <= ev$loss + control$armijo * slope &&
            (next_ev$loss < ev$loss || sum(candidate * next_ev$grad) - min(next_ev$grad) < gap)) {moved <- TRUE; break}
      }
      shrink <- 0.5
      if (direction == "lbfgs" && !length(memory) && !is.null(next_ev) && is.finite(slope) && slope < 0) {
        # A rejected cold probe also supplies a secant curvature estimate.
        # Use it to avoid many expensive halvings on a sharply curved face.
        curvature <- sum((next_ev$grad - ev$grad) * (candidate - a))
        if (is.finite(curvature) && curvature > 0)
          shrink <- max(0.05, min(0.5, -slope / curvature))
      }
      alpha <- alpha * shrink
    }
    timings["line_search"] <- timings["line_search"] + proc.time()[["elapsed"]] - t0
    if (!moved) {status <- "line_search_failed"; loss_resolution <- TRUE; break}
    improvement <- ev$loss - next_ev$loss
    a <- candidate; solver_weights <- candidate_solver_weights
    ev <- next_ev; ev_is_cold <- probe_is_cold; accepted <- accepted + 1L
    new_gap <- max(0, sum(a * ev$grad) - min(ev$grad))
    precision_floor <- max(loss_floor, control$weight_eps * scale)
    if (improvement <= loss_floor ||
        (improvement <= precision_floor && new_gap <= resolution_threshold(ev$loss))) {
      loss_resolution <- TRUE
      new_gap <- max(0, sum(a * ev$grad) - min(ev$grad))
      certified <- classify(new_gap, ev$loss, TRUE)
      status <- if (is.null(certified)) "stalled_resolution" else certified
      break
    }
  }
  # Only a cold eps-scaled solve may certify the returned point.
  tr <- proc.time()[["elapsed"]]
  if (!ev_is_cold) {ev <- evaluate(a); readback_evals <- readback_evals + 1L}
  timings["readback"] <- timings["readback"] + proc.time()[["elapsed"]] - tr
  gap <- max(0, sum(a * ev$grad) - min(ev$grad))
  loss_resolution <- loss_resolution || status %in% c("line_search_failed", "no_descent")
  certified <- classify(gap, ev$loss, loss_resolution)
  if (!is.null(certified)) status <- certified
  noise <- 10 * .Machine$double.eps * max(scale, abs(start_ev$loss))
  if (!cot_polish_is_converged(status) && start_ev$loss - ev$loss <= noise) {
    a <- initial; solver_weights <- initial_solver_weights
    ev <- start_ev; gap <- gap0; status <- "no_improvement"
  }
  elapsed <- proc.time()[["elapsed"]] - started
  timings <- c(timings, other = max(0, elapsed - sum(timings)))
  list(weights = a, solver_weights = solver_weights, loss = ev$loss, loss_start = start_ev$loss, fw_gap_start = gap0,
    fw_gap = gap, iter = iter, n_accepted = accepted, status = status,
    nnz = sum(a > control$zero_tol), max_shift = max_shift, direction = direction,
    divergence = ev$loss - y_self, divergence_start = start_ev$loss - y_self,
    y_self = y_self, fw_threshold = threshold(ev$loss), fw_threshold_start = threshold(start_ev$loss),
    resolution_threshold = resolution_threshold(ev$loss),
    resolution_threshold_start = resolution_threshold(start_ev$loss),
    loss_resolution = loss_resolution, weight_eps = control$weight_eps,
    weight_roundoff_floor = control$roundoff * control$weight_eps * scale,
    conv_scale = scale, roundoff_floor = loss_floor,
    direction_history = direction_history, free_history = free_history, binding_history = binding_history,
    line_search_evals = line_search_evals, first_step_norm = first_step_norm,
    factorizations = if (is.null(qp_state$n_factor)) 0L else qp_state$n_factor,
    time = elapsed, timings = timings, readback_evals = readback_evals, probe_timings = probe_timings, probe_iterations = probe_iterations)
}

# Rebuild in bounded blocks from the pre-cast points. Custom costs are checked
# on every block; the hook falls back as a unit if any block is unsupported.
cot_polish_data_cost <- function(C, x, y, block_size = 512L) {
  x <- cot_polish_tensor(x); y <- cot_polish_tensor(y); p <- C$p
  if (inherits(C, "costOnline") && identical(C$algorithm, "user"))
    stop("Custom online cost expressions are unsupported")
  custom <- identical(C$algorithm, "user")
  get <- function(i, j) {
    value <- if (custom) C$fun(x[i, , drop = FALSE], y[j, , drop = FALSE], p) else
      torch::torch_cdist(x[i, , drop = FALSE], y[j, , drop = FALSE], p = p)^p / cot_polish_tensor(p)
    if (!inherits(value, "torch_tensor") || value$dtype != torch::torch_float64())
      stop("User cost did not return float64 tensors")
    if (!identical(as.integer(dim(value)), c(length(i), length(j)))) stop("Invalid user cost dimensions")
    value$detach()
  }
  cost <- structure(list(nrow = nrow(x), ncol = nrow(y), get = get,
    block_size = block_size, dense = NULL,
    online = if (custom) NULL else list(x = x, y = y, p = as.numeric(p))), class = "cot_polish_cost")
  if (max(nrow(x), nrow(y)) <= block_size) cost$dense <- get(seq_len(nrow(x)), seq_len(nrow(y)))
  cost$diameter <- cot_polish_diameter(cost)
  cost
}

cot_polish_warm <- function(a, b, xy, xx, lambda, potentials, control = cot_polish_control()) {
  a <- cot_polish_tensor(a); a <- a / a$sum()
  b <- cot_polish_tensor(b); b <- b / b$sum()
  tol <- torch::jit_scalar(control$marginal_tol)
  residual <- torch::jit_scalar(max(control$sinkhorn_tol, 8 * .Machine$double.eps * lambda))
  accelerate <- torch::jit_scalar(lambda >= cot_polish_diameter(xx) / 4)
  eps <- torch::jit_scalar(lambda); maxit <- torch::jit_scalar(as.integer(control$sinkhorn_niter))
  if (!is.null(xy$dense) && !is.null(xx$dense)) {
    pot <- cot_polish_dense()$polish_warm_dense(a, b, xy$dense, xx$dense, eps,
      potentials$f_xy, potentials$g_yx, potentials$f_xx, tol, residual, maxit, accelerate)
  } else if (!is.null(xy$online) && !is.null(xx$online)) {
    pot <- cot_polish_dense()$polish_warm_online(a, b, xy$online$x, xy$online$y,
      torch::jit_scalar(xy$online$p), torch::jit_scalar(as.integer(control$block_size)), eps,
      potentials$f_xy, potentials$g_yx, potentials$f_xx, tol, residual, maxit, accelerate)
  } else {
    f <- potentials$f_xy; g <- potentials$g_yx; s <- potentials$f_xx
    la <- a$log(); lb <- b$log(); done <- FALSE
    for (it in seq_len(control$sinkhorn_niter)) {
      f <- cot_polish_softmin(xy, g + lambda * lb, lambda); f <- f - a$dot(f)
      g <- cot_polish_softmin(xy, f + lambda * la, lambda, TRUE)
      row <- cot_polish_softmin(xy, g + lambda * lb, lambda)
      col <- cot_polish_softmin(xy, f + lambda * la, lambda, TRUE)
      err <- ((a * ((f - row) / lambda)$exp() - a)$abs()$sum() +
        (b * ((g - col) / lambda)$exp() - b)$abs()$sum())$item()
      target <- cot_polish_softmin(xx, s + lambda * la, lambda)
      self_err <- (target - s)$abs()$max()$item()
      s <- if (lambda >= cot_polish_diameter(xx) / 4) target + a$dot(s - target) / 2 else (s + target) / 2
      if (err <= control$marginal_tol && self_err <= max(control$sinkhorn_tol, 8 * .Machine$double.eps * lambda)) {done <- TRUE; break}
    }
    if (!done) stop("Warm polish evaluator did not converge")
    pot <- list(f_xy = f, g_yx = g, f_xx = s, iterations = it)
  }
  grad <- pot$f_xy - pot$f_xx
  list(loss = (a$dot(grad) + b$dot(pot$g_yx))$item(), grad = as.numeric(grad),
    pot = pot[c("f_xy", "g_yx", "f_xx")], iterations = as.integer(pot$iterations))
}
