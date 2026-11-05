# Float64 CPU refinement of a single, unconstrained primal COT problem.
# Weights and costs are converted to float64 on the CPU once, and the polished
# weights are returned in float64. All tensors are detached; nothing here
# mutates the OT object.

# Internal constants; these are not user settings.
cot_polish_ctl <- list(maxit = 20L, tol = 1e-8, sinkhorn_tol = 1e-13,
  sinkhorn_niter = 10000L, max_free = 2000L, armijo = 1e-4, alpha_min = 1e-8,
  zero_tol = 1e-12, lbfgs_m = 10L, block_size = 512L, marginal_tol = 1e-13,
  warm_start = TRUE)

cot_polish_tensor <- function(x) {
  # MPS has no float64, so move to the CPU before changing dtype.
  if (inherits(x, "torch_tensor")) return(x$detach()$to(device = "cpu")$to(dtype = torch::torch_float64()))
  torch::torch_tensor(x, dtype = torch::torch_float64(), device = "cpu")
}

# Stored costs become one dense float64 CPU matrix. Online costs keep float64
# copies of their points and are evaluated in blocks. Self costs computed in
# float32 can be slightly asymmetric, which breaks the debiased gradient, so
# `symmetric = TRUE` averages a stored matrix with its transpose.
cot_polish_cost <- function(x, block_size = cot_polish_ctl$block_size, symmetric = FALSE) {
  if (inherits(x, "cot_polish_cost")) return(x)
  if (inherits(x, "costOnline")) {
    if (identical(x$algorithm, "user")) stop("Custom online cost expressions are not supported by polish")
    xx <- cot_polish_tensor(x$data$x); yy <- cot_polish_tensor(x$data$y); p <- as.numeric(x$p)
    get <- function(i, j) torch::torch_cdist(xx[i, , drop = FALSE], yy[j, , drop = FALSE], p = p)^p / p
    cost <- list(nrow = nrow(xx), ncol = nrow(yy), get = get, dense = NULL,
      online = list(x = xx, y = yy, p = p))
  } else {
    if (inherits(x, "costTensor")) x <- x$data
    dense <- cot_polish_tensor(x)
    if (length(dim(dense)) != 2L) stop("Polish costs must be matrices or cost objects")
    if (symmetric) dense <- (dense + dense$t()) / 2
    cost <- list(nrow = nrow(dense), ncol = ncol(dense),
      get = function(i, j) dense[i, j, drop = FALSE], dense = dense, online = NULL)
  }
  cost$block_size <- block_size
  cost <- structure(cost, class = "cot_polish_cost")
  cost$diameter <- cot_polish_diameter(cost)
  cost
}

cot_polish_diameter <- function(cost) {
  if (!is.null(cost$diameter)) return(cost$diameter)
  blocks <- if (is.null(cost$dense)) {
    bs <- cost$block_size
    rows <- lapply(seq.int(1L, cost$nrow, by = bs), function(lo) lo:min(cost$nrow, lo + bs - 1L))
    cols <- lapply(seq.int(1L, cost$ncol, by = bs), function(lo) lo:min(cost$ncol, lo + bs - 1L))
    unlist(lapply(rows, function(i) lapply(cols, function(j) list(i, j))), recursive = FALSE)
  } else list(NULL)
  ans <- 0
  for (ij in blocks) {
    C <- if (is.null(ij)) cost$dense else cost$get(ij[[1L]], ij[[2L]])
    if (!torch::torch_isfinite(C)$all()$item()) stop("Polish costs must be finite")
    ans <- max(ans, C$max()$item())
  }
  ans
}

# Compiled evaluators. Keeping the inner Sinkhorn loop in TorchScript avoids an
# R/C++ round trip for every update. Online costs are evaluated in blocks.
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

cot_polish_path <- function(lambda, diameter) {
  if (lambda >= diameter) return(lambda)
  unique(pmax(c(diameter * exp(-0.9 * (0:ceiling(log(diameter / lambda) / 0.9))), lambda), lambda))
}

cot_polish_sinkhorn_tol <- function(lambda, ctl) max(ctl$sinkhorn_tol, 8 * .Machine$double.eps * lambda)

cot_polish_value <- function(a, b, pot) {
  grad <- pot$f_xy - pot$f_xx
  list(loss = (a$dot(grad) + b$dot(pot$g_yx))$item(), grad = as.numeric(grad),
    pot = pot[c("f_xy", "g_yx", "f_xx")])
}

cot_polish_online <- function(xy, xx) {
  if (!is.null(xy$dense) && !is.null(xx$dense)) return(FALSE)
  if (!is.null(xy$online) && !is.null(xx$online)) return(TRUE)
  stop("Polish costs must be all stored or all online")
}

# The target-only constant is solved once per lambda, in the same float64
# arithmetic as the weight objective.
cot_polish_self <- function(b, yy, lambda, ctl = cot_polish_ctl) {
  b <- cot_polish_tensor(b); yy <- cot_polish_cost(yy, ctl$block_size)
  path <- cot_polish_tensor(cot_polish_path(lambda, yy$diameter))
  tol <- torch::jit_scalar(cot_polish_sinkhorn_tol(lambda, ctl))
  niter <- torch::jit_scalar(as.integer(ctl$sinkhorn_niter))
  f <- if (is.null(yy$dense)) {
    cot_polish_dense()$polish_self_online(b, yy$online$x, torch::jit_scalar(yy$online$p),
      path, tol, niter, torch::jit_scalar(as.integer(yy$block_size)))
  } else cot_polish_dense()$polish_self_dense(b, yy$dense, path, tol, niter)
  b$dot(f)$item()
}

# Cold evaluation with epsilon scaling from the cost diameter.
cot_polish_eval <- function(a, b, xy, xx, lambda, ctl = cot_polish_ctl) {
  a <- cot_polish_tensor(a); b <- cot_polish_tensor(b)
  a <- a / a$sum(); b <- b / b$sum()
  path <- cot_polish_tensor(cot_polish_path(lambda, max(xy$diameter, xx$diameter)))
  tol <- torch::jit_scalar(cot_polish_sinkhorn_tol(lambda, ctl))
  niter <- torch::jit_scalar(as.integer(ctl$sinkhorn_niter))
  pot <- if (cot_polish_online(xy, xx)) {
    cot_polish_dense()$polish_online(a, b, xy$online$x, xy$online$y,
      torch::jit_scalar(xy$online$p), path, tol, niter, torch::jit_scalar(as.integer(xy$block_size)))
  } else cot_polish_dense()$polish_dense(a, b, xy$dense, xx$dense, path, tol, niter)
  cot_polish_value(a, b, pot)
}

# Warm evaluation for line-search probes, started from the current potentials
# and stopped on marginal residuals.
cot_polish_warm <- function(a, b, xy, xx, lambda, pot, ctl = cot_polish_ctl) {
  a <- cot_polish_tensor(a); b <- cot_polish_tensor(b)
  a <- a / a$sum(); b <- b / b$sum()
  tol <- torch::jit_scalar(ctl$marginal_tol)
  residual <- torch::jit_scalar(cot_polish_sinkhorn_tol(lambda, ctl))
  accelerate <- torch::jit_scalar(lambda >= xx$diameter / 4)
  eps <- torch::jit_scalar(lambda); maxit <- torch::jit_scalar(as.integer(ctl$sinkhorn_niter))
  out <- if (cot_polish_online(xy, xx)) {
    cot_polish_dense()$polish_warm_online(a, b, xy$online$x, xy$online$y,
      torch::jit_scalar(xy$online$p), torch::jit_scalar(as.integer(xy$block_size)), eps,
      pot$f_xy, pot$g_yx, pot$f_xx, tol, residual, maxit, accelerate)
  } else cot_polish_dense()$polish_warm_dense(a, b, xy$dense, xx$dense, eps,
      pot$f_xy, pot$g_yx, pot$f_xx, tol, residual, maxit, accelerate)
  cot_polish_value(a, b, out)
}

# Exact candidate block. S includes every strictly positive mass so the block
# triangular identity remains exact even below the support selection tolerance.
cot_polish_hessian <- function(a, b, xy, xx, lambda, pot, candidates = seq_along(a), support = which(a > 0)) {
  block_size <- xy$block_size
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
  H <- (X + Y) * lambda; H <- (H + H$t()) / 2
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

# Active-set solve of the simplex QP min d'Hd/2 + g'd, sum(d) = 0, a + d >= 0.
# The binding set is retained across Levenberg shifts and outer steps.
cot_polish_exact <- function(H, g, a, zero_tol, state = new.env(parent = emptyenv())) {
  n <- length(a)
  binding <- if (length(state$active) == n) state$active else a <= zero_tol & g > sum(a * g) / sum(a)
  if (all(binding)) binding[which.min(g)] <- FALSE
  kkt_tol <- 128 * .Machine$double.eps * max(1, max(abs(g)))
  reduced_step <- function(F, rhs) {
    Z <- cot_polish_tensor(cot_polish_basis(length(F)))
    L <- torch::linalg_cholesky(Z$t()$matmul(H[F, F, drop = FALSE])$matmul(Z))
    as.numeric(Z$matmul(torch::torch_cholesky_solve(-Z$t()$matmul(cot_polish_tensor(rhs)$view(c(-1L, 1L))), L)))
  }
  seen <- character()
  for (attempt in seq_len(20L)) {
    key <- paste(which(binding), collapse = ",")
    if (key %in% seen) break
    seen <- c(seen, key)
    state$active <- binding
    F <- which(!binding)
    if (!length(F)) {binding[which.min(g)] <- FALSE; next}
    trial <- -a; trial[F] <- sum(a[binding]) / length(F)
    if (length(F) > 1L)
      trial[F] <- trial[F] + reduced_step(F, g[F] + as.numeric(H[F, , drop = FALSE]$matmul(cot_polish_tensor(trial))))
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
    if (length(F) > 1L) p[F] <- reduced_step(F, q[F])
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
  binding <- which(a <= epsilon & g > nu + 8 * ctl$sinkhorn_tol)
  zeros <- which(a <= ctl$zero_tol & g < nu)
  list(free = sort(unique(c(setdiff(which(a > 0), binding), zeros))), binding = binding)
}

# Final fit convergence uses the user's absolute-plus-relative tolerance.
cot_converged <- function(gap, D, atol, tol) {
  values <- list(gap, D, atol, tol)
  if (any(lengths(values) != 1L) || any(!is.finite(unlist(values)))) return(FALSE)
  gap <= atol + tol * max(0, D)
}

cot_fit_valid <- function(weights, loss) {
  if (inherits(weights, "torch_tensor")) weights <- as.numeric(weights$detach()$to(device = "cpu"))
  length(weights) > 0L && all(is.finite(weights)) &&
    all(weights >= -64 * .Machine$double.eps) && sum(weights) > 0 &&
    length(loss) == 1L && is.finite(as.numeric(loss))
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

# Newton direction on the free set, with a Levenberg shift when the reduced
# Hessian is not usable. Falls back to the base step if no shift works.
cot_polish_newton <- function(a, b, xy, xx, lambda, ev, base, F, state, ctl) {
  H <- tryCatch(cot_polish_hessian(a, b, xy, xx, lambda, ev$pot, F,
    support = intersect(which(a > 0), F)), error = function(e) NULL)
  if (is.null(H)) return(NULL)
  # The Hessian gauge is irrelevant on sum(d) = 0. Remove it before ridge
  # scaling and active-set multiplier tests, especially at large lambda.
  rows <- H$mean(dim = 2, keepdim = TRUE)
  H <- H - rows - rows$t() + H$mean()
  model_g <- ev$grad[F] + as.numeric(H$matmul(cot_polish_tensor(base[F])))
  centered_g <- model_g - mean(model_g)
  hscale <- max(abs(as.numeric(H$diag())), .Machine$double.eps)
  eye <- torch::torch_eye(length(F), dtype = torch::torch_float64())
  shift <- 0
  repeat {
    dc <- tryCatch(cot_polish_exact(H + eye * shift, centered_g, a[F] + base[F], ctl$zero_tol, state),
      error = function(e) NULL)
    if (!is.null(dc) && all(is.finite(dc))) {
      d <- base; d[F] <- d[F] + dc
      if (sum(ev$grad * d) < 0) return(d)
    }
    shift <- if (shift == 0) 1e-10 * hscale else shift * 10
    if (shift > hscale * 1.01) return(base)
  }
}

# Monotone Newton polish. Uses exact Newton steps on the free set, or L-BFGS
# when the free set exceeds max_free. Returns float64 weights; if the polish
# cannot improve an unconverged start, the start is returned unchanged.
cot_newton_polish <- function(a, b, C_xy, C_xx, C_yy, lambda, scale = NULL, cache = NULL, ctl = cot_polish_ctl) {
  a <- as.numeric(cot_polish_tensor(a)); b <- as.numeric(cot_polish_tensor(b))
  if (any(!is.finite(a)) || any(a < 0) || sum(a) <= 0 || any(!is.finite(b)) || any(b < 0) ||
      abs(sum(b) - 1) > 1e-6 || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0)
    stop("Polish requires probability weights and a finite positive lambda")
  a <- a / sum(a); b <- b / sum(b)
  xy <- cot_polish_cost(C_xy, ctl$block_size); xx <- cot_polish_cost(C_xx, ctl$block_size)
  yy <- cot_polish_cost(C_yy, ctl$block_size)
  if (xy$nrow != length(a) || xy$ncol != length(b) || xx$nrow != length(a) || xx$ncol != length(a) ||
      yy$nrow != length(b) || yy$ncol != length(b)) stop("Incompatible polish cost dimensions")
  eps <- .Machine$double.eps
  if (is.null(scale)) scale <- max(xy$diameter, xx$diameter, lambda * eps, eps)
  scale <- max(scale, yy$diameter)
  key <- as.character(lambda)
  y_self <- if (!is.null(cache) && !is.null(cache[[key]])) cache[[key]] else {
    value <- cot_polish_self(b, yy, lambda, ctl)
    if (!is.null(cache)) cache[[key]] <- value
    value
  }
  # Loss resolution, and the gap resolution: a potential error of size delta
  # moves the computed gap by up to 2 * (delta_xy + delta_xx), and each delta is
  # about the Sinkhorn tolerance over (1 - contraction rate). 16 covers rates up
  # to 3/4.
  floor <- 64 * eps * scale
  atol <- floor + 16 * cot_polish_sinkhorn_tol(lambda, ctl)
  done <- function(gap, loss) gap <= atol + ctl$tol * max(0, loss - y_self)
  fw_gap <- function(w, g) max(0, sum(w * g) - min(g))
  evaluate <- function(w) cot_polish_eval(w, b, xy, xx, lambda, ctl)

  initial <- a
  start <- ev <- evaluate(a); cold <- TRUE
  gap_start <- fw_gap(a, ev$grad)
  status <- "maxit"; accepted <- iter <- 0L; use_warm <- ctl$warm_start
  memory <- list(); prevF <- integer(); prev_a <- prev_g <- NULL; last_dir <- ""
  qp_state <- new.env(parent = emptyenv())
  for (k in seq_len(ctl$maxit)) {
    gap <- fw_gap(a, ev$grad)
    if (done(gap, ev$loss)) {
      # Only a cold evaluation may end the polish.
      if (!cold) {ev <- evaluate(a); cold <- TRUE; use_warm <- FALSE; gap <- fw_gap(a, ev$grad)}
      if (done(gap, ev$loss)) break
    }
    iter <- k
    sets <- cot_polish_sets(a, ev$grad, scale, ctl)
    F <- sets$free
    if (!identical(qp_state$indices, F)) {qp_state$active <- NULL; qp_state$indices <- F}
    dir <- if (length(F) > ctl$max_free) "lbfgs" else "newton"
    if (dir != last_dir) {memory <- list(); prevF <- integer(); prev_a <- prev_g <- NULL}
    last_dir <- dir
    base <- numeric(length(a))
    base[sets$binding] <- -a[sets$binding]
    base[F] <- -sum(base) / length(F)
    if (dir == "lbfgs") {
      if (!identical(F, prevF)) memory <- list() else if (!is.null(prev_a)) {
        s <- a[F] - prev_a[F]; y <- ev$grad[F] - prev_g[F]
        s <- s - mean(s); y <- y - mean(y); sy <- sum(s * y)
        if (sy > 1e-12 * sqrt(sum(s^2) * sum(y^2)))
          memory <- utils::tail(c(memory, list(list(s = s, y = y, sy = sy))), ctl$lbfgs_m)
      }
      d <- base
      d[F] <- d[F] + cot_polish_lbfgs(ev$grad[F], memory, mean(a))
      # The first step is scaled to probability-mass units.
      if (!length(memory) && max(abs(d)) > mean(a)) d <- d * mean(a) / max(abs(d))
      prevF <- F; prev_a <- a; prev_g <- ev$grad
    } else d <- cot_polish_newton(a, b, xy, xx, lambda, ev, base, F, qp_state, ctl)
    if (is.null(d) || !all(is.finite(d)) || sum(ev$grad * d) >= 0) {status <- "no_descent"; break}

    alpha <- 1; moved <- FALSE
    while (alpha >= ctl$alpha_min) {
      next_ev <- NULL; slope <- NA_real_
      cand <- pmax(a + alpha * d, 0)
      cand[cand < ctl$zero_tol] <- 0
      if (sum(cand) > 0) {
        cand <- cand / sum(cand)
        # Centered: cand - a sums to zero only up to rounding, which an
        # uncentered gradient would turn into a slope error of about 1e-16.
        slope <- sum((ev$grad - sum(a * ev$grad)) * (cand - a))
        if (slope < 0) {
          probe_cold <- !use_warm
          if (use_warm) next_ev <- tryCatch(cot_polish_warm(cand, b, xy, xx, lambda, ev$pot, ctl), error = function(e) NULL)
          if (is.null(next_ev)) {
            next_ev <- evaluate(cand); probe_cold <- TRUE
          } else if (abs(ev$loss - next_ev$loss) <= 1e-11 * scale) {
            # Close to resolution, marginal residual noise must not choose the
            # support: switch to cold evaluations for the rest of the polish.
            next_ev <- evaluate(cand); probe_cold <- TRUE; use_warm <- FALSE
          }
          if (is.finite(next_ev$loss) && next_ev$loss <= ev$loss &&
              next_ev$loss <= ev$loss + ctl$armijo * slope &&
              (next_ev$loss < ev$loss || fw_gap(cand, next_ev$grad) < gap)) {moved <- TRUE; break}
          # Below loss resolution the loss cannot rank a Newton step; a cold
          # probe that at least halves the gap is accepted on the gap alone.
          if (probe_cold && is.finite(next_ev$loss) && next_ev$loss <= ev$loss + floor &&
              fw_gap(cand, next_ev$grad) < 0.5 * gap) {moved <- TRUE; break}
        }
      }
      shrink <- 0.5
      if (dir == "lbfgs" && !length(memory) && !is.null(next_ev) && slope < 0) {
        # A rejected probe also supplies a secant curvature estimate, which
        # avoids many halvings on a sharply curved face.
        curvature <- sum((next_ev$grad - ev$grad) * (cand - a))
        if (is.finite(curvature) && curvature > 0) shrink <- max(0.05, min(0.5, -slope / curvature))
      }
      alpha <- alpha * shrink
    }
    if (!moved) {status <- "line_search_failed"; break}
    improvement <- ev$loss - next_ev$loss
    a <- cand; ev <- next_ev; cold <- probe_cold; accepted <- accepted + 1L
    # Stalled only when neither the loss nor the gap made real progress.
    if (improvement <= floor && !(fw_gap(a, ev$grad) < 0.5 * gap)) {status <- "stalled"; break}
  }
  if (!cold) ev <- evaluate(a)
  gap <- fw_gap(a, ev$grad)
  if (done(gap, ev$loss)) {
    status <- if (accepted > 0L) "refined" else "already_optimal"
  } else if (start$loss - ev$loss <= 10 * eps * max(scale, abs(start$loss))) {
    a <- initial; ev <- start; gap <- gap_start; status <- "no_improvement"
  }
  list(weights = a, loss = ev$loss, divergence = ev$loss - y_self, gap = gap,
    changed = !identical(a, initial), status = status, iterations = iter, atol = atol)
}
