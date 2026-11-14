# "auto" uses treatment-prediction RBF; "bootstrap" uses the legacy selector
cot_lambda_selector_choices <- function() {
  c("auto", "bootstrap")
}

cot_normalize_lambda_selector <- function(lambda.selector) {
  stopifnot(
    is.character(lambda.selector),
    length(lambda.selector) == 1L,
    !is.na(lambda.selector)
  )

  allowed <- cot_lambda_selector_choices()
  if (!(lambda.selector %in% allowed)) {
    stop(
      sprintf(
        "lambda.selector must be one of %s.",
        paste(sprintf("'%s'", allowed), collapse = ", ")
      )
    )
  }
  lambda.selector
}

cot_normalize_lambda_selector_control <- function(control = list()) {
  if (is.null(control)) {
    control <- list()
  }
  if (!is.list(control)) {
    stop("lambda.selector.control must be a list.")
  }
  if (length(control) > 0L && is.null(names(control))) {
    stop("lambda.selector.control entries must be named.")
  }
  if (anyNA(names(control)) || any(names(control) == "")) {
    stop("lambda.selector.control entries must be named.")
  }

  unknown <- setdiff(names(control), "details")
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown lambda.selector.control field(s): %s",
        paste(unknown, collapse = ", ")
      )
    )
  }

  out <- list(details = FALSE)
  if ("details" %in% names(control)) {
    details <- control[["details"]]
    if (!is.logical(details) || length(details) != 1L || is.na(details)) {
      stop("lambda.selector.control$details must be TRUE or FALSE.")
    }
    out$details <- details
  }

  out
}

cot_to_numeric_matrix <- function(x) {
  if (inherits(x, "costTensor")) {
    x <- x$data
  }
  if (inherits(x, "torch_tensor")) {
    return(as.matrix(x$to(device = "cpu")))
  }
  as.matrix(x)
}

# treatment prediction with solver cost and unpenalized intercept
# deletion identity uses the full-sample ridge penalty alpha * n
cot_trt_rbf_score <- function(C, z, t, alphas = 10^seq(-3, 1, 0.5)) {
  n <- length(z)
  K <- exp(-C / t)
  Kc <- sweep(sweep(K, 1L, rowMeans(K)), 2L, colMeans(K)) + mean(K)
  eig <- eigen((Kc + t(Kc)) / 2, symmetric = TRUE)
  tol <- 1e-8 * max(1, max(abs(eig$values)))
  if (min(eig$values) < -tol) {
    stop("Treatment kernel is materially indefinite.")
  }
  ev <- pmax(eig$values, 0)
  zc <- z - mean(z)
  proj <- as.vector(crossprod(eig$vectors, zc))
  scores <- vapply(
    alphas,
    function(alpha) {
      s <- ev / (ev + alpha * n)
      fitted <- as.vector(eig$vectors %*% (s * proj))
      h <- as.vector((eig$vectors^2) %*% s) + 1 / n
      if (any(!is.finite(h)) || any(1 - h <= 1e-10)) {
        return(Inf)
      }
      mean(((zc - fitted) / (1 - h))^2)
    },
    numeric(1)
  )
  scores[!is.finite(scores)] <- Inf
  i <- which.min(scores)
  if (!length(i)) {
    return(list(score = Inf, alpha = NA_real_))
  }
  list(score = scores[i], alpha = alphas[i])
}

# bounded search over scaled costs; return the best evaluated point
cot_trt_rbf_search <- function(
  score,
  max_eval = 50L,
  max_expand = 2L,
  tolerance = 1e-3,
  tie_tolerance = 1e-10
) {
  if (
    length(max_eval) != 1L ||
      !is.finite(max_eval) ||
      max_eval < 1 ||
      max_eval != floor(max_eval)
  ) {
    stop("max_eval must be a positive integer.")
  }
  if (
    length(max_expand) != 1L ||
      !is.finite(max_expand) ||
      max_expand < 0 ||
      max_expand != floor(max_expand)
  ) {
    stop("max_expand must be a nonnegative integer.")
  }
  if (
    length(tolerance) != 1L ||
      !is.finite(tolerance) ||
      tolerance <= 0 ||
      length(tie_tolerance) != 1L ||
      !is.finite(tie_tolerance) ||
      tie_tolerance < 0
  ) {
    stop(
      "Search tolerances must be finite, with tolerance > 0 and tie_tolerance >= 0."
    )
  }
  trace <- data.frame(u = numeric(), score = numeric(), alpha = numeric())
  evaluate <- function(u) {
    old <- which(abs(trace$u - u) < 1e-12)
    if (length(old)) {
      return(
        if (is.finite(trace$score[old[1]])) {
          trace$score[old[1]]
        } else {
          .Machine$double.xmax / 100
        }
      )
    }
    if (nrow(trace) >= max_eval) {
      stop(structure(
        list(
          message = "Lambda search evaluation budget exhausted.",
          call = NULL
        ),
        class = c("cot_budget", "error", "condition")
      ))
    }
    value <- score(exp(u))
    if (length(value$score) != 1L || !is.finite(value$score)) {
      value$score <- Inf
    }
    trace[nrow(trace) + 1L, ] <<- c(u, value$score, value$alpha)
    if (is.finite(value$score)) value$score else .Machine$double.xmax / 100
  }
  bounds <- log(2) + log(10) * c(-1.5, 2.5)
  initial <- sort(unique(c(seq(bounds[1], bounds[2], length.out = 7), log(2))))
  # always score the starting bandwidth
  initial <- unique(c(log(2), bounds, initial))
  for (u in utils::head(initial, max_eval)) {
    evaluate(u)
  }
  expansions <- 0L
  # expand when the endpoint beats its neighbor
  repeat {
    o <- order(trace$u)
    a <- trace[o, ]
    if (nrow(a) < 2L) {
      break
    }
    improving <- c(
      a$score[1] < a$score[2],
      utils::tail(a$score, 1) < utils::tail(a$score, 2)[1]
    )
    edges <- which(improving)
    direction <- if (length(edges)) {
      edge_scores <- c(a$score[1], tail(a$score, 1))
      c(-1L, 1L)[edges[which.min(edge_scores[edges])]]
    } else {
      0L
    }
    if (!direction || expansions >= max_expand || nrow(trace) >= max_eval) {
      break
    }
    edge <- if (direction < 0) bounds[1] else bounds[2]
    next_edge <- max(-30, min(30, edge + direction * log(10)))
    if (next_edge == edge) {
      break
    }
    bounds[if (direction < 0) 1 else 2] <- next_edge
    evaluate(next_edge)
    expansions <- expansions + 1L
  }
  # refine sampled local minima by score, sharing the budget across brackets
  a <- trace[order(trace$u), ]
  n <- nrow(a)
  candidates <- which(
    c(TRUE, a$score[-1] <= a$score[-n]) & c(a$score[-n] <= a$score[-1], TRUE)
  )
  candidates <- candidates[order(a$score[candidates])]
  budget_hit <- FALSE
  for (i in candidates) {
    lo <- a$u[max(1, i - 1L)]
    hi <- a$u[min(n, i + 1L)]
    if (hi <= lo) {
      next
    }
    tryCatch(
      stats::optimize(evaluate, c(lo, hi), tol = log1p(tolerance)),
      cot_budget = function(e) {
        budget_hit <<- TRUE
        NULL
      }
    )
    if (budget_hit) break
  }
  valid <- which(is.finite(trace$score))
  if (!length(valid)) {
    stop("No valid treatment-kernel lambda candidate was found.")
  }
  best_score <- min(trace$score[valid])
  tied <- valid[
    trace$score[valid] <= best_score + tie_tolerance * max(1, abs(best_score))
  ]
  best <- tied[which.max(trace$u[tied])]
  ordered <- trace[order(trace$u), , drop = FALSE]
  nn <- nrow(ordered)
  minima <- which(
    c(TRUE, ordered$score[-1] <= ordered$score[-nn]) &
      c(ordered$score[-nn] <= ordered$score[-1], TRUE) &
      is.finite(ordered$score)
  )
  competing <- minima[abs(ordered$u[minima] - trace$u[best]) > log1p(tolerance)]
  gap <- if (length(competing)) {
    min(ordered$score[competing]) - trace$score[best]
  } else {
    NA_real_
  }
  list(
    multiplier = exp(trace$u[best]),
    alpha = trace$alpha[best],
    score = trace$score[best],
    status = if (budget_hit || nrow(trace) >= max_eval) {
      "budget_exhausted"
    } else {
      "search_complete"
    },
    bounds = exp(bounds),
    boundary = c(
      lower = abs(trace$u[best] - bounds[1]) < log1p(tolerance),
      upper = abs(trace$u[best] - bounds[2]) < log1p(tolerance)
    ),
    expansions = expansions,
    score_gap = gap,
    n_evaluations = nrow(trace),
    relative_tolerance = tolerance,
    tie_tolerance = tie_tolerance,
    evaluations = transform(trace, multiplier = exp(u))
  )
}

# numeric tuning data; IDs identify observations, not equal covariate rows
cot_trt_rbf_context <- function(x, z, unit.id = NULL) {
  x <- as.matrix(x)
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop("lambda selection by 'auto' requires finite numeric covariates.")
  }
  if (
    (!is.numeric(z) && !is.logical(z)) ||
      length(z) != nrow(x) ||
      anyNA(z) ||
      any(!z %in% c(0, 1)) ||
      any(tabulate(as.integer(z) + 1L, 2L) < 2L)
  ) {
    stop(
      "lambda selection by 'auto' requires at least two observations in each binary class."
    )
  }
  if (
    !is.null(unit.id) &&
      (length(unit.id) != nrow(x) || anyNA(unit.id) || anyDuplicated(unit.id))
  ) {
    stop("lambda.selection.data unit.id must be unique and nonmissing.")
  }
  list(
    x = x + 0,
    z = as.integer(z),
    unit.id = if (is.null(unit.id)) NULL else as.character(unit.id)
  )
}

cot_trt_rbf_prepare <- function(context, max_sample = 2000L, seed = 2026L) {
  x <- context$x
  z <- context$z
  n <- length(z)
  if (
    length(max_sample) != 1L ||
      !is.finite(max_sample) ||
      max_sample < 4 ||
      max_sample != floor(max_sample)
  ) {
    stop("max_sample must be an integer of at least four.")
  }
  indices <- seq_len(n)
  if (n > max_sample) {
    existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (existed) {
      saved <- get(".Random.seed", envir = .GlobalEnv)
    }
    on.exit(
      if (existed) {
        assign(".Random.seed", saved, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      },
      add = TRUE
    )
    set.seed(seed)
    counts <- tabulate(z + 1L, 2L)
    n0 <- max(
      2L,
      min(counts[1], max_sample - 2L, round(max_sample * counts[1] / n))
    )
    n0 <- max(n0, max_sample - counts[2])
    indices <- sort(c(
      sample(which(z == 0L), n0),
      sample(which(z == 1L), max_sample - n0)
    ))
  }
  xs <- x[indices, , drop = FALSE]
  # use CPU float64 costs for spectral regression to avoid RBF roundoff errors
  # the transport fit uses the requested device and dtype
  xs <- torch::torch_tensor(
    xs,
    dtype = torch::torch_double(),
    device = torch::torch_device("cpu")
  )$contiguous()
  C <- cot_to_numeric_matrix(cost(xs, xs, p = 2, tensorized = TRUE))
  if (any(!is.finite(C))) {
    stop("lambda selection by 'auto' requires finite solver costs.")
  }
  m <- stats::median(C[upper.tri(C)])
  if (!is.finite(m) || m <= 0) {
    stop(
      "lambda selection by 'auto' requires a positive median off-diagonal cost, including duplicate observations."
    )
  }
  list(
    C = C / m,
    z = z[indices],
    median_cost = m,
    sample_size = length(indices),
    sample_indices = indices,
    sample_seed = seed,
    median_is_subsample = length(indices) < n,
    criterion = "centered kernel ridge, unpenalized intercept, LOO with full-sample alpha*n penalty"
  )
}

cot_trt_rbf_grid <- function(
  prepared,
  lambda,
  eligible = rep(TRUE, length(lambda)),
  details = FALSE,
  cache = NULL
) {
  reasons <- rep(NA_character_, length(lambda))
  reasons[
    !is.finite(lambda) | lambda <= 0
  ] <- "lambda is not finite and positive"
  reasons[!eligible] <- "numerically failed fit"
  eligible <- eligible & is.finite(lambda) & lambda > 0
  scores <- rep(Inf, length(lambda))
  alpha <- rep(NA_real_, length(lambda))
  for (value in unique(lambda[eligible])) {
    i <- which(eligible & lambda == value)
    key <- paste(
      sprintf("%02x", as.integer(writeBin(value, raw(), size = 8L))),
      collapse = ""
    )
    result <- if (
      !is.null(cache) && exists(key, envir = cache, inherits = FALSE)
    ) {
      get(key, envir = cache)
    } else {
      value_score <- cot_trt_rbf_score(
        prepared$C,
        prepared$z,
        value / prepared$median_cost
      )
      if (!is.null(cache)) {
        assign(key, value_score, envir = cache)
      }
      value_score
    }
    scores[i] <- result$score
    alpha[i] <- result$alpha
  }
  reasons[eligible & !is.finite(scores)] <- "non-finite prediction score"
  valid <- which(is.finite(scores))
  if (!length(valid)) {
    stop("No eligible finite lambda scores remain.")
  }
  best <- min(scores[valid])
  tied <- valid[scores[valid] <= best + 1e-10 * max(1, abs(best))]
  selected <- tied[which.max(lambda[tied])]
  list(
    scores = scores,
    alpha = alpha,
    eligible = eligible,
    reasons = reasons,
    selected_index = selected,
    lambda = lambda[selected],
    score = scores[selected],
    selected_alpha = alpha[selected],
    multiplier = lambda / prepared$median_cost,
    boundary = lambda[selected] %in% range(lambda[eligible]),
    mode = "grid",
    method = "auto",
    median_cost = prepared$median_cost,
    sample_size = prepared$sample_size,
    sample_indices = prepared$sample_indices,
    sample_seed = prepared$sample_seed,
    median_is_subsample = prepared$median_is_subsample,
    criterion = prepared$criterion
  )
}

cot_select_trt_rbf <- function(
  data,
  options,
  max_sample = 2000L,
  seed = 2026L,
  max_eval = 50L,
  max_expand = 2L
) {
  if (!isTRUE(options$p == 2) || !is.null(options$cost.function)) {
    stop(
      "lambda = 'auto' requires the standard p = 2 cost, without a custom cost.function."
    )
  }
  n <- length(data@z)
  if (
    any(!is.finite(data@weights)) || max(abs(data@weights - 1 / n)) > 1e-10 / n
  ) {
    stop(
      "lambda selection by 'auto' requires uniform base observation weights."
    )
  }
  prepared <- cot_trt_rbf_prepare(
    cot_trt_rbf_context(data@x, data@z),
    max_sample,
    seed
  )
  m <- prepared$median_cost
  result <- cot_trt_rbf_search(
    function(t) cot_trt_rbf_score(prepared$C, prepared$z, t),
    max_eval,
    max_expand
  )
  result$method <- "auto"
  result$median_cost <- m
  result$lambda <- m * result$multiplier
  result$criterion <- "centered kernel ridge, unpenalized intercept, LOO with full-sample alpha*n penalty"
  result$sample_size <- prepared$sample_size
  result$sample_indices <- prepared$sample_indices
  result$sample_seed <- prepared$sample_seed
  result$median_is_subsample <- prepared$median_is_subsample
  result$mode <- "adaptive"
  if (!isTRUE(options$lambda.selector.control$details)) {
    result$evaluations <- NULL
  }
  result
}
