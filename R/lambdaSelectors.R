cot_selector_loo_convention <- "self_and_cross_overlap_renormalized_empty_minus_inf"
cot_selector_marginal_loo_convention <- "cross_overlap_renormalized_empty_minus_inf"

cot_lambda_selector_choices <- function(include_unsupported = TRUE) {
  out <- c(
    "profile_symmetric_loo",
    "profile_marginal_loo",
    "symmetric_loo",
    "marginal_loo",
    "mass_bootstrap_symmetric",
    "bootstrap_legacy",
    "sym_evidence",
    "legacy_bootstrap_balance"
  )
  if (isTRUE(include_unsupported)) {
    out <- c(out, "kfold_symmetric")
  }
  out
}

cot_normalize_lambda_selector <- function(lambda.selector) {
  stopifnot(
    is.character(lambda.selector),
    length(lambda.selector) == 1L,
    !is.na(lambda.selector)
  )

  allowed <- cot_lambda_selector_choices(include_unsupported = TRUE)
  if (!(lambda.selector %in% allowed)) {
    stop(
      sprintf(
        "lambda.selector must be one of %s.",
        paste(sprintf("'%s'", allowed), collapse = ", ")
      )
    )
  }

  if (identical(lambda.selector, "kfold_symmetric")) {
    stop(
      "lambda.selector = 'kfold_symmetric' is not implemented yet. It requires explicit fold identities, held-out scoring rules, and refit plumbing."
    )
  }

  switch(
    lambda.selector,
    "sym_evidence" = "profile_symmetric_loo",
    "legacy_bootstrap_balance" = "bootstrap_legacy",
    lambda.selector
  )
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

  allowed <- c(
    "pilot.weights",
    "details",
    "mass.bootstrap.summary",
    "mass.bootstrap.quantile",
    "mass.bootstrap.kappa"
  )
  unknown <- setdiff(names(control), allowed)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown lambda.selector.control field(s): %s",
        paste(unknown, collapse = ", ")
      )
    )
  }

  out <- list(
    pilot.weights = NULL,
    details = FALSE,
    mass.bootstrap.summary = "mean",
    mass.bootstrap.quantile = 0.1,
    mass.bootstrap.kappa = 1
  )

  if ("pilot.weights" %in% names(control)) {
    if (!is.null(control[["pilot.weights"]])) {
      pilot.weights <- as.numeric(control[["pilot.weights"]])
      if (
        length(pilot.weights) == 0L ||
          any(!is.finite(pilot.weights)) ||
          any(pilot.weights < 0)
      ) {
        stop(
          "lambda.selector.control$pilot.weights must be a nonnegative numeric vector with positive length."
        )
      }
      out$pilot.weights <- pilot.weights
    }
  }

  if ("details" %in% names(control)) {
    details <- control[["details"]]
    if (!is.logical(details) || length(details) != 1L || is.na(details)) {
      stop("lambda.selector.control$details must be TRUE or FALSE.")
    }
    out$details <- details
  }

  if ("mass.bootstrap.summary" %in% names(control)) {
    out$mass.bootstrap.summary <- match.arg(
      control[["mass.bootstrap.summary"]],
      c("mean", "quantile", "mean_minus_sd")
    )
  }

  if ("mass.bootstrap.quantile" %in% names(control)) {
    mass.bootstrap.quantile <- as.numeric(control[["mass.bootstrap.quantile"]])
    if (
      length(mass.bootstrap.quantile) != 1L ||
        !is.finite(mass.bootstrap.quantile) ||
        mass.bootstrap.quantile < 0 ||
        mass.bootstrap.quantile > 1
    ) {
      stop(
        "lambda.selector.control$mass.bootstrap.quantile must be a single number in [0, 1]."
      )
    }
    out$mass.bootstrap.quantile <- mass.bootstrap.quantile
  }

  if ("mass.bootstrap.kappa" %in% names(control)) {
    mass.bootstrap.kappa <- as.numeric(control[["mass.bootstrap.kappa"]])
    if (
      length(mass.bootstrap.kappa) != 1L ||
        !is.finite(mass.bootstrap.kappa) ||
        mass.bootstrap.kappa < 0
    ) {
      stop(
        "lambda.selector.control$mass.bootstrap.kappa must be a single nonnegative number."
      )
    }
    out$mass.bootstrap.kappa <- mass.bootstrap.kappa
  }

  out
}

cot_resolve_pilot_weights <- function(pilot.weights, n) {
  stopifnot(length(n) == 1L, is.numeric(n), is.finite(n), n >= 1)
  n <- as.integer(n)

  if (is.null(pilot.weights)) {
    return(rep(1.0 / n, n))
  }

  pilot.weights <- as.numeric(pilot.weights)
  if (length(pilot.weights) != n) {
    stop(
      sprintf(
        "lambda.selector.control$pilot.weights must have length %s.",
        n
      )
    )
  }
  if (any(!is.finite(pilot.weights)) || any(pilot.weights < 0)) {
    stop(
      "lambda.selector.control$pilot.weights must be a nonnegative numeric vector."
    )
  }

  total_mass <- sum(pilot.weights)
  if (total_mass <= sqrt(.Machine$double.eps)) {
    stop(
      "lambda.selector.control$pilot.weights must have positive total mass."
    )
  }

  pilot.weights / total_mass
}

cot_safe_log_weights <- function(weights) {
  log_weights <- log(as.numeric(weights))
  log_weights[!is.finite(log_weights)] <- -Inf
  log_weights
}

cot_to_numeric_vector <- function(x) {
  if (inherits(x, "torch_tensor")) {
    return(as.numeric(x$to(device = "cpu")))
  }
  as.numeric(x)
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

cot_log_mixture <- function(cost, log_weights, lambda) {
  cost <- cot_to_numeric_matrix(cost)
  log_weights <- cot_to_numeric_vector(log_weights)
  stopifnot(is.matrix(cost))
  stopifnot(length(log_weights) == ncol(cost))
  stopifnot(
    is.numeric(lambda),
    length(lambda) == 1L,
    is.finite(lambda),
    lambda > 0
  )

  matrixStats::rowLogSumExps(
    sweep(-cost / lambda, 2, log_weights, `+`)
  )
}

cot_loo_log_mixture <- function(cost, log_weights, weights, lambda) {
  cost <- cot_to_numeric_matrix(cost)
  log_weights <- cot_to_numeric_vector(log_weights)
  weights <- cot_to_numeric_vector(weights)
  stopifnot(is.matrix(cost))
  stopifnot(nrow(cost) == ncol(cost))
  stopifnot(length(log_weights) == ncol(cost))
  stopifnot(length(weights) == nrow(cost))
  stopifnot(
    is.numeric(lambda),
    length(lambda) == 1L,
    is.finite(lambda),
    lambda > 0
  )

  loo_mix <- sweep(-cost / lambda, 2, log_weights, `+`)
  diag(loo_mix) <- -Inf
  loo_mix <- matrixStats::rowLogSumExps(loo_mix)

  remaining_mass <- pmax(1.0 - weights, 0.0)
  valid <- remaining_mass > sqrt(.Machine$double.eps)
  loo_mix[valid] <- loo_mix[valid] - log(remaining_mass[valid])
  loo_mix[!valid] <- -Inf

  loo_mix
}

cot_cross_loo_log_mixture <- function(cost, log_weights, weights, lambda) {
  cost <- cot_to_numeric_matrix(cost)
  log_weights <- cot_to_numeric_vector(log_weights)
  weights <- cot_to_numeric_vector(weights)
  stopifnot(is.matrix(cost))
  stopifnot(nrow(cost) == ncol(cost))
  stopifnot(length(log_weights) == nrow(cost))
  stopifnot(length(weights) == nrow(cost))
  stopifnot(
    is.numeric(lambda),
    length(lambda) == 1L,
    is.finite(lambda),
    lambda > 0
  )

  loo_mix <- sweep(-cost / lambda, 1, log_weights, `+`)
  diag(loo_mix) <- -Inf
  loo_mix <- matrixStats::colLogSumExps(loo_mix)

  remaining_mass <- pmax(1.0 - weights, 0.0)
  valid <- remaining_mass > sqrt(.Machine$double.eps)
  loo_mix[valid] <- loo_mix[valid] - log(remaining_mass[valid])
  loo_mix[!valid] <- -Inf

  loo_mix
}

cot_cross_masked_log_mixture <- function(
  cost,
  log_weights,
  weights,
  lambda,
  exclude
) {
  cost <- cot_to_numeric_matrix(cost)
  log_weights <- cot_to_numeric_vector(log_weights)
  weights <- cot_to_numeric_vector(weights)
  exclude <- as.matrix(exclude)
  stopifnot(is.matrix(cost))
  stopifnot(all(dim(exclude) == dim(cost)))
  stopifnot(length(log_weights) == nrow(cost))
  stopifnot(length(weights) == nrow(cost))
  stopifnot(
    is.numeric(lambda),
    length(lambda) == 1L,
    is.finite(lambda),
    lambda > 0
  )

  masked_mix <- sweep(-cost / lambda, 1, log_weights, `+`)
  masked_mix[exclude] <- -Inf
  masked_mix <- matrixStats::colLogSumExps(masked_mix)

  excluded_mass <- as.numeric(crossprod(weights, exclude))
  remaining_mass <- pmax(1.0 - excluded_mass, 0.0)
  valid <- remaining_mass > sqrt(.Machine$double.eps)
  masked_mix[valid] <- masked_mix[valid] - log(remaining_mass[valid])
  masked_mix[!valid] <- -Inf

  masked_mix
}

cot_overlap_exclusion_mask <- function(
  source_x,
  target_x,
  tol = sqrt(.Machine$double.eps)
) {
  # Conservative fallback when explicit unit IDs are unavailable. This avoids
  # excluding every duplicate-covariate pairing by matching exact covariate
  # duplicates one-to-one instead.
  source_x <- cot_to_numeric_matrix(source_x)
  target_x <- cot_to_numeric_matrix(target_x)
  stopifnot(ncol(source_x) == ncol(target_x))

  out <- matrix(FALSE, nrow(source_x), nrow(target_x))
  if (nrow(source_x) == 0L || nrow(target_x) == 0L) {
    return(out)
  }

  target_available <- rep(TRUE, nrow(target_x))
  for (i in seq_len(nrow(source_x))) {
    diffs <- abs(sweep(target_x, 2, source_x[i, ], `-`))
    matches <- which(target_available & (rowSums(diffs > tol) == 0L))
    if (length(matches) > 0L) {
      match_idx <- matches[[1L]]
      out[i, match_idx] <- TRUE
      target_available[[match_idx]] <- FALSE
    }
  }

  out
}

cot_unit_exclusion_mask <- function(source_unit_id, target_unit_id) {
  source_unit_id <- as.character(source_unit_id)
  target_unit_id <- as.character(target_unit_id)

  out <- matrix(FALSE, nrow = length(source_unit_id), ncol = length(target_unit_id))
  if (length(source_unit_id) == 0L || length(target_unit_id) == 0L) {
    return(out)
  }

  valid_target <- !is.na(target_unit_id)
  for (i in seq_along(source_unit_id)) {
    if (is.na(source_unit_id[[i]])) {
      next
    }
    out[i, ] <- valid_target & (target_unit_id == source_unit_id[[i]])
  }

  out
}

cot_cross_log_mixture <- function(cross_cost, w, lambda, cross_exclude = NULL) {
  if (!is.null(cross_exclude)) {
    return(cot_cross_masked_log_mixture(
      cross_cost,
      log_weights = log_weights(w),
      weights = w,
      lambda = lambda,
      exclude = cross_exclude
    ))
  }

  cot_log_mixture(
    t(cross_cost),
    log_weights = log_weights(w),
    lambda = lambda
  )
}

cot_marginal_evidence <- function(
  cross_cost,
  cross_exclude = NULL,
  w,
  a,
  lambda
) {
  cross_logmix <- cot_cross_log_mixture(
    cross_cost = cross_cost,
    w = w,
    lambda = lambda,
    cross_exclude = cross_exclude
  )
  cross_idx <- a > 0
  sum(a[cross_idx] * cross_logmix[cross_idx])
}

cot_sym_evidence <- function(
  cross_cost,
  self_a_cost,
  self_w_cost,
  cross_exclude = NULL,
  w,
  a,
  lambda
) {
  cross_logmix <- cot_cross_log_mixture(
    cross_cost = cross_cost,
    w = w,
    lambda = lambda,
    cross_exclude = cross_exclude
  )
  self_a_logmix <- cot_loo_log_mixture(
    self_a_cost,
    log_weights = log_weights(a),
    weights = a,
    lambda = lambda
  )
  self_w_logmix <- cot_loo_log_mixture(
    self_w_cost,
    log_weights = log_weights(w),
    weights = w,
    lambda = lambda
  )

  cross_idx <- a > 0
  self_a_idx <- a > 0
  self_w_idx <- w > 0

  sum(a[cross_idx] * cross_logmix[cross_idx]) -
    0.5 * sum(a[self_a_idx] * self_a_logmix[self_a_idx]) -
    0.5 * sum(w[self_w_idx] * self_w_logmix[self_w_idx])
}

cot_draw_bootstrap_mass <- function(prob, size = length(prob)) {
  prob <- cot_to_numeric_vector(prob)
  stopifnot(length(size) == 1L, is.numeric(size), is.finite(size), size >= 1)
  size <- as.integer(size)

  total_mass <- sum(prob)
  if (any(prob < 0) || total_mass <= sqrt(.Machine$double.eps)) {
    stop("Bootstrap masses must be sampled from nonnegative probabilities.")
  }

  counts <- stats::rmultinom(1L, size = size, prob = prob / total_mass)[, 1L]
  counts / sum(counts)
}

cot_summarize_bootstrap_curve <- function(
  score_matrix,
  summary = c("mean", "quantile", "mean_minus_sd"),
  quantile = 0.1,
  kappa = 1
) {
  summary <- match.arg(summary)
  score_matrix <- as.matrix(score_matrix)

  summarize_one <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(NA_real_)
    }

    switch(
      summary,
      "mean" = mean(x),
      "quantile" = as.numeric(stats::quantile(
        x,
        probs = quantile,
        names = FALSE,
        type = 7
      )),
      "mean_minus_sd" = {
        penalty_sd <- if (length(x) > 1L) stats::sd(x) else 0.0
        mean(x) - kappa * penalty_sd
      }
    )
  }

  apply(score_matrix, 2, summarize_one)
}
