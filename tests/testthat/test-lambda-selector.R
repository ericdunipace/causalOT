sq_cost <- function(x, y) {
  out <- matrix(0.0, nrow(x), nrow(y))
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(y))) {
      out[i, j] <- sum((x[i, ] - y[j, ])^2) / 2.0
    }
  }
  out
}

force_converged_flags <- function(cot) {
  priv <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
  priv$converged_run <- lapply(priv$converged_run, function(x) {
    rep(TRUE, length(x))
  })
  priv
}

mock_detachable <- function(x) {
  structure(
    list(detach = function() x),
    class = "mock_detachable"
  )
}

make_weight_env <- function(weight, key = "source", metric = NULL) {
  out <- rlang::env()
  class(out) <- c("weightEnv", class(out))
  out[[key]] <- mock_detachable(weight)
  if (!is.null(metric)) {
    out$metric <- metric
  }
  out
}

selector_state_fixture <- function(cross_exclude = NULL) {
  list(
    supported = TRUE,
    source_add = "source",
    target_add = "target",
    source_n = 3L,
    target_n = 2L,
    target_weights = c(0.4, 0.6),
    cross_exclude = cross_exclude,
    cross_cost = matrix(
      c(0.2, 1.1, 0.8, 0.3, 1.4, 0.4),
      nrow = 3,
      byrow = TRUE
    ),
    self_w_cost = matrix(
      c(0.0, 0.6, 1.1, 0.6, 0.0, 0.4, 1.1, 0.4, 0.0),
      nrow = 3,
      byrow = TRUE
    ),
    self_a_cost = matrix(
      c(0.0, 0.7, 0.7, 0.0),
      nrow = 2,
      byrow = TRUE
    )
  )
}

selector_test_problem <- function() {
  causalOT:::torch_check()
  causalOT:::OTProblem(
    causalOT:::Measure(matrix(0, 1, 1), adapt = "weights"),
    causalOT:::Measure(matrix(1, 1, 1), adapt = "none")
  )
}

replace_private_binding <- function(env, name, value) {
  was_locked <- bindingIsLocked(name, env)
  old_value <- env[[name]]
  if (was_locked) {
    unlockBinding(name, env)
  }
  env[[name]] <- value
  if (was_locked) {
    lockBinding(name, env)
  }
  old_value
}

testthat::test_that("Measure ignores automatic row names for default unit ids", {
  df <- data.frame(x = 1:2)
  testthat::expect_null(causalOT:::Measure(df)$unit_id)
})

testthat::test_that("Measure preserves explicit non-automatic row names as unit ids", {
  df <- data.frame(x = 1:2)
  rownames(df) <- c("u1", "u2")
  testthat::expect_equal(causalOT:::Measure(df)$unit_id, c("u1", "u2"))
})

testthat::test_that("log-mixture helpers match hand calculations", {
  lambda <- 2.5
  weights <- c(0.2, 0.5, 0.3)
  log_weights <- log(weights)
  cross_cost <- matrix(
    c(0.0, 1.0, 2.0, 0.5, 0.25, 1.5),
    nrow = 2,
    byrow = TRUE
  )
  self_cost <- matrix(
    c(0.0, 1.0, 2.0, 1.0, 0.0, 0.5, 2.0, 0.5, 0.0),
    nrow = 3,
    byrow = TRUE
  )

  expected_cross <- apply(cross_cost, 1, function(row) {
    log(sum(weights * exp(-row / lambda)))
  })
  expected_loo <- vapply(
    seq_len(nrow(self_cost)),
    function(i) {
      keep <- setdiff(seq_len(nrow(self_cost)), i)
      renorm <- weights[keep] / sum(weights[keep])
      log(sum(renorm * exp(-self_cost[i, keep] / lambda)))
    },
    numeric(1)
  )

  testthat::expect_equal(
    causalOT:::cot_log_mixture(cross_cost, log_weights, lambda),
    expected_cross,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    causalOT:::cot_loo_log_mixture(self_cost, log_weights, weights, lambda),
    expected_loo,
    tolerance = 1e-12
  )

  expected_cross_loo <- vapply(
    seq_len(ncol(self_cost)),
    function(j) {
      keep <- setdiff(seq_len(nrow(self_cost)), j)
      renorm <- weights[keep] / sum(weights[keep])
      log(sum(renorm * exp(-self_cost[keep, j] / lambda)))
    },
    numeric(1)
  )
  testthat::expect_equal(
    causalOT:::cot_cross_loo_log_mixture(
      self_cost,
      log_weights,
      weights,
      lambda
    ),
    expected_cross_loo,
    tolerance = 1e-12
  )

  exclude <- matrix(
    c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    nrow = 3,
    byrow = TRUE
  )
  expected_masked <- vapply(
    seq_len(ncol(self_cost)),
    function(j) {
      keep <- which(!exclude[, j])
      renorm <- weights[keep] / sum(weights[keep])
      log(sum(renorm * exp(-self_cost[keep, j] / lambda)))
    },
    numeric(1)
  )
  testthat::expect_equal(
    causalOT:::cot_cross_masked_log_mixture(
      self_cost,
      log_weights,
      weights,
      lambda,
      exclude
    ),
    expected_masked,
    tolerance = 1e-12
  )
})

testthat::test_that("symmetric evidence cancels a shared volume term", {
  lambda <- 1.7
  d <- 2L
  x <- rbind(c(0, 0), c(1, 0), c(0, 1))
  y <- rbind(c(0.25, 0.1), c(1.1, -0.2), c(-0.1, 1.0))
  w <- c(0.25, 0.5, 0.25)
  a <- c(0.2, 0.3, 0.5)

  cross_cost <- sq_cost(x, y)
  self_w_cost <- sq_cost(x, x)
  self_a_cost <- sq_cost(y, y)

  raw <- causalOT:::cot_sym_evidence(
    cross_cost = cross_cost,
    self_a_cost = self_a_cost,
    self_w_cost = self_w_cost,
    w = w,
    a = a,
    lambda = lambda
  )

  norm_const <- -0.5 * d * log(2 * pi * lambda)
  cross_logmix <- causalOT:::cot_log_mixture(
    t(cross_cost),
    log(w),
    lambda
  ) +
    norm_const
  self_a_logmix <- causalOT:::cot_loo_log_mixture(
    self_a_cost,
    log(a),
    a,
    lambda
  ) +
    norm_const
  self_w_logmix <- causalOT:::cot_loo_log_mixture(
    self_w_cost,
    log(w),
    w,
    lambda
  ) +
    norm_const
  normalized <- sum(a * cross_logmix) -
    0.5 * sum(a * self_a_logmix) -
    0.5 * sum(w * self_w_logmix)

  testthat::expect_equal(raw, normalized, tolerance = 1e-12)
})

testthat::test_that("symmetric evidence cross term uses target likelihood under donor mixture", {
  lambda <- 1.3
  w <- c(0.4, 0.6)
  a <- c(0.3, 0.7)
  cross_cost <- matrix(c(0.2, 1.0, 1.5, 0.4), nrow = 2, byrow = TRUE)
  self_w_cost <- matrix(c(0.0, 0.7, 0.7, 0.0), nrow = 2, byrow = TRUE)
  self_a_cost <- matrix(c(0.0, 0.5, 0.5, 0.0), nrow = 2, byrow = TRUE)

  manual_cross <- vapply(
    seq_len(ncol(cross_cost)),
    function(j) {
      log(sum(w * exp(-cross_cost[, j] / lambda)))
    },
    numeric(1)
  )
  expected <- sum(a * manual_cross) -
    0.5 *
      sum(a * causalOT:::cot_loo_log_mixture(self_a_cost, log(a), a, lambda)) -
    0.5 *
      sum(w * causalOT:::cot_loo_log_mixture(self_w_cost, log(w), w, lambda))

  observed <- causalOT:::cot_sym_evidence(
    cross_cost = cross_cost,
    self_a_cost = self_a_cost,
    self_w_cost = self_w_cost,
    w = w,
    a = a,
    lambda = lambda
  )

  testthat::expect_equal(observed, expected, tolerance = 1e-12)
})

testthat::test_that("symmetric evidence applies cross LOO when donor and target supports coincide", {
  lambda <- 0.9
  w <- c(0.2, 0.3, 0.5)
  a <- c(0.1, 0.4, 0.5)
  shared_cost <- matrix(
    c(0.0, 0.6, 1.2, 0.6, 0.0, 0.3, 1.2, 0.3, 0.0),
    nrow = 3,
    byrow = TRUE
  )

  manual_cross_loo <- vapply(
    seq_len(ncol(shared_cost)),
    function(j) {
      keep <- setdiff(seq_len(nrow(shared_cost)), j)
      renorm <- w[keep] / sum(w[keep])
      log(sum(renorm * exp(-shared_cost[keep, j] / lambda)))
    },
    numeric(1)
  )
  expected <- sum(a * manual_cross_loo) -
    0.5 *
      sum(a * causalOT:::cot_loo_log_mixture(shared_cost, log(a), a, lambda)) -
    0.5 *
      sum(w * causalOT:::cot_loo_log_mixture(shared_cost, log(w), w, lambda))

  observed <- causalOT:::cot_sym_evidence(
    cross_cost = shared_cost,
    self_a_cost = shared_cost,
    self_w_cost = shared_cost,
    cross_exclude = diag(TRUE, nrow(shared_cost)),
    w = w,
    a = a,
    lambda = lambda
  )

  testthat::expect_equal(observed, expected, tolerance = 1e-12)
})

testthat::test_that("unit exclusion mask leaves duplicate covariates alone when ids differ", {
  source_x <- rbind(c(0, 0), c(0, 0))
  target_x <- rbind(c(0, 0), c(0, 0))
  mask <- causalOT:::cot_unit_exclusion_mask(
    source_unit_id = c("s1", "s2"),
    target_unit_id = c("t1", "t2")
  )

  testthat::expect_equal(dim(mask), c(nrow(source_x), nrow(target_x)))
  testthat::expect_equal(mask, matrix(FALSE, nrow(source_x), nrow(target_x)))
})

testthat::test_that("unit exclusion mask only excludes the shared unit id", {
  mask <- causalOT:::cot_unit_exclusion_mask(
    source_unit_id = c("shared", "donor_only"),
    target_unit_id = c("target_only", "shared", "other")
  )

  expected <- matrix(
    c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    nrow = 2,
    byrow = TRUE
  )
  testthat::expect_equal(mask, expected)
})

testthat::test_that("covariate overlap exclusion mask fallback catches donor subset inside larger target sample", {
  source_x <- rbind(c(0, 0), c(1, 0))
  target_x <- rbind(c(-1, -1), c(0, 0), c(2, 2), c(1, 0))
  mask <- causalOT:::cot_overlap_exclusion_mask(source_x, target_x)

  expected <- matrix(
    c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    nrow = 2,
    byrow = TRUE
  )
  testthat::expect_equal(mask, expected)
})

testthat::test_that("covariate overlap exclusion mask fallback matches duplicate covariates one-to-one", {
  source_x <- rbind(c(0, 0), c(0, 0))
  target_x <- rbind(c(0, 0), c(0, 0), c(1, 1))
  mask <- causalOT:::cot_overlap_exclusion_mask(source_x, target_x)

  testthat::expect_equal(dim(mask), c(2L, 3L))
  testthat::expect_equal(rowSums(mask), c(1, 1))
  testthat::expect_equal(colSums(mask), c(1, 1, 0))
})

testthat::test_that("sym evidence state uses exact unit ids before covariate fallback", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(419)

  x <- rbind(c(0, 0), c(0, 0))
  y <- rbind(c(0, 0), c(0, 0), c(1, 0))
  m1 <- causalOT:::Measure(
    x = x,
    adapt = "weights",
    unit.id = c("source_a", "shared")
  )
  m2 <- causalOT:::Measure(
    x = y,
    adapt = "none",
    unit.id = c("shared", "target_b", "target_c")
  )
  ot <- causalOT:::OTProblem(m1, m2)
  ot$setup_arguments(
    lambda = 5,
    debias = TRUE,
    cost.online = "tensorized"
  )
  ot$solve(
    niter = 3L,
    tol = 1e-2,
    torch_optim = torch::optim_rmsprop,
    torch_args = list(lr = 1e-2)
  )

  state <- ot$.__enclos_env__$private$get_sym_evidence_state()

  testthat::expect_equal(state$source_unit_id, c("source_a", "shared"))
  testthat::expect_equal(
    state$target_unit_id,
    c("shared", "target_b", "target_c")
  )
  testthat::expect_equal(
    state$cross_exclude,
    matrix(
      c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
      nrow = 2,
      byrow = TRUE
    )
  )
})

testthat::test_that("symmetric evidence applies cross exclusion for donor subset of target sample", {
  lambda <- 1.1
  w <- c(0.25, 0.75)
  a <- c(0.2, 0.3, 0.1, 0.4)
  cross_cost <- matrix(
    c(3.0, 0.0, 2.0, 1.0, 4.0, 1.5, 2.5, 0.0),
    nrow = 2,
    byrow = TRUE
  )
  self_w_cost <- matrix(c(0.0, 0.6, 0.6, 0.0), nrow = 2, byrow = TRUE)
  self_a_cost <- matrix(
    c(
      0.0,
      3.0,
      1.0,
      4.0,
      3.0,
      0.0,
      2.0,
      1.5,
      1.0,
      2.0,
      0.0,
      2.5,
      4.0,
      1.5,
      2.5,
      0.0
    ),
    nrow = 4,
    byrow = TRUE
  )
  exclude <- matrix(
    c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    nrow = 2,
    byrow = TRUE
  )

  manual_cross <- vapply(
    seq_len(ncol(cross_cost)),
    function(j) {
      keep <- which(!exclude[, j])
      renorm <- w[keep] / sum(w[keep])
      log(sum(renorm * exp(-cross_cost[keep, j] / lambda)))
    },
    numeric(1)
  )
  expected <- sum(a * manual_cross) -
    0.5 *
      sum(a * causalOT:::cot_loo_log_mixture(self_a_cost, log(a), a, lambda)) -
    0.5 *
      sum(w * causalOT:::cot_loo_log_mixture(self_w_cost, log(w), w, lambda))

  observed <- causalOT:::cot_sym_evidence(
    cross_cost = cross_cost,
    self_a_cost = self_a_cost,
    self_w_cost = self_w_cost,
    cross_exclude = exclude,
    w = w,
    a = a,
    lambda = lambda
  )

  testthat::expect_equal(observed, expected, tolerance = 1e-12)
})

testthat::test_that("mean_minus_sd summary is stable with one bootstrap draw", {
  score_matrix <- matrix(c(1.5, -0.25), nrow = 1L)

  testthat::expect_equal(
    causalOT:::cot_summarize_bootstrap_curve(
      score_matrix,
      summary = "mean_minus_sd",
      kappa = 3
    ),
    c(1.5, -0.25),
    tolerance = 1e-12
  )
})

testthat::test_that("helper selectors reject invalid lambdas and treat zero weights as zero mass", {
  cost <- matrix(c(0.0, 1.0, 1.0, 0.0), nrow = 2, byrow = TRUE)
  weights <- c(0.75, 0.25)
  log_weights <- log(weights)
  invalid_lambda <- list(0, -1, Inf, NA_real_, NaN)

  for (lambda in invalid_lambda) {
    testthat::expect_error(
      causalOT:::cot_log_mixture(cost, log_weights, lambda)
    )
    testthat::expect_error(
      causalOT:::cot_loo_log_mixture(cost, log_weights, weights, lambda)
    )
    testthat::expect_error(
      causalOT:::cot_cross_loo_log_mixture(cost, log_weights, weights, lambda)
    )
    testthat::expect_error(
      causalOT:::cot_cross_masked_log_mixture(
        cost,
        log_weights,
        weights,
        lambda,
        diag(TRUE, 2L)
      )
    )
  }

  weights_zero <- c(0.5, 0.0, 0.5)
  log_weights_zero <- causalOT:::cot_safe_log_weights(weights_zero)
  cost_zero <- matrix(
    c(0.1, 1.0, 0.3, 0.4, 0.2, 0.6),
    nrow = 2,
    byrow = TRUE
  )
  expected_zero <- apply(cost_zero, 1, function(row) {
    keep <- weights_zero > 0
    log(sum(weights_zero[keep] * exp(-row[keep] / 1.3)))
  })
  testthat::expect_equal(
    causalOT:::cot_log_mixture(cost_zero, log_weights_zero, 1.3),
    expected_zero,
    tolerance = 1e-12
  )
})

testthat::test_that("helper selectors return -Inf when no positive support remains", {
  cost <- matrix(
    c(0.0, 0.5, 1.0, 0.5, 0.0, 0.25, 1.0, 0.25, 0.0),
    nrow = 3,
    byrow = TRUE
  )
  weights <- c(1.0, 0.0, 0.0)
  log_weights <- causalOT:::cot_safe_log_weights(weights)

  loo_self <- causalOT:::cot_loo_log_mixture(cost, log_weights, weights, 0.9)
  loo_cross <- causalOT:::cot_cross_loo_log_mixture(cost, log_weights, weights, 0.9)
  testthat::expect_equal(loo_self[[1L]], -Inf)
  testthat::expect_equal(loo_cross[[1L]], -Inf)
  testthat::expect_true(all(is.finite(loo_self[-1L])))
  testthat::expect_true(all(is.finite(loo_cross[-1L])))
  testthat::expect_equal(
    causalOT:::cot_cross_masked_log_mixture(
      cost,
      log_weights,
      weights,
      0.9,
      matrix(TRUE, nrow = 3L, ncol = 3L)
    ),
    rep(-Inf, 3L)
  )
})

testthat::test_that("pilot weights normalize internally and invalid pilot weights fail clearly", {
  testthat::expect_equal(
    causalOT:::cot_resolve_pilot_weights(c(2, 2, 6), 3L),
    c(0.2, 0.2, 0.6),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    causalOT:::cot_resolve_pilot_weights(NULL, 3L),
    rep(1 / 3, 3L),
    tolerance = 1e-12
  )

  testthat::expect_error(
    causalOT:::cot_resolve_pilot_weights(c(1, 2), 3L),
    "must have length 3"
  )
  testthat::expect_error(
    causalOT:::cot_resolve_pilot_weights(c(1, -1, 2), 3L),
    "nonnegative"
  )
  testthat::expect_error(
    causalOT:::cot_resolve_pilot_weights(c(0, 0, 0), 3L),
    "positive total mass"
  )
})

testthat::test_that("bootstrap_legacy routes through the legacy selector", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(410)

  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(20), 10, 2)
  m1 <- causalOT:::Measure(x = x, adapt = "weights")
  m2 <- causalOT:::Measure(x = y, adapt = "none")
  ot <- causalOT:::OTProblem(m1, m2)
  ot$setup_arguments(
    lambda = c(10, 1),
    debias = TRUE,
    cost.online = "tensorized"
  )
  ot$solve(
    niter = 3L,
    tol = 1e-2,
    torch_optim = torch::optim_rmsprop,
    torch_args = list(lr = 1e-2)
  )

  priv <- ot$.__enclos_env__$private
  lambda_input <- priv$penalty_input_list$lambda

  set.seed(777)
  ot$choose_hyperparameters(
    n_boot_lambda = 1L,
    lambda.selector = "bootstrap_legacy"
  )

  metrics <- ot$info()$hyperparam.metrics$lambda
  testthat::expect_equal(
    ot$selected_lambda,
    lambda_input[[which.min(metrics)[1L]]]
  )
  testthat::expect_equal(
    ot$info()$hyperparam.metrics$lambda.selector,
    "bootstrap_legacy"
  )
})

testthat::test_that("cotOptions normalizes selector aliases and rejects unsupported selectors", {
  opt_sym <- causalOT:::cotOptions(lambda.selector = "sym_evidence")
  opt_boot <- causalOT:::cotOptions(
    lambda.selector = "legacy_bootstrap_balance"
  )

  testthat::expect_equal(opt_sym$lambda.selector, "profile_symmetric_loo")
  testthat::expect_equal(opt_boot$lambda.selector, "bootstrap_legacy")
  testthat::expect_error(
    causalOT:::cotOptions(lambda.selector = "kfold_symmetric"),
    "not implemented yet"
  )
})

testthat::test_that("fixed-pilot marginal selector ignores fitted path weights and scores only interior lambdas", {
  ot <- selector_test_problem()
  priv <- ot$.__enclos_env__$private
  state <- selector_state_fixture()
  old_state <- replace_private_binding(
    priv,
    "get_sym_evidence_state",
    function() state
  )
  on.exit({
    replace_private_binding(priv, "get_sym_evidence_state", old_state)
  }, add = TRUE)

  lambda_grid <- c(Inf, 2, 0.5, 0)
  weights_list_a <- list(
    make_weight_env(c(0.99, 0.005, 0.005)),
    make_weight_env(c(0.85, 0.10, 0.05)),
    make_weight_env(c(0.10, 0.10, 0.80)),
    make_weight_env(c(0.05, 0.05, 0.90))
  )
  weights_list_b <- list(
    make_weight_env(c(0.05, 0.90, 0.05)),
    make_weight_env(c(0.20, 0.20, 0.60)),
    make_weight_env(c(0.70, 0.20, 0.10)),
    make_weight_env(c(0.34, 0.33, 0.33))
  )

  res_default <- priv$select_lambda_metric(
    weights_list = weights_list_a,
    weights_converged = rep(TRUE, length(lambda_grid)),
    lambda_values = lambda_grid,
    lambda_input_values = lambda_grid,
    lambda_selector = "marginal_loo",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(
      list(details = TRUE)
    ),
    n_boot_lambda = 0L,
    lambda_bootstrap = 1
  )
  res_explicit_uniform <- priv$select_lambda_metric(
    weights_list = weights_list_b,
    weights_converged = rep(TRUE, length(lambda_grid)),
    lambda_values = lambda_grid,
    lambda_input_values = lambda_grid,
    lambda_selector = "marginal_loo",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(
      list(
        pilot.weights = rep(1, state$source_n),
        details = TRUE
      )
    ),
    n_boot_lambda = 0L,
    lambda_bootstrap = 1
  )

  pilot_uniform <- rep(1 / state$source_n, state$source_n)
  expected_curve <- c(
    NA_real_,
    causalOT:::cot_marginal_evidence(
      cross_cost = state$cross_cost,
      cross_exclude = state$cross_exclude,
      w = pilot_uniform,
      a = state$target_weights,
      lambda = 2
    ),
    causalOT:::cot_marginal_evidence(
      cross_cost = state$cross_cost,
      cross_exclude = state$cross_exclude,
      w = pilot_uniform,
      a = state$target_weights,
      lambda = 0.5
    ),
    NA_real_
  )

  testthat::expect_equal(res_default$metrics, expected_curve, tolerance = 1e-12)
  testthat::expect_equal(
    res_default$metrics,
    res_explicit_uniform$metrics,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    res_default$details$pilot.weights,
    pilot_uniform,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    res_default$selected_lambda,
    lambda_grid[[which.max(expected_curve)]]
  )
})

testthat::test_that("fixed-pilot symmetric selector uses one lambda across all evidence terms", {
  ot <- selector_test_problem()
  priv <- ot$.__enclos_env__$private
  state <- selector_state_fixture(
    cross_exclude = matrix(
      c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
      nrow = 3,
      byrow = TRUE
    )
  )
  old_state <- replace_private_binding(
    priv,
    "get_sym_evidence_state",
    function() state
  )
  on.exit({
    replace_private_binding(priv, "get_sym_evidence_state", old_state)
  }, add = TRUE)

  lambda_grid <- c(Inf, 1.5, 0.4)
  pilot <- c(2, 3, 5)
  weights_list <- list(
    make_weight_env(c(0.1, 0.2, 0.7)),
    make_weight_env(c(0.7, 0.2, 0.1)),
    make_weight_env(c(0.3, 0.3, 0.4))
  )
  res <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = rep(TRUE, length(lambda_grid)),
    lambda_values = lambda_grid,
    lambda_input_values = lambda_grid,
    lambda_selector = "symmetric_loo",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(
      list(pilot.weights = pilot)
    ),
    n_boot_lambda = 0L,
    lambda_bootstrap = 1
  )

  pilot_norm <- pilot / sum(pilot)
  expected_curve <- c(
    NA_real_,
    causalOT:::cot_sym_evidence(
      cross_cost = state$cross_cost,
      self_a_cost = state$self_a_cost,
      self_w_cost = state$self_w_cost,
      cross_exclude = state$cross_exclude,
      w = pilot_norm,
      a = state$target_weights,
      lambda = 1.5
    ),
    causalOT:::cot_sym_evidence(
      cross_cost = state$cross_cost,
      self_a_cost = state$self_a_cost,
      self_w_cost = state$self_w_cost,
      cross_exclude = state$cross_exclude,
      w = pilot_norm,
      a = state$target_weights,
      lambda = 0.4
    )
  )

  testthat::expect_equal(res$metrics, expected_curve, tolerance = 1e-12)
  testthat::expect_equal(
    res$selected_lambda,
    lambda_grid[[which.max(expected_curve)]]
  )
})

testthat::test_that("profile marginal selector uses the matching fitted path weights for each lambda", {
  ot <- selector_test_problem()
  priv <- ot$.__enclos_env__$private
  state <- selector_state_fixture()
  old_state <- replace_private_binding(
    priv,
    "get_sym_evidence_state",
    function() state
  )
  on.exit({
    replace_private_binding(priv, "get_sym_evidence_state", old_state)
  }, add = TRUE)

  lambda_grid <- c(Inf, 2, 0.5)
  donor_weights <- list(
    c(0.1, 0.1, 0.8),
    c(0.8, 0.1, 0.1),
    c(0.2, 0.6, 0.2)
  )
  res <- priv$select_lambda_metric(
    weights_list = lapply(donor_weights, make_weight_env),
    weights_converged = rep(TRUE, length(lambda_grid)),
    lambda_values = lambda_grid,
    lambda_input_values = lambda_grid,
    lambda_selector = "profile_marginal_loo",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(),
    n_boot_lambda = 0L,
    lambda_bootstrap = 1
  )

  expected_curve <- c(
    NA_real_,
    causalOT:::cot_marginal_evidence(
      cross_cost = state$cross_cost,
      cross_exclude = state$cross_exclude,
      w = donor_weights[[2L]],
      a = state$target_weights,
      lambda = 2
    ),
    causalOT:::cot_marginal_evidence(
      cross_cost = state$cross_cost,
      cross_exclude = state$cross_exclude,
      w = donor_weights[[3L]],
      a = state$target_weights,
      lambda = 0.5
    )
  )

  testthat::expect_equal(res$metrics, expected_curve, tolerance = 1e-12)
  testthat::expect_equal(
    res$selected_lambda,
    lambda_grid[[which.max(expected_curve)]]
  )
})

testthat::test_that("profile selectors fall back explicitly when eligible path points are not converged", {
  ot <- selector_test_problem()
  priv <- ot$.__enclos_env__$private
  state <- selector_state_fixture()
  old_state <- replace_private_binding(
    priv,
    "get_sym_evidence_state",
    function() state
  )
  old_draw_boot_measure <- replace_private_binding(
    priv,
    "draw_boot_measure",
    function() NULL
  )
  old_set_lambda <- replace_private_binding(
    priv,
    "set_lambda",
    function(lambda) {
      invisible(lambda)
    }
  )
  old_eval_lambda <- replace_private_binding(
    priv,
    "eval_lambda",
    function(w, boot) {
      w$metric
    }
  )
  on.exit({
    replace_private_binding(priv, "get_sym_evidence_state", old_state)
    replace_private_binding(priv, "draw_boot_measure", old_draw_boot_measure)
    replace_private_binding(priv, "eval_lambda", old_eval_lambda)
    replace_private_binding(priv, "set_lambda", old_set_lambda)
  }, add = TRUE)

  lambda_grid <- c(2, 0.5)
  weights_list <- list(
    make_weight_env(c(0.8, 0.1, 0.1), metric = 3),
    make_weight_env(c(0.2, 0.2, 0.6), metric = 1)
  )

  res <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, FALSE),
    lambda_values = lambda_grid,
    lambda_input_values = lambda_grid,
    lambda_selector = "profile_symmetric_loo",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(),
    n_boot_lambda = 1L,
    lambda_bootstrap = 0.25
  )

  testthat::expect_equal(res$selector, "bootstrap_legacy")
  testthat::expect_equal(res$metrics, c(3, 1))
  testthat::expect_equal(res$selected_lambda, 0.5)
  testthat::expect_match(res$note, "not converged")
})

testthat::test_that("profile symmetric selector reports explicit unsupported fallback reasons", {
  ot <- selector_test_problem()
  priv <- ot$.__enclos_env__$private
  old_state <- replace_private_binding(
    priv,
    "get_sym_evidence_state",
    function() {
      list(supported = FALSE, reason = "online")
    }
  )
  old_draw_boot_measure <- replace_private_binding(
    priv,
    "draw_boot_measure",
    function() NULL
  )
  old_set_lambda <- replace_private_binding(
    priv,
    "set_lambda",
    function(lambda) {
      invisible(lambda)
    }
  )
  old_eval_lambda <- replace_private_binding(
    priv,
    "eval_lambda",
    function(w, boot) {
      w$metric
    }
  )
  on.exit({
    replace_private_binding(priv, "get_sym_evidence_state", old_state)
    replace_private_binding(priv, "draw_boot_measure", old_draw_boot_measure)
    replace_private_binding(priv, "eval_lambda", old_eval_lambda)
    replace_private_binding(priv, "set_lambda", old_set_lambda)
  }, add = TRUE)

  weights_list <- list(
    make_weight_env(c(0.8, 0.1, 0.1), metric = 1),
    make_weight_env(c(0.2, 0.2, 0.6), metric = 2)
  )

  res_online <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = "profile_symmetric_loo",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(),
    n_boot_lambda = 1L,
    lambda_bootstrap = 0.25
  )
  testthat::expect_equal(res_online$selector, "bootstrap_legacy")
  testthat::expect_match(res_online$note, "online costs")

  replace_private_binding(
    priv,
    "get_sym_evidence_state",
    function() {
      list(supported = FALSE, reason = "structure")
    }
  )
  res_structure <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = "profile_symmetric_loo",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(),
    n_boot_lambda = 1L,
    lambda_bootstrap = 0.25
  )
  testthat::expect_equal(res_structure$selector, "bootstrap_legacy")
  testthat::expect_match(
    res_structure$note,
    "single-problem weight-fitting tensorized COT paths"
  )
})

testthat::test_that("mass bootstrap selector is reproducible and avoids legacy bootstrap machinery", {
  ot <- selector_test_problem()
  priv <- ot$.__enclos_env__$private
  state <- selector_state_fixture()
  old_state <- replace_private_binding(
    priv,
    "get_sym_evidence_state",
    function() state
  )
  old_draw_boot_measure <- replace_private_binding(
    priv,
    "draw_boot_measure",
    function() {
      stop("legacy bootstrap draw should not be used")
    }
  )
  old_set_lambda <- replace_private_binding(
    priv,
    "set_lambda",
    function(lambda) {
      stop("lambda refits should not be triggered")
    }
  )
  old_eval_lambda <- replace_private_binding(
    priv,
    "eval_lambda",
    function(w, boot) {
      stop("legacy bootstrap scoring should not be used")
    }
  )
  on.exit({
    replace_private_binding(priv, "get_sym_evidence_state", old_state)
    replace_private_binding(priv, "draw_boot_measure", old_draw_boot_measure)
    replace_private_binding(priv, "eval_lambda", old_eval_lambda)
    replace_private_binding(priv, "set_lambda", old_set_lambda)
  }, add = TRUE)

  weights_list <- list(
    make_weight_env(c(0.7, 0.2, 0.1)),
    make_weight_env(c(0.2, 0.3, 0.5))
  )

  set.seed(901)
  res_a <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = "mass_bootstrap_symmetric",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(
      list(details = TRUE)
    ),
    n_boot_lambda = 4L,
    lambda_bootstrap = 0.25
  )
  set.seed(901)
  res_b <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = "mass_bootstrap_symmetric",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(
      list(details = TRUE)
    ),
    n_boot_lambda = 4L,
    lambda_bootstrap = 0.25
  )

  testthat::expect_equal(res_a$metrics, res_b$metrics, tolerance = 1e-12)
  testthat::expect_equal(
    res_a$details$bootstrap.scores,
    res_b$details$bootstrap.scores,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    weights_list[[1L]][["source"]]$detach(),
    c(0.7, 0.2, 0.1),
    tolerance = 1e-12
  )
})

testthat::test_that("mass bootstrap summaries select from the requested bootstrap functional", {
  ot <- selector_test_problem()
  priv <- ot$.__enclos_env__$private
  state <- selector_state_fixture()
  old_state <- replace_private_binding(
    priv,
    "get_sym_evidence_state",
    function() state
  )
  on.exit({
    replace_private_binding(priv, "get_sym_evidence_state", old_state)
  }, add = TRUE)

  weights_list <- list(
    make_weight_env(c(0.7, 0.2, 0.1)),
    make_weight_env(c(0.2, 0.3, 0.5))
  )

  set.seed(902)
  res_mean <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = "mass_bootstrap_symmetric",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(
      list(details = TRUE, mass.bootstrap.summary = "mean")
    ),
    n_boot_lambda = 6L,
    lambda_bootstrap = 0.25
  )
  expected_mean <- causalOT:::cot_summarize_bootstrap_curve(
    res_mean$details$bootstrap.scores,
    summary = "mean"
  )
  expected_mean_idx <- if (any(is.finite(expected_mean))) {
    which.max(expected_mean)
  } else {
    1L
  }
  testthat::expect_equal(
    res_mean$selected_lambda,
    c(2, 0.5)[[expected_mean_idx]]
  )

  set.seed(902)
  res_quantile <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = "mass_bootstrap_symmetric",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(
      list(
        details = TRUE,
        mass.bootstrap.summary = "quantile",
        mass.bootstrap.quantile = 0.25
      )
    ),
    n_boot_lambda = 6L,
    lambda_bootstrap = 0.25
  )
  expected_quantile <- causalOT:::cot_summarize_bootstrap_curve(
    res_quantile$details$bootstrap.scores,
    summary = "quantile",
    quantile = 0.25
  )
  expected_quantile_idx <- if (any(is.finite(expected_quantile))) {
    which.max(expected_quantile)
  } else {
    1L
  }
  testthat::expect_equal(
    res_quantile$selected_lambda,
    c(2, 0.5)[[expected_quantile_idx]]
  )

  set.seed(902)
  res_mean_minus_sd <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = "mass_bootstrap_symmetric",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(
      list(
        details = TRUE,
        mass.bootstrap.summary = "mean_minus_sd",
        mass.bootstrap.kappa = 1.5
      )
    ),
    n_boot_lambda = 6L,
    lambda_bootstrap = 0.25
  )
  expected_mean_minus_sd <- causalOT:::cot_summarize_bootstrap_curve(
    res_mean_minus_sd$details$bootstrap.scores,
    summary = "mean_minus_sd",
    kappa = 1.5
  )
  expected_mean_minus_sd_idx <- if (any(is.finite(expected_mean_minus_sd))) {
    which.max(expected_mean_minus_sd)
  } else {
    1L
  }
  testthat::expect_equal(
    res_mean_minus_sd$selected_lambda,
    c(2, 0.5)[[expected_mean_minus_sd_idx]]
  )

  set.seed(902)
  res_no_details <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = "mass_bootstrap_symmetric",
    lambda_selector_control = causalOT:::cot_normalize_lambda_selector_control(),
    n_boot_lambda = 6L,
    lambda_bootstrap = 0.25
  )
  testthat::expect_null(res_no_details$details)
})

testthat::test_that("legacy bootstrap alias matches canonical selector and honors lambda.bootstrap", {
  ot <- selector_test_problem()
  priv <- ot$.__enclos_env__$private
  old_set_lambda <- replace_private_binding(
    priv,
    "set_lambda",
    function(lambda) {
      set_lambda_calls <<- c(set_lambda_calls, lambda)
      invisible(lambda)
    }
  )
  old_draw_boot_measure <- replace_private_binding(
    priv,
    "draw_boot_measure",
    function() NULL
  )
  old_eval_lambda <- replace_private_binding(
    priv,
    "eval_lambda",
    function(w, boot) {
      w$metric
    }
  )
  on.exit({
    replace_private_binding(priv, "set_lambda", old_set_lambda)
    replace_private_binding(priv, "draw_boot_measure", old_draw_boot_measure)
    replace_private_binding(priv, "eval_lambda", old_eval_lambda)
  }, add = TRUE)
  set_lambda_calls <- numeric()

  weights_list <- list(
    make_weight_env(c(0.8, 0.1, 0.1), metric = 2),
    make_weight_env(c(0.2, 0.2, 0.6), metric = 1)
  )
  control <- causalOT:::cot_normalize_lambda_selector_control()

  res_canonical <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = "bootstrap_legacy",
    lambda_selector_control = control,
    n_boot_lambda = 2L,
    lambda_bootstrap = 0.75
  )
  set.seed(904)
  res_alias <- priv$select_lambda_metric(
    weights_list = weights_list,
    weights_converged = c(TRUE, TRUE),
    lambda_values = c(2, 0.5),
    lambda_input_values = c(2, 0.5),
    lambda_selector = causalOT:::cot_normalize_lambda_selector(
      "legacy_bootstrap_balance"
    ),
    lambda_selector_control = control,
    n_boot_lambda = 2L,
    lambda_bootstrap = 0.75
  )

  testthat::expect_equal(res_canonical$selector, "bootstrap_legacy")
  testthat::expect_equal(res_alias$selector, "bootstrap_legacy")
  testthat::expect_equal(res_canonical$metrics, c(2, 1))
  testthat::expect_equal(res_alias$metrics, c(2, 1))
  testthat::expect_equal(res_canonical$selected_lambda, 0.5)
  testthat::expect_true(all(set_lambda_calls == 0.75))
})

testthat::test_that("kfold selector errors through choose_hyperparameters as not implemented", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(903)

  x <- matrix(stats::rnorm(18), 9, 2)
  y <- matrix(stats::rnorm(16), 8, 2)
  ot <- causalOT:::OTProblem(
    causalOT:::Measure(x = x, adapt = "weights"),
    causalOT:::Measure(x = y, adapt = "none")
  )
  ot$setup_arguments(lambda = c(5, 1), debias = TRUE, cost.online = "tensorized")
  ot$solve(
    niter = 3L,
    tol = 1e-2,
    torch_optim = torch::optim_rmsprop,
    torch_args = list(lr = 1e-2)
  )

  testthat::expect_error(
    ot$choose_hyperparameters(
      n_boot_lambda = 1L,
      lambda.selector = "kfold_symmetric"
    ),
    "not implemented yet"
  )
})

testthat::test_that("sym_evidence selector code stays outcome-free", {
  x <- matrix(0, 4, 2)
  y <- matrix(1, 5, 2)
  ot <- causalOT:::OTProblem(
    causalOT:::Measure(x = x, adapt = "weights"),
    causalOT:::Measure(x = y, adapt = "none")
  )
  priv <- ot$.__enclos_env__$private
  selector_code <- paste(
    c(
      deparse(body(causalOT:::cot_sym_evidence)),
      deparse(body(causalOT:::cot_log_mixture)),
      deparse(body(causalOT:::cot_loo_log_mixture)),
      deparse(body(priv$select_lambda_metric)),
      deparse(body(priv$get_sym_evidence_state)),
      deparse(body(priv$eval_lambda_sym_evidence))
    ),
    collapse = "\n"
  )

  forbidden <- c("Y", "outcome.model", "effect.estimate", "residual")
  for (token in forbidden) {
    testthat::expect_false(grepl(token, selector_code, fixed = TRUE))
  }
})

testthat::test_that("sym_evidence integrates with primal finite-lambda COT and excludes boundary lambdas", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(412)

  x <- matrix(stats::rnorm(28), 14, 2)
  y <- matrix(stats::rnorm(24), 12, 2)
  opt <- causalOT:::cotOptions(
    lambda = c(Inf, 25, 5),
    opt.direction = "primal",
    solver = "lbfgs",
    cost.online = "tensorized",
    niter = 10L,
    nboot = 1L,
    tol = 1e-2,
    lambda.selector = "sym_evidence",
    line_search_fn = "strong_wolfe",
    history_size = 5L
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  cot$solve()
  force_converged_flags(cot)

  res <- cot$grid_search()

  testthat::expect_equal(
    res$metric$lambda.selector,
    "profile_symmetric_loo"
  )
  testthat::expect_equal(
    res$metric$lambda.loo,
    causalOT:::cot_selector_loo_convention
  )
  testthat::expect_true(is.na(res$metric$lambda[[1L]]))
  testthat::expect_true(is.finite(res$penalty[["lambda"]]))
  testthat::expect_true(res$penalty[["lambda"]] %in% c(25, 5))
})

testthat::test_that("sym_evidence integrates with dual finite-lambda COT", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(413)

  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(20), 10, 2)
  opt <- causalOT:::cotOptions(
    lambda = c(20, 5),
    opt.direction = "dual",
    cost.online = "tensorized",
    niter = 10L,
    nboot = 1L,
    tol = 1e-2,
    lambda.selector = "sym_evidence",
    lr = 1e-2
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  cot$solve()
  force_converged_flags(cot)

  res <- cot$grid_search()

  testthat::expect_equal(
    res$metric$lambda.selector,
    "profile_symmetric_loo"
  )
  testthat::expect_true(all(is.finite(res$metric$lambda)))
  testthat::expect_true(res$penalty[["lambda"]] %in% c(20, 5))
})

testthat::test_that("profile_marginal_loo integrates with dual finite-lambda COT", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(415)

  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(20), 10, 2)
  opt <- causalOT:::cotOptions(
    lambda = c(20, 5),
    opt.direction = "dual",
    cost.online = "tensorized",
    niter = 10L,
    nboot = 1L,
    tol = 1e-2,
    lambda.selector = "profile_marginal_loo",
    lr = 1e-2
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  cot$solve()
  force_converged_flags(cot)

  res <- cot$grid_search()

  testthat::expect_equal(res$metric$lambda.selector, "profile_marginal_loo")
  testthat::expect_equal(
    res$metric$lambda.loo,
    causalOT:::cot_selector_marginal_loo_convention
  )
  testthat::expect_true(all(is.finite(res$metric$lambda)))
  testthat::expect_true(res$penalty[["lambda"]] %in% c(20, 5))
})

testthat::test_that("symmetric_loo stores normalized pilot masses in optional details", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(416)

  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(20), 10, 2)
  pilot <- seq_len(nrow(x))
  pilot <- pilot / sum(pilot)
  opt <- causalOT:::cotOptions(
    lambda = c(20, 5),
    opt.direction = "dual",
    cost.online = "tensorized",
    niter = 10L,
    nboot = 1L,
    tol = 1e-2,
    lambda.selector = "symmetric_loo",
    lambda.selector.control = list(
      pilot.weights = pilot,
      details = TRUE
    ),
    lr = 1e-2
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  cot$solve()
  force_converged_flags(cot)

  res <- cot$grid_search()

  testthat::expect_equal(res$metric$lambda.selector, "symmetric_loo")
  testthat::expect_equal(
    res$metric$lambda.details$pilot.weights,
    pilot,
    tolerance = 1e-12
  )
})

testthat::test_that("mass_bootstrap_symmetric stores bootstrap details when requested", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(417)

  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(20), 10, 2)
  opt <- causalOT:::cotOptions(
    lambda = c(20, 5),
    opt.direction = "dual",
    cost.online = "tensorized",
    niter = 10L,
    nboot = 2L,
    tol = 1e-2,
    lambda.selector = "mass_bootstrap_symmetric",
    lambda.selector.control = list(
      details = TRUE,
      mass.bootstrap.summary = "quantile",
      mass.bootstrap.quantile = 0.25
    ),
    lr = 1e-2
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  cot$solve()
  force_converged_flags(cot)

  res <- cot$grid_search()

  testthat::expect_equal(
    res$metric$lambda.selector,
    "mass_bootstrap_symmetric"
  )
  testthat::expect_equal(
    res$metric$lambda.details$mass.bootstrap.summary,
    "quantile"
  )
  testthat::expect_equal(
    res$metric$lambda.details$mass.bootstrap.quantile,
    0.25
  )
  testthat::expect_equal(
    dim(res$metric$lambda.details$bootstrap.scores),
    c(2L, 2L)
  )
})

testthat::test_that("mass_bootstrap_symmetric handles mean_minus_sd with one bootstrap draw", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(418)

  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(20), 10, 2)
  opt <- causalOT:::cotOptions(
    lambda = c(20, 5),
    opt.direction = "dual",
    cost.online = "tensorized",
    niter = 10L,
    nboot = 1L,
    tol = 1e-2,
    lambda.selector = "mass_bootstrap_symmetric",
    lambda.selector.control = list(
      mass.bootstrap.summary = "mean_minus_sd",
      mass.bootstrap.kappa = 1
    ),
    lr = 1e-2
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  cot$solve()
  force_converged_flags(cot)

  res <- cot$grid_search()

  testthat::expect_equal(
    res$metric$lambda.selector,
    "mass_bootstrap_symmetric"
  )
  testthat::expect_true(all(is.finite(res$metric$lambda)))
  testthat::expect_true(res$penalty[["lambda"]] %in% c(20, 5))
})

testthat::test_that("selector aliases match canonical curves in integration and selected weights stay on the chosen path", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(420)

  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(20), 10, 2)

  fit_once <- function(selector) {
    opt <- causalOT:::cotOptions(
      lambda = c(20, 5),
      opt.direction = "dual",
      cost.online = "tensorized",
      niter = 10L,
      nboot = 2L,
      tol = 1e-2,
      lambda.selector = selector,
      lr = 1e-2
    )
    cot <- causalOT:::COT$new(source = x, target = y, options = opt)
    cot$solve()
    force_converged_flags(cot)
    cot
  }

  cot_alias <- fit_once("sym_evidence")
  cot_canonical <- fit_once("profile_symmetric_loo")
  res_alias <- cot_alias$grid_search()
  res_canonical <- cot_canonical$grid_search()

  testthat::expect_equal(
    res_alias$metric$lambda,
    res_canonical$metric$lambda,
    tolerance = 1e-10
  )
  testthat::expect_equal(
    res_alias$penalty[["lambda"]],
    res_canonical$penalty[["lambda"]]
  )

  ot <- causalOT:::OTProblem(
    causalOT:::Measure(x = x, adapt = "weights"),
    causalOT:::Measure(x = y, adapt = "none")
  )
  ot$setup_arguments(lambda = c(20, 5), debias = TRUE, cost.online = "tensorized")
  ot$solve(
    niter = 10L,
    tol = 1e-2,
    torch_optim = torch::optim_rmsprop,
    torch_args = list(lr = 1e-2)
  )
  priv <- ot$.__enclos_env__$private
  priv$converged_run <- lapply(priv$converged_run, function(x) {
    rep(TRUE, length(x))
  })
  ot$choose_hyperparameters(
    n_boot_lambda = 1L,
    lambda.selector = "profile_symmetric_loo"
  )

  selected_env <- priv$weights_list[[as.character(ot$selected_lambda)]][[1L]]
  selected_key <- ls(selected_env)[[1L]]
  testthat::expect_equal(
    causalOT:::cot_to_numeric_vector(
      priv$measures[[selected_key]]$weights$detach()
    ),
    causalOT:::cot_to_numeric_vector(selected_env[[selected_key]]$detach()),
    tolerance = 1e-10
  )
})

testthat::test_that("calc_weight COT stores lambda selector metadata in causalWeights info", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(421)

  x <- matrix(stats::rnorm(36), 18, 2)
  z <- rep(c(0L, 1L), each = 9L)
  cw <- causalOT:::calc_weight(
    x = x,
    z = z,
    estimand = "ATT",
    method = "COT",
    options = list(
      lambda = c(10, 2),
      opt.direction = "dual",
      cost.online = "tensorized",
      niter = 10L,
      nboot = 1L,
      tol = 1e-2,
      lambda.selector = "marginal_loo",
      lr = 1e-2
    )
  )

  testthat::expect_s4_class(cw, "causalWeights")
  testthat::expect_true("metric" %in% names(cw@info))
  testthat::expect_equal(cw@info$metric$lambda.selector, "marginal_loo")
  testthat::expect_equal(cw@info$metric$lambda.grid, c(10, 2))
  testthat::expect_equal(
    cw@penalty[["lambda"]],
    cw@info$metric$lambda.selected
  )
  testthat::expect_true(all(is.finite(cw@info$metric$lambda)))
})

testthat::test_that("calc_weight ATE stores separate control and treated lambda selector metadata", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(422)

  x <- matrix(stats::rnorm(40), 20, 2)
  z <- rep(c(0L, 1L), each = 10L)
  cw <- causalOT:::calc_weight(
    x = x,
    z = z,
    estimand = "ATE",
    method = "COT",
    options = list(
      lambda = c(12, 3),
      opt.direction = "dual",
      cost.online = "tensorized",
      niter = 10L,
      nboot = 1L,
      tol = 1e-2,
      lambda.selector = "symmetric_loo",
      lr = 1e-2
    )
  )

  testthat::expect_s4_class(cw, "causalWeights")
  testthat::expect_equal(cw@estimand, "ATE")
  testthat::expect_true(all(c("w0", "w1") %in% names(cw@info)))

  for (arm in c("w0", "w1")) {
    metric <- cw@info[[arm]]$metric
    testthat::expect_equal(metric$lambda.selector, "symmetric_loo")
    testthat::expect_equal(metric$lambda.grid, c(12, 3))
    testthat::expect_equal(
      cw@penalty[[arm]][["lambda"]],
      metric$lambda.selected
    )
    testthat::expect_true(all(is.finite(metric$lambda)))
  }
})

testthat::test_that("fixed-lambda CCP behavior is not broken by selector changes", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(414)

  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(20), 10, 2)
  opt <- causalOT:::cotOptions(
    lambda = 50,
    opt.direction = "primal",
    solver = "mirror",
    cost.online = "tensorized",
    niter = 3L,
    nboot = 1L,
    ot_niter = 50L,
    ot_tol = 1e-5,
    tol = 1e-5,
    lambda.selector = "sym_evidence"
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)

  testthat::expect_silent(cot$solve())
  testthat::expect_silent(res <- cot$grid_search())
  testthat::expect_equal(res$penalty[["lambda"]], 50)
})
