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
