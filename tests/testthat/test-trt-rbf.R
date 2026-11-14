test_that("centered kernel score equals explicit deleted fits", {
  x <- matrix(c(-2, -1, -.3, .1, .4, 1, 2, 3), ncol = 1)
  z <- c(0, 0, 0, 0, 0, 1, 1, 1)
  C <- as.matrix(dist(x))^2 / 2
  t <- 1.3; alpha <- .07; n <- length(z)
  score <- cot_trt_rbf_score(C, z, t, alpha)$score
  K <- exp(-C / t)
  pred <- vapply(seq_len(n), function(i) {
    keep <- setdiff(seq_len(n), i)
    # Block system for an unpenalized intercept and fixed full-data penalty.
    A <- rbind(c(0, rep(1, n - 1)), cbind(1, K[keep, keep] + diag(alpha * n, n - 1)))
    coef <- solve(A, c(0, z[keep]))
    coef[1] + sum(K[i, keep] * coef[-1])
  }, numeric(1))
  expect_equal(score, mean((z - pred)^2), tolerance = 1e-9)
  S <- K %*% solve(K + diag(alpha * n, n))
  prototype <- mean(((z - mean(z) - S %*% (z - mean(z))) / (1 - diag(S) - 1/n))^2)
  expect_gt(abs(score - prototype), 1e-6)
})

test_that("treatment selection respects solver costs and scale", {
  x <- matrix(seq(-2, 2, length.out = 12), ncol = 1)
  z <- rep(c(0, 1), 6)
  d <- dataHolder(x, z)
  o <- cotOptions(lambda = "auto", lambda.selector.control = list(details = TRUE))
  a <- cot_select_trt_rbf(d, o)
  b <- cot_select_trt_rbf(dataHolder(x * 3, z), o)
  C <- cot_to_numeric_matrix(cost(x, x, p = 2))
  expect_equal(a$median_cost, median((as.matrix(dist(x))^2)[upper.tri(C)]) / 2)
  expect_equal(exp(-C/a$lambda), exp(-(C/a$median_cost)/a$multiplier))
  expect_equal(a$multiplier, b$multiplier, tolerance = 1e-6)
  expect_equal(b$lambda, 9 * a$lambda, tolerance = 1e-6)
  expect_equal(a$lambda, a$median_cost * a$multiplier)
})

test_that("search reports budget and handles flat and multimodal objectives", {
  f <- function(t) list(score = 1, alpha = 1)
  a <- cot_trt_rbf_search(f, max_eval = 50)
  expect_equal(a$multiplier, max(a$evaluations$multiplier))
  expect_equal(cot_trt_rbf_search(f, max_eval = 8)$status, "budget_exhausted")
  f2 <- function(t) list(score = min((log(t) + 1)^2 + .1, (log(t) - 2)^2), alpha = 1)
  a <- cot_trt_rbf_search(f2)
  ref <- min(vapply(exp(seq(-3, 6, length.out = 2000)), function(t) f2(t)$score, numeric(1)))
  expect_lt(a$score - ref, 1e-4)
})

test_that("input validation and deterministic stratified sampling work", {
  x <- matrix(seq_len(80), ncol = 2); z <- rep(0:1, each = 20)
  o <- cotOptions(lambda = "auto")
  expect_error(cot_select_trt_rbf(dataHolder(x, z), cotOptions(lambda = "auto", p = 1)), "p = 2")
  expect_error(cot_select_trt_rbf(dataHolder(x * 0, z), o), "positive median")
  expect_error(cot_select_trt_rbf(dataHolder(x, z, weights = seq_len(40)), o), "uniform")
  set.seed(33); state <- .Random.seed
  a <- cot_select_trt_rbf(dataHolder(x, z), o, max_sample = 12, max_eval = 10)
  expect_identical(.Random.seed, state)
  b <- cot_select_trt_rbf(dataHolder(x, z), o, max_sample = 12, max_eval = 10)
  expect_identical(a$sample_indices, b$sample_indices)
  expect_true(a$median_is_subsample)
  expect_true(all(table(z[a$sample_indices]) >= 2))
})

test_that("selected lambda is shared and passed directly to ATE fits", {
  set.seed(4)
  x <- matrix(rnorm(24), ncol = 1); z <- rep(0:1, each = 12)
  w <- suppressWarnings(calc_weight(x, z, estimand = "ATE", method = "COT",
    options = list(lambda = "auto", niter = 3, nboot = 0, polish = "none")))
  expect_identical(unname(w@penalty$w0[1]), w@info$lambda_selection$lambda)
  expect_identical(unname(w@penalty$w1[1]), w@info$lambda_selection$lambda)
  expect_null(w@info$lambda_selection$evaluations)
  # without pooled data from calc_weight, "auto" fits the default grid
  cot <- COT$new(x[z == 0, , drop = FALSE], x[z == 1, , drop = FALSE],
                 options = list(lambda = "auto"))
  expect_length(cot$penalty$lambda, 8L)
})

test_that("search expands improving boundaries and rejects invalid scores", {
  f <- function(t) list(score = -log(t), alpha = 1)
  a <- cot_trt_rbf_search(f, max_expand = 2)
  expect_equal(a$expansions, 2L)
  expect_true(a$boundary[["upper"]])
  expect_error(cot_trt_rbf_search(function(t) list(score = Inf, alpha = 1)), "No valid")
  # non-binary labels cannot be predicted, so "auto" uses the bootstrap
  w <- suppressWarnings(calc_weight(matrix(1:8, ncol = 1), c(0, 0, 0, .5, 1, 1, 1, 1),
    estimand = "ATT", method = "COT", options = list(niter = 2, nboot = 1, grid.length = 2)))
  expect_identical(w@info$metric$lambda.selector, "bootstrap")
  expect_null(w@info$lambda_selection)
})

test_that("automatic mode does not break omitted lambda in direct solver setup", {
  ot <- OTProblem(Measure(matrix(1:4, ncol = 1), adapt = "weights"),
                  Measure(matrix(5:8, ncol = 1)))
  expect_no_error(ot$setup_arguments())
  expect_warning(ot$setup_arguments(lambda = "auto"), "already set up")
  expect_true(is.finite(ot$penalty$lambda))
})

test_that("integer p and single precision fits can select a stable bandwidth", {
  set.seed(1)
  d <- dataHolder(matrix(rnorm(300), 100, 3), rep(0:1, each = 50))
  a <- cot_select_trt_rbf(d, cotOptions(lambda = "auto", p = 2L,
                         dtype = torch::torch_float()), max_eval = 10)
  b <- cot_select_trt_rbf(d, cotOptions(lambda = "auto",
                         dtype = torch::torch_double()), max_eval = 10)
  expect_true(is.finite(a$lambda) && a$lambda > 0)
  expect_equal(a$lambda, b$lambda)
})

test_that("small search budgets return their best candidate without overrunning", {
  f <- function(t) list(score = log(t)^2, alpha = 1)
  for (budget in c(1, 3, 7)) {
    a <- cot_trt_rbf_search(f, max_eval = budget)
    expect_equal(a$status, "budget_exhausted")
    expect_equal(a$n_evaluations, budget)
    expect_equal(a$score, min(a$evaluations$score))
  }
  expect_error(cot_trt_rbf_search(f, max_eval = 0), "positive integer")
  expect_error(cot_trt_rbf_search(f, max_expand = -1), "nonnegative integer")
})

test_that("boundary expansion prioritizes the better of two improving edges", {
  # Both tails improve, but the upper tail has the better objective.
  f <- function(t) list(score = -log(t)^2, alpha = 1)
  a <- cot_trt_rbf_search(f, max_expand = 1)
  expect_equal(a$bounds[2], 2 * 10^3.5)
  expect_true(a$boundary[["upper"]])
})


test_that("grid selection scores solver-unit candidates directly", {
  x <- matrix(c(-3, -2, -1, 0, .3, 1, 2, 4), ncol = 1)
  z <- rep(0:1, each = 4)
  prep <- cot_trt_rbf_prepare(cot_trt_rbf_context(x, z))
  grid <- c(10, 0, 1, Inf, 10, 2)
  res <- cot_trt_rbf_grid(prep, grid)
  expected <- vapply(grid[c(1, 3, 5, 6)], function(v)
    cot_trt_rbf_score(prep$C, prep$z, v / prep$median_cost)$score, numeric(1))
  expect_equal(res$scores[c(1, 3, 5, 6)], expected)
  expect_true(all(is.infinite(res$scores[c(2, 4)])))
  expect_equal(res$scores[1], res$scores[5])
  expect_equal(res$lambda, grid[res$selected_index])
  expect_error(cot_trt_rbf_grid(prep, c(0, Inf)), "No eligible")
  scaled <- prep; scaled$median_cost <- 9 * prep$median_cost
  res2 <- cot_trt_rbf_grid(scaled, grid * 9)
  expect_equal(res2$scores, res$scores)
  expect_equal(res2$selected_index, res$selected_index)
  expect_equal(cot_trt_rbf_score(prep$C, 1 - prep$z, 1)$score,
               cot_trt_rbf_score(prep$C, prep$z, 1)$score)
})

test_that("OTProblem selects and reselects fitted RBF grids", {
  x <- matrix(seq_len(8), ncol = 1)
  ot <- OTProblem(Measure(x[1:4, , drop = FALSE], adapt = "weights"),
                  Measure(x[5:8, , drop = FALSE]))
  ot$setup_arguments(lambda = c(1, 10, 2))
  suppressWarnings(ot$solve(niter = 2))
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
  info <- ot$info()$hyperparam.metrics
  expect_identical(info$lambda.selector, "auto")
  expect_identical(info$lambda.details$original_index, c(2L, 3L, 1L))
  expect_equal(info$lambda.selected, info$lambda.grid[info$lambda.details$selected_index])
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
  expect_equal(ot$info()$hyperparam.metrics$lambda, info$lambda)
})

test_that("RBF contexts validate identities and unsupported geometry", {
  x <- matrix(seq_len(8), ncol = 1)
  expect_error(cot_trt_rbf_context(x, rep(0:1, each = 4), rep(1, 8)), "unique")
  a <- Measure(x[1:4, , drop = FALSE], unit.id = 1:4, adapt = "weights")
  b <- Measure(x[5:8, , drop = FALSE], unit.id = 4:7)
  ot <- OTProblem(a, b)
  expect_error(ot$setup_arguments(lambda = "auto"), "Overlapping")
  expect_no_error(ot$setup_arguments(lambda = "auto",
    lambda.selection.data = list(x = x, z = rep(0:1, each = 4), unit.id = 1:8)))
})


test_that("calc_weight uses pooled data for single-arm RBF grids", {
  set.seed(11)
  x <- matrix(rnorm(20), ncol = 1); z <- rep(0:1, each = 10)
  options <- list(lambda = c(1, 2), lambda.selector = "auto", niter = 2, nboot = 0)
  w <- suppressWarnings(calc_weight(x, z, estimand = "ATT", method = "COT", options = options))
  expect_identical(w@info$metric$lambda.selector, "auto")
  expect_equal(w@info$metric$lambda.details$sample_size, length(z))
  ate <- suppressWarnings(calc_weight(x, z, estimand = "ATE", method = "COT", options = options))
  expect_equal(unname(ate@penalty$w0[1]), unname(ate@penalty$w1[1]))
})


test_that("dual selectors restore their exact stored fitted state", {
  set.seed(18)
  x <- matrix(rnorm(16), ncol = 1)
  for (cls in list(cotDualLegacyTrain, cotDualTrain, cotDualEntropyFixedPointTrain)) {
    ot <- cls$new(Measure(x[1:8, , drop = FALSE], adapt = "weights"),
                  Measure(x[9:16, , drop = FALSE]))
    ot$setup_arguments(lambda = c(2, 1, 2))
    suppressWarnings(ot$solve(niter = 2, torch_optim = torch::optim_rmsprop))
    suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
    p <- ot$.__enclos_env__$private
    info <- ot$info()$hyperparam.metrics$lambda.details
    saved <- p$weights_list[[info$selected_index]][[1L]]
    expected <- as.numeric(saved[[ls(saved)[1L]]])
    expect_equal(as.numeric(ot$weights), expected, tolerance = 1e-6)
    expect_identical(info$original_index[info$selected_index],
                     p$rbf_input_order[info$selected_index])
    suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
    expect_equal(as.numeric(ot$weights), expected, tolerance = 1e-6)
  }
})

test_that("shared RBF eligibility excludes candidates failed in either ATE arm", {
  set.seed(19)
  x <- matrix(rnorm(20), ncol = 1); z <- rep(0:1, each = 10)
  opts <- cotOptions(lambda = c(.1, 1, 10), lambda.selector = "auto", niter = 2, nboot = 0)
  opts$lambda.selection.data <- list(x = x, z = z)
  pair <- gridSearch(dataHolder(x, z), "ATE", "COT", opts)
  a <- pair@control@solver; b <- pair@treated@solver
  suppressWarnings(a$solve()); suppressWarnings(b$solve())
  bp <- b$.__enclos_env__$private$optimizer$.__enclos_env__$private
  w <- bp$weights_list[[1L]][[1L]]
  for (key in ls(w)) w[[key]] <- torch::torch_full_like(w[[key]], NA_real_)
  p0 <- a$rbf_path_eligibility(); p1 <- b$rbf_path_eligibility()
  expect_false(p1$eligible[1L])
  common <- p0$eligible & p1$eligible
  shared <- a$rbf_shared_scores()
  r0 <- suppressWarnings(a$grid_search(.rbf.eligible = common, .rbf.shared = shared))
  r1 <- suppressWarnings(b$grid_search(.rbf.eligible = common, .rbf.shared = shared))
  expect_equal(r0$penalty[[1L]], r1$penalty[[1L]])
  expect_false(r0$metric$lambda.details$eligible[1L])
  expect_equal(r0$metric$lambda, r1$metric$lambda)
})

test_that("adaptive OTProblem and pooled selector agree", {
  x <- matrix(c(-3, -2, -1, -.5, .2, .7, 2, 4), ncol = 1)
  z <- rep(0:1, each = 4)
  direct <- cot_select_trt_rbf(dataHolder(x, z), cotOptions(lambda = "auto"))
  make <- function(context = NULL) {
    ot <- OTProblem(Measure(x[1:4, , drop = FALSE], adapt = "weights"),
                    Measure(x[5:8, , drop = FALSE]))
    ot$setup_arguments(lambda = "auto", lambda.selection.data = context,
                       lambda.selector.control = list(details = TRUE))
    ot
  }
  inferred <- make(); explicit <- make(list(x = x, z = z))
  expect_equal(inferred$penalty$lambda, direct$lambda, tolerance = 1e-7)
  expect_equal(explicit$penalty$lambda, direct$lambda, tolerance = 1e-7)
  suppressWarnings(inferred$solve(niter = 2))
  suppressWarnings(inferred$choose_hyperparameters())
  diagnostics <- inferred$info()$hyperparam.metrics
  expect_identical(diagnostics$lambda.selector, "auto")
  expect_identical(diagnostics$lambda.details$mode, "adaptive")
  expect_equal(diagnostics$lambda.details$score, direct$score, tolerance = 1e-7)
  expect_true(nrow(diagnostics$lambda.details$evaluations) > 1L)
})

test_that("adaptive snapshots retain original uniform masses after fitting", {
  x <- matrix(seq_len(10), ncol = 1)
  a <- Measure(x[1:5, , drop = FALSE], adapt = "weights")
  b <- Measure(x[6:10, , drop = FALSE])
  ot <- OTProblem(a, b)
  ot$setup_arguments(lambda = c(1, 2))
  suppressWarnings(ot$solve(niter = 2))
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
  before <- ot$info()$hyperparam.metrics$lambda
  a$weights <- torch::torch_tensor(c(.4, .3, .1, .1, .1), dtype = a$dtype, device = a$device)
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
  expect_equal(ot$info()$hyperparam.metrics$lambda, before)
  expect_warning(ot$setup_arguments(lambda = 1), "already set up")
  suppressWarnings(ot$solve(niter = 2))
  expect_no_error(suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto")))
})

test_that("replacing the RBF context invalidates cached scores", {
  x <- matrix(c(-3, -2, -1, 0, .1, 1, 2, 4), ncol = 1)
  z <- rep(0:1, each = 4)
  ot <- OTProblem(Measure(x[1:4, , drop = FALSE], adapt = "weights"),
                  Measure(x[5:8, , drop = FALSE]))
  ot$setup_arguments(lambda = c(1, 3), lambda.selection.data = list(x = x, z = z))
  suppressWarnings(ot$solve(niter = 2))
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
  before <- ot$info()$hyperparam.metrics$lambda
  z2 <- rep(0:1, 4)
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto",
    lambda.selection.data = list(x = x, z = z2)))
  after <- ot$info()$hyperparam.metrics$lambda
  expect_gt(max(abs(before - after)), 1e-6)
  expect_error(cot_trt_rbf_context(x, factor(z)), "binary")
})

test_that("composite objectives use one explicit pooled criterion", {
  x <- matrix(seq_len(8), ncol = 1); z <- rep(0:1, each = 4)
  a <- Measure(x[1:4, , drop = FALSE], adapt = "weights")
  b <- Measure(x[5:8, , drop = FALSE])
  term <- OTProblem(a, b)
  composite <- term + term
  expect_error(composite$setup_arguments(lambda = "auto"), "Composite")
  suppressWarnings(composite$setup_arguments(lambda = c(1, 2),
    lambda.selection.data = list(x = x, z = z)))
  suppressWarnings(composite$solve(niter = 2))
  suppressWarnings(composite$choose_hyperparameters(lambda.selector = "auto"))
  observed <- composite$info()$hyperparam.metrics$lambda
  prepared <- cot_trt_rbf_prepare(cot_trt_rbf_context(x, z))
  expect_equal(observed, cot_trt_rbf_grid(prepared, c(2, 1))$scores)
})

test_that("zero and infinite input candidates cannot be selected", {
  x <- matrix(seq_len(8), ncol = 1)
  ot <- OTProblem(Measure(x[1:4, , drop = FALSE], adapt = "weights"),
                  Measure(x[5:8, , drop = FALSE]))
  ot$setup_arguments(lambda = c(0, 2, Inf, 2))
  suppressWarnings(ot$solve(niter = 2))
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
  info <- ot$info()$hyperparam.metrics
  expect_equal(info$lambda.selected, 2)
  expect_equal(info$lambda.details$eligible, c(FALSE, TRUE, TRUE, FALSE))
  expect_true(all(is.infinite(info$lambda[c(1, 4)])))
})

test_that("grid selection does not snap the adaptive optimum to the grid", {
  x <- matrix(c(-.591566350158614, .722962029282451, -.364993894031002,
    -.245858948803953, -.161587704274692, .71444891843344, .0653653255014037,
    1.41584111899645, .984983285227803, .290975219761947, -1.11889606371966,
    .933836353985498), ncol = 1)
  z <- c(1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1)
  prep <- cot_trt_rbf_prepare(cot_trt_rbf_context(x, z))
  adaptive <- cot_trt_rbf_search(function(t) cot_trt_rbf_score(prep$C, prep$z, t))
  grid <- prep$median_cost * c(.01, .1, 1, 10, 100)
  result <- cot_trt_rbf_grid(prep, grid)
  expect_false(result$selected_index == which.min(abs(grid - adaptive$multiplier * prep$median_cost)))
  expect_equal(result$selected_index, 1L)
})

test_that("mirror RBF selection preserves the fitted potentials", {
  x <- matrix(c(-2, -1, -.5, 0, .2, .6, 1, 3), ncol = 1)
  ot <- cotCCPTrain$new(Measure(x[1:4, , drop = FALSE], adapt = "weights"),
                        Measure(x[5:8, , drop = FALSE]))
  ot$setup_arguments(lambda = "auto")
  suppressWarnings(ot$solve(niter = 2))
  p <- ot$.__enclos_env__$private
  expect_length(p$weights_list, 1L)
  expect_length(p$state_run, 1L)
  if (length(p$state_run) != 1L) return(invisible(NULL))
  # No selection call may launch another transport fit.
  unlockBinding("optimization_loop", p)
  p$optimization_loop <- function(...) stop("Unexpected transport refit")
  lockBinding("optimization_loop", p)
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
  # one adaptive lambda, so there is nothing to choose between
  expect_identical(ot$info()$hyperparam.metrics$lambda.details$mode, "adaptive")
  i <- 1L
  saved <- p$state_run[[i]][[1L]]
  for (key in names(saved$potentials)) {
    actual <- p$ot_objects[[key]]$potentials
    expect_equal(lapply(actual, as.numeric), lapply(saved$potentials[[key]], as.numeric))
  }
  expect_length(ot$info()$hyperparam.metrics$lambda, 1L)
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "bootstrap"))
  suppressWarnings(ot$choose_hyperparameters(lambda.selector = "auto"))
  expect_identical(ot$info()$hyperparam.metrics$lambda.details$mode, "adaptive")
})
