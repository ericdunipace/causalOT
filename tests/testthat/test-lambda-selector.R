testthat::test_that("Measure ignores automatic row names for default unit ids", {
  causalOT:::torch_check()
  df <- data.frame(x = 1:2)
  testthat::expect_null(causalOT:::Measure(df)$unit_id)
})

testthat::test_that("Measure preserves explicit non-automatic row names as unit ids", {
  causalOT:::torch_check()
  df <- data.frame(x = 1:2)
  rownames(df) <- c("u1", "u2")
  testthat::expect_equal(causalOT:::Measure(df)$unit_id, c("u1", "u2"))
})

testthat::test_that("cotOptions defaults to auto lambda and auto selector", {
  causalOT:::torch_check()
  opt <- causalOT:::cotOptions()
  testthat::expect_equal(opt$lambda, "auto")
  testthat::expect_equal(opt$lambda.selector, "auto")
  testthat::expect_equal(opt$lambda.selector.control, list(details = FALSE))
  testthat::expect_equal(
    causalOT:::cotOptions(lambda.selector = "bootstrap")$lambda.selector,
    "bootstrap"
  )
  testthat::expect_null(causalOT:::cotOptions(lambda = NULL)$lambda)
  testthat::expect_setequal(causalOT:::cotOptions(lambda = c(1, 10))$lambda, c(10, 1))
})

testthat::test_that("cotOptions rejects unknown lambda selectors and controls", {
  causalOT:::torch_check()
  testthat::expect_error(
    causalOT:::cotOptions(lambda.selector = "other"),
    "must be one of 'auto', 'bootstrap'"
  )
  testthat::expect_error(causalOT:::cotOptions(lambda = "other"), "'auto'")
  testthat::expect_error(
    causalOT:::cotOptions(lambda.selector.control = list(other = 1)),
    "Unknown lambda.selector.control"
  )
  testthat::expect_error(
    causalOT:::cotOptions(lambda.selector.control = list(details = "yes")),
    "TRUE or FALSE"
  )
})

testthat::test_that("auto falls back to the bootstrap selector when the RBF rule cannot be used", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(415)

  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(20), 10, 2)
  ot <- causalOT:::OTProblem(
    causalOT:::Measure(x = x, adapt = "weights"),
    causalOT:::Measure(x = y, adapt = "none")
  )
  ot$setup_arguments(lambda = c(10, 1), p = 1, cost.online = "tensorized")
  ot$solve(niter = 3L, tol = 1e-2, torch_optim = torch::optim_rmsprop,
           torch_args = list(lr = 1e-2))
  # three iterations leave the fits unconverged, which warns
  suppressWarnings(ot$choose_hyperparameters(n_boot_lambda = 2L))

  info <- ot$info()$hyperparam.metrics
  testthat::expect_equal(info$lambda.selector, "bootstrap")
  testthat::expect_match(info$lambda.note, "auto used the bootstrap selector")
  testthat::expect_match(info$lambda.note, "p = 2")
  testthat::expect_true(all(is.finite(info$lambda)))
})

testthat::test_that("calc_weight uses the bootstrap selector for non-uniform weights", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(416)

  x <- matrix(stats::rnorm(36), 18, 2)
  z <- rep(c(0L, 1L), each = 9L)
  cw <- causalOT:::calc_weight(
    x = x, z = z, estimand = "ATT", method = "COT",
    weights = stats::runif(18),
    options = list(niter = 5L, nboot = 2L, grid.length = 3L)
  )

  testthat::expect_equal(cw@info$metric$lambda.selector, "bootstrap")
  testthat::expect_length(cw@info$metric$lambda.grid, 4L)
  testthat::expect_null(cw@info$lambda_selection)
})

testthat::test_that("calc_weight with default options selects lambda before fitting", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(417)

  x <- matrix(stats::rnorm(36), 18, 2)
  z <- rep(c(0L, 1L), each = 9L)
  cw <- causalOT:::calc_weight(x = x, z = z, estimand = "ATT", method = "COT")

  testthat::expect_equal(cw@info$lambda_selection$mode, "adaptive")
  testthat::expect_equal(cw@info$lambda_selection$method, "auto")
  testthat::expect_equal(unname(cw@penalty[["lambda"]]), cw@info$lambda_selection$lambda)
  testthat::expect_length(cw@info$metric$lambda.grid, 1L)
})

testthat::test_that("bootstrap selector picks the smallest bootstrap score", {
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
    cost.online = "tensorized",
    polish = "auto"
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
    lambda.selector = "bootstrap"
  )

  metrics <- ot$info()$hyperparam.metrics$lambda
  testthat::expect_equal(
    ot$selected_lambda,
    lambda_input[[which.min(metrics)[1L]]]
  )
  testthat::expect_equal(
    ot$info()$hyperparam.metrics$lambda.selector,
    "bootstrap"
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
      lambda.selector = "bootstrap",
      lr = 1e-2
    )
  )

  testthat::expect_s4_class(cw, "causalWeights")
  testthat::expect_true("metric" %in% names(cw@info))
  testthat::expect_equal(cw@info$metric$lambda.selector, "bootstrap")
  # dual solvers fit (and report) the grid from small to large lambda
  testthat::expect_setequal(cw@info$metric$lambda.grid, c(10, 2))
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
      lr = 1e-2
    )
  )

  testthat::expect_s4_class(cw, "causalWeights")
  testthat::expect_equal(cw@estimand, "ATE")
  testthat::expect_true(all(c("w0", "w1") %in% names(cw@info)))

  for (arm in c("w0", "w1")) {
    metric <- cw@info[[arm]]$metric
    testthat::expect_equal(metric$lambda.selector, "auto")
    testthat::expect_setequal(metric$lambda.grid, c(12, 3))
    testthat::expect_equal(
      cw@penalty[[arm]],
      metric$lambda.selected
    )
    testthat::expect_true(all(is.finite(metric$lambda)))
  }
  # both arms score the same pooled data, so they choose the same lambda
  testthat::expect_equal(cw@penalty$w0, cw@penalty$w1)
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
    polish = "auto",
    cost.online = "tensorized",
    niter = 3L,
    nboot = 1L,
    ot_niter = 50L,
    ot_tol = 1e-5,
    tol = 1e-5
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)

  testthat::expect_silent(cot$solve())
  testthat::expect_silent(res <- cot$grid_search())
  testthat::expect_equal(res$penalty[["lambda"]], 50)
})

testthat::test_that("bootstrap selector scores the candidate weights, not the fitted ones", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(1)
  x <- matrix(stats::rnorm(60), 30, 2)
  y <- matrix(stats::rnorm(80, 0.5), 40, 2)
  opts <- causalOT:::cotOptions(
    lambda = c(10, 1),
    nboot = 0L,
    lambda.selector = "bootstrap",
    cost.online = "tensorized",
    niter = 20L,
    solver = "lbfgs",
    device = torch::torch_device("cpu")
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opts)
  cot$solve()
  priv <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
  adds <- ls(priv$measures)
  add <- adds[vapply(adds, function(k) priv$measures[[k]]$adapt == "weights", logical(1))]
  dtype <- priv$measures[[add]]$init_weights$dtype

  uniform <- rlang::env()
  class(uniform) <- c("weightEnv", class(uniform))
  uniform[[add]] <- torch::torch_tensor(rep(1 / 30, 30), dtype = dtype)
  spiky <- rlang::env()
  class(spiky) <- c("weightEnv", class(spiky))
  spiky_w <- c(rep(1, 3), rep(1e-6, 27))
  spiky[[add]] <- torch::torch_tensor(spiky_w / sum(spiky_w), dtype = dtype)

  for (lambda_boot in c(Inf, 1)) {
    priv$set_lambda(lambda_boot)
    set.seed(2)
    boot <- priv$draw_boot_measure()
    testthat::expect_false(isTRUE(all.equal(
      priv$eval_lambda(uniform, boot),
      priv$eval_lambda(spiky, boot)
    )))
  }
})

testthat::test_that("bootstrap selector renormalizes the bootstrapped candidate weights", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(1)
  x <- matrix(stats::rnorm(60), 30, 2)
  y <- matrix(stats::rnorm(80, 0.5), 40, 2)
  opts <- causalOT:::cotOptions(
    lambda = c(10, 1),
    p = 1,
    nboot = 0L,
    lambda.selector = "bootstrap",
    cost.online = "tensorized",
    niter = 20L,
    solver = "lbfgs",
    device = torch::torch_device("cpu"),
    dtype = torch::torch_float64()
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opts)
  cot$solve()
  priv <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
  adds <- ls(priv$measures)
  add <- adds[vapply(adds, function(k) priv$measures[[k]]$adapt == "weights", logical(1))]
  add_y <- setdiff(adds, add)

  spiky <- rlang::env()
  class(spiky) <- c("weightEnv", class(spiky))
  spiky_w <- c(rep(1, 3), rep(1e-3, 27))
  spiky_w <- spiky_w / sum(spiky_w)
  spiky[[add]] <- torch::torch_tensor(spiky_w, dtype = torch::torch_float64())

  priv$set_lambda(Inf)
  set.seed(2)
  boot <- priv$draw_boot_measure()
  a <- spiky_w * as.numeric(boot[[add]])
  a <- a / sum(a)
  b <- as.numeric(boot[[add_y]]) / 40

  # p = 1 is the L1 cost, and the Inf-penalty divergence is half the energy distance
  l1 <- function(u, v) as.matrix(stats::dist(rbind(u, v), method = "manhattan"))[
    seq_len(nrow(u)), nrow(u) + seq_len(nrow(v)), drop = FALSE]
  energy <- 2 * sum(a * (l1(x, y) %*% b)) - sum(a * (l1(x, x) %*% a)) -
    sum(b * (l1(y, y) %*% b))
  testthat::expect_equal(priv$eval_lambda(spiky, boot), energy / 2, tolerance = 1e-8)
})
