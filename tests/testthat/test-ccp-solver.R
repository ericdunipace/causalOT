testthat::test_that("semi-relaxed initializer is log-domain column softmax", {
  causalOT:::torch_check()
  C <- torch::torch_tensor(
    matrix(
      c(0.0, 2.0, 1.0, 3.0, 1.0, 4.0, 2.0, 0.5, 1.5, 1.0, 3.0, 0.0),
      nrow = 4,
      byrow = TRUE
    ),
    dtype = torch::torch_double()
  )
  b <- torch::torch_tensor(c(0.2, 0.5, 0.3), dtype = torch::torch_double())
  lambda <- 10.0

  w <- causalOT:::cot_init_semirelaxed(C, b, lambda)
  eta <- -C / lambda
  expected <- (eta -
    eta$logsumexp(1)$view(c(1L, -1L)) +
    log(b)$view(c(1L, -1L)))$logsumexp(2)$exp()

  testthat::expect_true(all(as.logical((w >= 0)$to(device = "cpu"))))
  testthat::expect_equal(w$sum()$item(), 1.0, tolerance = 1e-10)
  testthat::expect_equal(as.numeric(w), as.numeric(expected), tolerance = 1e-10)
})

testthat::test_that("tilted softmax is invariant to additive tilt gauges", {
  causalOT:::torch_check()
  C <- torch::torch_tensor(
    matrix(stats::runif(20), nrow = 5),
    dtype = torch::torch_double()
  )
  b <- torch::torch_tensor(rep(0.25, 4), dtype = torch::torch_double())
  tilt <- torch::torch_tensor(stats::rnorm(5), dtype = torch::torch_double())

  w1 <- causalOT:::cot_tilted_softmax_weights(C, b, 100.0, tilt)
  w2 <- causalOT:::cot_tilted_softmax_weights(C, b, 100.0, tilt + 25.0)

  testthat::expect_equal(as.numeric(w1), as.numeric(w2), tolerance = 1e-10)
})

testthat::test_that("ccp solver runs fixed lambda primal tensorized COT", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(2141)
  n <- 15
  m <- 13
  d <- 2
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)

  opt <- cotOptions(
    lambda = 100.0,
    opt.direction = "primal",
    solver = "mirror",
    niter = 10L,
    ot_niter = 500L,
    ot_tol = 1e-8,
    tol = 1e-10
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)

  testthat::expect_silent(cot$solve())
  testthat::expect_equal(sum(cot$weights), 1.0, tolerance = 1e-6)
  testthat::expect_true(all(cot$weights >= 0))

  ots <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  ot <- ots[[ls(ots)[[1L]]]]
  diag <- causalOT:::cot_ccp_diagnostics(ot)

  testthat::expect_true(is.finite(diag$objective))
  testthat::expect_true(diag$objective >= -1e-7)
  testthat::expect_lt(diag$row_marginal_max_error, 1e-5)
  testthat::expect_lt(diag$column_marginal_max_error, 1e-5)
  testthat::expect_equal(
    cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ccp_kkt_resid,
    diag$kkt_max_deviation,
    tolerance = 1e-12
  )
})

testthat::test_that("ccp option rejects unsupported first-slice configurations", {
  causalOT:::torch_check()
  testthat::expect_error(
    cotOptions(lambda = 100.0, opt.direction = "dual", solver = "mirror"),
    "opt.direction"
  )
  testthat::expect_error(
    cotOptions(
      lambda = c(100.0, 10.0),
      opt.direction = "primal",
      solver = "mirror"
    ),
    "single finite positive lambda"
  )
  testthat::expect_error(
    cotOptions(
      lambda = 100.0,
      opt.direction = "primal",
      solver = "mirror",
      balance.formula = "~X1"
    ),
    "balance constraints"
  )
})

testthat::test_that("solver option routes to the expected optimizer class", {
  causalOT:::torch_check()
  set.seed(619)
  x <- matrix(stats::rnorm(20), 10, 2)
  y <- matrix(stats::rnorm(18), 9, 2)

  opt_ccp <- cotOptions(
    lambda = 100.0,
    opt.direction = "primal",
    solver = "mirror",
    cost.online = "tensorized"
  )
  cot_ccp <- causalOT:::COT$new(source = x, target = y, options = opt_ccp)
  testthat::expect_s3_class(
    cot_ccp$.__enclos_env__$private$optimizer,
    "cotCCPTrain"
  )

  opt_lbfgs <- cotOptions(
    lambda = 100.0,
    opt.direction = "primal",
    solver = "lbfgs",
    cost.online = "tensorized"
  )
  cot_lbfgs <- causalOT:::COT$new(source = x, target = y, options = opt_lbfgs)
  testthat::expect_s3_class(
    cot_lbfgs$.__enclos_env__$private$optimizer,
    "OTProblem"
  )
  testthat::expect_false(inherits(
    cot_lbfgs$.__enclos_env__$private$optimizer,
    "cotCCPTrain"
  ))
})

testthat::test_that("ccp solve carries updated weights across iterations", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(991)
  n <- 15
  m <- 13
  d <- 2
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)

  run_ccp <- function(niter) {
    opt <- cotOptions(
      lambda = 100.0,
      opt.direction = "primal",
      solver = "mirror",
      cost.online = "tensorized",
      niter = niter,
      ot_niter = 100L,
      ot_tol = 1e-6,
      tol = 1e-12
    )
    cot <- causalOT:::COT$new(source = x, target = y, options = opt)
    cot$solve()
    cot$weights
  }

  w1 <- run_ccp(1L)
  w2 <- run_ccp(2L)

  testthat::expect_equal(sum(w1), 1.0, tolerance = 1e-6)
  testthat::expect_equal(sum(w2), 1.0, tolerance = 1e-6)
  testthat::expect_gt(max(abs(w1 - w2)), 1e-8)
})

testthat::test_that("ccp default fast path tries one tau per iteration", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(416)
  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(22), 11, 2)

  opt <- cotOptions(
    lambda = 100.0,
    opt.direction = "primal",
    solver = "mirror",
    cost.online = "tensorized",
    niter = 3L,
    ot_niter = 50L,
    ot_tol = 1e-5,
    tol = 0
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  cot$solve()

  priv <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
  testthat::expect_false(priv$ccp_backtracking)
  testthat::expect_equal(priv$ccp_candidate_count, priv$ccp_iter)
})

testthat::test_that("ccp sparse KKT diagnostic matches post-hoc diagnostic", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(417)
  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(22), 11, 2)

  opt <- cotOptions(
    lambda = 100.0,
    opt.direction = "primal",
    solver = "mirror",
    cost.online = "tensorized",
    niter = 4L,
    ot_niter = 50L,
    ot_tol = 1e-5,
    tol = 0
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ccp_diagnostics_every <- 100L
  cot$solve()

  priv <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
  ot <- priv$ot_objects[[ls(priv$ot_objects)[[1L]]]]
  diag <- causalOT:::cot_ccp_diagnostics(ot)
  testthat::expect_equal(
    priv$ccp_kkt_resid,
    diag$kkt_max_deviation,
    tolerance = 1e-12
  )
})

testthat::test_that("ccp fast path does not mutate the fixed y self potential", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(418)
  x <- matrix(stats::rnorm(24), 12, 2)
  y <- matrix(stats::rnorm(22), 11, 2)

  opt <- cotOptions(
    lambda = 100.0,
    opt.direction = "primal",
    solver = "mirror",
    cost.online = "tensorized",
    niter = 2L,
    ot_niter = 50L,
    ot_tol = 1e-5,
    tol = 0
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  priv <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
  ot <- priv$ot_objects[[ls(priv$ot_objects)[[1L]]]]
  ot$sinkhorn_cot(niter = 50L, tol = 1e-5)
  g_yy_before <- as.numeric(ot$potentials$g_yy$to(device = "cpu"))
  cot$solve()

  testthat::expect_equal(
    as.numeric(ot$potentials$g_yy$to(device = "cpu")),
    g_yy_before,
    tolerance = 1e-12
  )
})

testthat::test_that("dual solve writes solved weights back to COT object", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(971)
  n <- 15
  m <- 13
  d <- 2
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)

  opt <- cotOptions(
    lambda = 100.0,
    opt.direction = "dual",
    cost.online = "tensorized",
    niter = 5L,
    ot_niter = 50L,
    ot_tol = 1e-5,
    tol = 1e-5
  )
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  start_weights <- cot$weights

  testthat::expect_silent(cot$solve())
  direct_weights <- as.numeric(
    cot$.__enclos_env__$private$optimizer$weights$to(device = "cpu")
  )

  testthat::expect_equal(cot$weights, direct_weights, tolerance = 1e-7)
  testthat::expect_equal(sum(cot$weights), 1.0, tolerance = 1e-7)
  testthat::expect_gt(max(abs(cot$weights - start_weights)), 1e-6)
})

testthat::test_that("ccp tilted update is scale invariant", {
  causalOT:::torch_check()
  C <- torch::torch_tensor(
    matrix(stats::runif(24), nrow = 6),
    dtype = torch::torch_double()
  )
  b <- torch::torch_tensor(rep(0.25, 4), dtype = torch::torch_double())
  tilt <- torch::torch_tensor(stats::rnorm(6), dtype = torch::torch_double())

  w1 <- causalOT:::cot_tilted_softmax_weights(C, b, 100.0, tilt)
  w2 <- causalOT:::cot_tilted_softmax_weights(C * 7.0, b, 700.0, tilt * 7.0)

  testthat::expect_equal(as.numeric(w1), as.numeric(w2), tolerance = 1e-10)
})
