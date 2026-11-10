testthat::test_that("COT objects form", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)
  a <- rep(1 / n, n)
  b <- rep(1 / m, m)
  delta <- 0.5

  #cot defaults
  testthat::expect_silent(causalOT:::COT$new(source = x, target = y))

  #setup
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y))
  testthat::expect_silent(cot$weights <- c(1, rep(0, n - 1)))
  testthat::expect_true(
    length(cot$.__enclos_env__$private$optimizer$penalty$lambda) == 8
  )
  testthat::expect_true(is.na(
    cot$.__enclos_env__$private$optimizer$penalty$delta
  ))
  testthat::expect_equal(as_numeric(cot$weights), c(1, rep(0, n - 1)))
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(balance.formula = "~X1+X2")
    )
  )
  testthat::expect_true(inherits(
    cot$.__enclos_env__$private$optimizer,
    "OTProblem"
  ))

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(lambda = 0, debias = TRUE)
    )
  )

  # entropy not debiased
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = FALSE
      )
    )
  )
  testthat::expect_true(all(
    class(cot$.__enclos_env__$private$optimizer) != "cotDualTrain"
  ))
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = TRUE,
        balance.formula = "~X1*X3"
      )
    )
  )
  testthat::expect_true(
    length(cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$target_objects) == 1L
  )
  testthat::expect_true(inherits(
    cot$.__enclos_env__$private$optimizer,
    "OTProblem"
  ))
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = TRUE,
        opt.direction = "dual",
        dual.solver = "joint"
      )
    )
  )
  testthat::expect_true(inherits(
    cot$.__enclos_env__$private$optimizer,
    "cotDualTrain"
  ))
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = TRUE,
        opt.direction = "dual",
        dual.solver = "entropy"
      )
    )
  )
  testthat::expect_true(inherits(
    cot$.__enclos_env__$private$optimizer,
    "cotDualLegacyTrain"
  ))
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = TRUE,
        opt.direction = "dual",
        dual.solver = "fixed_point_joint",
        fixed_point_eta = 1,
        fixed_point_beta = 0.25,
        fixed_point_weight_floor = 1e-6,
        fixed_point_eta_anneal_every = 5L,
        fixed_point_eta_min = 0.1
      )
    )
  )
  testthat::expect_true(inherits(
    cot$.__enclos_env__$private$optimizer,
    "cotDualFixedPointJointTrain"
  ))
  testthat::expect_equal(
    cot$.__enclos_env__$private$solve_args$fixed_point_eta,
    1
  )
  testthat::expect_equal(
    cot$.__enclos_env__$private$solve_args$fixed_point_beta,
    0.25
  )
  testthat::expect_equal(
    cot$.__enclos_env__$private$solve_args$fixed_point_eta_anneal_every,
    5L
  )
  testthat::expect_equal(
    cot$.__enclos_env__$private$solve_args$fixed_point_eta_min,
    0.1
  )

  # debias without torch optim should throw error
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = TRUE,
        torch.optimizer = NULL
      )
    )
  )

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = TRUE,
        opt.direction = "primal"
      )
    )
  )
  testthat::expect_true(inherits(
    cot$.__enclos_env__$private$optimizer,
    "OTProblem"
  ))
})

testthat::test_that("ent works", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()

  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)
  a <- rep(1 / n, n)
  b <- rep(1 / m, m)
  delta <- 0.5
  lambda <- 10
  wlist <- list(a, a)

  #debias = FALSE
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(debias = FALSE, niter = 4, opt.direction = "primal")
    )
  )
  ot <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  testthat::expect_true(inherits(ot[[ls(ot)]]$C_xy, "costTensor"))
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$torch_optim, "optim_lbfgs")
  )
  testthat::expect_no_warning(cot$solve())
  cot <- causalOT:::COT$new(
    source = x,
    target = y,
    options = list(
      debias = FALSE,
      niter = 1,
      opt.direction = "primal",
      line_search_fn = "strong_wolfe"
    )
  )
  testthat::expect_equivalent(
    cot$.__enclos_env__$private$torch_optim_args,
    list(line_search_fn = "strong_wolfe")
  )
  testthat::expect_silent(cot$solve())

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(niter = 1, balance.formula = "~X1", lambda = c(1, 10), debias = FALSE)
    )
  )

  osqpout <- testthat::capture_output(testthat::expect_no_warning(cot$solve()))

  testthat::expect_true(
    as_numeric((cot$.__enclos_env__$private$source$balance_functions$transpose(
      2,
      1
    )$matmul(cot$.__enclos_env__$private$source$weights) -
      cot$.__enclos_env__$private$source$balance_target)$abs()$max()$to(
      device = "cpu"
    )$item()) <
      1.5
  )

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = FALSE,
        balance.formula = "~X1+X2",
        ot_niter = 50L,
        niter = 1,
        verbose = FALSE,
        delta = 0.1
      )
    )
  )
  osqpout <- testthat::capture_output(testthat::expect_no_warning(cot$solve()))
  testthat::expect_true(
    as_numeric((cot$.__enclos_env__$private$source$balance_functions$transpose(
      2,
      1
    )$matmul(cot$.__enclos_env__$private$source$weights) -
      cot$.__enclos_env__$private$source$balance_target)$abs()$max()$item()) <
      0.9
  )
})

testthat::test_that("ent debiased works, online", {
  causalOT:::torch_check()

  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)
  a <- rep(1 / n, n)
  b <- rep(1 / m, m)
  delta <- 0.5
  lambda <- 10
  wlist <- list(a, a)

  #debias = FALSE
  mess <- testthat::capture_output(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = FALSE,
        niter = 4,
        opt.direction = "primal",
        cost.online = "online"
      )
    )
  )
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$torch_optim, "optim_lbfgs")
  )
  testthat::expect_no_warning(cot$solve())
  cot <- causalOT:::COT$new(
    source = x,
    target = y,
    options = list(
      debias = FALSE,
      niter = 1,
      opt.direction = "primal",
      line_search_fn = "strong_wolfe",
      cost.online = "online"
    )
  )
  testthat::expect_equivalent(
    cot$.__enclos_env__$private$torch_optim_args,
    list(line_search_fn = "strong_wolfe")
  )
  testthat::expect_silent(cot$solve())

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        niter = 1,
        balance.formula = "~X1",
        lambda = c(1, 10),
        debias = FALSE,
        cost.online = "online"
      )
    )
  )

  osqpout <- testthat::capture_output(testthat::expect_no_warning(cot$solve()))

  testthat::expect_true(
    as_numeric(max(abs(
      cot$.__enclos_env__$private$source$balance_functions$transpose(
        2,
        1
      )$matmul(cot$.__enclos_env__$private$source$weights) -
        cot$.__enclos_env__$private$source$balance_target
    ))) <
      1.5
  )

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = FALSE,
        balance.formula = "~X1+X2",
        lambda = c(1, 10),
        niter = 1,
        verbose = FALSE,
        delta = 0.1,
        cost.online = "online"
      )
    )
  )
  osqpout <- testthat::capture_output(testthat::expect_no_warning(cot$solve()))
  testthat::expect_true(
    as_numeric(max(abs(
      cot$.__enclos_env__$private$source$balance_functions$transpose(
        2,
        1
      )$matmul(cot$.__enclos_env__$private$source$weights) -
        cot$.__enclos_env__$private$source$balance_target
    ))) <
      0.9
  )

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        debias = TRUE,
        balance.formula = "~X1+X2",
        lambda = c(10, Inf),
        niter = 1,
        verbose = FALSE,
        delta = 0.1,
        cost.online = "online"
      )
    )
  )

  osqpout <- testthat::capture_output(
    cot$solve()
  )

  testthat::expect_true(
    as_numeric(max(abs(
      cot$.__enclos_env__$private$source$balance_functions$transpose(
        2,
        1
      )$matmul(cot$.__enclos_env__$private$source$weights) -
        cot$.__enclos_env__$private$source$balance_target
    ))) <
      0.9
  )
})

testthat::test_that("grid_search function works", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()

  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)
  a <- rep(1 / n, n)
  b <- rep(1 / m, m)
  delta <- 0.5

  #setup
  torch::torch_manual_seed(123123)
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        verbose = FALSE,
        opt.direction = "dual",
        # torch.scheduler = NULL,
        niter = 1L,
        nboot = 2L
      )
    )
  )
  cot$solve()
  output <- cot$grid_search()
  testthat::expect_true(output$penalty[1] < 2704243.27850342 + 1)
  testthat::expect_warning(cot$grid_search())
})

testthat::test_that("NNM works", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()

  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)
  a <- rep(1 / n, n)
  b <- rep(1 / m, m)
  delta <- 0.5
  lambda <- 10
  wlist <- list(a, a)

  #### tensor ####
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        lambda = 0,
        debias = FALSE,
        niter = 4,
        opt.direction = "primal"
      )
    )
  )
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$optimizer, "NNM")
  )

  ot <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  testthat::expect_true(inherits(ot[[ls(ot)]]$C_xy, "costTensor"))
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$torch_optim, "optim_lbfgs")
  )
  testthat::expect_true(isFALSE(ot[[ls(ot)]]$debias))
  testthat::expect_invisible(cot$solve())
  res <- cot$grid_search()
  testthat::expect_true(all(
    names(res) %in% c("weight", "penalty", "metric", "penalty.grid")
  ))

  cot <- causalOT:::COT$new(
    source = x,
    target = y,
    options = list(debias = TRUE, lambda = 0)
  )
  ot <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  testthat::expect_true(isFALSE(ot[[ls(ot)]]$debias))

  cot <- causalOT:::COT$new(
    source = x,
    target = y,
    options = list(balance.formula = "~.", lambda = 0)
  )
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$optimizer, "cotDualFixedPointJointTrain")
  )

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(niter = 1, balance.formula = "~X1", lambda = c(1, 10), debias = FALSE)
    )
  )
  testthat::expect_equal(
    class(cot$.__enclos_env__$private$optimizer),
    c("OTProblem", "R6")
  )
  osqpout <- testthat::capture_output(testthat::expect_no_warning(cot$solve()))

  #### online ####

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        lambda = 0,
        debias = FALSE,
        niter = 4,
        cost.online = "online",
        opt.direction = "primal"
      )
    )
  )
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$optimizer, "NNM")
  )

  ot <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  testthat::expect_true(inherits(ot[[ls(ot)]]$C_xy, "costOnline"))
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$torch_optim, "optim_lbfgs")
  )
  testthat::expect_true(isFALSE(ot[[ls(ot)]]$debias))
  testthat::expect_invisible(cot$solve())
  res <- cot$grid_search()
  testthat::expect_true(all(
    names(res) %in% c("weight", "penalty", "metric", "penalty.grid")
  ))

  cot <- causalOT:::COT$new(
    source = x,
    target = y,
    options = list(debias = TRUE, cost.online = "online", lambda = 0)
  )
  ot <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  testthat::expect_true(isFALSE(ot[[ls(ot)]]$debias))

  cot <- causalOT:::COT$new(
    source = x,
    target = y,
    options = list(balance.formula = "~.", cost.online = "online", lambda = 0)
  )
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$optimizer, "cotDualFixedPointJointTrain")
  )

  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        niter = 1,
        cost.online = "online",
        balance.formula = "~X1",
        lambda = c(1, 10),
        debias = FALSE
      )
    )
  )
  testthat::expect_equal(
    class(cot$.__enclos_env__$private$optimizer),
    c("OTProblem", "R6")
  )
  osqpout <- testthat::capture_output(testthat::expect_no_warning(cot$solve()))
})

testthat::test_that("weights function works", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()

  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)
  a <- rep(1 / n, n)
  b <- rep(1 / m, m)
  delta <- 0.5

  #setup
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        verbose = FALSE,
        opt.direction = "dual",
        # torch.scheduler = NULL,
        niter = 1L,
        nboot = 10L
      )
    )
  )
  cot$solve()
  testthat::expect_silent(w <- cot$weights)
  # testthat::expect_true(all(w[1] == w))

  # run online version
  testthat::expect_silent(
    cot <- causalOT:::COT$new(
      source = x,
      target = y,
      options = list(
        verbose = FALSE,
        opt.direction = "dual",
        cost.online = "online",
        # torch.scheduler = NULL,
        niter = 1L,
        nboot = 10L
      )
    )
  )
  cot$solve()
  testthat::expect_silent(w <- cot$weights)
  # testthat::expect_true(all(w[1] == w))
})

testthat::test_that("cotOptions error checking works", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n * d), n, d)
  y <- matrix(stats::rnorm(m * d), m, d)
  a <- rep(1 / n, n)
  b <- rep(1 / m, m)
  delta <- 0.5

  testthat::expect_silent(opt <- cotOptions())

  testthat::expect_warning(
    opt <- cotOptions(
      opt.direction = "dual",
      torch.optimizer = torch::optim_lbfgs
    )
  )
  testthat::expect_true(inherits(opt$torch.optimizer, "optim_rmsprop"))
  opt <- cotOptions(opt.direction = "dual", dual.solver = "entropy")
  testthat::expect_true(inherits(opt$torch.optimizer, "optim_rmsprop"))
  testthat::expect_error(opt <- cotOptions(torch.optimizer = sum))
  testthat::expect_error(opt <- cotOptions(torch.scheduler = sum))
  opt <- cotOptions(debias = FALSE, opt.direction = "dual")
  testthat::expect_true(opt$opt.direction == "primal")
})

testthat::test_that("cotOptions warns about unused arguments", {
  testthat::skip_on_cran()
  causalOT:::torch_check()

  testthat::expect_warning(
    opt <- cotOptions(lambda = 1, dual.sover = "mirror"),
    "dual.sover"
  )
  testthat::expect_equal(opt$dual.solver, "fixed_point_joint")

  # arguments used by the torch optimizer are passed through quietly
  testthat::expect_silent(opt <- cotOptions(lambda = 1, lr = 0.01))
  testthat::expect_equal(opt$solver.options$lr, 0.01)

  # inner Sinkhorn controls reach the optimizer setup
  testthat::expect_silent(
    opt <- cotOptions(lambda = 1, ot_niter = 50L, ot_tol = 1e-5)
  )
  testthat::expect_equal(opt$ot_niter, 50L)
  testthat::expect_equal(opt$ot_tol, 1e-5)
})

testthat::test_that("cotOptions defaults to the bootstrap lambda selector", {
  testthat::skip_on_cran()
  causalOT:::torch_check()

  testthat::expect_equal(cotOptions()$lambda.selector, "bootstrap_legacy")
})
