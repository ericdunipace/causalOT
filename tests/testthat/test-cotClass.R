testthat::test_that("COT objects form", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  
  #cot defaults
  testthat::expect_silent(causalOT:::COT$new(source = x, target = y))
  
  #setup
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y))
  testthat::expect_silent(cot$weights <- c(1, rep(0,n-1)))
  testthat::expect_true(length(cot$.__enclos_env__$private$optimizer$penalty$lambda) == 8)
  testthat::expect_true(is.na(cot$.__enclos_env__$private$optimizer$penalty$delta))
  testthat::expect_equal(as_numeric(cot$weights), c(1, rep(0,n-1)))
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(balance.formula = "~X1+X2")))
  testthat::expect_true(inherits(cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$nn_holder, "cotDualBfOpt"))
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(lambda = 0,
                                                                   debias = TRUE)))
  
  # entropy not debiased
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(
                                                                   debias = FALSE)))
  testthat::expect_true(all(class(cot$.__enclos_env__$private$optimizer) !=
                                 "cotDualTrain"))
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(
                                                                   debias = TRUE,
                                                                   balance.formula = "~X1*X3")))
  testthat::expect_true("bf" %in% names(cot$.__enclos_env__$private$optimizer$.__enclos_env__$private))
  testthat::expect_true(inherits(cot$.__enclos_env__$private$optimizer,
                              "cotDualTrain"))
  
  # debias without torch optim should throw error
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(
                                                                   debias = TRUE,
                                                                   torch.optimizer = NULL)
                                                    ))
  
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(
                                                      debias = TRUE,
                                                      opt.direction = "primal")
  ))
  testthat::expect_true(inherits(cot$.__enclos_env__$private$optimizer, "OTProblem"))
  
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
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  lambda <- 10
  wlist <- list(a,a)
  
  #debias = FALSE
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(debias = FALSE, 
                                                                   niter = 4,
                                                                   opt.direction = "primal")))
  ot <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  testthat::expect_true(inherits(ot[[ls(ot)]]$C_xy, "costTensor"))
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$torch_optim, "optim_lbfgs")
  )
  testthat::expect_warning(cot$solve())
  cot <- causalOT:::COT$new(source = x, target = y,
                            options = list(debias = FALSE, 
                                           niter = 1,
                                           opt.direction = "primal",
                                           line_search_fn = "strong_wolfe"))
  testthat::expect_equivalent(cot$.__enclos_env__$private$torch_optim_args, list(line_search_fn = "strong_wolfe"))
  testthat::expect_silent(cot$solve())
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(niter = 1,
                                                                   balance.formula = "~X1",
                                                                   debias = FALSE)))
  
  osqpout <- testthat::capture_output(testthat::expect_warning(cot$solve()))
  
  testthat::expect_true(as_numeric((cot$.__enclos_env__$private$source$balance_functions$transpose(2,1)$matmul(cot$.__enclos_env__$private$source$weights) - cot$.__enclos_env__$private$source$balance_target)$abs()$max()$to(device = "cpu")$item()) < 1.5)
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(debias = FALSE,
                                                                   balance.formula = "~X1+X2", niter = 1, verbose = FALSE, delta = 0.1)))
  osqpout <- testthat::capture_output(testthat::expect_warning(cot$solve() ))
  testthat::expect_true(as_numeric((cot$.__enclos_env__$private$source$balance_functions$transpose(2,1)$matmul(cot$.__enclos_env__$private$source$weights) - cot$.__enclos_env__$private$source$balance_target)$abs()$max()$item()) < 0.9)
  
})

testthat::test_that("ent debiased works, online", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  causalOT:::rkeops_check()
  
  testthat::skip_if_not_installed("rkeops")
  
  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  lambda <- 10
  wlist <- list(a,a)
  
  #debias = FALSE
  mess <- testthat::capture_output(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(debias = FALSE, 
                                                                   niter = 4,
                                                                   opt.direction = "primal",
                                                                   cost.online = "online")))
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$torch_optim, "optim_lbfgs")
  )
  testthat::expect_warning(cot$solve())
  cot <- causalOT:::COT$new(source = x, target = y,
                            options = list(debias = FALSE, 
                                           niter = 1,
                                           opt.direction = "primal",
                                           line_search_fn = "strong_wolfe",
                                           cost.online = "online"))
  testthat::expect_equivalent(cot$.__enclos_env__$private$torch_optim_args, list(line_search_fn = "strong_wolfe"))
  testthat::expect_silent(cot$solve())
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(niter = 1,
                                                                   balance.formula = "~X1",
                                                                   debias = FALSE,
                                                                   cost.online = "online")))
  
  osqpout <- testthat::capture_output(testthat::expect_warning(cot$solve()))
  
  testthat::expect_true(as_numeric(max(abs(cot$.__enclos_env__$private$source$balance_functions$transpose(2,1)$matmul(cot$.__enclos_env__$private$source$weights) - cot$.__enclos_env__$private$source$balance_target))) < 1.5)
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(debias = FALSE,
                                                                   balance.formula = "~X1+X2", niter = 1, verbose = FALSE, delta = 0.1,
                                                                   cost.online = "online")))
  osqpout <- testthat::capture_output(testthat::expect_warning(cot$solve() ))
  testthat::expect_true(as_numeric(max(abs(cot$.__enclos_env__$private$source$balance_functions$transpose(2,1)$matmul(cot$.__enclos_env__$private$source$weights) - cot$.__enclos_env__$private$source$balance_target))) < 0.9)
  
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(debias = TRUE,
                                                                   balance.formula = "~X1+X2", niter = 1, verbose = FALSE, delta = 0.1,
                                                                   cost.online = "online")))
  
  osqpout <- testthat::capture_output(
    cot$solve() )

  testthat::expect_true(as_numeric(max(abs(cot$.__enclos_env__$private$source$balance_functions$transpose(2,1)$matmul(cot$.__enclos_env__$private$source$weights) - cot$.__enclos_env__$private$source$balance_target))) < 0.9)
  
})

testthat::test_that("grid_search function works",{
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  
  #setup
  torch::torch_manual_seed(123123)
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(verbose = FALSE, 
                                                                   opt.direction = "dual",
                                                                   # torch.scheduler = NULL,
                                                                   niter = 1L,
                                                                   nboot = 10L
                                                                   )))
  cot$solve()
  output <- cot$grid_search()
  testthat::expect_true(output$penalty[1] < 2704243.27850342 + 1 ) 
  testthat::expect_warning(cot$grid_search())
})

testthat::test_that("NNM works",{
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  
  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  lambda <- 10
  wlist <- list(a,a)
  
  #### tensor ####
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(lambda = 0,
                                                                   debias = FALSE, 
                                                                   niter = 4,
                                                                   opt.direction = "primal")))
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
  testthat::expect_true(all(names(res) %in% c("weight", "penalty", "metric","penalty.grid")))
  
  cot <- causalOT:::COT$new(source = x, target = y,
                            options = list(debias = TRUE, lambda = 0))
  ot <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  testthat::expect_true(isFALSE(ot[[ls(ot)]]$debias))
  
  cot <- causalOT:::COT$new(source = x, target = y,
                            options = list(balance.formula = "~.", lambda = 0))
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$optimizer, "cotDualTrain")
  )
  
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(niter = 1,
                                                                   balance.formula = "~X1",
                                                                   debias = FALSE)))
  testthat::expect_equal(class(cot$.__enclos_env__$private$optimizer), c("OTProblem", "R6"))
  osqpout <- testthat::capture_output(testthat::expect_warning(cot$solve()))
  
  #### keops ####
  causalOT:::rkeops_check()
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(lambda = 0,
                                                                   debias = FALSE, 
                                                                   niter = 4,
                                                                   cost.online = "online",
                                                                   opt.direction = "primal")))
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
  testthat::expect_true(all(names(res) %in% c("weight", "penalty", "metric","penalty.grid")))
  
  cot <- causalOT:::COT$new(source = x, target = y,
                            options = list(debias = TRUE, 
                                           cost.online = "online",lambda = 0))
  ot <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  testthat::expect_true(isFALSE(ot[[ls(ot)]]$debias))
  
  cot <- causalOT:::COT$new(source = x, target = y,
                            options = list(balance.formula = "~.", 
                                           cost.online = "online",lambda = 0))
  testthat::expect_true(
    inherits(cot$.__enclos_env__$private$optimizer, "cotDualTrain")
  )
  
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(niter = 1,
                                                                   cost.online = "online",
                                                                   balance.formula = "~X1",
                                                                   debias = FALSE)))
  testthat::expect_equal(class(cot$.__enclos_env__$private$optimizer), c("OTProblem", "R6"))
  osqpout <- testthat::capture_output(testthat::expect_warning(cot$solve()))
})

testthat::test_that("weights function works",{
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  
  #setup
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(verbose = FALSE, 
                                                                   opt.direction = "dual",
                                                                   # torch.scheduler = NULL,
                                                                   niter = 1L,
                                                                   nboot = 10L
                                                    )))
  cot$solve()
  testthat::expect_silent(w <- cot$weights)
  testthat::expect_true(all(w[1] == w))
  
  # run rkeops version
  causalOT:::rkeops_check() #skips if rkeops fails or is not installed
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(verbose = FALSE, 
                                                                   opt.direction = "dual",
                                                                   cost.online = "online",
                                                                   # torch.scheduler = NULL,
                                                                   niter = 1L,
                                                                   nboot = 10L
                                                    )))
  cot$solve()
  testthat::expect_silent(w <- cot$weights)
  testthat::expect_true(all(w[1] == w))
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
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  
  testthat::expect_silent(opt <- cotOptions())
  
  testthat::expect_warning(opt <- cotOptions(opt.direction = "dual", torch.optimizer = torch::optim_lbfgs))
  testthat::expect_true(inherits(opt$torch.optimizer, "optim_rmsprop"))
  testthat::expect_error(opt <- cotOptions(torch.optimizer = sum))
  testthat::expect_error(opt <- cotOptions(torch.scheduler = sum))
  opt <- cotOptions(debias = FALSE, opt.direction = "dual")
  testthat::expect_true(opt$opt.direction == "primal")
})
