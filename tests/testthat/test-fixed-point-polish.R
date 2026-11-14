testthat::test_that("fixed-point polish certifies and preserves exact sparse weights", {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(42)
  x <- matrix(rnorm(60), 30)
  y <- matrix(rnorm(60, .5), 30)
  fits <- lapply(c("none", "auto"), function(polish) {
    opts <- causalOT::cotOptions(lambda = 2, polish = polish, ot_tol = 1e-8,
      dtype = torch::torch_double(), device = torch::torch_device("cpu"),
      cost.online = "tensorized")
    fit <- causalOT:::COT$new(x, y, options = opts)
    fit$solve()
    fit
  })
  opt <- lapply(fits, function(fit) fit$.__enclos_env__$private$optimizer)
  info <- lapply(opt, function(o) o$info())
  testthat::expect_identical(info[[1]]$polish[[1]][[1]]$status, "none")
  testthat::expect_true(info[[2]]$polish[[1]][[1]]$status %in% c("refined", "already_optimal"))
  testthat::expect_true(info[[2]]$converged[[1]][[1]])
  p <- lapply(opt, function(o) o$.__enclos_env__$private)
  testthat::expect_lte(p[[2]]$polish_gap[[1]][[1]][["gap"]],
    p[[1]]$polish_gap[[1]][[1]][["gap"]] + 1e-10)
  polished <- opt[[2]]$polished_weights
  testthat::expect_true(any(polished == 0))
  testthat::expect_equal(sum(polished), 1, tolerance = 1e-12)
  testthat::expect_equal(fits[[2]]$weights, polished, tolerance = 1e-14)
  fits[[2]]$grid_search()
  testthat::expect_identical(fits[[2]]$weights, polished)
  testthat::expect_equal(as.numeric(opt[[2]]$weights), polished, tolerance = 1e-14)
})

testthat::test_that("fixed-point polish respects the free-set limit", {
  causalOT:::torch_check()
  set.seed(19)
  x <- matrix(rnorm(60), 30)
  y <- matrix(rnorm(60, .5), 30)
  ctl <- utils::modifyList(causalOT:::cot_polish_ctl, list(max_free = 1L))
  testthat::local_mocked_bindings(cot_polish_ctl = ctl, .package = "causalOT")
  opts <- causalOT::cotOptions(lambda = 2, niter = 1L,
    dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  fit <- causalOT:::COT$new(x, y, options = opts)
  fit$solve()
  info <- fit$.__enclos_env__$private$optimizer$info()
  testthat::expect_identical(info$polish[[1]][[1]]$status,
    "skipped: free set larger than 1")
})

testthat::test_that("legacy dual solvers remain excluded from polish", {
  causalOT:::torch_check()
  set.seed(29)
  x <- matrix(rnorm(24), 12)
  y <- matrix(rnorm(24, .5), 12)
  for (solver in c("joint", "entropy", "fixed_point_joint")) {
    opts <- suppressWarnings(causalOT::cotOptions(lambda = 2, niter = 1L,
      opt.direction = "dual", dual.solver = solver, polish = "auto",
      dtype = torch::torch_double(), device = torch::torch_device("cpu")))
    fit <- causalOT:::COT$new(x, y, options = opts)
    fit$solve()
    info <- fit$.__enclos_env__$private$optimizer$info()
    testthat::expect_match(info$polish[[1]][[1]]$status,
      "^skipped: dual solvers", info = solver)
  }
})

testthat::test_that("fixed-point returns float64 polish after float32 fitting", {
  causalOT:::torch_check()
  set.seed(42)
  x <- matrix(rnorm(60), 30)
  y <- matrix(rnorm(60, .5), 30)
  opts <- causalOT::cotOptions(lambda = 2, niter = 50L,
    dtype = torch::torch_float(), device = torch::torch_device("cpu"),
    cost.online = "tensorized")
  fit <- causalOT:::COT$new(x, y, options = opts)
  fit$solve()
  opt <- fit$.__enclos_env__$private$optimizer
  polished <- opt$polished_weights
  testthat::expect_true(any(polished == 0))
  testthat::expect_type(polished, "double")
  testthat::expect_equal(sum(polished), 1, tolerance = 1e-12)
  fit$grid_search()
  testthat::expect_identical(fit$weights, polished)
})

testthat::test_that("a finished polish does not rebuild the fixed-point potential", {
  causalOT:::torch_check()
  set.seed(42)
  x <- matrix(rnorm(60), 30)
  y <- matrix(rnorm(60, .5), 30)
  opts <- causalOT::cotOptions(lambda = 2, niter = 50L,
    dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  fit <- causalOT:::COT$new(x, y, options = opts)
  p <- fit$.__enclos_env__$private$optimizer$.__enclos_env__$private
  original_build <- p$fp_build_potential
  original_polish <- causalOT:::cot_newton_polish
  polish_done <- FALSE
  builds_after_polish <- 0L
  unlockBinding("fp_build_potential", p)
  p$fp_build_potential <- function(w) {
    if (polish_done) builds_after_polish <<- builds_after_polish + 1L
    original_build(w)
  }
  lockBinding("fp_build_potential", p)
  testthat::local_mocked_bindings(cot_newton_polish = function(...) {
    result <- original_polish(...)
    polish_done <<- TRUE
    result
  }, .package = "causalOT")
  fit$solve()
  fit$grid_search()
  testthat::expect_true(polish_done)
  testthat::expect_identical(builds_after_polish, 0L)
  testthat::expect_equal(sum(fit$weights), 1, tolerance = 1e-12)
})
