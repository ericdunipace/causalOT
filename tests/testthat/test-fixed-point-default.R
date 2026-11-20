testthat::test_that("fixed-point routing and compatibility defaults are explicit", {
  default <- causalOT::cotOptions()
  testthat::expect_identical(default$solver, "fixed_point")
  testthat::expect_identical(default$opt.direction, "dual")
  testthat::expect_identical(default$dual.solver, "fixed_point")
  testthat::expect_identical(default$niter, 500L)
  explicit <- causalOT::cotOptions(solver = "fixed_point")
  testthat::expect_identical(explicit$opt.direction, "dual")
  testthat::expect_identical(explicit$dual.solver, "fixed_point")
  for (args in list(list(debias = FALSE), list(delta = .1),
                    list(balance.formula = ~ x))) {
    testthat::expect_identical(do.call(causalOT::cotOptions, args)$solver, "lbfgs")
    testthat::expect_error(do.call(causalOT::cotOptions,
      c(list(solver = "fixed_point"), args)), "requires debias|balance constraints")
  }
  legacy <- suppressWarnings(causalOT::cotOptions(opt.direction = "dual", dual.solver = "joint"))
  testthat::expect_identical(legacy$dual.solver, "joint")
  testthat::expect_identical(legacy$opt.direction, "dual")
})

testthat::test_that("default cap is quiet and returns valid fixed-point weights", {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(172)
  x <- matrix(rnorm(40), 20)
  y <- matrix(rnorm(40, .5), 20)
  options <- causalOT::cotOptions(lambda = 1, polish = "none", tol = 0, ot_tol = 1e-14,
    dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  fit <- causalOT:::COT$new(x, y, options = options)
  testthat::expect_no_warning(fit$solve())
  p <- fit$.__enclos_env__$private$optimizer$.__enclos_env__$private
  testthat::expect_lte(max(unlist(p$iterations_run)), 500)
  testthat::expect_true(all(is.finite(fit$weights) & fit$weights >= 0))
  testthat::expect_equal(sum(fit$weights), 1, tolerance = 1e-10)
  testthat::expect_true(any(unlist(p$reason_run) == "niter"))
  testthat::expect_no_warning(fit$grid_search())
})

testthat::test_that("explicit 2000 iterations retain legacy fixed-point routing", {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(173)
  x <- matrix(rnorm(24), 12)
  y <- matrix(rnorm(24, .5), 12)
  args <- list(lambda = 1, niter = 2000L, tol = 1e-4, ot_tol = 1e-6,
    dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  current <- causalOT:::COT$new(x, y, options = do.call(causalOT::cotOptions,
    c(args, list(solver = "fixed_point"))))
  legacy <- causalOT:::COT$new(x, y, options = suppressWarnings(do.call(causalOT::cotOptions,
    c(args, list(opt.direction = "dual", dual.solver = "fixed_point")))))
  current$solve()
  legacy$solve()
  testthat::expect_equal(current$weights, legacy$weights, tolerance = 1e-10)
})
