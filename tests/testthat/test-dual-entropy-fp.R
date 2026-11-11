# Entropy fixed-point tests: scalar fits, sparse starts, and lambda paths.
defp_fixture <- function(n = 100L, seed = 1L, mult = 1, online = "tensorized") {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(seed)
  x <- matrix(rnorm(n * 2L), n)
  y <- matrix(rnorm(n * 2L, .5), n)
  C <- pmax(0, outer(rowSums(x^2), rowSums(y^2), "+") / 2 - tcrossprod(x, y))
  lambda <- mult * median(C)
  source <- causalOT:::Measure(x, adapt = "weights", dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  target <- causalOT:::Measure(y, adapt = "none", dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  fit <- causalOT:::cotDualEntropyFixedPointTrain$new(source, target)
  fit$setup_arguments(lambda = lambda, cost.online = online, ot_tol = 1e-8)
  list(fit = fit, private = fit$.__enclos_env__$private, lambda = lambda, x = x, y = y)
}

testthat::test_that("entropy start is the self potential in the correct form", {
  z <- defp_fixture()
  u <- z$private$nn_holder$gamma$detach()$clone()
  f0 <- z$fit$ot$potentials$f_xx / z$lambda
  testthat::expect_equal(as.numeric(u), as.numeric(f0 - f0$max()), tolerance = 1e-12)
  a <- z$fit$weights
  implied <- z$lambda * (u - a$log())
  z$fit$ot$a <- a
  z$fit$ot$sinkhorn_opt(niter = 10000L, tol = 1e-11)
  actual <- z$fit$ot$potentials$f_xx
  testthat::expect_equal(as.numeric(implied - implied$mean()),
    as.numeric(actual - actual$mean()), tolerance = 1e-7)
})

testthat::test_that("sparse large-lambda updates and underflowed residuals stay finite", {
  z <- defp_fixture(mult = 10)
  initial <- z$private$fp_checkpoint(z$fit$weights)$fw_gap
  result <- z$private$optimization_loop(500L, 1e-4)
  res <- z$private$fp_forward()
  testthat::expect_true(all(is.finite(as.numeric(res$d))))
  testthat::expect_true(all(is.finite(as.numeric(z$private$nn_holder$gamma))))
  testthat::expect_true(all(is.finite(as.numeric(res$w2))))
  testthat::expect_equal(sum(as.numeric(res$w2)), 1, tolerance = 1e-12)
  testthat::expect_lt(result$fw_gap, initial)
  torch::with_no_grad(z$private$nn_holder$gamma[1:80]$sub_(1e4))
  res <- z$private$fp_forward()
  testthat::expect_true(any(as.numeric(res$w2) == 0))
  testthat::expect_true(all(is.finite(as.numeric(res$d))))
})

testthat::test_that("public routing preserves honest capped-fit diagnostics and options", {
  z <- defp_fixture(n = 20L)
  options <- suppressWarnings(causalOT::cotOptions(lambda = z$lambda * c(1, 2),
    opt.direction = "dual", dual.solver = "fixed_point", niter = 2L,
    fixed_point_eta = .7, fixed_point_update = "momentum", anderson_m = 3L,
    fixed_point_lambda_order = "descending",
    fixed_point_warm_start_mix = .1,
    check_every = 1L, dtype = torch::torch_double(), device = torch::torch_device("cpu"), ot_tol = 1e-8))
  fit <- causalOT:::COT$new(z$x, z$y, options = options)
  fit$solve()
  opt <- fit$.__enclos_env__$private$optimizer
  p <- opt$.__enclos_env__$private
  testthat::expect_s3_class(opt, "cotDualEntropyFixedPointTrain")
  testthat::expect_false(any(unlist(p$converged_run)))
  testthat::expect_equal(unname(unlist(p$iterations_run)), c(2, 2))
  testthat::expect_equal(p$fixed_point_eta, .7)
  testthat::expect_equal(p$fixed_point_update, "momentum")
  testthat::expect_equal(p$anderson_m, 3L)
  testthat::expect_equal(p$check_every, 1L)
  testthat::expect_identical(p$lambda_order, "descending")
  testthat::expect_equal(p$warm_start_mix, .1)
  testthat::expect_equal(p$penalty_list$lambda, z$lambda * c(2, 1))
  testthat::expect_false(inherits(p$opt, "optim_optimizer"))
  testthat::expect_true(all(unlist(p$reason_run) == "niter"))
})

testthat::test_that("returned state has the smallest checked divergence", {
  z <- defp_fixture(seed = 4L)
  # Force the final proposal uphill so returning the last iterate would fail.
  original <- z$private$fp_forward
  forwards <- 0L
  unlockBinding("fp_forward", z$private)
  z$private$fp_forward <- function() {
    forwards <<- forwards + 1L
    if (forwards == 101L) z$private$nn_holder$gamma[1]$add_(100)
    original()
  }
  lockBinding("fp_forward", z$private)
  result <- z$private$optimization_loop(100L, 1e-8)
  checked <- vapply(z$private$check_history, `[[`, numeric(1), "divergence")
  testthat::expect_lt(which.min(checked), length(checked))
  testthat::expect_gt(z$private$memory_resets, 0L)
  testthat::expect_equal(result$divergence, min(checked), tolerance = 1e-7)
  testthat::expect_equal(as.numeric(z$fit$weights), as.numeric(z$private$fp_forward()$w2), tolerance = 1e-12)
  testthat::expect_identical(result$converged,
    causalOT:::cot_converged(result$fw_gap, result$divergence, 1e-8, 1e-8))
})

testthat::test_that("online and tensorized entropy updates agree", {
  dense <- defp_fixture(n = 20L, mult = .1)
  online <- defp_fixture(n = 20L, mult = .1, online = "online")
  a <- dense$private$optimization_loop(100L, 1e-6)
  b <- online$private$optimization_loop(100L, 1e-6)
  testthat::expect_equal(as.numeric(dense$fit$weights), as.numeric(online$fit$weights), tolerance = 1e-5)
  testthat::expect_equal(a$divergence, b$divergence, tolerance = 1e-6)
  torch::with_no_grad(online$private$nn_holder$gamma[1:15]$sub_(1e4))
  testthat::expect_true(all(is.finite(as.numeric(online$private$fp_forward()$d))))
})

testthat::test_that("collapse case improves on joint solver and matches polished reference", {
  z <- defp_fixture()
  options <- suppressWarnings(causalOT::cotOptions(lambda = z$lambda,
    opt.direction = "dual", dual.solver = "fixed_point_joint", niter = 2000L,
    ot_tol = 1e-8, tol = 1e-4, dtype = torch::torch_double(), device = torch::torch_device("cpu")))
  legacy <- causalOT:::COT$new(z$x, z$y, options = options)
  legacy$solve()
  testthat::expect_gt(max(legacy$weights), .5)
  reference <- causalOT:::COT$new(z$x, z$y,
    options = causalOT::cotOptions(lambda = z$lambda, solver = "mirror",
      polish = "auto", niter = 2000L, ot_tol = 1e-8, tol = 1e-6,
      dtype = torch::torch_double(), device = torch::torch_device("cpu")))
  reference$solve()
  ref_ot <- reference$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
  ref_ot <- ref_ot[[ls(ref_ot)[1L]]]
  ref_ot$a <- torch::torch_tensor(reference$weights, dtype = torch::torch_double())
  ref_ot$sinkhorn_opt(niter = 10000L, tol = 1e-11)
  S_star <- causalOT:::oop_loss_select(ref_ot)$item()
  testthat::expect_true(causalOT:::cot_converged(causalOT:::cot_ccp_kkt(ref_ot)$fw_gap, S_star, 1e-8, 1e-6))
  result <- z$private$optimization_loop(10000L, 1e-4)
  testthat::expect_lt(max(as.numeric(z$fit$weights)), .5)
  testthat::expect_lte(result$divergence, S_star + 1e-4)
  z$fit$ot$sinkhorn_opt(niter = 10000L, tol = 1e-11)
  testthat::expect_identical(result$converged,
    causalOT:::cot_converged(causalOT:::cot_ccp_kkt(z$fit$ot)$fw_gap,
      causalOT:::oop_loss_select(z$fit$ot)$item(), 1e-8, 1e-4))
})

testthat::test_that("balance functions are rejected clearly", {
  source <- causalOT:::Measure(matrix(1:12, 6), adapt = "weights",
    target.values = c(2, 8), dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  target <- causalOT:::Measure(matrix(2:13, 6), adapt = "none", dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  fit <- causalOT:::cotDualEntropyFixedPointTrain$new(source, target)
  testthat::expect_error(fit$setup_arguments(lambda = 1), "does not support balance functions")
})

testthat::test_that("Anderson agrees with an independent base-R update", {
  z <- defp_fixture(n = 20L, mult = .1)
  u <- as.numeric(z$private$nn_holder$gamma$detach())
  xy <- -as.matrix(z$fit$ot$C_xy$data) / z$lambda
  xx <- -as.matrix(z$fit$ot$C_xx$data) / z$lambda
  lb <- as.numeric(z$private$b_log)
  lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
  U <- R <- NULL
  for (i in 1:10) {
    g <- lb - apply(xy + u, 2, lse)
    s <- apply(sweep(xx, 2, u, "+"), 1, lse)
    Z <- lse(u + s)
    d <- apply(sweep(xy, 2, g, "+"), 1, lse) - s + Z
    w <- exp(u + s - Z)
    U <- cbind(U, u + .5 * d)
    R <- cbind(R, sqrt(w) * .5 * d)
    if (ncol(U) > 5) { U <- U[, -1]; R <- R[, -1] }
    G <- crossprod(R)
    alpha <- solve(G + 1e-10 * max(diag(G)) * diag(ncol(G)), rep(1, ncol(G)))
    u <- drop(U %*% (alpha / sum(alpha)))
    u <- u - max(u)
  }
  z$private$optimization_loop(10L, 0)
  testthat::expect_equal(as.numeric(z$private$nn_holder$gamma$detach()), u, tolerance = 1e-6)
})

testthat::test_that("entropy forwards retain legacy weights and gradients", {
  z <- defp_fixture(n = 20L)
  h <- z$private$nn_holder
  torch::with_no_grad(h$gamma$copy_(torch::torch_linspace(-4, 2, 20, dtype = torch::torch_double())))
  res <- z$private$fp_forward()
  w1 <- h$calc_w1(h$gamma$detach(), z$fit$ot$C_xy$data,
    z$private$a_log, z$private$b_log, torch::jit_scalar(z$lambda), torch::jit_scalar(20L))
  w2 <- z$fit$weights
  testthat::expect_equal(as.numeric(res$w1), as.numeric(w1), tolerance = 1e-12)
  testthat::expect_equal(as.numeric(res$w2), as.numeric(w2), tolerance = 1e-12)
  testthat::expect_equal(as.numeric(res$d), as.numeric(w1$log() - w2$log()), tolerance = 1e-12)
  res$loss$backward()
  expected <- -(w1 - w2)
  testthat::expect_equal(as.numeric(h$gamma$grad), as.numeric(expected - expected$mean()), tolerance = 1e-12)
})

testthat::test_that("stall stops stay distinct from certified public convergence", {
  stalled <- defp_fixture(n = 20L)
  stalled$private$fixed_point_eta <- 1e-14
  result <- stalled$private$optimization_loop(100L, 1e-4)
  testthat::expect_false(result$converged)
  testthat::expect_true(result$stopped_by_rule)
  testthat::expect_identical(result$reason, "stalled")
  testthat::expect_equal(result$iter, 1L)
  z <- defp_fixture(mult = .1)
  z$fit$solve(niter = 2000L, tol = 1e-4)
  testthat::expect_true(all(unlist(z$private$converged_run)))
  z$fit$ot$sinkhorn_opt(niter = 10000L, tol = 1e-11)
  testthat::expect_true(causalOT:::cot_converged(
    causalOT:::cot_ccp_kkt(z$fit$ot)$fw_gap,
    causalOT:::oop_loss_select(z$fit$ot)$item(), 1e-8, 1e-4))
})

testthat::test_that("refitting a certified lambda stops at iteration zero without mixing", {
  z <- defp_fixture(mult = .1)
  z$fit$solve(niter = 2000L, tol = 1e-4)
  testthat::expect_true(all(unlist(z$private$converged_run)))
  before <- as.numeric(z$fit$weights)
  z$fit$solve(niter = 2000L, tol = 1e-4)
  testthat::expect_equal(unname(unlist(z$private$iterations_run)), 0)
  testthat::expect_true(all(unlist(z$private$converged_run)))
  testthat::expect_identical(as.numeric(z$fit$weights), before)
  testthat::expect_equal(z$private$check_history[[1L]]$iter, 0L)
  testthat::expect_equal(z$private$memory_resets, 0L)
})

testthat::test_that("a sparse certified warm start is checked and returned exactly", {
  x <- matrix(seq_len(8), ncol = 1)
  a <- c(0, .1, 0, .2, 0, .3, 0, .4)
  source <- causalOT:::Measure(x, weights = a, adapt = "weights", dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  target <- causalOT:::Measure(x, weights = a, adapt = "none", dtype = torch::torch_double(), device = torch::torch_device("cpu"))
  fit <- causalOT:::cotDualEntropyFixedPointTrain$new(source, target)
  fit$setup_arguments(lambda = c(1, 2), cost.online = "tensorized", ot_tol = 1e-10)
  p <- fit$.__enclos_env__$private
  before <- as.numeric(p$pending_weights)
  # S = 0 here, so the gap threshold is ot_tol alone, and a check whose Sinkhorn
  # starts from the cold-restart (uniform) potentials misses it by Sinkhorn error.
  # On a real path the previous fit leaves potentials at the incoming weights.
  fit$ot$a <- torch::torch_tensor(a, dtype = torch::torch_double())
  fit$ot$sinkhorn_opt(niter = 1e5L, tol = 1e-13)
  # Reconstruction at either lambda would mix away the exact zeros.
  unlockBinding("fp_build_potential", p)
  p$fp_build_potential <- function(...) stop("unexpected mixing of a certified start")
  lockBinding("fp_build_potential", p)
  fit$solve(niter = 100L, tol = 1e-4)
  testthat::expect_true(all(unlist(p$converged_run)))
  testthat::expect_equal(unname(unlist(p$iterations_run)), c(0, 0))
  testthat::expect_equal(as.numeric(fit$weights), a, tolerance = 1e-12)
  testthat::expect_identical(as.numeric(fit$weights), before)
  testthat::expect_identical(as.numeric(fit$weights) == 0, a == 0)
})

testthat::test_that("the weight setter reconstructs a finite potential from exact zeros", {
  for (online in c("tensorized", "online")) {
    z <- defp_fixture(n = 20L, mult = 10, online = online)
    a <- c(rep(0, 18), .2, .8)
    # A small mix keeps the exact zeros nearly in place, which is the hard case.
    z$private$warm_start_mix <- 1e-3
    z$private$parameters_get_set(list(a))
    mixed <- (1 - z$private$warm_start_mix) * a + z$private$warm_start_mix / length(a)
    testthat::expect_true(all(is.finite(as.numeric(z$private$nn_holder$gamma))))
    testthat::expect_equal(as.numeric(z$fit$weights), mixed, tolerance = 1e-7)
    w <- z$private$parameters_get_set(clone = TRUE)
    testthat::expect_equal(as.numeric(w[[ls(w)[1L]]]), mixed, tolerance = 1e-7)
    z$private$set_lambda(z$lambda / 10)
    testthat::expect_identical(as.numeric(z$private$pending_weights), a)
    result <- z$private$optimization_loop(2L, 1e-6)
    testthat::expect_false(result$converged)
    testthat::expect_true(all(is.finite(as.numeric(z$private$nn_holder$gamma))))
    testthat::expect_true(all(is.finite(as.numeric(z$private$fp_forward()$d))))
  }
})

testthat::test_that("an iteration error restores phase_start before the next lambda", {
  z <- defp_fixture(n = 20L)
  p <- z$private
  p$penalty_list$lambda <- p$penalty_input_list$lambda <- z$lambda * c(1, 2, 3)
  starts <- list()
  checkpoint <- p$fp_checkpoint
  unlockBinding("fp_checkpoint", p)
  p$fp_checkpoint <- function(w) {
    key <- as.character(z$fit$ot$penalty)
    if (is.null(starts[[key]])) starts[[key]] <<- as.numeric(w)
    checkpoint(w)
  }
  lockBinding("fp_checkpoint", p)
  forward <- p$fp_forward
  failed <- FALSE
  unlockBinding("fp_forward", p)
  p$fp_forward <- function() {
    if (!failed && abs(z$fit$ot$penalty - 2 * z$lambda) < 1e-10) {
      failed <<- TRUE
      p$nn_holder$gamma$fill_(NaN)
      stop("forced iteration failure")
    }
    forward()
  }
  lockBinding("fp_forward", p)
  z$fit$solve(niter = 2L, tol = 1e-6, fixed_point_lambda_order = "ascending")
  keys <- as.character(z$lambda * c(1, 2, 3))
  testthat::expect_match(p$reason_run[[keys[2]]][[1]], "error: forced iteration failure")
  failed_weights <- p$weights_list[[keys[2]]][[1]]
  testthat::expect_true(all(is.na(as.numeric(failed_weights[[ls(failed_weights)[1]]]))))
  testthat::expect_identical(starts[[keys[3]]], starts[[keys[2]]])
  testthat::expect_identical(p$reason_run[[keys[3]]][[1]], "niter")
  testthat::expect_true(all(is.finite(as.numeric(z$fit$weights))))
  testthat::expect_true(all(is.finite(as.numeric(p$nn_holder$gamma))))
})

testthat::test_that("native three-lambda paths agree with cold fits at the default and a partial mix", {
  testthat::expect_equal(causalOT:::cotDualEntropyFixedPointTrain$private_fields$warm_start_mix, 1)
  x <- matrix(c(0, 1, 2), ncol = 1)
  # A 1e-3 mix lost the middle point on this support-changing path.
  y <- matrix(c(.2, 1.8), ncol = 1)
  lambdas <- c(.05, .1, .2)
  make <- function(lambda) causalOT:::COT$new(x, y,
    options = suppressWarnings(causalOT::cotOptions(lambda = lambda,
      opt.direction = "dual", dual.solver = "fixed_point", niter = 2000L,
      ot_tol = 1e-10, tol = 1e-6, dtype = torch::torch_double(), device = torch::torch_device("cpu"))))
  cold <- lapply(lambdas, function(l) {
    f <- make(l); f$solve()
    p <- f$.__enclos_env__$private$optimizer$.__enclos_env__$private
    testthat::expect_true(all(unlist(p$converged_run)))
    f$weights
  })
  for (mix in c(1, .5)) for (order in c("ascending", "descending")) {
    f <- make(lambdas)
    f$solve(fixed_point_lambda_order = order, fixed_point_warm_start_mix = mix)
    p <- f$.__enclos_env__$private$optimizer$.__enclos_env__$private
    testthat::expect_equal(p$warm_start_mix, mix)
    testthat::expect_equal(p$penalty_list$lambda, sort(lambdas, decreasing = order == "descending"))
    testthat::expect_true(all(unlist(p$converged_run)))
    for (i in seq_along(lambdas)) {
      w <- p$weights_list[[as.character(lambdas[i])]][[1]]
      testthat::expect_equal(as.numeric(w[[ls(w)[1]]]), cold[[i]], tolerance = 1e-5)
    }
  }
})

testthat::test_that("a full cold restart rebuilds uniform weights after the incoming check", {
  z <- defp_fixture(n = 20L)
  a <- c(rep(0, 18), .2, .8)
  z$private$parameters_get_set(list(a))
  incoming <- rebuilt <- NULL
  checkpoint <- z$private$fp_checkpoint
  build <- z$private$fp_build_potential
  unlockBinding("fp_checkpoint", z$private)
  z$private$fp_checkpoint <- function(w) {
    if (is.null(incoming)) incoming <<- as.numeric(w)
    checkpoint(w)
  }
  lockBinding("fp_checkpoint", z$private)
  unlockBinding("fp_build_potential", z$private)
  z$private$fp_build_potential <- function(w) {
    build(w)
    rebuilt <<- as.numeric(z$fit$weights)
  }
  lockBinding("fp_build_potential", z$private)
  z$fit$solve(niter = 2L, tol = 1e-6, fixed_point_warm_start_mix = 1)
  testthat::expect_identical(incoming, a)
  testthat::expect_equal(rebuilt, rep(1 / 20, 20), tolerance = 1e-7)
  testthat::expect_equal(z$private$warm_start_mix, 1)
  # Rebuilding from unrelated cached potentials gives exactly the same start.
  z$private$fp_build_potential(z$fit$weights)
  cold_u <- as.numeric(z$private$nn_holder$gamma)
  z$private$fp_checkpoint(torch::torch_tensor(a, dtype = torch::torch_double()))
  z$private$fp_build_potential(torch::torch_tensor(a, dtype = torch::torch_double()))
  testthat::expect_identical(as.numeric(z$private$nn_holder$gamma), cold_u)
})

testthat::test_that("warm-start mixing requires a positive fraction at most one", {
  z <- defp_fixture(n = 20L)
  for (mix in list(0, -.1, 1.1, NA_real_, Inf, numeric(), c(.1, .5))) {
    testthat::expect_error(z$fit$solve(fixed_point_warm_start_mix = mix),
      "fixed_point_warm_start_mix must be in")
  }
})
