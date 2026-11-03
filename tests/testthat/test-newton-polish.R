polish_fixture <- function() {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(42)
  x <- matrix(rnorm(24), 12, 2); y <- matrix(rnorm(30), 15, 2)
  cost <- function(x, y) pmax(outer(rowSums(x^2), rowSums(y^2), '+') / 2 - tcrossprod(x, y), 0)
  list(x = x, y = y, xy = cost(x, y), xx = cost(x, x), yy = cost(y, y),
       a = rep(1 / 12, 12), b = rep(1 / 15, 15), lambda = 2)
}
polish_run <- function(f, a = f$a, direction = 'exact', ...) {
  causalOT:::cot_newton_polish(a, f$b, f$xy, f$xx, f$lambda,
    direction = direction, control = causalOT:::cot_polish_control(...), C_yy = f$yy)
}

test_that('polish Hessian matches tangent differences, including zero candidates', {
  f <- polish_fixture()
  for (sparse in c(FALSE, TRUE)) {
    a <- f$a
    if (sparse) { a[10:12] <- 0; a <- a / sum(a) }
    ev <- causalOT:::cot_polish_eval(a, f$b, f$xy, f$xx, f$lambda)
    C <- if (sparse) 1:11 else 1:12
    H <- as.matrix(causalOT:::cot_polish_hessian(a, f$b, f$xy, f$xx, f$lambda, ev$pot, C))
    # At a boundary use tangent directions within S for central differences;
    # one-sided release derivatives are checked separately below.
    S <- which(a[C] > 0)
    Z <- matrix(0, length(C), length(S) - 1)
    Z[S, ] <- causalOT:::cot_polish_basis(length(S))
    step <- 1e-5
    fd <- sapply(seq_len(ncol(Z)), function(j) {
      d <- numeric(length(a)); d[C] <- Z[, j]
      (causalOT:::cot_polish_eval(a + step * d, f$b, f$xy, f$xx, f$lambda)$grad[C] -
       causalOT:::cot_polish_eval(a - step * d, f$b, f$xy, f$xx, f$lambda)$grad[C]) / (2 * step)
    })
    expect_lt(norm(crossprod(Z, H %*% Z - fd), 'F') / norm(crossprod(Z, fd), 'F'), 1e-5)
    if (sparse) {
      full <- as.matrix(causalOT:::cot_polish_hessian(a, f$b, f$xy, f$xx, f$lambda, ev$pot))
      expect_equal(H, full[C, C], tolerance = 1e-10)
      # Independent full formula, rather than comparing two calls of the block algorithm.
      K <- exp((outer(as.numeric(ev$pot$f_xy), as.numeric(ev$pot$g_yx), '+') - f$xy) / f$lambda)
      Kxx <- exp((outer(as.numeric(ev$pot$f_xx), as.numeric(ev$pot$f_xx), '+') - f$xx) / f$lambda)
      R <- K %*% (f$b * t(K)); A <- diag(a); I <- diag(length(a))
      ref <- f$lambda * (solve(I - R %*% A + outer(rep(1, length(a)), a), R) + solve(I + Kxx %*% A, Kxx))
      ref <- (ref + t(ref)) / 2
      expect_equal(H, ref[C, C], tolerance = 1e-10)
      d <- numeric(length(a)); d[10] <- 1; d[1] <- -1
      ep <- 1e-5
      g1 <- causalOT:::cot_polish_eval(a + ep*d, f$b, f$xy, f$xx, f$lambda)$grad
      g2 <- causalOT:::cot_polish_eval(a + 2*ep*d, f$b, f$xy, f$xx, f$lambda)$grad
      derivative <- (-3 * ev$grad + 4 * g1 - g2) / (2 * ep)
      P <- diag(length(a)) - 1 / length(a)
      expect_lt(max(abs(P %*% (derivative - full %*% d))), 1e-5)
    }
  }
})

test_that('Newton polish is monotone, sparse, and releases missing support', {
  f <- polish_fixture()
  ans <- polish_run(f, tol = 1e-11)
  expect_lte(ans$fw_gap, 1e-8 * max(f$xy))
  expect_lte(ans$loss, ans$loss_start)
  expect_equal(sum(ans$weights), 1, tolerance = 1e-14)
  expect_true(all(ans$weights >= 0))
  expect_lt(ans$nnz, length(f$a))
  damaged <- ans$weights; damaged[head(which(damaged > 0), 2)] <- 0
  restored <- polish_run(f, a = damaged / sum(damaged), tol = 1e-11)
  expect_identical(restored$weights > 1e-8, ans$weights > 1e-8)
  expect_equal(restored$loss, ans$loss, tolerance = 1e-10)
  again <- polish_run(f, a = ans$weights, tol = 1e-8)
  expect_true(again$status %in% c('certified', 'certified_resolution'))
  expect_lte(again$loss, ans$loss)
  expect_equal(again$weights, ans$weights, tolerance = 1e-7)
})

test_that('a failed direction returns the cold-evaluated phase-1 weights', {
  f <- polish_fixture()
  testthat::local_mocked_bindings(cot_polish_exact = function(H, g, a, zero_tol, state) g,
    .package = 'causalOT')
  ans <- polish_run(f, active_tol = 0)
  expect_identical(ans$status, 'no_improvement')
  expect_equal(ans$weights, f$a)
  expect_identical(ans$loss, ans$loss_start)
})

test_that('OSQP, L-BFGS and cap fallback agree with exact Newton', {
  f <- polish_fixture()
  exact <- polish_run(f, tol = 1e-11)
  qp <- polish_run(f, direction = 'osqp', tol = 1e-11)
  expect_equal(qp$loss, exact$loss, tolerance = 1e-10)
  expect_identical(qp$weights > 1e-8, exact$weights > 1e-8)
  limited <- polish_run(f, max_free = 2L, maxit = 80L, tol = 1e-9)
  expect_identical(limited$direction, 'lbfgs')
  expect_lte(limited$loss, limited$loss_start)
  expect_lt(abs(limited$loss - exact$loss), 1e-6 * abs(exact$loss))
})

test_that('online cost blocks match dense evaluation without KeOps', {
  f <- polish_fixture()
  xy <- causalOT:::costOnline$new(causalOT:::cot_polish_tensor(f$x), causalOT:::cot_polish_tensor(f$y), p = 2)
  xx <- causalOT:::costOnline$new(causalOT:::cot_polish_tensor(f$x), causalOT:::cot_polish_tensor(f$x), p = 2)
  dense <- polish_run(f, direction = 'lbfgs', maxit = 80L, tol = 1e-9)
  online <- causalOT:::cot_newton_polish(f$a, f$b, xy, xx, f$lambda, direction = 'lbfgs',
    control = causalOT:::cot_polish_control(block_size = 512L, maxit = 80L, tol = 1e-9), C_yy = f$yy)
  expect_equal(online$loss, dense$loss, tolerance = 1e-10)
  expect_equal(online$weights, dense$weights, tolerance = 1e-5)
  expect_lte(online$loss, online$loss_start)
})

test_that('mirror and two-lambda primal paths save polished weights and diagnostics', {
  skip_on_cran()
  f <- polish_fixture()
  for (solver in c('mirror', 'lbfgs')) {
    opts <- cotOptions(lambda = if (solver == 'mirror') 2 else c(4, 2),
      solver = solver, polish = 'exact', niter = 50L, nboot = 0,
      cost.online = 'tensorized', tol = 1e-10, ot_tol = 1e-10,
      ot_niter = 1000L, line_search_fn = 'strong_wolfe',
      dtype = torch::torch_double(), device = torch::torch_device('cpu'))
    cot <- causalOT:::COT$new(source = f$x, target = f$y, options = opts)
    expect_silent(cot$solve())
    optimizer <- cot$.__enclos_env__$private$optimizer
    info <- optimizer$info()$polish
    expect_length(info, if (solver == 'mirror') 1L else 2L)
    expect_true(all(is.finite(cot$weights)))
    expect_true(any(cot$weights == 0))
    for (l in names(info)) {
      expect_lte(info[[l]][[1]]$loss, info[[l]][[1]]$loss_start)
      expect_true(optimizer$info()$converged[[l]][[1]])
    }
    expect_equal(cot$weights, polish_run(f, tol = 1e-11)$weights, tolerance = 1e-6)
    if (solver == 'lbfgs') {
      saved <- optimizer$.__enclos_env__$private$weights_list
      nnz <- vapply(saved, function(e) {w <- e[[1L]]; sum(as.numeric(w[[ls(w)[1L]]]) > 0)}, numeric(1))
      expect_gt(nnz[['2']], nnz[['4']])
    }
  }
})

test_that('none preserves phase-1 results and unsupported refinement records a reason', {
  skip_on_cran()
  f <- polish_fixture()
  run <- function(polish, debias = TRUE, bypass = FALSE) {
    opts <- cotOptions(lambda = 2, solver = 'lbfgs', polish = polish,
      debias = debias, niter = 5, nboot = 0, cost.online = 'tensorized',
      dtype = torch::torch_double(), device = torch::torch_device('cpu'), line_search_fn = 'strong_wolfe')
    cot <- causalOT:::COT$new(source = f$x, target = f$y, options = opts)
    if (bypass) {
      priv <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
      unlockBinding('polish_result', priv)
      priv$polish_result <- function(res, tol, delta) res
      lockBinding('polish_result', priv)
    }
    cot$solve(); cot
  }
  first <- run('none'); second <- run('none')
  expect_identical(first$weights, second$weights)
  expect_identical(first$weights, run('none', bypass = TRUE)$weights)
  # Bypass the hook entirely to compare with the original phase-1 behavior.
  # The disabled hook is also compared directly with a no-op on a fresh fit below.
  unsupported <- run('auto', FALSE)
  info <- unsupported$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]
  expect_identical(info$status, 'skipped')
  expect_identical(info$reason, 'gap unavailable')
  expect_match(info$skip_reason, 'debiased')
})

test_that('polish options validate controls', {
  expect_identical(cotOptions()$polish, 'auto')
  expect_error(cotOptions(polish = 'invalid'), 'arg')
  expect_error(cotOptions(polish.control = list(maxit = 0)), 'maxit')
  expect_error(cotOptions(polish.control = list(unknown = 1)), 'Unknown')
})

test_that('polish has its own divergence-relative default tolerance', {
  f <- polish_fixture()
  opts <- cotOptions(tol = 1e-2)
  expect_identical(opts$polish.control$tol, 1e-8)
  ans <- polish_run(f)
  expect_equal(ans$divergence, ans$loss - ans$y_self, tolerance = 1e-14)
  expect_equal(ans$fw_threshold, 1e-8 * max(0, ans$divergence))
  if (ans$status == 'certified') expect_lte(ans$fw_gap, ans$fw_threshold)
  # Near-zero divergence uses arithmetic precision, never the solver tolerance.
  near <- causalOT:::cot_newton_polish(f$a, f$a, f$xx, f$xx, f$lambda, C_yy = f$xx)
  expect_true(near$status %in% c('certified', 'certified_resolution'))
  expect_lte(near$fw_gap, near$roundoff_floor)
  expect_identical(near$n_accepted, 0L)
  expect_equal(near$weights, f$a, tolerance = 1e-15)
  expect_error(causalOT:::cot_newton_polish(f$a, f$b, f$xy, f$xx, f$lambda), 'C_yy')
})

test_that('epsilon-active positive weights do not fill the Hessian or trigger the cap', {
  f <- polish_fixture()
  # Four nearly absent, distant atoms have strictly positive mirror-like mass.
  a <- c(rep((1 - 4e-9) / 8, 8), rep(1e-9, 4))
  f$x[9:12, ] <- f$x[9:12, ] + 20
  cost <- function(x,y) pmax(outer(rowSums(x^2),rowSums(y^2),'+')/2-tcrossprod(x,y),0)
  f$xy <- cost(f$x,f$y); f$xx <- cost(f$x,f$x)
  ans <- polish_run(f, a = a, max_free = 8L, maxit = 1L)
  expect_true(all(a > 0))
  expect_identical(ans$direction_history[[1]], 'exact')
  expect_lte(ans$free_history[[1]], 8L)
  expect_gte(ans$binding_history[[1]], 4L)
  expect_lte(ans$loss, ans$loss_start)
})

test_that('cap fallback can return to Newton and explicit L-BFGS remains explicit', {
  f <- polish_fixture()
  switched <- polish_run(f, max_free = 10L, active_tol = 0, maxit = 60L)
  expect_identical(switched$direction_history[[1]], 'lbfgs')
  expect_true('exact' %in% switched$direction_history)
  only <- polish_run(f, direction = 'lbfgs', max_free = 10L, maxit = 3L)
  expect_true(all(only$direction_history == 'lbfgs'))
})

test_that('the first L-BFGS step uses probability-mass units', {
  f <- polish_fixture()
  ans <- polish_run(f, direction = 'lbfgs', active_tol = 0, maxit = 1L)
  expect_lte(ans$first_step_norm, mean(f$a) * (1 + 1e-14))
  expect_lte(ans$line_search_evals, 4L)
  expect_lt(ans$loss, ans$loss_start)
})

test_that('default mirror polishing does not certify a diameter-sized FW gap', {
  skip_on_cran()
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(1); n <- 256L
  x <- matrix(rnorm(n * 3), n, 3)
  z <- rbinom(n, 1, plogis(.8*x[,1] - .5*x[,2] + .3*x[,3]))
  scale <- median(as.numeric(dist(x))^2)
  for (dtype in list(torch::torch_double(), torch::torch_float32())) for (mult in c(1, 5, 25)) {
    # Deliberately omit both tol and polish.control to exercise the defaults.
    opts <- cotOptions(lambda = mult * scale, solver = 'mirror', niter = 100L,
      nboot = 0, cost.online = 'tensorized', dtype = dtype, device = torch::torch_device('cpu'))
    cot <- causalOT:::COT$new(source = x[z == 0, , drop = FALSE], target = x, options = opts)
    cot$solve()
    optimizer <- cot$.__enclos_env__$private$optimizer
    info <- optimizer$info()$polish[[1]][[1]]
    expect_identical(opts$tol, 1e-4)
    expect_identical(opts$polish.control$tol, 1e-8)
    expect_lte(info$loss, info$loss_start)
    expect_equal(info$fw_threshold, 1e-8 * max(0, info$divergence))
    expect_identical(optimizer$.__enclos_env__$private$converged_run[[1]][[1]], causalOT:::cot_converged(info$convergence$gap, info$convergence$D, info$convergence$atol, info$convergence$tol))
    expect_true(optimizer$.__enclos_env__$private$converged_run[[1]][[1]])
    expect_lte(info$fw_gap, info$convergence$atol + opts$tol * max(0, info$divergence))
    expect_identical(info$nnz, c('1'=39L,'5'=34L,'25'=31L)[[as.character(mult)]])
    if (mult == 25 && dtype == torch::torch_float32()) expect_lt(info$iter, opts$polish.control$maxit)
    ots <- optimizer$.__enclos_env__$private$ot_objects
    ot <- ots[[ls(ots)[[1L]]]]
    returned <- cot$weights / sum(cot$weights)
    C <- optimizer$.__enclos_env__$private$polish_cache$costs
    actual <- causalOT:::cot_polish_eval(returned, ot$b, C$xy, C$xx, ot$penalty)
    returned_gap <- max(0, sum(returned * actual$grad) - min(actual$grad))
    expect_equal(info$fw_gap, returned_gap, tolerance = 1e-11)
    expect_lte(returned_gap, info$convergence$atol + info$convergence$tol * max(0, info$convergence$D))

  }
})

test_that('resolution convergence requires a stalled solve and a bounded gap', {
  ctl <- causalOT:::cot_polish_control()
  classify <- causalOT:::cot_polish_convergence
  expect_null(classify(2e-11, 1e-6, 0, 31, ctl, FALSE))
  expect_identical(classify(2e-11, 1e-6, 0, 31, ctl, TRUE), 'certified_resolution')
  expect_null(classify(1e-8, 1e-6, 0, 31, ctl, TRUE))
  expect_null(classify(1e-8, 0, 0, 31, ctl, TRUE))
  ctl$weight_eps <- 2^-23
  expect_identical(classify(1e-8, 0, 0, 31, ctl, TRUE), 'certified_resolution')
  expect_null(classify(1e-5, 0, 0, 31, ctl, TRUE))
})

test_that('epsilon Hessians keep current mass and use only the free operator', {
  f <- polish_fixture()
  a <- c(rep((1 - 4e-9) / 8, 8), rep(1e-9, 4))
  f$x[9:12, ] <- f$x[9:12, ] + 20
  cost <- function(x,y) pmax(outer(rowSums(x^2),rowSums(y^2),'+')/2-tcrossprod(x,y),0)
  f$xy <- cost(f$x,f$y); f$xx <- cost(f$x,f$x)
  original <- causalOT:::cot_polish_hessian
  calls <- list()
  testthat::local_mocked_bindings(cot_polish_hessian = function(a,b,C_xy,C_xx,lambda,pot,candidates,block_size,support) {
    calls[[length(calls)+1L]] <<- list(a = a, support = support, candidates = candidates)
    original(a,b,C_xy,C_xx,lambda,pot,candidates,block_size,support)
  }, .package = 'causalOT')
  ans <- polish_run(f, a = a, max_free = 8L, maxit = 1L)
  expect_length(calls, 1L)
  expect_true(all(calls[[1]]$a > 0))
  expect_equal(calls[[1]]$a, a)
  expect_lte(length(calls[[1]]$support), 8L)
  expect_identical(ans$free_history[[1]], length(calls[[1]]$candidates))
})

test_that('resolution-converged grids keep the requested path selector', {
  skip_on_cran()
  f <- polish_fixture()
  for (dtype in list(torch::torch_double(), torch::torch_float32())) {
    opts <- cotOptions(lambda = c(4, 2), solver = 'lbfgs', niter = 50L, nboot = 1L,
      cost.online = 'tensorized', lambda.selector = 'profile_symmetric_loo',
      line_search_fn = 'strong_wolfe', dtype = dtype, device = torch::torch_device('cpu'))
    cot <- causalOT:::COT$new(source = f$x, target = f$y, options = opts)
    cot$solve()
    priv <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
    expect_true(all(unlist(priv$converged_run)))
    result <- cot$grid_search()
    expect_identical(result$metric$lambda.selector, 'profile_symmetric_loo')
    expect_false(grepl('not converged', result$metric$lambda.note, fixed = TRUE))
  }
})

test_that('failed float32 installation restores the original parameterization exactly', {
  skip_on_cran()
  f <- polish_fixture()
  run <- function(polish) {
    opts <- cotOptions(lambda = 2, solver = 'lbfgs', polish = polish, niter = 5L,
      nboot = 0, cost.online = 'tensorized', dtype = torch::torch_float32(),
      device = torch::torch_device('cpu'), line_search_fn = 'strong_wolfe')
    cot <- causalOT:::COT$new(source = f$x, target = f$y, options = opts)
    cot$solve(); cot
  }
  baseline <- run('none')
  original <- causalOT:::cot_newton_polish
  testthat::local_mocked_bindings(cot_newton_polish = function(a,b,C_xy,C_xx,lambda,...) {
    result <- original(a,b,C_xy,C_xx,lambda,...)
    bad <- numeric(length(a)); bad[which.max(rowMeans(as.matrix(C_xy$get(seq_len(C_xy$nrow), seq_len(C_xy$ncol)))))] <- 1
    result$solver_weights <- result$weights <- bad
    result
  }, .package = 'causalOT')
  failed <- run('auto')
  expect_identical(failed$weights, baseline$weights)
  info <- failed$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]
  expect_identical(info$status, 'no_improvement')
  expect_match(info$rejection_reason, 'No improvement')
  expect_identical(info$convergence$phase, 'polish')
})
