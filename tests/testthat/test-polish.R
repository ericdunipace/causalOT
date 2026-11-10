polish_fixture <- function() {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(42)
  x <- matrix(rnorm(24), 12, 2)
  y <- matrix(rnorm(30), 15, 2)
  cost <- function(x, y) {
    pmax(outer(rowSums(x^2), rowSums(y^2), '+') / 2 - tcrossprod(x, y), 0)
  }
  f <- list(
    x = x,
    y = y,
    xy = cost(x, y),
    xx = cost(x, x),
    yy = cost(y, y),
    a = rep(1 / 12, 12),
    b = rep(1 / 15, 15),
    lambda = 2
  )
  f
}
polish_twin <- function(f, a = f$a, mode = 'tensorized') {
  causalOT:::OT$new(
    x = f$x,
    y = f$y,
    a = a,
    b = f$b,
    penalty = f$lambda,
    p = 2,
    debias = TRUE,
    tensorized = mode,
    device = torch::torch_device('cpu'),
    dtype = torch::torch_double()
  )
}
polish_run <- function(f, a = f$a, mode = 'tensorized', ...) {
  ctl <- utils::modifyList(causalOT:::cot_polish_ctl, list(...))
  causalOT:::cot_newton_polish(polish_twin(f, a, mode), ctl = ctl)
}
polish_loss <- function(f, a) {
  causalOT:::cot_polish_eval(polish_twin(f, a), a)$loss
}

fit_cot <- function(
  polish = 'auto',
  lambda = 2,
  dtype = torch::torch_double(),
  debias = TRUE,
  niter = 2L,
  solver = 'lbfgs',
  source = NULL,
  target = NULL
) {
  f <- polish_fixture()
  opts <- cotOptions(
    lambda = lambda,
    solver = solver,
    line_search_fn = 'strong_wolfe',
    polish = polish,
    debias = debias,
    niter = niter,
    nboot = 0,
    cost.online = 'tensorized',
    dtype = dtype,
    device = torch::torch_device('cpu')
  )
  causalOT:::COT$new(
    source = if (is.null(source)) f$x else source,
    target = if (is.null(target)) f$y else target,
    options = opts
  )
}
fit_optimizer <- function(cot) cot$.__enclos_env__$private$optimizer
fit_private <- function(cot) fit_optimizer(cot)$.__enclos_env__$private
replace_private <- function(p, name, value) {
  unlockBinding(name, p)
  assign(name, value, envir = p)
  lockBinding(name, p)
}
fake_phase1 <- function(...) {
  out <- list(
    loss = 1,
    iter = 1L,
    converged = TRUE,
    stopped_by_rule = TRUE,
    fw_gap = 0,
    divergence = 1
  )
  utils::modifyList(out, list(...))
}

# Polish algorithm ------------------------------------------------------------

test_that('one tolerance test handles absolute and relative accuracy and finite inputs', {
  check <- causalOT:::cot_converged
  expect_true(check(1e-3, 0, 1e-3, 1e-4))
  expect_false(check(1.01e-3, 0, 1e-3, 1e-4))
  expect_true(check(.101, 1000, .001, 1e-4))
  expect_false(check(.102, 1000, .001, 1e-4))
  expect_true(check(.001, -10, .001, 1e-4))
  expect_false(check(.002, -10, .001, 1e-4))
  for (i in 1:4) {
    for (bad in c(NA_real_, NaN, Inf, -Inf)) {
      args <- list(.001, 1, .001, 1e-4)
      args[[i]] <- bad
      expect_false(do.call(check, args))
    }
  }
  expect_false(check(NULL, 1, .001, 1e-4))
})

test_that('weight Hessian matches finite differences and the full formula', {
  f <- polish_fixture()
  eval_grad <- function(a) causalOT:::cot_polish_eval(polish_twin(f, a), a)$grad
  for (sparse in c(FALSE, TRUE)) {
    a <- f$a
    if (sparse) {
      a[10:12] <- 0
      a <- a / sum(a)
    }
    twin <- polish_twin(f, a)
    causalOT:::cot_polish_eval(twin, a)
    C <- if (sparse) 1:11 else 1:12
    H <- as.matrix(twin$hessian_weights(C))
    S <- which(a[C] > 0)
    Z <- matrix(0, length(C), length(S) - 1)
    Z[S, ] <- stats::contr.sum(length(S))
    step <- 1e-5
    fd <- sapply(seq_len(ncol(Z)), function(j) {
      d <- numeric(length(a))
      d[C] <- Z[, j]
      (eval_grad(a + step * d)[C] - eval_grad(a - step * d)[C]) / (2 * step)
    })
    expect_lt(
      norm(crossprod(Z, H %*% Z - fd), 'F') / norm(crossprod(Z, fd), 'F'),
      1e-5
    )
    if (sparse) {
      p <- twin$.__enclos_env__$private
      K <- as.matrix(p$kernel('xy'))
      Kxx <- as.matrix(p$kernel('xx'))
      R <- K %*% (f$b * t(K))
      A <- diag(a)
      I <- diag(length(a))
      ref <- f$lambda *
        (solve(I - R %*% A + outer(rep(1, length(a)), a), R) +
          solve(I + Kxx %*% A, Kxx))
      expect_equal(H, ((ref + t(ref)) / 2)[C, C], tolerance = 1e-10)
    }
  }
})

test_that('Newton polish converges monotonically to a sparse optimum and recovers removed support', {
  f <- polish_fixture()
  ans <- polish_run(f)
  expect_identical(ans$status, 'refined')
  expect_true(causalOT:::cot_converged(ans$gap, ans$divergence, ans$atol, 1e-8))
  expect_lte(ans$loss, polish_loss(f, f$a))
  expect_equal(sum(ans$weights), 1, tolerance = 1e-14)
  expect_true(all(ans$weights >= 0))
  expect_true(any(ans$weights == 0))
  damaged <- ans$weights
  damaged[head(which(damaged > 0), 2)] <- 0
  restored <- polish_run(f, a = damaged / sum(damaged))
  expect_identical(restored$weights > 1e-8, ans$weights > 1e-8)
  expect_equal(restored$loss, ans$loss, tolerance = 1e-10)
  again <- polish_run(f, a = ans$weights)
  expect_identical(again$status, 'already_optimal')
  expect_equal(again$weights, ans$weights, tolerance = 1e-15)
})

test_that('a polish that cannot descend returns the start unchanged', {
  f <- polish_fixture()
  testthat::local_mocked_bindings(
    cot_polish_exact = function(H, g, a, ctl, state) g,
    .package = 'causalOT'
  )
  ans <- polish_run(f)
  expect_identical(ans$status, 'no_improvement')
  expect_false(ans$changed)
  expect_identical(ans$weights, f$a)
})

test_that('L-BFGS and online costs achieve the Newton answer', {
  f <- polish_fixture()
  newton <- polish_run(f)
  lbfgs <- polish_run(f, max_free = 2L, lbfgs_maxit = 80L)
  expect_lte(lbfgs$loss, polish_loss(f, f$a))
  expect_lt(abs(lbfgs$loss - newton$loss), 1e-6 * abs(newton$loss))
  online <- polish_run(f, mode = 'online')
  expect_equal(online$loss, newton$loss, tolerance = 1e-10)
  expect_equal(online$weights, newton$weights, tolerance = 1e-6)
})

# Fits ------------------------------------------------------------------------

test_that('polished fits converge, keep exact zeros, and return float64 weights', {
  skip_on_cran()
  f <- polish_fixture()
  reference <- polish_run(f)$weights
  for (solver in c('mirror', 'lbfgs')) {
    cot <- fit_cot(
      solver = solver,
      lambda = if (solver == 'mirror') 2 else c(4, 2),
      niter = 50L
    )
    expect_silent(cot$solve())
    if (solver == 'lbfgs') {
      # two lambdas and nboot = 0, so nothing is selected
      expect_warning(cot$grid_search(), 'nboot >= 1')
    } else {
      cot$grid_search()
    }
    info <- fit_optimizer(cot)$info()
    expect_true(all(unlist(info$converged)))
    expect_true(all(unlist(info$reason) == 'ok'))
    entry <- info$polish[['2']][[1]]
    expect_named(entry, c('status', 'iterations'))
    expect_true(entry$status %in% c('refined', 'already_optimal'))
    if (fit_optimizer(cot)$selected_lambda == 2) {
      expect_equal(cot$weights, reference, tolerance = 1e-6)
    }
    expect_type(cot$weights, 'double')
    expect_true(any(cot$weights == 0))
  }
})

test_that('float32 and float64 fits agree on the support and the weights', {
  causalOT:::torch_check()
  skip_on_cran()
  w <- lapply(c(torch::torch_double(), torch::torch_float()), function(dtype) {
    cot <- fit_cot(dtype = dtype, niter = 50L)
    cot$solve()
    cot$grid_search()
    expect_true(fit_optimizer(cot)$info()$converged[[1]][[1]])
    cot$weights
  })
  expect_equal(sum(w[[2]]), 1, tolerance = 1e-12)
  expect_identical(w[[1]] > 0, w[[2]] > 0)
  expect_lt(max(abs(w[[1]] - w[[2]])), 1e-4)
})

test_that('default accelerator fits are polished and return float64 weights', {
  causalOT:::torch_check()
  skip_if_not(
    torch::backends_mps_is_available() || torch::cuda_is_available(),
    'MPS/CUDA unavailable'
  )
  skip_if(
    causalOT:::cuda_device_check(NULL)$type == 'cpu',
    'Default device is CPU'
  )
  f <- polish_fixture()
  opts <- cotOptions(
    lambda = 2,
    solver = 'mirror',
    polish = 'auto',
    niter = 30L,
    nboot = 0,
    cost.online = 'tensorized'
  )
  cot <- causalOT:::COT$new(source = f$x, target = f$y, options = opts)
  cot$solve()
  cot$grid_search()
  expect_true(
    fit_optimizer(cot)$info()$polish[[1]][[1]]$status %in%
      c('refined', 'already_optimal')
  )
  expect_equal(sum(cot$weights), 1, tolerance = 1e-12)
  expect_equal(cot$weights, polish_run(f)$weights, tolerance = 1e-4)
})

test_that('a source containing the target converges at zero divergence through atol', {
  causalOT:::torch_check()
  cot <- fit_cot(
    lambda = .5,
    niter = 50L,
    source = matrix(c(-1, 0, 1), 3, 1),
    target = matrix(0, 1, 1)
  )
  cot$solve()
  info <- fit_optimizer(cot)$info()
  expect_true(info$converged[[1]][[1]])
  expect_identical(info$reason[[1]][[1]], 'ok')
  # A perturbed feasible point must fail the same test.
  p <- fit_private(cot)
  ot <- p$ot_objects[[ls(p$ot_objects)[[1]]]]
  a <- c(.01, .98, .01)
  ev <- causalOT:::cot_polish_eval(p$polish_cache$twin, a)
  D <- ev$loss - ev$y_self
  ctl <- causalOT:::cot_polish_ctl
  atol <- ctl$round_factor * .Machine$double.eps * causalOT:::cot_conv_scale(ot) +
    ctl$sinkhorn_factor * ctl$sinkhorn_tol
  expect_false(causalOT:::cot_converged(
    sum(a * ev$grad) - min(ev$grad),
    D,
    atol,
    1e-4
  ))
})

test_that('polished convergence uses the user tol and an unchanged polish keeps phase-1 weights', {
  baseline <- fit_cot('none')
  baseline$solve()
  for (changed in c(TRUE, FALSE)) {
    for (passes in c(TRUE, FALSE)) {
      cot <- fit_cot()
      local({
        testthat::local_mocked_bindings(
          cot_newton_polish = function(twin, ...) {
            a <- as.numeric(twin$a)
            a <- a / sum(a)
            threshold <- 1e-12 + 1e-4 * 2
            list(
              weights = a,
              loss = 2,
              divergence = 2,
              gap = if (passes) threshold / 2 else 2 * threshold,
              changed = changed,
              status = 'maxit',
              iterations = 3L,
              atol = 1e-12
            )
          },
          .package = 'causalOT'
        )
        cot$solve()
      })
      info <- fit_optimizer(cot)$info()
      expect_identical(info$converged[[1]][[1]], passes)
      expect_identical(
        info$reason[[1]][[1]],
        if (passes) 'ok' else 'gap above tolerance'
      )
      expect_identical(
        info$polish[[1]][[1]],
        list(status = 'maxit', iterations = 3L)
      )
      if (!changed) expect_identical(cot$weights, baseline$weights)
    }
  }
})

# Fit outcome rules -----------------------------------------------------------

test_that('unpolished convergence requires normal stopping and uses ot_tol', {
  for (stopped in c(FALSE, TRUE)) {
    cot <- fit_cot('none')
    replace_private(
      fit_private(cot),
      'optimization_loop',
      function(niter, tol) {
        fake_phase1(
          stopped_by_rule = stopped,
          fw_gap = 1e-10,
          divergence = 1e-5
        )
      }
    )
    cot$tol <- 1e-8
    cot$solve()
    info <- fit_optimizer(cot)$info()
    expect_identical(info$converged[[1]][[1]], stopped)
    expect_identical(info$reason[[1]][[1]], if (stopped) 'ok' else 'niter')
  }
  cot <- fit_cot('none')
  cot$solve()
  expect_identical(fit_optimizer(cot)$info()$reason[[1]][[1]], 'niter')
})

test_that('polish errors fall back to the phase-1 rule', {
  for (passes in c(TRUE, FALSE)) {
    cot <- fit_cot()
    replace_private(fit_private(cot), 'optimization_loop', function(...) {
      fake_phase1(fw_gap = if (passes) 1e-4 else .1)
    })
    local({
      testthat::local_mocked_bindings(
        cot_newton_polish = function(...) stop('test failure'),
        .package = 'causalOT'
      )
      cot$solve()
    })
    info <- fit_optimizer(cot)$info()
    expect_identical(info$converged[[1]][[1]], passes)
    expect_identical(
      info$reason[[1]][[1]],
      if (passes) 'ok' else 'polish error: test failure'
    )
    expect_identical(info$polish[[1]][[1]]$status, 'error')
  }
})

test_that('out-of-scope fits are skipped and retain the solver flag', {
  for (stopped in c(FALSE, TRUE)) {
    cot <- fit_cot(debias = FALSE)
    replace_private(fit_private(cot), 'optimization_loop', function(...) {
      list(loss = 1, iter = 1L, converged = stopped)
    })
    cot$solve()
    info <- fit_optimizer(cot)$info()
    expect_identical(info$converged[[1]][[1]], stopped)
    expect_identical(info$reason[[1]][[1]], 'gap unavailable')
    expect_match(info$polish[[1]][[1]]$status, '^skipped: .*debiased')
    expect_silent(fit_optimizer(cot)$choose_hyperparameters(
      n_boot_lambda = 0,
      n_boot_delta = 0
    ))
  }
  cot <- fit_cot('none')
  replace_private(fit_private(cot), 'optimization_loop', function(...) {
    fake_phase1(fw_gap = NA_real_)
  })
  cot$solve()
  expect_identical(
    fit_optimizer(cot)$info()$reason[[1]][[1]],
    'gap unavailable'
  )
})

test_that('numerical failures store NA, are excluded from selection, and all-failed paths error', {
  cot <- fit_cot('none', lambda = c(4, 2))
  p <- fit_private(cot)
  original <- p$optimization_loop
  replace_private(p, 'optimization_loop', function(niter, tol) {
    result <- original(niter, tol)
    if (p$ot_objects[[ls(p$ot_objects)[[1L]]]]$penalty == 4) {
      torch::with_no_grad(
        for (k in ls(p$parameters)) {
          p$parameters[[k]]$fill_(NaN)
        }
      )
      result$loss <- NaN
    }
    result
  })
  cot$solve()
  expect_match(p$reason_run[['4']][[1]], '^error: ')
  saved <- p$weights_list[['4']][[1]]
  expect_true(all(is.na(as.numeric(saved[[ls(saved)[[1]]]]))))
  optimizer <- fit_optimizer(cot)
  for (selector in c(
    'bootstrap_legacy',
    'profile_symmetric_loo',
    'mass_bootstrap_symmetric'
  )) {
    suppressWarnings(optimizer$choose_hyperparameters(
      n_boot_lambda = 1,
      n_boot_delta = 0,
      lambda.selector = selector
    ))
    expect_equal(optimizer$selected_lambda, 2)
  }
  expect_match(optimizer$info()$hyperparam.metrics$lambda.note, 'Excluded.*4')
  for (lambda in list(c(4, 2), 2)) {
    bad <- fit_cot(lambda = lambda)
    replace_private(
      fit_private(bad),
      'optimization_loop',
      function(niter, tol) list(loss = NaN, iter = 1L, converged = FALSE)
    )
    expect_error(bad$solve(), 'All lambda fits failed')
  }
  broken <- fit_cot()
  replace_private(
    fit_private(broken),
    'optimization_loop',
    function(niter, tol) stop('programming error')
  )
  expect_error(broken$solve(), 'programming error')
})

test_that('unpolished selection warns only when phase 1 ran out of iterations', {
  # gap above tolerance after phase 1 stopped by its rule: recorded, no warning
  cot <- fit_cot('none', lambda = c(4, 2))
  p <- fit_private(cot)
  replace_private(p, 'optimization_loop', function(niter, tol) {
    lambda <- p$ot_objects[[ls(p$ot_objects)[[1L]]]]$penalty
    fake_phase1(loss = 0, fw_gap = if (lambda == 4) .1 else 0)
  })
  cot$solve()
  expect_identical(unname(unlist(p$reason_run)), c('gap above tolerance', 'ok'))
  calls <- 0L
  replace_private(p, 'draw_boot_measure', function() NULL)
  replace_private(p, 'eval_lambda', function(w, boot) {
    calls <<- calls + 1L
    c(0, 1)[[calls]]
  })
  expect_no_warning(
    fit_optimizer(cot)$choose_hyperparameters(
      n_boot_lambda = 1,
      n_boot_delta = 0,
      lambda.selector = 'bootstrap_legacy'
    )
  )
  expect_equal(fit_optimizer(cot)$selected_lambda, 4)
  expect_false(fit_optimizer(cot)$info()$converged[['4']][[1]])

  # phase 1 hit niter: warn once
  cot <- fit_cot('none', lambda = c(4, 2))
  p <- fit_private(cot)
  replace_private(p, 'optimization_loop', function(niter, tol) {
    lambda <- p$ot_objects[[ls(p$ot_objects)[[1L]]]]$penalty
    fake_phase1(loss = 0, fw_gap = .1, stopped_by_rule = lambda != 4)
  })
  cot$solve()
  expect_identical(
    unname(unlist(p$reason_run)),
    c('niter', 'gap above tolerance')
  )
  calls <- 0L
  replace_private(p, 'draw_boot_measure', function() NULL)
  replace_private(p, 'eval_lambda', function(w, boot) {
    calls <<- calls + 1L
    c(0, 1)[[calls]]
  })
  warnings <- character()
  withCallingHandlers(
    fit_optimizer(cot)$choose_hyperparameters(
      n_boot_lambda = 1,
      n_boot_delta = 0,
      lambda.selector = 'bootstrap_legacy'
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart('muffleWarning')
    }
  )
  expect_equal(fit_optimizer(cot)$selected_lambda, 4)
  expect_length(warnings, 1L)
  expect_match(warnings, 'not converged \\(niter.*gap/D')
})

test_that('only auto and none are polish options', {
  causalOT:::torch_check()
  expect_error(cotOptions(polish = 'osqp'))
  expect_identical(cotOptions(polish = 'none')$polish, 'none')
})
