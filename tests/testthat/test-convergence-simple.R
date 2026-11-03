convergence_fit <- function(polish = 'none', niter = 2L, debias = TRUE) {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(42)
  opts <- cotOptions(lambda = 2, solver = 'lbfgs', line_search_fn = 'strong_wolfe',
    polish = polish, debias = debias, niter = niter, nboot = 0, cost.online = 'tensorized',
    dtype = torch::torch_double(), device = torch::torch_device('cpu'))
  causalOT:::COT$new(source = matrix(rnorm(24),12,2), target = matrix(rnorm(30),15,2), options = opts)
}
convergence_optimizer <- function(cot) cot$.__enclos_env__$private$optimizer
convergence_private <- function(cot) convergence_optimizer(cot)$.__enclos_env__$private

test_that('one tolerance test handles absolute and relative accuracy and finite inputs', {
  check <- causalOT:::cot_converged
  expect_true(check(1e-3, 0, 1e-3, 1e-4))
  expect_false(check(1.01e-3, 0, 1e-3, 1e-4))
  expect_true(check(.101, 1000, .001, 1e-4))
  expect_false(check(.102, 1000, .001, 1e-4))
  expect_true(check(.001, -10, .001, 1e-4))
  expect_false(check(.002, -10, .001, 1e-4))
  for (i in 1:4) for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    args <- list(.001, 1, .001, 1e-4); args[[i]] <- bad
    expect_false(do.call(check, args))
  }
  expect_false(check(NULL, 1, .001, 1e-4))
})

test_that('iteration cap fails even when the phase-1 gap passes', {
  cot <- convergence_fit(); p <- convergence_private(cot)
  original <- p$optimization_loop
  unlockBinding('optimization_loop', p)
  p$optimization_loop <- function(niter, tol) {
    res <- original(niter, tol)
    expect_identical(res$iter, 2L)
    expect_false(res$stopped_by_rule)
    res$fw_gap <- 0; res$divergence <- 0
    res
  }
  lockBinding('optimization_loop', p)
  cot$solve(); info <- convergence_optimizer(cot)$info()
  expect_false(info$converged[[1]][[1]])
  expect_identical(info$reason[[1]][[1]], 'niter')
  expect_null(info$fit_status)
  expect_null(info$polish[[1]][[1]]$certification)
})

test_that('a source containing the target converges at zero divergence through atol', {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  opts <- cotOptions(lambda = .5, solver = 'lbfgs', line_search_fn = 'strong_wolfe',
    polish = 'auto', niter = 50L, nboot = 0, cost.online = 'tensorized',
    dtype = torch::torch_double(), device = torch::torch_device('cpu'))
  cot <- causalOT:::COT$new(source = matrix(c(-1,0,1),3,1), target = matrix(0,1,1), options = opts)
  cot$solve(); info <- convergence_optimizer(cot)$info(); diag <- info$polish[[1]][[1]]$convergence
  expect_true(info$converged[[1]][[1]])
  expect_identical(info$reason[[1]][[1]], 'ok')
  expect_identical(diag$phase, 'polish')
  expect_lt(abs(diag$D), 1e-12)
  expect_gt(diag$atol, 0)
  expect_lte(diag$gap, diag$atol + diag$tol * max(0, diag$D))
  # A cold evaluation of a perturbed feasible point must fail that same test.
  p <- convergence_private(cot); ot <- p$ot_objects[[ls(p$ot_objects)[[1]]]]; C <- p$polish_cache$costs
  ev <- causalOT:::cot_polish_eval(c(.01,.98,.01), ot$b, C$xy, C$xx, ot$penalty)
  yy <- causalOT:::cot_polish_self(ot$b, C$yy, ot$penalty, causalOT:::cot_polish_control())
  gap <- sum(c(.01,.98,.01) * ev$grad) - min(ev$grad)
  expect_false(causalOT:::cot_converged(gap, ev$loss - yy, diag$atol, diag$tol))
})

test_that('final polished convergence uses user tol and rejects a larger gap', {
  original <- causalOT:::cot_newton_polish
  for (passes in c(TRUE, FALSE)) {
    cot <- convergence_fit('auto')
    local({
      testthat::local_mocked_bindings(cot_newton_polish = function(...) {
        res <- original(...)
        threshold <- 64 * .Machine$double.eps * res$conv_scale + 1e-4 * max(0, res$divergence)
        res$fw_gap <- if (passes) threshold / 2 else 2 * threshold
        res
      }, .package = 'causalOT')
      cot$solve()
    })
    info <- convergence_optimizer(cot)$info(); diag <- info$polish[[1]][[1]]$convergence
    expect_identical(diag$phase, 'polish')
    expect_identical(diag$tol, 1e-4)
    if (passes) expect_gt(diag$gap, 1e-8 * max(0, diag$D))
    expect_identical(info$converged[[1]][[1]], passes)
    expect_identical(info$reason[[1]][[1]], if (passes) 'ok' else 'gap above tolerance')
  }
})

test_that('rejected polish judges the cold starting point regardless of phase-1 cap', {
  for (passes in c(TRUE, FALSE)) {
    cot <- convergence_fit('auto')
    local({
      testthat::local_mocked_bindings(cot_newton_polish = function(a, ...) {
        list(weights = as.numeric(a), solver_weights = as.numeric(a),
          status = 'maxit', n_accepted = 0L, loss = 2, loss_start = 1,
          fw_gap = 1, fw_gap_start = if (passes) 1e-6 else 1e-2,
          divergence = 2, divergence_start = 1, conv_scale = 10, timings = c(other = 0))
      }, .package = 'causalOT')
      cot$solve()
    })
    info <- convergence_optimizer(cot)$info(); diag <- info$polish[[1]][[1]]
    expect_identical(diag$status, 'maxit')
    expect_identical(diag$convergence$phase, 'polish')
    expect_equal(diag$convergence$gap, if (passes) 1e-6 else 1e-2)
    expect_equal(diag$convergence$D, 1)
    expect_identical(info$converged[[1]][[1]], passes)
    expect_identical(info$reason[[1]][[1]], if (passes) 'ok' else 'gap above tolerance')
  }
})

test_that('polish errors use phase-1 convergence and only report an error if it fails', {
  for (passes in c(TRUE,FALSE)) {
    cot <- convergence_fit('auto'); p <- convergence_private(cot)
    unlockBinding('optimization_loop',p)
    p$optimization_loop <- function(...) list(loss=1,iter=1L,converged=TRUE,
      stopped_by_rule=TRUE,fw_gap=if(passes) 1e-4 else .1,divergence=1)
    lockBinding('optimization_loop',p)
    local({
      testthat::local_mocked_bindings(cot_newton_polish=function(...) stop('test failure'),.package='causalOT')
      cot$solve()
    })
    info <- convergence_optimizer(cot)$info()
    expect_identical(info$converged[[1]][[1]],passes)
    expect_identical(info$reason[[1]][[1]],if(passes) 'ok' else 'polish error: test failure')
    expect_identical(info$polish[[1]][[1]]$convergence$phase,'phase1')
  }
})


test_that('out-of-scope fits retain either solver flag with gap unavailable', {
  for(stopped in c(FALSE,TRUE)) {
    cot <- convergence_fit('auto',debias=FALSE); p <- convergence_private(cot)
    unlockBinding('optimization_loop',p)
    p$optimization_loop <- function(...) list(loss=1,iter=1L,converged=stopped)
    lockBinding('optimization_loop',p)
    cot$solve(); info <- convergence_optimizer(cot)$info()
    expect_identical(info$converged[[1]][[1]],stopped)
    expect_identical(info$reason[[1]][[1]],'gap unavailable')
    expect_true(is.na(info$polish[[1]][[1]]$convergence$gap))
    expect_match(info$polish[[1]][[1]]$skip_reason,'debiased')
    expect_silent(convergence_optimizer(cot)$choose_hyperparameters(n_boot_lambda=0,n_boot_delta=0))
  }
})

test_that('a normally stopped phase-1 fit reports unavailable non-finite gap', {
  cot <- convergence_fit(); p <- convergence_private(cot)
  unlockBinding('optimization_loop',p)
  p$optimization_loop <- function(...) list(loss=1,iter=1L,converged=TRUE,
    stopped_by_rule=TRUE,fw_gap=NA_real_,divergence=1)
  lockBinding('optimization_loop',p)
  cot$solve(); info <- convergence_optimizer(cot)$info()
  expect_false(info$converged[[1]][[1]])
  expect_identical(info$reason[[1]][[1]],'gap unavailable')
})

test_that('a zero-step optimal polish retains all refinement diagnostics', {
  cot <- convergence_fit('auto')
  testthat::local_mocked_bindings(cot_newton_polish=function(a,...) {
    list(weights=as.numeric(a),solver_weights=as.numeric(a),status='certified',
      n_accepted=0L,nnz=length(a),iter=0L,loss=1,loss_start=1,fw_gap=0,fw_gap_start=0,
      divergence=1,divergence_start=1,conv_scale=10,timings=c(start=.01,other=0),
      line_search_evals=0L,direction='exact')
  },.package='causalOT')
  cot$solve(); info <- convergence_optimizer(cot)$info(); diag <- info$polish[[1]][[1]]
  expect_identical(diag$status,'already_optimal')
  expect_identical(diag$n_accepted,0L)
  expect_identical(diag$iter,0L)
  expect_identical(diag$nnz,12L)
  expect_identical(diag$line_search_evals,0L)
  expect_identical(diag$direction,'exact')
  expect_equal(diag$timings[['start']],.01)
  expect_true(all(c('dtype_recheck','device_refresh','other') %in% names(diag$timings)))
  expect_true(info$converged[[1]][[1]])
  expect_identical(info$reason[[1]][[1]],'ok')
  expect_identical(diag$convergence$phase,'polish')
})

test_that('rejected float32 dtype round-trip retains diagnostics and is not a polish error', {
  for(passes in c(FALSE,TRUE)) {
    causalOT:::torch_check()
    # Use the actual float32 weight installation and rounding rejection path.
    cot <- causalOT:::COT$new(source=matrix(seq_len(24)/10,12,2),target=matrix(seq_len(30)/10,15,2),
      options=cotOptions(lambda=2,solver='lbfgs',line_search_fn='strong_wolfe',polish='auto',
        niter=2L,nboot=0,cost.online='tensorized',dtype=torch::torch_float32(),device=torch::torch_device('cpu')))
    p <- convergence_private(cot)
    baseline <- as.numeric(p$ot_objects[[ls(p$ot_objects)[[1]]]]$a)
    unlockBinding('optimization_loop',p)
    p$optimization_loop <- function(...) list(loss=1,iter=2L,converged=FALSE,
      stopped_by_rule=FALSE,fw_gap=.5,divergence=1)
    lockBinding('optimization_loop',p)
    local({
      testthat::local_mocked_bindings(cot_newton_polish=function(a,...) {
        list(weights=as.numeric(a),solver_weights=as.numeric(a),status='maxit',
          n_accepted=1L,nnz=length(a),iter=1L,loss=.5,loss_start=1,fw_gap=0,
          fw_gap_start=if(passes) 1e-6 else .01,divergence=.5,divergence_start=1,
          y_self=0,conv_scale=10,timings=c(start=.01,other=0),line_search_evals=2L)
      },cot_polish_eval=function(...) list(loss=2,grad=rep(0,12)),.package='causalOT')
      cot$solve()
    })
    info <- convergence_optimizer(cot)$info(); diag <- info$polish[[1]][[1]]
    expect_equal(cot$weights,baseline,tolerance=0)
    expect_identical(diag$status,'no_improvement')
    expect_match(diag$rejection_reason,'No improvement after restoring the solver dtype',fixed=TRUE)
    expect_null(diag$polish_error)
    expect_identical(diag$n_accepted,1L)
    expect_identical(diag$nnz,12L)
    expect_identical(diag$line_search_evals,2L)
    expect_equal(diag$timings[['start']],.01)
    expect_equal(diag$fw_gap,if(passes) 1e-6 else .01)
    expect_equal(diag$divergence,1)
    expect_identical(diag$convergence$phase,'polish')
    expect_identical(info$converged[[1]][[1]],passes)
    expect_identical(info$reason[[1]][[1]],if(passes) 'ok' else 'gap above tolerance')
  }
})

test_that('a real dtype evaluator error still falls back to the phase-1 rule', {
  for(passes in c(FALSE,TRUE)) {
    causalOT:::torch_check()
    cot <- causalOT:::COT$new(source=matrix(seq_len(24)/10,12,2),target=matrix(seq_len(30)/10,15,2),
      options=cotOptions(lambda=2,solver='lbfgs',line_search_fn='strong_wolfe',polish='auto',niter=2L,nboot=0,
        cost.online='tensorized',dtype=torch::torch_float32(),device=torch::torch_device('cpu')))
    p <- convergence_private(cot); before <- as.numeric(cot$weights)
    unlockBinding('optimization_loop',p)
    p$optimization_loop <- function(...) list(loss=1,iter=1L,converged=TRUE,
      stopped_by_rule=TRUE,fw_gap=if(passes) 1e-4 else .1,divergence=1)
    lockBinding('optimization_loop',p)
    local({
      testthat::local_mocked_bindings(cot_newton_polish=function(a,...) {
        list(weights=as.numeric(a),solver_weights=as.numeric(a),status='maxit',n_accepted=1L,
          nnz=length(a),loss=.5,loss_start=1,fw_gap=0,fw_gap_start=0,divergence=1,
          divergence_start=1,conv_scale=10,timings=c(other=0))
      },cot_polish_eval=function(...) stop('dtype evaluator failed'),.package='causalOT')
      cot$solve()
    })
    info <- convergence_optimizer(cot)$info();diag <- info$polish[[1]][[1]]
    expect_equal(cot$weights,before,tolerance=0)
    expect_identical(info$converged[[1]][[1]],passes)
    expect_identical(info$reason[[1]][[1]],if(passes) 'ok' else 'polish error: dtype evaluator failed')
    expect_identical(diag$convergence$phase,'phase1')
    expect_identical(diag$polish_error,'dtype evaluator failed')
    expect_null(diag$rejection_reason)
    expect_identical(diag$n_accepted,1L)
  }
})
