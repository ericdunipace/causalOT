followup_fixture <- function() {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(42)
  x <- matrix(rnorm(24), 12, 2); y <- matrix(rnorm(30), 15, 2)
  C <- function(x, y) pmax(outer(rowSums(x^2), rowSums(y^2), '+') / 2 - tcrossprod(x, y), 0)
  list(x = x, y = y, xy = C(x,y), xx = C(x,x), yy = C(y,y), a = rep(1/12,12), b = rep(1/15,15))
}
followup_fit <- function(polish = 'none', lambda = 2, dtype = torch::torch_double(), debias = TRUE) {
  f <- followup_fixture()
  opts <- cotOptions(lambda = lambda, solver = 'lbfgs', line_search_fn = 'strong_wolfe', polish = polish, debias = debias,
    niter = 2L, nboot = 0,
    cost.online = 'tensorized', dtype = dtype, device = torch::torch_device('cpu'))
  causalOT:::COT$new(source = f$x, target = f$y, options = opts)
}
followup_private <- function(cot) cot$.__enclos_env__$private$optimizer$.__enclos_env__$private

test_that('shared certificate and validity distinguish all numerical cases', {
  expect_identical(causalOT:::cot_certify(1e-12, 1e-3, 10, .Machine$double.eps, 1e-8), 'certified')
  expect_identical(causalOT:::cot_certify(1e-9, 1e-3, 10, .Machine$double.eps, 1e-8, TRUE), 'certified_resolution')
  expect_identical(causalOT:::cot_certify(1e-3, 1e-3, 10, .Machine$double.eps, 1e-8, TRUE), 'uncertified')
  expect_identical(causalOT:::cot_certify(NaN, 1, 1, .Machine$double.eps, 1e-8), 'uncertified')
  for (w in list(c(NaN,1), c(-1,2), c(0,0))) expect_false(causalOT:::cot_fit_valid(w, 0))
  expect_false(causalOT:::cot_fit_valid(c(.5,.5), Inf))
  expect_true(causalOT:::cot_fit_valid(c(.5,.5), 0))
})

test_that('warm probes agree with cold evaluations across lambda scales and supports', {
  f <- followup_fixture()
  xy <- causalOT:::cot_polish_cost(f$xy); xx <- causalOT:::cot_polish_cost(f$xx)
  for (mult in c(.05, .25, 1, 5, 25)) {
    lambda <- mult * median(f$xy)
    initial <- causalOT:::cot_polish_eval(f$a, f$b, xy, xx, lambda)
    a <- f$a; a[c(1,4,7)] <- 0; a <- a/sum(a)
    warm <- causalOT:::cot_polish_warm(a, f$b, xy, xx, lambda, initial$pot)
    cold <- causalOT:::cot_polish_eval(a, f$b, xy, xx, lambda)
    expect_lte(abs(warm$loss-cold$loss), 1e-12 * max(1,abs(cold$loss)))
    expect_equal(warm$grad, cold$grad, tolerance = 1e-10)
    run <- function(warm) causalOT:::cot_newton_polish(f$a,f$b,xy,xx,lambda,C_yy=f$yy,
      control=causalOT:::cot_polish_control(warm_start=warm))
    w <- run(TRUE); c <- run(FALSE)
    expect_identical(w$weights > 1e-12, c$weights > 1e-12)
    expect_lte(abs(w$loss-c$loss), 1e-12 * max(1,abs(c$loss)))
    expect_lte(w$timings[['other']], .05 * w$time + .02)
    expect_equal(length(w$probe_timings), w$line_search_evals)
  }
})

test_that('data costs preserve float64 inputs and bound memory', {
  f <- followup_fixture()
  measure <- causalOT:::Measure(f$x, dtype = torch::torch_float32(), device = torch::torch_device('cpu'))
  expect_equal(as.matrix(measure$x_float64), f$x, tolerance = 0)
  x <- causalOT:::cot_polish_tensor(f$x); y <- causalOT:::cot_polish_tensor(f$y)
  C <- causalOT:::costTensor$new(x$float(), y$float(), p=2)
  rebuilt <- causalOT:::cot_polish_data_cost(C,x,y,block_size=4L)
  expect_null(rebuilt$dense)
  expect_equal(as.matrix(rebuilt$get(1:12,1:15)), f$xy, tolerance=1e-12)
  l1 <- function(x,y,p) torch::torch_cdist(x,y,p=1)
  user <- causalOT:::costTensor$new(x$float(),y$float(),p=1,cost_function=l1)
  rebuilt <- causalOT:::cot_polish_data_cost(user,x,y,block_size=4L)
  expect_equal(as.matrix(rebuilt$get(1:12,1:15)), as.matrix(l1(x,y,1)), tolerance=1e-14)
  user$fun <- function(x,y,p) torch::torch_cdist(x,y,p=1)$float()
  expect_error(causalOT:::cot_polish_data_cost(user,x,y), 'did not return float64')
  user$fun <- function(x,y,p) stop('intentional float64 failure')
  expect_error(causalOT:::cot_polish_data_cost(user,x,y), 'intentional')
})

test_that('none is certified without changing phase-1 stopping and selected uncertified fits warn once', {
  cot <- followup_fit(lambda=c(4,2)); cot$solve()
  priv <- followup_private(cot)
  expect_true(all(unlist(priv$fit_status_run) == 'uncertified'))
  expect_false(any(unlist(priv$converged_run)))
  optimizer <- cot$.__enclos_env__$private$optimizer
  warnings <- character()
  withCallingHandlers(optimizer$choose_hyperparameters(n_boot_lambda=0,n_boot_delta=0,
    lambda.selector='profile_symmetric_loo'), warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')
    })
  expect_length(warnings,1L)
  expect_match(warnings,'niter')
  expect_match(optimizer$info()$hyperparam.metrics$lambda.note,'certification.*lambda')
  expect_identical(optimizer$info()$hyperparam.metrics$lambda.selector,'bootstrap_legacy')
  # Identical source/target measures have a zero FW gap and certify without refinement.
  f <- followup_fixture()
  opts <- cotOptions(lambda=2,solver='mirror',polish='none',niter=2L,nboot=0,
    cost.online='tensorized',dtype=torch::torch_double(),device=torch::torch_device('cpu'))
  same <- causalOT:::COT$new(source=f$x,target=f$x,options=opts)
  # Supply the known optimum directly: phase 1's loose Sinkhorn gradient can
  # otherwise move equal empirical measures away from their optimum.
  p <- followup_private(same)
  unlockBinding('optimization_loop',p)
  p$optimization_loop <- function(niter,tol) list(loss=0,iter=1L,converged=TRUE,stopped_by_rule=TRUE,fw_gap=0,divergence=0)
  lockBinding('optimization_loop',p)
  same$solve()
  expect_true(all(unlist(followup_private(same)$fit_status_run) == 'certified'))
})

test_that('numerical failures store NA and are excluded; optimization errors propagate', {
  cot <- followup_fit(lambda=c(4,2)); priv <- followup_private(cot)
  original <- priv$optimization_loop
  unlockBinding('optimization_loop',priv)
  priv$optimization_loop <- function(niter,tol) {
    result <- original(niter,tol)
    if (priv$ot_objects[[ls(priv$ot_objects)[[1L]]]]$penalty == 4) {
      torch::with_no_grad(for(k in ls(priv$parameters)) priv$parameters[[k]]$fill_(NaN))
      result$loss <- NaN
    }
    result
  }
  lockBinding('optimization_loop',priv)
  cot$solve()
  expect_identical(priv$fit_status_run[['4']][[1]],'failed')
  saved <- priv$weights_list[['4']][[1]]
  expect_true(all(is.na(as.numeric(saved[[ls(saved)[[1]]]]))))
  optimizer <- cot$.__enclos_env__$private$optimizer
  suppressWarnings(optimizer$choose_hyperparameters(n_boot_lambda=0,n_boot_delta=0))
  expect_equal(optimizer$selected_lambda,2)
  expect_match(optimizer$info()$hyperparam.metrics$lambda.note,'Excluded.*4')
  for (selector in c('profile_symmetric_loo','profile_marginal_loo',
      'symmetric_loo','marginal_loo','mass_bootstrap_symmetric')) {
    suppressWarnings(optimizer$choose_hyperparameters(n_boot_lambda=1,n_boot_delta=0,lambda.selector=selector))
    expect_equal(optimizer$selected_lambda,2)
  }
  bad <- followup_fit(lambda=c(4,2)); p <- followup_private(bad)
  unlockBinding('optimization_loop',p)
  p$optimization_loop <- function(niter,tol) list(loss=NaN,iter=1L,converged=FALSE)
  lockBinding('optimization_loop',p)
  expect_error(bad$solve(),'All lambda fits failed.*Phase-1')
  broken <- followup_fit(); p <- followup_private(broken)
  unlockBinding('optimization_loop',p)
  p$optimization_loop <- function(niter,tol) stop('programming error')
  lockBinding('optimization_loop',p)
  expect_error(broken$solve(),'programming error')
  single <- followup_fit(); p <- followup_private(single)
  unlockBinding('optimization_loop',p)
  p$optimization_loop <- function(niter,tol) list(loss=NaN,iter=1L,converged=FALSE)
  lockBinding('optimization_loop',p)
  expect_error(single$solve(),'All lambda fits failed')
})

test_that('ineligible configurations retain the stopping flag but numerical failures do not', {
  cot <- followup_fit(polish='auto',debias=FALSE); p <- followup_private(cot)
  original <- p$optimization_loop
  unlockBinding('optimization_loop',p)
  p$optimization_loop <- function(niter,tol) {r <- original(niter,tol); r$converged <- TRUE; r}
  lockBinding('optimization_loop',p)
  cot$solve()
  expect_identical(p$fit_status_run[[1]][[1]],'not_applicable')
  expect_true(p$converged_run[[1]][[1]])
  unlockBinding('optimization_loop',p)
  p$optimization_loop <- function(niter,tol) list(loss=NaN,iter=1L,converged=TRUE)
  lockBinding('optimization_loop',p)
  expect_error(cot$solve(),'All lambda fits failed')
})

test_that('hook records data costs, custom-cost fallback and redundant-solve removal', {
  cot <- followup_fit(polish='auto',dtype=torch::torch_float32()); cot$solve()
  p <- followup_private(cot); info <- cot$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]
  expect_identical(info$cost_source,'data_float64')
  expect_true(is.finite(info$loss_device))
  expect_gte(info$divergence,-1e-12)
  double <- followup_fit(polish='auto'); double$solve()
  info <- double$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]
  expect_true(info$dtype_recheck_skipped)
  expect_equal(info$timings[['dtype_recheck']],0,tolerance=.01)
  expect_lte(info$timings[['other']], .05 * info$time + .02)
  custom <- followup_fit(polish='auto',dtype=torch::torch_float32()); p <- followup_private(custom)
  ot <- p$ot_objects[[ls(p$ot_objects)[[1]]]]
  for (C in list(ot$C_xy,ot$C_xx,ot$C_yy)) {
    C$algorithm <- 'user'
    C$fun <- function(x,y,p) torch::torch_cdist(x,y,p=p)^p/p
  }
  custom$solve()
  expect_identical(custom$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]$cost_source,'data_float64')
  custom <- followup_fit(polish='auto',dtype=torch::torch_float32()); p <- followup_private(custom)
  ot <- p$ot_objects[[ls(p$ot_objects)[[1]]]]
  ot$C_xy$algorithm <- 'user'; ot$C_xy$fun <- function(x,y,p) stop('float64 unsupported')
  custom$solve()
  info <- custom$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]
  expect_identical(info$cost_source,'solver_rounded')
  expect_match(info$cost_fallback_reason,'float64 unsupported')
})


test_that('target self terms are cached per lambda without changing the answer', {
  f <- followup_fixture(); cache <- new.env(parent=emptyenv()); calls <- 0L
  original <- causalOT:::cot_polish_self
  testthat::local_mocked_bindings(cot_polish_self=function(...) {
    calls <<- calls+1L; original(...)
  }, .package='causalOT')
  run <- function(lambda) causalOT:::cot_newton_polish(f$a,f$b,f$xy,f$xx,lambda,C_yy=f$yy,cache=cache)
  first <- run(2); second <- run(2); third <- run(4)
  expect_equal(calls,2L)
  expect_equal(second$weights,first$weights,tolerance=0)
  expect_identical(second$loss,first$loss)
  expect_equal(length(ls(cache)),2L)
})

test_that('polish evaluator failure preserves valid phase-1 fits', {
  baseline <- followup_fit(); baseline$solve()
  cot <- followup_fit(polish='auto')
  testthat::local_mocked_bindings(cot_polish_eval=function(...) stop('cold evaluator failed'),.package='causalOT')
  cot$solve()
  p <- followup_private(cot)
  expect_identical(p$fit_status_run[[1]][[1]],'uncertified')
  expect_identical(cot$weights,baseline$weights)
  saved <- p$weights_list[[1]][[1]]
  expect_true(all(is.finite(as.numeric(saved[[ls(saved)[[1]]]]))))
  expect_match(cot$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]$reason,'cold evaluator failed')
})

test_that('rebuilt block and dense costs agree for default p costs', {
  f <- followup_fixture(); x <- causalOT:::cot_polish_tensor(f$x); y <- causalOT:::cot_polish_tensor(f$y)
  for(p in c(1,2,3)) {
    Cxy <- causalOT:::costTensor$new(x,y,p=p); Cxx <- causalOT:::costTensor$new(x,x,p=p)
    xy <- causalOT:::cot_polish_data_cost(Cxy,x,y,4L)
    xx <- causalOT:::cot_polish_data_cost(Cxx,x,x,4L)
    block <- causalOT:::cot_polish_eval(f$a,f$b,xy,xx,2)
    dense <- causalOT:::cot_polish_eval(f$a,f$b,Cxy,Cxx,2)
    expect_lte(abs(block$loss-dense$loss),1e-12*max(1,abs(dense$loss)))
    expect_equal(block$grad,dense$grad,tolerance=1e-12)
    warm <- causalOT:::cot_polish_warm(f$a,f$b,xy,xx,2,dense$pot)
    expect_lte(abs(warm$loss-dense$loss),1e-12*max(1,abs(dense$loss)))
  }
})

test_that('cold probes reuse the final evaluation', {
  f <- followup_fixture()
  ans <- causalOT:::cot_newton_polish(f$a,f$b,f$xy,f$xx,2,C_yy=f$yy,
    control=causalOT:::cot_polish_control(maxit=1L,warm_start=FALSE))
  expect_identical(ans$readback_evals,0L)
  expect_equal(ans$timings[['readback']],0,tolerance=.01)
})

test_that('unpolished certification performs no Sinkhorn solves after phase 1', {
  cot <- followup_fit(); p <- followup_private(cot)
  in_hook <- FALSE; calls <- c(eval=0L,self=0L,y_self=0L); phase_time <- 0
  original_hook <- p$polish_result; original_loop <- p$optimization_loop
  unlockBinding('polish_result',p); unlockBinding('optimization_loop',p)
  p$optimization_loop <- function(...) {t <- proc.time()[['elapsed']]; out <- original_loop(...); phase_time <<- proc.time()[['elapsed']]-t; out}
  p$polish_result <- function(...) {in_hook <<- TRUE; on.exit(in_hook <<- FALSE); original_hook(...)}
  lockBinding('polish_result',p); lockBinding('optimization_loop',p)
  ev <- causalOT:::cot_polish_eval; self <- causalOT:::cot_polish_self; yy <- causalOT:::cot_ccp_y_self_term
  testthat::local_mocked_bindings(
    cot_polish_eval=function(...) {if(in_hook) calls[['eval']] <<- calls[['eval']]+1L; ev(...)},
    cot_polish_self=function(...) {if(in_hook) calls[['self']] <<- calls[['self']]+1L; self(...)},
    cot_ccp_y_self_term=function(...) {if(in_hook) calls[['y_self']] <<- calls[['y_self']]+1L; yy(...)}, .package='causalOT')
  cot$solve()
  expect_equal(unname(calls),c(0L,0L,0L))
  info <- cot$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]
  expect_identical(info$certification$source,'phase1')
  expect_lte(info$time,.05*phase_time+.02)
})

test_that('phase-1 certification requires its stop rule and has no loose resolution tier', {
  for(stopped in c(FALSE,TRUE)) {
    cot <- followup_fit(); p <- followup_private(cot)
    unlockBinding('optimization_loop',p)
    p$optimization_loop <- function(niter,tol) list(loss=0,iter=1L,converged=TRUE,
      stopped_by_rule=stopped,fw_gap=1e-10,divergence=1e-5)
    lockBinding('optimization_loop',p)
    # Set a strict solver tolerance while leaving polish disabled.
    cot$tol <- 1e-8
    cot$solve()
    expect_identical(p$fit_status_run[[1]][[1]],'uncertified')
    diag <- cot$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]
    expect_match(diag$reason,if(stopped) 'gap above' else 'niter')
  }
  expect_error(cotOptions(polish.control=list(certify='float64')),'Unknown polish control')
})

test_that('bootstrap pool includes uncertified fits alongside certified fits', {
  cot <- followup_fit(lambda=c(4,2)); p <- followup_private(cot)
  unlockBinding('optimization_loop',p)
  p$optimization_loop <- function(niter,tol) {
    lambda <- p$ot_objects[[ls(p$ot_objects)[[1L]]]]$penalty
    list(loss=0,iter=1L,converged=TRUE,stopped_by_rule=TRUE,
      fw_gap=if(lambda==4) .1 else 0,divergence=1)
  }
  lockBinding('optimization_loop',p)
  cot$solve()
  expect_identical(unname(unlist(p$fit_status_run)),c('uncertified','certified'))
  calls <- 0L
  unlockBinding('draw_boot_measure',p); unlockBinding('eval_lambda',p)
  p$draw_boot_measure <- function() NULL
  p$eval_lambda <- function(w,boot) {calls <<- calls+1L; c(0,1)[[calls]]}
  lockBinding('draw_boot_measure',p); lockBinding('eval_lambda',p)
  optimizer <- cot$.__enclos_env__$private$optimizer
  warnings <- character()
  withCallingHandlers(optimizer$choose_hyperparameters(n_boot_lambda=1,n_boot_delta=0,
    lambda.selector='bootstrap_legacy'),warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  expect_equal(calls,2L)
  expect_equal(optimizer$selected_lambda,4)
  expect_length(warnings,1L)
  expect_match(warnings,'uncertified')
})

test_that('default accelerator fits actually refine using float64 data costs', {
  skip_if_not(torch::backends_mps_is_available() || torch::cuda_is_available(), 'MPS/CUDA unavailable')
  skip_if(causalOT:::cuda_device_check(NULL)$type == 'cpu','Default device is CPU')
  f <- followup_fixture()
  opts <- cotOptions(lambda=2,solver='mirror',niter=30L,nboot=0,cost.online='tensorized')
  cot <- causalOT:::COT$new(source=f$x,target=f$y,options=opts)
  cot$solve()
  info <- cot$.__enclos_env__$private$optimizer$info()$polish[[1]][[1]]
  expect_identical(info$cost_source,'data_float64')
  expect_false(info$status %in% c('skipped','failed'))
  expect_true(causalOT:::cot_polish_tensor(torch::torch_tensor(1:3,device=causalOT:::cuda_device_check(NULL),dtype=torch::torch_float32()))$dtype == torch::torch_float64())
})
