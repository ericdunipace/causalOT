# One test file at a time; run outside the sandbox to exercise MPS.
# Rscript dev/convergence_simple_validate.R tests default [test-file-filter ...]
# Rscript dev/convergence_simple_validate.R regression [package-root] [output-stem]
# Rscript dev/convergence_simple_validate.R phase1 [package-root] [output-stem]
# Rscript dev/convergence_simple_validate.R compare [base-stem] [branch-stem]
# Rscript dev/convergence_simple_validate.R equivalence [base-snapshot-root]
# Regression is intentionally opt-in: ask the user before running it.
args <- commandArgs(TRUE)
mode <- if(length(args)) args[[1]] else 'tests'
if(mode == 'equivalence') {
  # The base snapshot must be extracted from 74fa9c9, never the user's live tree.
  base <- if(length(args)>1) args[[2]] else '/private/tmp/causalOT-convergence-base'
  stages <- list(
    base_phase1=c('phase1',base,'dev/convergence_simple_base_phase1'),
    branch_phase1=c('phase1','.','dev/convergence_simple_branch_phase1'),
    compare_phase1=c('compare','dev/convergence_simple_base_phase1','dev/convergence_simple_branch_phase1'),
    base_regression=c('regression',base,'dev/convergence_simple_base_regression'),
    branch_regression=c('regression','.','dev/convergence_simple_branch_regression'),
    compare_regression=c('compare','dev/convergence_simple_base_regression','dev/convergence_simple_branch_regression'))
  for(stage in names(stages)) {
    log <- paste0('dev/convergence_simple_',stage,'.log')
    cat('Starting',stage,'; log:',log,'\n'); flush.console()
    status <- system2(file.path(R.home('bin'),'Rscript'),
      c('dev/convergence_simple_validate.R',shQuote(stages[[stage]])),stdout=log,stderr=log)
    cat('Finished',stage,'; exit:',status,'\n'); flush.console()
    if(status!=0) {cat(tail(readLines(log),60),sep='\n'); stop('Validation failed: ',stage)}
  }
  quit(save='no')
}
if (mode == 'compare') {
  base <- readRDS(paste0(args[[2]], '.rds')); branch <- readRDS(paste0(args[[3]], '.rds'))
  stopifnot(length(base) == length(branch))
  rows <- lapply(seq_along(base), function(i) {
    fields <- intersect(c('configuration','profile','method','mult','polish'), names(base[[i]]))
    stopifnot(identical(base[[i]][fields], branch[[i]][fields]))
    if ('phase1' %in% names(base[[i]])) {
      phase1_equal <- identical(base[[i]]$phase1,branch[[i]]$phase1)
      weight_equal <- identical(base[[i]]$weights,branch[[i]]$weights)
      support_equal <- identical(base[[i]]$weights>0,branch[[i]]$weights>0)
      weight_difference <- max(abs(base[[i]]$weights-branch[[i]]$weights))
      polished <- identical(base[[i]]$polish,'auto')
      loss_difference <- if(polished) abs(base[[i]]$info$loss-branch[[i]]$info$loss) else 0
      loss_floor <- if(polished) 64*.Machine$double.eps*base[[i]]$info$conv_scale else 0
      # Phase 1 and unpolished weights must be bit-for-bit identical. Polished
      # points must retain support and agree in loss to evaluator roundoff.
      stopifnot(phase1_equal,support_equal,polished || weight_equal,loss_difference<=loss_floor)
      data.frame(base[[i]][fields],phase1_identical=phase1_equal,weights_identical=weight_equal,
        support_identical=support_equal,max_weight_difference=weight_difference,
        polished_loss_difference=loss_difference,polished_loss_floor=loss_floor)
    } else {
      keys <- c('iterations','loss','weights','evaluation','potentials')
      equal <- vapply(keys,function(k) identical(base[[i]][[k]],branch[[i]][[k]]),logical(1))
      stopifnot(all(equal))
      data.frame(base[[i]][fields], identical = all(equal))
    }
  })
  print(do.call(rbind, rows))
  write.csv(do.call(rbind, rows), paste0(args[[3]], '_equivalence.csv'),row.names=FALSE)
  quit(save='no')
}
root <- if(mode %in% c('regression','phase1') && length(args)>=2) args[[2]] else '.'
stem <- if(length(args)>=3 && mode != 'tests') args[[3]] else paste0('dev/convergence_simple_',mode)
suppressMessages(pkgload::load_all(root, quiet=TRUE))
torch::torch_set_num_threads(1L)
Sys.setenv(NOT_CRAN='true')
cat('Package:', normalizePath(root), '\nMPS:', torch::backends_mps_is_available(),
    ' Default:', as.character(causalOT:::cuda_device_check(NULL)), '\n')
replace <- function(env,name,value) {
  unlockBinding(name,env); env[[name]] <- value; lockBinding(name,env)
}
priv <- function(cot) cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
ot_get <- function(p) p$ot_objects[[ls(p$ot_objects)[[1]]]]
if(mode == 'tests') {
  profile <- if(length(args)>=2) args[[2]] else 'default'
  run <- function() {
    stopifnot(profile == 'default')
    suites <- if(length(args)>2) args[-c(1,2)] else
      c('convergence-simple','newton-followup','newton-polish','lambda-selector','ccp-solver')
    summary_file <- paste0('dev/convergence_simple_tests_',profile,'.csv')
    rows <- if(file.exists(summary_file)) {
      previous <- read.csv(summary_file)
      split(previous,previous$suite)
    } else list()
    for(suite in suites) {
      ans <- devtools::test(filter=paste0('^',suite,'$'),reporter='summary',stop_on_failure=FALSE)
      d <- as.data.frame(ans)
      rows[[suite]] <- data.frame(profile=profile,suite=suite,failed=sum(d$failed),
        errors=sum(d$error),warnings=sum(d$warning),skipped=sum(d$skipped),passed=sum(d$passed))
      write.csv(do.call(rbind,rows),paste0('dev/convergence_simple_tests_',profile,'.csv'),row.names=FALSE)
      saveRDS(ans,paste0('dev/convergence_simple_test_',profile,'_',suite,'.rds'))
    }
    print(do.call(rbind,rows))
    if(any(vapply(rows,function(x) x$failed+x$errors>0,logical(1)))) stop('Focused tests failed')
  }
  run()
} else if(mode == 'phase1') {
  # Ported from 22b2231: compare all legacy phase-1 outputs with the hook bypassed.
  rows <- list()
  for(profile in c('default','cpu_float64')) for(method in c('mirror','lbfgs','ccp_accelerated')) {
    set.seed(804); x <- matrix(rnorm(96),32,3); y <- matrix(rnorm(144),48,3)
    opts <- cotOptions(lambda=5,solver=if(method=='lbfgs') 'lbfgs' else 'mirror',polish='none',niter=100L,
      tol=1e-4,nboot=0,cost.online='tensorized',line_search_fn='strong_wolfe',
      dtype=if(profile=='default') NULL else torch::torch_double(),
      device=if(profile=='default') NULL else torch::torch_device('cpu'))
    cot <- COT$new(source=x,target=y,options=opts); p <- priv(cot)
    replace(p,'polish_result',function(res,tol,delta) res)
    if(method=='ccp_accelerated') p$ccp_step_rule <- 'accelerated'
    cot$solve(); ot <- ot_get(p)
    ev <- ot$sinkhorn_ccp_eval(niter=11L,tol=1e-9)
    rows[[length(rows)+1L]] <- list(profile=profile,method=method,
      iterations=p$iterations_run,loss=p$final_loss,weights=as.numeric(cot$weights),
      evaluation=lapply(ev,as.numeric),potentials=lapply(ot$potentials,as.numeric))
  }
  saveRDS(rows,paste0(stem,'.rds'))
  print(lapply(rows,function(x) x[c('profile','method','iterations','loss')]))
} else if(mode == 'regression') {
  set.seed(1); n <- 256L; x <- matrix(rnorm(n*3),n,3)
  z <- rbinom(n,1,plogis(.8*x[,1]-.5*x[,2]+.3*x[,3])); scale <- median(as.numeric(dist(x))^2)
  rows <- list(); diagnostics <- list()
  configurations <- if(length(args)>=4) args[[4]] else c('cpu_float64','cpu_float32','default_float32')
  multipliers <- if(length(args)>=5) as.numeric(args[[5]]) else c(1,5,25)
  polishes <- if(length(args)>=6) args[[6]] else c('none','auto')
  for(configuration in configurations) for(mult in multipliers) for(polish in polishes) {
    if(configuration=='default_float32' && !torch::backends_mps_is_available() && !torch::cuda_is_available())
      stop('Accelerator unavailable: run this harness outside the sandbox')
    opts <- cotOptions(lambda=mult*scale,solver='mirror',polish=polish,niter=100L,nboot=0,
      cost.online='tensorized',device=if(configuration=='default_float32') NULL else torch::torch_device('cpu'),
      dtype=if(configuration=='cpu_float64') torch::torch_double() else torch::torch_float32())
    cot <- COT$new(source=x[z==0,,drop=FALSE],target=x,options=opts); p <- priv(cot)
    loop_res <- NULL; loop <- p$optimization_loop
    replace(p,'optimization_loop',function(...) {loop_res <<- loop(...); loop_res})
    t <- proc.time()[['elapsed']]; cot$solve(); elapsed <- proc.time()[['elapsed']]-t
    outcome <- cot$.__enclos_env__$private$optimizer$info()
    info <- outcome$polish[[1]][[1]]
    simple <- !is.null(info$convergence)
    diag <- if(simple) info$convergence else list(phase=info$certification$source,
      gap=info$certification$fw_gap,D=info$certification$divergence,atol=NA_real_,tol=opts$tol)
    converged <- outcome$converged[[1]][[1]]
    expected <- if(polish=='none') 134L else c('1'=39L,'5'=34L,'25'=31L)[[as.character(mult)]]
    stopifnot(sum(cot$weights>0)==expected)
    if(simple) {
      stopifnot(identical(converged, polish=='auto'))
      stopifnot(identical(converged, causalOT:::cot_converged(diag$gap,diag$D,diag$atol,diag$tol)))
    }
    rows[[length(rows)+1L]] <- data.frame(configuration=configuration,device=as.character(ot_get(p)$device),
      mult=mult,polish=polish,iterations=loop_res$iter,phase1_loss=loop_res$loss,
      stopped=isTRUE(loop_res$stopped_by_rule),time=elapsed,converged=converged,
      reason=if(simple) outcome$reason[[1]][[1]] else outcome$fit_status[[1]][[1]],
      phase=diag$phase,gap=diag$gap,D=diag$D,atol=diag$atol,tol=diag$tol,
      threshold=diag$atol+diag$tol*max(0,diag$D),nnz=sum(cot$weights>0))
    diagnostics[[length(diagnostics)+1L]] <- list(configuration=configuration,mult=mult,polish=polish,
      phase1=loop_res,weights=as.numeric(cot$weights),info=info)
    print(tail(rows,1)[[1]])
    write.csv(do.call(rbind,rows),paste0(stem,'.csv'),row.names=FALSE)
    saveRDS(diagnostics,paste0(stem,'.rds'))
  }
} else stop('Unknown mode')
