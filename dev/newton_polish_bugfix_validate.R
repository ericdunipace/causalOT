# Run from the package root:
# NOT_CRAN=true Rscript dev/newton_polish_bugfix_validate.R default
# NOT_CRAN=true Rscript dev/newton_polish_bugfix_validate.R cpu_float64
# Rscript dev/newton_polish_bugfix_validate.R regression
args <- commandArgs(TRUE)
profile <- args[[1L]]
package_root <- if(profile=='phase1' && length(args)>=2L) args[[2L]] else '.'
suppressMessages(pkgload::load_all(package_root,quiet=TRUE))
torch::torch_set_num_threads(1L)
cat('MPS available:',torch::backends_mps_is_available(),' CUDA available:',torch::cuda_is_available(),'\n')
cat('Actual default device:',as.character(causalOT:::cuda_device_check(NULL)),'\n')
run_suites <- function(profile) {
  Sys.setenv(NOT_CRAN="true")
  if(profile=='cpu_float64') {
    old_device <- causalOT:::cuda_device_check; old_dtype <- causalOT:::cuda_dtype_check
    testthat::local_mocked_bindings(
      cuda_device_check=function(device=NULL) old_device(if(is.null(device)) torch::torch_device('cpu') else device),
      cuda_dtype_check=function(dtype=NULL,device=NULL) old_dtype(if(is.null(dtype)) torch::torch_double() else dtype,device),
      .package='causalOT')
  }
  suites <- c('newton-followup','newton-polish','lambda-selector','ccp-solver','cotOOP',
    'cotProblem','cotClass','gridSearch','costMethod','dual_opts')
  if(length(args)>1L) suites <- args[-1L]
  results <- list(); summary <- list(); issues <- list()
  for(suite in suites) {
    result <- testthat::test_file(paste0('tests/testthat/test-',suite,'.R'),reporter='summary')
    results[[suite]] <- result
    d <- as.data.frame(result)
    summary[[suite]] <- data.frame(profile=profile,suite=suite,failed=sum(d$failed),
      errors=sum(d$error),warnings=sum(d$warning),skipped=sum(d$skipped),passed=sum(d$passed))
    for(i in seq_len(nrow(d))) for(expectation in d$result[[i]]) {
      type <- if(inherits(expectation,'expectation_failure')) 'failure' else
        if(inherits(expectation,'expectation_error')) 'error' else
        if(inherits(expectation,'expectation_warning')) 'warning' else NULL
      if(!is.null(type)) issues[[length(issues)+1L]] <- data.frame(profile=profile,suite=suite,
        test=d$test[[i]],type=type,message=conditionMessage(expectation))
    }
    saveRDS(results,paste0('dev/newton_polish_bugfix_tests_',profile,'.rds'))
    write.csv(do.call(rbind,summary),paste0('dev/newton_polish_bugfix_tests_',profile,'.csv'),row.names=FALSE)
    if(length(issues)) write.csv(do.call(rbind,issues),paste0('dev/newton_polish_bugfix_issues_',profile,'.csv'),row.names=FALSE)
  }
  print(do.call(rbind,summary))
}
if(profile=='phase1') {
rows <- list()
for(profile in c('default','cpu_float64')) for(method in c('mirror','lbfgs','ccp_accelerated')) {
  set.seed(804)
  x <- matrix(rnorm(96),32,3); y <- matrix(rnorm(144),48,3)
  opts <- cotOptions(lambda=5,solver=if(method=='lbfgs') 'lbfgs' else 'mirror',
    polish='none',niter=100L,tol=1e-4,nboot=0,cost.online='tensorized',
    dtype=if(profile=='default') NULL else torch::torch_double(),
    device=if(profile=='default') NULL else torch::torch_device('cpu'),
    line_search_fn='strong_wolfe')
  cot <- COT$new(source=x,target=y,options=opts)
  p <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
  unlockBinding('polish_result',p)
  p$polish_result <- function(res,tol,delta) res
  lockBinding('polish_result',p)
  if(method=='ccp_accelerated') p$ccp_step_rule <- 'accelerated'
  cot$solve()
  rows[[length(rows)+1L]] <- data.frame(profile=profile,method=method,
    iterations=as.numeric(p$iterations_run[[1]][[1]]),loss=as.numeric(p$final_loss[[1]][[1]]))
}
  out<-do.call(rbind,rows);print(out)
  saveRDS(out,if(length(args)>=3L) args[[3L]] else 'dev/newton_polish_bugfix_phase1.rds')
} else if(profile=='regression') {
  set.seed(1); n<-256L; x<-matrix(rnorm(n*3),n,3)
  z<-rbinom(n,1,plogis(.8*x[,1]-.5*x[,2]+.3*x[,3])); scale<-median(as.numeric(dist(x))^2)
  rows<-list()
  for(configuration in c('cpu_float64','cpu_float32','default_float32')) for(mult in c(1,5,25)) for(mode in c('none','auto')) {
    opts<-cotOptions(lambda=mult*scale,solver='mirror',polish=mode,niter=100L,nboot=0,
      cost.online='tensorized',device=if(configuration=='default_float32') NULL else torch::torch_device('cpu'),
      dtype=if(configuration=='cpu_float64') torch::torch_double() else torch::torch_float32())
    cot<-COT$new(source=x[z==0,,drop=FALSE],target=x,options=opts);cot$solve()
    optimizer<-cot$.__enclos_env__$private$optimizer; info<-optimizer$info(); p<-info$polish[[1]][[1]]
    status<-p$certification$status
    expected<-if(mode=='none') 134L else c('1'=39L,'5'=34L,'25'=31L)[[as.character(mult)]]
    stopifnot(sum(cot$weights>0)==expected)
    if(mode=='none') stopifnot(status=='uncertified',p$certification$source=='phase1') else
      stopifnot(status %in% c('certified','certified_resolution'),p$divergence>=0,
        configuration=='cpu_float64' || p$cost_source=='data_float64')
    rows[[length(rows)+1L]]<-data.frame(configuration=configuration,
      device=as.character(optimizer$device),mult=mult,polish=mode,status=status,
      refinement=p$status,gap=p$fw_gap,divergence=p$divergence,nnz=sum(cot$weights>0),
      cost_source=p$cost_source,certificate_source=p$certification$source)
  }
  out<-do.call(rbind,rows);print(out)
  write.csv(out,'dev/newton_polish_bugfix_regression.csv',row.names=FALSE)
} else run_suites(profile)
