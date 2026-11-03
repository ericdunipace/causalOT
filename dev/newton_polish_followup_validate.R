# Focused, paired regression and performance validation. Run from package root.
# POLISH_FOLLOWUP_NS=256,1000,2000,5000 enables the deferred scaling sizes.
suppressMessages(pkgload::load_all('.', quiet=TRUE))
torch::torch_set_num_threads(1L)
rows <- list()
set.seed(1); n <- 256L
x <- matrix(rnorm(n*3),n,3)
z <- rbinom(n,1,plogis(.8*x[,1]-.5*x[,2]+.3*x[,3]))
scale <- median(as.numeric(dist(x))^2)
for (dtype in c('float64','float32')) for (mult in c(1,5,25)) for (mode in c('none','auto')) {
  opts <- cotOptions(lambda=mult*scale,solver='mirror',polish=mode,niter=100L,nboot=0,
    cost.online='tensorized',device=torch::torch_device('cpu'),
    dtype=if(dtype=='float64') torch::torch_double() else torch::torch_float32())
  cot <- COT$new(source=x[z==0,,drop=FALSE],target=x,options=opts); cot$solve()
  optimizer <- cot$.__enclos_env__$private$optimizer
  info <- optimizer$info(); p <- info$polish[[1]][[1]]
  status <- info$fit_status[[1]][[1]]
  if (mode=='auto') {
    stopifnot(status %in% c('certified','certified_resolution'),p$divergence>=0,
      sum(cot$weights>0)==c('1'=39L,'5'=34L,'25'=31L)[[as.character(mult)]])
  } else stopifnot(status=='uncertified',sum(cot$weights>0)==134L)
  rows[[length(rows)+1L]] <- data.frame(dtype=dtype,mult=mult,polish=mode,status=status,
    gap=p$fw_gap,divergence=p$divergence,nnz=sum(cot$weights>0),cost_source=p$cost_source)
}
regression <- do.call(rbind,rows)
write.csv(regression,'dev/newton_polish_followup_regression.csv',row.names=FALSE)
print(regression)
ns <- as.integer(strsplit(Sys.getenv('POLISH_FOLLOWUP_NS','256,1000,2000'),',',fixed=TRUE)[[1]])
rows <- list()
for(n in ns) {
  set.seed(501+n)
  x <- matrix(rnorm(n*3),n,3); y <- matrix(rnorm(n*3),n,3)
  lambda <- 5*median(as.numeric(dist(x[seq_len(min(n,256L)),,drop=FALSE]))^2)
  opts <- cotOptions(lambda=lambda,solver='mirror',polish='none',niter=50L,tol=1e-10,nboot=0,
    cost.online='tensorized',device=torch::torch_device('cpu'),dtype=torch::torch_float32())
  cot <- COT$new(source=x,target=y,options=opts); cot$solve()
  p <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private
  ot <- p$ot_objects[[ls(p$ot_objects)[[1L]]]]
  addresses <- p$problems[[ls(p$ot_objects)[[1L]]]]
  xx <- p$measures[[addresses[[1]]]]$x_float64; yy <- p$measures[[addresses[[2]]]]$x_float64
  Cxy <- cot_polish_data_cost(ot$C_xy,xx,yy); Cxx <- cot_polish_data_cost(ot$C_xx,xx,xx)
  Cyy <- cot_polish_data_cost(ot$C_yy,yy,yy)
  # Compile before either timed arm.
  invisible(cot_polish_eval(rep(1/4,4),rep(1/4,4),matrix(0,4,4),matrix(0,4,4),1))
  results <- list()
  for(warm in c(FALSE,TRUE)) {
    ans <- cot_newton_polish(cot$weights,ot$b,Cxy,Cxx,lambda,C_yy=Cyy,
      control=cot_polish_control(maxit=1L,warm_start=warm))
    results[[as.character(warm)]] <- ans
    row <- data.frame(n=n,warm=warm,time=ans$time,first_free=ans$free_history[[1]],
      loss=ans$loss,nnz=ans$nnz,line_search_evals=ans$line_search_evals,
      probe_iterations=paste(ans$probe_iterations,collapse=','),as.list(ans$timings))
    rows[[length(rows)+1L]] <- row
    print(row)
    write.csv(do.call(rbind,rows),'dev/newton_polish_followup_scaling.csv',row.names=FALSE)
  }
  stopifnot(abs(results[['TRUE']]$loss-results[['FALSE']]$loss)<=1e-12*max(1,abs(results[['FALSE']]$loss)),
    identical(results[['TRUE']]$weights>1e-12,results[['FALSE']]$weights>1e-12))
  rm(cot,p,ot,Cxy,Cxx,Cyy,xx,yy,results); gc()
}
