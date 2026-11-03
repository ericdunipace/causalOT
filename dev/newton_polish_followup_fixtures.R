suppressMessages(pkgload::load_all('.',quiet=TRUE))
torch::torch_set_num_threads(1L)
set.seed(1); n<-256L; x<-matrix(rnorm(n*3),n,3)
z<-rbinom(n,1,plogis(.8*x[,1]-.5*x[,2]+.3*x[,3]))
x<-x[z==0,,drop=FALSE]; set.seed(1); y<-matrix(rnorm(n*3),n,3)
C<-function(x,y) pmax(outer(rowSums(x*x),rowSums(y*y),'+')/2-tcrossprod(x,y),0)
xy<-C(x,y); xx<-C(x,x); yy<-C(y,y); b<-rep(1/n,n)
rows<-list()
for(mult in c(.05,.25,1,5,25)) {
  fixture<-readRDS(sprintf('dev/polish_proto/warmstart_mult%s.rds',mult))
  for(name in names(fixture$runs)) {
    cat("fixture",mult,name,"\n")
    a<-fixture$runs[[name]]$w
    run<-function(warm) cot_newton_polish(a,b,xy,xx,fixture$eps,C_yy=yy,
      control=cot_polish_control(warm_start=warm,sinkhorn_niter=100000L))
    cold<-run(FALSE); warm<-run(TRUE)
    difference<-abs(warm$loss-cold$loss)
    support_match<-identical(warm$weights>1e-12,cold$weights>1e-12)
    stopifnot(difference<=1e-12*max(1,abs(cold$loss)),support_match)
    rows[[length(rows)+1L]]<-data.frame(mult=mult,fixture=name,difference=difference,
      support_match=support_match,cold_iterations=sum(cold$probe_iterations,na.rm=TRUE),
      warm_iterations=sum(warm$probe_iterations,na.rm=TRUE),cold_time=cold$time,warm_time=warm$time)
    print(tail(rows,1)[[1]])
    write.csv(do.call(rbind,rows),'dev/newton_polish_followup_fixtures.csv',row.names=FALSE)
  }
}
