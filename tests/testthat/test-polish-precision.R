precision_fixture <- function() {
  causalOT:::torch_check()
  torch::torch_set_num_threads(1L)
  set.seed(42)
  x <- matrix(rnorm(24),12,2); y <- matrix(rnorm(30),15,2)
  cost <- function(x,y) pmax(outer(rowSums(x^2),rowSums(y^2),'+')/2-tcrossprod(x,y),0)
  list(xy=cost(x,y),xx=cost(x,x),yy=cost(y,y),a=rep(1/12,12),b=rep(1/15,15),lambda=.22916026113137714)
}
precision_residual <- function(a,b,f,xy,xx,lambda) {
  plan <- outer(a,b)*exp((outer(as.numeric(f$f_xy),as.numeric(f$g_yx),'+')-xy)/lambda)
  self <- -lambda * matrixStats::rowLogSumExps(sweep(-xx,2,as.numeric(f$f_xx)+lambda*log(a),'+')/lambda)
  c(marginal=sum(abs(rowSums(plan)-a))+sum(abs(colSums(plan)-b)),
    self=max(abs(self-as.numeric(f$f_xx))))
}

test_that('blocked softmin retains the scalar precision of a non-float32 lambda', {
  f <- precision_fixture(); h <- seq(-.5,.5,length.out=length(f$b))
  cost <- causalOT:::cot_polish_cost(f$xy,block_size=4L)
  expect_null(cost$dense)
  actual <- as.numeric(causalOT:::cot_polish_softmin(cost,causalOT:::cot_polish_tensor(h),f$lambda))
  expected <- -f$lambda*matrixStats::rowLogSumExps(sweep(-f$xy,2,h,'+')/f$lambda)
  expect_lt(max(abs(actual-expected)),1e-14)
})

test_that('blocked cold, warm and target self agree with float64 dense evaluations', {
  f <- precision_fixture()
  blocked <- causalOT:::cot_polish_control(block_size=4L,sinkhorn_niter=1000L)
  dense <- causalOT:::cot_polish_control(sinkhorn_niter=1000L)
  xy <- causalOT:::cot_polish_cost(f$xy,4L); xx <- causalOT:::cot_polish_cost(f$xx,4L)
  cold <- causalOT:::cot_polish_eval(f$a,f$b,xy,xx,f$lambda,blocked)
  ref <- causalOT:::cot_polish_eval(f$a,f$b,f$xy,f$xx,f$lambda,dense)
  expect_equal(cold$loss,ref$loss,tolerance=1e-12)
  expect_equal(cold$grad,ref$grad,tolerance=1e-11)
  expect_equal(causalOT:::cot_polish_self(f$b,f$yy,f$lambda,blocked),
    causalOT:::cot_polish_self(f$b,f$yy,f$lambda,dense),tolerance=1e-12)
  unchanged <- causalOT:::cot_polish_warm(f$a,f$b,xy,xx,f$lambda,cold$pot,blocked)
  expect_lte(unchanged$iterations,5L)
  expect_lte(precision_residual(f$a,f$b,unchanged$pot,f$xy,f$xx,f$lambda)[['marginal']],1.2e-13)
  expect_lte(precision_residual(f$a,f$b,unchanged$pot,f$xy,f$xx,f$lambda)[['self']],1.2e-13)
  a <- f$a; a[c(1,4,7)] <- 0; a <- a/sum(a)
  changed <- causalOT:::cot_polish_eval(a,f$b,xy,xx,f$lambda,blocked)
  for(start in list(cold$pot,changed$pot)) {
    warm <- causalOT:::cot_polish_warm(a,f$b,xy,xx,f$lambda,start,blocked)
    expect_lt(warm$iterations,blocked$sinkhorn_niter)
    expect_equal(warm$loss,changed$loss,tolerance=1e-12)
    expect_equal(warm$grad,changed$grad,tolerance=1e-11)
    residual <- precision_residual(a,f$b,warm$pot,f$xy,f$xx,f$lambda)
    expect_lte(residual[['marginal']],1.2e-13)
    expect_lte(residual[['self']],1.2e-13)
  }
  # Independent matrix formula also checks scalar precision in the Hessian.
  K <- exp((outer(as.numeric(ref$pot$f_xy),as.numeric(ref$pot$g_yx),'+')-f$xy)/f$lambda)
  Kxx <- exp((outer(as.numeric(ref$pot$f_xx),as.numeric(ref$pot$f_xx),'+')-f$xx)/f$lambda)
  R <- K %*% (f$b*t(K)); A <- diag(f$a); I <- diag(length(f$a))
  H <- f$lambda*(solve(I-R%*%A+outer(rep(1,length(f$a)),f$a),R)+solve(I+Kxx%*%A,Kxx))
  H <- (H+t(H))/2
  actual <- as.matrix(causalOT:::cot_polish_hessian(f$a,f$b,xy,xx,f$lambda,ref$pot,block_size=4L))
  expect_lt(max(abs(actual-H)),1e-11)
})

test_that('complete blocked polish agrees with cold-only and dense polishing', {
  f <- precision_fixture()
  run <- function(block,warm) causalOT:::cot_newton_polish(f$a,f$b,f$xy,f$xx,f$lambda,C_yy=f$yy,
    control=causalOT:::cot_polish_control(block_size=block,warm_start=warm,sinkhorn_niter=1000L,maxit=40L))
  warm <- run(4L,TRUE); cold <- run(4L,FALSE); dense <- run(512L,FALSE)
  for(ans in list(warm,cold)) {
    expect_true(causalOT:::cot_polish_is_converged(ans$status))
    expect_lte(ans$loss,ans$loss_start)
    expect_identical(ans$weights>1e-8,dense$weights>1e-8)
    expect_equal(ans$loss,dense$loss,tolerance=1e-12)
    expect_equal(ans$weights,dense$weights,tolerance=1e-6)
  }
  expect_equal(warm$loss,cold$loss,tolerance=1e-12)
})
