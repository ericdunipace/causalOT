#check large ot
Rcpp::sourceCpp('src/large_matrix_entry.cpp')
n <- 1e4
if(n %% 1e4) stop("n must be multiple of 10,000")
p <- 5

x0 <- matrix(rnorm(n*p), p, n) #observations are column-wise
x1 <- matrix(rnorm(n,p), p,n)
iter_seq <- 0:(n/1e4 - 1)

cost <- matrix(0, n,n)
tracemem(cost)
for (i in iter_seq) {
  entry(cost , 
   causalOT::cost_calc_lp(x0, x1[,(i * 1e4) + 1:1e4, drop = FALSE], 2, "colwise"),
   (i * 1e4) + 1)
}

microbenchmark::microbenchmark(transport::transport(a = rep(1/n,n), b = rep(1/n,n), costm = cost, p = 1),
                               times = 10)

set.seed(234897)
n <- 1e7
max <- 5e5
p <- 8
x0 <- matrix(rnorm(n*p), p, n) #observations are column-wise
x1 <- matrix(rnorm(n,p), p, n)

if(n <= max) tplan <- approxOT::transport_plan(x0, x1, 2, 2, "colwise", "hilbert")$tplan
if(n > max) {
  tplan <- list()
  tplan$from <-  approxOT::hilbert_proj_(x1) + 1
  tplan$to <- approxOT::hilbert_proj_(x0) + 1
  tplan$mass <- rep(1/n,n)
}
hilbert.dist <- sqrt(weighted.mean(colSums((x0[,tplan$to] - x1[,tplan$from])^2), w = tplan$mass))


# swap.dist <- approxOT::transport_plan(x0, x1, 2, 2, "colwise", "swapping") # too slow


set.seed(234897)
n <- 1e7
max <- 5e5
p <- 8
x0 <- matrix(rnorm(n*p), p, n) #observations are column-wise
x1 <- matrix(rnorm(n,p), p, n)
tplan_full <- list()
tplan_full$z <- c(rep(0, n), rep(1,n))
tplan_full$idx <- approxOT::hilbert_proj_(cbind(x0,x1)) + 1

