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
