testthat::test_that("multiplication works", {
  n0 <- 100
  n1 <- 55
  d <- 5
  
  x1 <- matrix(rnorm(n1*d), n1, d)
  x0 <- matrix(rnorm(n0*d), n0, d)
  
  cov_mat <- 0.5*(cov(x1) + cov(x0))
  inv_cov <- solve(cov_mat)
  U <- chol(inv_cov)
  
  mhdefault  <- matrix(NA, n0,n1)
  for(i in 1:n0){
    for(j in 1:n1) {
      mhdefault[i,j] <- sqrt(sum((U %*% (x0[i,]- x1[j,]))^2))
    }
  }
  mhown <- cost_mahalanobis(x0,x1, 2, "rowwise")
  testthat::expect_equal(mhdefault, mhown)
})
