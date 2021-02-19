testthat::test_that("cost mahalanobis 2.0", {
  n0 <- 100
  n1 <- 55
  n <- n0 + n1
  d <- 5
  
  x1 <- matrix(rnorm(n1*d), n1, d)
  x0 <- matrix(rnorm(n0*d), n0, d)
  
  cov_mat <- cov(rbind(x1, x0)) #(n1 / n * cov(x1) + n0 / n * cov(x0))
  U <- solve(chol(cov_mat))
  
  mhdefault  <- matrix(NA, n0,n1)
  for(i in 1:n0){
    for(j in 1:n1) {
      mhdefault[i,j] <- sqrt(sum(((x0[i,]- x1[j,])%*% U )^2))
    }
  }
  mhown <- cost_mahalanobis(x0,x1, 2, "rowwise", estimand = "ATT")
  testthat::expect_equal(mhdefault, mhown)
})
testthat::test_that("cost mahalanobis 1.0", {
  n0 <- 100
  n1 <- 55
  d <- 5
  n <- n0 + n1
  
  x1 <- matrix(rnorm(n1*d), n1, d)
  x0 <- matrix(rnorm(n0*d), n0, d)
  
  cov_mat <- cov(rbind(x1, x0)) #(n1/n * cov(x1) + n0/n *cov(x0))
  U <- inv_sqrt_mat(cov_mat)
  
  mhdefault  <- matrix(NA, n0,n1)
  for(i in 1:n0){
    for(j in 1:n1) {
      mhdefault[i,j] <- sum(abs((x0[i,]- x1[j,])%*% U))
    }
  }
  mhown <- cost_mahalanobis(x0,x1, 1, "rowwise", estimand = "ATT")
  testthat::expect_equal(mhdefault, mhown)
})
testthat::test_that("cost mahalanobis 3.0", {
  n0 <- 100
  n1 <- 55
  n <- n0 + n1
  d <- 5
  
  x1 <- matrix(rnorm(n1*d), n1, d)
  x0 <- matrix(rnorm(n0*d), n0, d)
  
  cov_mat <- cov(rbind(x1, x0)) #(n1 / n * cov(x1) + n0 / n * cov(x0))
  U <- inv_sqrt_mat(cov_mat)
  
  mhdefault  <- matrix(NA, n0,n1)
  for(i in 1:n0){
    for(j in 1:n1) {
      mhdefault[i,j] <- sum(abs((x0[i,]- x1[j,])%*% U)^3)^(1/3)
    }
  }
  mhown <- cost_mahalanobis(x0,x1, 3, "rowwise", estimand = "ATT")
  testthat::expect_equal(mhdefault, mhown)
})
