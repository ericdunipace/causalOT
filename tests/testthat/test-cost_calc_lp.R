testthat::test_that("cost lp 2.0", {
  n0 <- 100
  n1 <- 55
  d <- 5
  
  x1 <- matrix(rnorm(n1*d), n1, d)
  x0 <- matrix(rnorm(n0*d), n0, d)
  
  mhdefault  <- matrix(NA, n0,n1)
  for(i in 1:n0){
    for(j in 1:n1) {
      mhdefault[i,j] <- sqrt(sum(((x0[i,]- x1[j,]) )^2))
    }
  }
  mhown <- causalOT:::cost_calc_lp(x0,x1, 2, "rowwise")
  testthat::expect_equal(mhdefault, mhown)
})
testthat::test_that("cost lp 1.0", {
  n0 <- 100
  n1 <- 55
  d <- 5
  
  x1 <- matrix(rnorm(n1*d), n1, d)
  x0 <- matrix(rnorm(n0*d), n0, d)
  
  mhdefault  <- matrix(NA, n0,n1)
  for(i in 1:n0){
    for(j in 1:n1) {
      mhdefault[i,j] <- sum(abs((x0[i,]- x1[j,])))
    }
  }
  mhown <- causalOT:::cost_calc_lp(x0,x1, 1, "rowwise")
  testthat::expect_equal(mhdefault, mhown)
})
testthat::test_that("cost lp 3.0", {
  n0 <- 100
  n1 <- 55
  d <- 5
  
  x1 <- matrix(rnorm(n1*d), n1, d)
  x0 <- matrix(rnorm(n0*d), n0, d)
  
  mhdefault  <- matrix(NA, n0,n1)
  for(i in 1:n0){
    for(j in 1:n1) {
      mhdefault[i,j] <- sum(abs((x0[i,]- x1[j,]))^3)^(1/3)
    }
  }
  mhown <- causalOT:::cost_calc_lp(x0,x1, 3, "rowwise")
  testthat::expect_equal(mhdefault, mhown)
})