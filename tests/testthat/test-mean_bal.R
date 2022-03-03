testthat::test_that("mean_bal works", {
  set.seed(23324)
  library(causalOT)
  n <- 100
  p <- 6
  x0 <- matrix(rnorm(2*n * p), 2*n, p)
  x1 <- matrix(rnorm(n * p), n, p)
  weights <- list(w0 = rep(1/(2*n), 2 * n), w1 = rep(1/n, n))
  
  data <- cbind(rbind(x0,x1), z = c(rep(0,2*n), rep(1, n)))
  colnames(data) <- c(paste0("x", 1:p), "z")
  
  v1 <- matrixStats::colVars(x0)
  v2 <- matrixStats::colVars(x1)
  
  poolsd <- sqrt(2*n/(3*n) *v1 + n/(3*n) * v2)
  
  mb <- mean_bal(data, weights, balance.covariates = paste0("x", 1:p), 
           treatment.indicator = "z")
  
  wc <- weights
  class(wc) <- "causalWeights"
  
  mb.c <- mean_bal(data, wc, balance.covariates = paste0("x", 1:p), 
                 treatment.indicator = "z")
  
  testthat::expect_equivalent(mb, abs(c(crossprod(x1, weights$w1) - crossprod(x0, weights$w0) )) / poolsd)
  testthat::expect_equal(mb, mb.c)
})
