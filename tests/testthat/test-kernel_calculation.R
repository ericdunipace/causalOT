testthat::test_that("confirm kernel functions correct, LP", {
  set.seed(239402398)
  p <- 1.0
  theta <- c(1.0,1.0)
  gamma <- c(1.0,1.0)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  
  K <- causalOT:::kernel_calculation(X = x, z = z, p = p, theta = theta, kernel = "polynomial",
                                    gamma = gamma, metric = "Lp", is.dose = TRUE)
  
  mu_x <- colMeans(x)
  mu_z <- mean(z)
  
  a <- x - matrix(mu_x, n, d, byrow=TRUE)
  b <- as.matrix(z - mu_z)
  
  x_mat <- tcrossprod(a)
  z_mat <- tcrossprod(b)
  
  K_temp <- gamma[1] * (1 + theta[1] * x_mat)^p * gamma[2] * (1 + theta[2] * z_mat)^p
  
  testthat::expect_equal(K$cov_kernel, K_temp)
  
})

testthat::test_that("confirm kernel functions correct, mahalanobis", {
  set.seed(0877)
  p <- 1.0
  theta <- c(1.0,1.0)
  gamma <- c(1.0,1.0)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  
  K <- causalOT:::kernel_calculation(X = x, z = z, p = p, theta = theta, 
                                    kernel = "polynomial", gamma = gamma, metric = "mahalanobis", is.dose = TRUE)
  
  mu_x <- colMeans(x)
  mu_z <- mean(z)
  
  a <- backsolve(chol(cov(x)), t(x - matrix(mu_x, n, d, byrow=TRUE)) , transpose = TRUE)
  b <- 1/sd(z) * as.matrix(z - mu_z)
  
  x_mat <- crossprod(a)
  z_mat <- tcrossprod(b)
  
  K_temp <- gamma[1] * (1 + theta[1] * x_mat)^p * gamma[2] * (1 + theta[2] * z_mat)^p
  
  testthat::expect_equal(K$cov_kernel, K_temp)
  
})

testthat::test_that("different values of parameters, Lp", {
  set.seed(0988)
  vals <- seq(1,10, length.out = 10)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  
  mu_x <- colMeans(x)
  mu_z <- mean(z)
  
  a <- x - matrix(mu_x, n, d, byrow=TRUE)
  b <- as.matrix(z - mu_z)
  
  x_mat <- tcrossprod(a)
  z_mat <- tcrossprod(b)
  
  
  for( v in vals) {
    p <- v
    theta <- c(v,v)
    gamma <- c(v,v)
    
    K <- causalOT:::kernel_calculation(X = x, z = z, p = p, kernel = "polynomial", theta = theta, gamma = gamma, metric = "Lp", is.dose = TRUE)
    
    K_temp <- gamma[1] * (1 + theta[1] * x_mat)^p * gamma[2] * (1 + theta[2] * z_mat)^p
    
    testthat::expect_equal(K$cov_kernel, K_temp)
  }
  
})

testthat::test_that("different values of parameters, mahalanobis", {
  set.seed(9808)
  vals <- seq(1,10, length.out = 10)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  mu_x <- colMeans(x)
  mu_z <- mean(z)
  
  a <- backsolve(chol(cov(x)), t(x - matrix(mu_x, n, d, byrow=TRUE)) , transpose = TRUE)
  b <- 1/sd(z) * as.matrix(z - mu_z)
  
  x_mat <- crossprod(a)
  z_mat <- tcrossprod(b)
  
  
  for( v in vals) {
    p <- v
    theta <- c(v,v)
    gamma <- c(v,v)
  
  
    K <- causalOT:::kernel_calculation(X = x, z = z, p = p, theta = theta, kernel = "polynomial",
                                      gamma = gamma, metric = "mahalanobis", is.dose = TRUE)

    K_temp <- gamma[1] * (1 + theta[1] * x_mat)^p * gamma[2] * (1 + theta[2] * z_mat)^p
    
    testthat::expect_equal(K$cov_kernel, K_temp)
  }
  
})
