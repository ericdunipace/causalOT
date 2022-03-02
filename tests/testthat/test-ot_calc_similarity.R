testthat::test_that("confirm ot kernel functions correct dose, LP", {
  set.seed(239402398)
  p <- 1.0
  theta <- c(1.0,1.0)
  gamma <- c(1.0,1.0)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  
  K <- causalOT:::ot_kernel_calculation(X = x, z = z, p = p, 
                                       theta = theta, gamma = gamma, 
                                       kernel = "polynomial",
                                       metric = "Lp", is.dose = TRUE)
  
  mu_x <- colMeans(x)
  mu_z <- mean(z)
  
  a <- x - matrix(mu_x, n, d, byrow=TRUE)
  b <- as.matrix(z - mu_z)
  
  x_mat <- tcrossprod(a)
  z_mat <- tcrossprod(b)
  
  K_temp <- gamma[1] * (1 + theta[1] * x_mat)^p * gamma[2] * (1 + theta[2] * z_mat)^p
  
  testthat::expect_equal(K, K_temp)
  
})

testthat::test_that("confirm ot kernel functions correct dose, mahalanobis", {
  set.seed(0877)
  p <- 1.0
  theta <- c(1.0,1.0)
  gamma <- c(1.0,1.0)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  
  K <- causalOT:::kernel_calculation(X = x, z = z, p = p, theta = theta, gamma = gamma, 
                                    kernel = "polynomial",
                                    metric = "mahalanobis", is.dose = TRUE)
  
  mu_x <- colMeans(x)
  mu_z <- mean(z)
  
  a <- backsolve(chol(cov(x)), t(x - matrix(mu_x, n, d, byrow=TRUE)) , transpose = TRUE)
  b <- 1/sd(z) * as.matrix(z - mu_z)
  
  x_mat <- crossprod(a)
  z_mat <- tcrossprod(b)
  
  K_temp <- gamma[1] * (1 + theta[1] * x_mat)^p * gamma[2] * (1 + theta[2] * z_mat)^p
  
  testthat::expect_equal(K$cov_kernel, K_temp)
  
})

testthat::test_that("confirm ot kernel functions correct, LP", {
  set.seed(239402398)
  p <- 1.0
  theta <- c(1.0,1.0)
  gamma <- c(1.0,1.0)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  
  K <- causalOT:::ot_kernel_calculation(X = x, z = z, p = p, theta = theta, gamma = gamma,
                                       kernel = "polynomial",
                                       metric = "Lp", is.dose = FALSE)
  
  orders <- order(z)
  x <- x[orders,]
  z <- z[orders]
  
  # mu_x <- colMeans(x)
  
  a <- x #- matrix(mu_x, n, d, byrow=TRUE)
  
  x_mat <- tcrossprod(a)
  
  K_temp <- list(gamma[1] * (1 + theta[1] * x_mat[z==0,])^p ,
                 gamma[2] * (1 + theta[2] * x_mat[z==1,])^p)
  
  testthat::expect_equal(K, K_temp)
  
  #ATT
  K <- causalOT:::ot_kernel_calculation(estimand = "ATT", X = x, z = z, p = p, theta = theta, 
                                       gamma = gamma, metric = "Lp", 
                                       kernel = "polynomial",
                                       is.dose = FALSE)
  
  orders <- order(z)
  x <- x[orders,]
  z <- z[orders]
  
  mu_x <- colMeans(x[z==1,])
  
  a <- x# - matrix(mu_x, n, d, byrow=TRUE)
  
  x_mat <- tcrossprod(a)
  
  K_temp <- list(gamma[2] * (1 + theta[2] * x_mat[z==0,z==1])^p ,
                 NULL)
  
  testthat::expect_equal(K, K_temp)
  
  
  #ATC
  K <- causalOT:::ot_kernel_calculation(estimand = "ATC", X = x, z = z, p = p, theta = theta, 
                                       gamma = gamma, metric = "Lp", 
                                       kernel = "polynomial",
                                       is.dose = FALSE)
  
  orders <- order(z)
  x <- x[orders,]
  z <- z[orders]
  
  mu_x <- colMeans(x[z==0,])
  
  a <- x #- matrix(mu_x, n, d, byrow=TRUE)
  
  x_mat <- tcrossprod(a)
  
  K_temp <- list(gamma[1] * (1 + theta[1] * x_mat[z==0,z==1])^p ,
                 NULL)
  
  testthat::expect_equal(K, K_temp)
  
})

testthat::test_that("confirm ot kernel functions correct, mahalanobis", {
  set.seed(239402398)
  p <- 1.0
  theta <- c(1.0,1.0)
  gamma <- c(1.0,1.0)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  
  K <- causalOT:::ot_kernel_calculation(estimand = "ATE", X = x, z = z, p = p, theta = theta, 
                                       gamma = gamma, metric = "mahalanobis", 
                                       kernel = "polynomial",
                                       is.dose = FALSE)
  
  orders <- order(z)
  x <- x[orders,]
  z <- z[orders]
  
  mu_x <- colMeans(x)
  
  a <- x - matrix(mu_x, n, d, byrow=TRUE)
  a <- a %*% solve(chol(cov(x)))
  
  x_mat <- tcrossprod(a)
  
  K_temp <- list(gamma[1] * (1 + theta[1] * x_mat[z==0,])^p ,
                 gamma[2] * (1 + theta[2] * x_mat[z==1,])^p)
  
  testthat::expect_equal(K, K_temp)
  
  #ATT
  K <- causalOT:::ot_kernel_calculation(estimand = "ATT", X = x, z = z, p = p, theta = theta, 
                                       gamma = gamma, metric = "mahalanobis", 
                                       kernel = "polynomial",
                                       is.dose = FALSE)
  
  orders <- order(z)
  x <- x[orders,]
  z <- z[orders]
  
  mu_x <- colMeans(x)
  
  a <- x - matrix(mu_x, n, d, byrow=TRUE)
  a <- a %*% solve(chol(cov(x)))
  
  x_mat <- tcrossprod(a)
  
  K_temp <- list(gamma[2] * (1 + theta[2] * x_mat[z==0,z==1])^p ,
                 NULL)
  
  testthat::expect_equal(K, K_temp)
  
  
  #ATC
  K <- causalOT:::ot_kernel_calculation(estimand = "ATC", X = x, z = z, p = p, theta = theta, 
                                       gamma = gamma, metric = "mahalanobis", 
                                       kernel = "polynomial",
                                       is.dose = FALSE)
  
  orders <- order(z)
  x <- x[orders,]
  z <- z[orders]
  
  mu_x <- colMeans(x)
  
  a <- x - matrix(mu_x, n, d, byrow=TRUE)
  a <- a %*% solve(chol(cov(x)))
  
  x_mat <- tcrossprod(a)
  
  K_temp <- list(gamma[1] * (1 + theta[1] * x_mat[z==0,z==1])^p ,
                 NULL)
  
  testthat::expect_equal(K, K_temp)
  
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
    
    K <- causalOT:::ot_kernel_calculation(X = x, z = z, p = p, theta = theta, gamma = gamma, 
                                         metric = "Lp", is.dose = TRUE,
                                         kernel = "polynomial")
    
    K_temp <- gamma[1] * (1 + theta[1] * x_mat)^p * gamma[2] * (1 + theta[2] * z_mat)^p
    
    testthat::expect_equal(K, K_temp)
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
    
    
    K <- causalOT:::ot_kernel_calculation(X = x, z = z, p = p, theta = theta, gamma = gamma, 
                                         metric = "mahalanobis", is.dose = TRUE,
                                         kernel = "polynomial")
    
    K_temp <- gamma[1] * (1 + theta[1] * x_mat)^p * gamma[2] * (1 + theta[2] * z_mat)^p
    
    testthat::expect_equal(K, 
                           K_temp)
  }
  
})

testthat::test_that("different values of parameters not dose, LP", {
  set.seed(239402398)
  vals <- seq(1,10, length.out = 10)
  p <- 1.0
  theta <- c(1.0,1.0)
  gamma <- c(1.0,1.0)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  
  for( v in vals) {
    p <- v
    theta <- c(v,v)
    gamma <- c(v,v)
    
    
    K <- causalOT:::ot_kernel_calculation(estimand = "ATE", X = x, z = z, p = p, theta = theta, 
                                         gamma = gamma, metric = "Lp", is.dose = FALSE,
                                         kernel = "polynomial")
    
    orders <- order(z)
    x <- x[orders,]
    z <- z[orders]
    
    mu_x <- colMeans(x)
    
    a <- x #- matrix(mu_x, n, d, byrow=TRUE)

    x_mat <- tcrossprod(a)
    
    K_temp <- list(gamma[1] * (1 + theta[1] * x_mat[z==0,])^p ,
                   gamma[2] * (1 + theta[2] * x_mat[z==1,])^p)
    
    testthat::expect_equal(K, K_temp)
  }
  
  
  #ATT
  for( v in vals) {
    p <- v
    theta <- c(v,v)
    gamma <- c(v,v)
    
    
    K <- causalOT:::ot_kernel_calculation(estimand = "ATT", X = x, z = z, p = p, theta = theta, 
                                         gamma = gamma, metric = "Lp", is.dose = FALSE,
                                         kernel = "polynomial")
    
    orders <- order(z)
    x <- x[orders,]
    z <- z[orders]
    
    mu_x <- colMeans(x[z==1,])
    
    a <- x #- matrix(mu_x, n, d, byrow=TRUE)

    x_mat <- tcrossprod(a)
    
    K_temp <- list(gamma[2] * (1 + theta[2] * x_mat[z==0,z==1])^p ,
                   NULL)
    
    testthat::expect_equal(K, K_temp)
  }
  
  
  
  #ATC
  for( v in vals) {
    p <- v
    theta <- c(v,v)
    gamma <- c(v,v)
    K <- causalOT:::ot_kernel_calculation(estimand = "ATC", X = x, z = z, p = p, theta = theta, 
                                         gamma = gamma, metric = "Lp", is.dose = FALSE,
                                         kernel = "polynomial")
    
    orders <- order(z)
    x <- x[orders,]
    z <- z[orders]
    
    mu_x <- colMeans(x[z==0,])
    
    a <- x #- matrix(mu_x, n, d, byrow=TRUE)

    x_mat <- tcrossprod(a)
    
    K_temp <- list(gamma[1] * (1 + theta[1] * x_mat[z==0,z==1])^p ,
                   NULL)
    
    testthat::expect_equal(K, K_temp)
  }
  
})

testthat::test_that("different values of parameters not dose, mahalanobis", {
  set.seed(239402398)
  vals <- seq(1,10, length.out = 10)
  p <- 1.0
  theta <- c(1.0,1.0)
  gamma <- c(1.0,1.0)
  n <- 100
  d <- 10
  
  x <- matrix(rnorm(n*d),n,d)
  z <- rbinom(n, 1, 0.5)
  
  for( v in vals) {
    p <- v
    theta <- c(v,v)
    gamma <- c(v,v)
    
    
    K <- causalOT:::ot_kernel_calculation(estimand = "ATE", X = x, z = z, p = p, theta = theta, 
                                         gamma = gamma, metric = "mahalanobis", is.dose = FALSE,
                                         kernel = "polynomial")
    
    orders <- order(z)
    x <- x[orders,]
    z <- z[orders]
    
    mu_x <- colMeans(x)
    
    a <- x - matrix(mu_x, n, d, byrow=TRUE)
    a <- a %*% solve(chol(cov(x)))
    
    x_mat <- tcrossprod(a)
    
    K_temp <- list(gamma[1] * (1 + theta[1] * x_mat[z==0,])^p ,
                   gamma[2] * (1 + theta[2] * x_mat[z==1,])^p)
    
    testthat::expect_equal(K, K_temp)
  }
  
  
  #ATT
  for( v in vals) {
    p <- v
    theta <- c(v,v)
    gamma <- c(v,v)
    
    
    K <- causalOT:::ot_kernel_calculation(estimand = "ATT", X = x, z = z, p = p, theta = theta, 
                                         gamma = gamma, metric = "mahalanobis", is.dose = FALSE,
                                         kernel = "polynomial")
    
    orders <- order(z)
    x <- x[orders,]
    z <- z[orders]
    
    mu_x <- colMeans(x)
    
    a <- x - matrix(mu_x, n, d, byrow=TRUE)
    a <- a %*% solve(chol(cov(x)))
    
    x_mat <- tcrossprod(a)
    
    K_temp <- list(gamma[2] * (1 + theta[2] * x_mat[z==0,z==1])^p ,
                   NULL)
    
    testthat::expect_equal(K, K_temp)
  }
  
  
  
  #ATC
  for( v in vals) {
    p <- v
    theta <- c(v,v)
    gamma <- c(v,v)
    K <- causalOT:::ot_kernel_calculation(estimand = "ATC", X = x, z = z, p = p, theta = theta, 
                                         gamma = gamma, metric = "mahalanobis", is.dose = FALSE,
                                         kernel = "polynomial")
    
    orders <- order(z)
    x <- x[orders,]
    z <- z[orders]
    
    mu_x <- colMeans(x)
    
    a <- x - matrix(mu_x, n, d, byrow=TRUE)
    a <- a %*% solve(chol(cov(x)))
    
    x_mat <- tcrossprod(a)
    
    K_temp <- list(gamma[1] * (1 + theta[1] * x_mat[z==0,z==1])^p ,
                   NULL)
    
    testthat::expect_equal(K, K_temp)
  }
  
})
