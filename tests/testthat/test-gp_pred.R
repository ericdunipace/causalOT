gauss_pred_r <- function(x,y,z,param,
             estimand = c("ATE","ATT","ATC")){
  n <- nrow(x)
  n1 <- sum(z)
  n0 <- n - n1
  estimand <- match.arg(estimand)
  
  if(param$metric == "mahalanobis") {
    cov <- cov(x)
    A <- scale(x, center = TRUE, scale = FALSE) %*% solve(chol(cov))
    
  } else {
    A <- x
  }
  A0 <- A[z==0,]
  A1 <- A[z==1,]
  
  if(param$kernel == "polynomial") {
    kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + 
      param$gamma[1] * (1+ param$theta[1] * tcrossprod(A0))^param$p
    kernel_cross0 <- param$gamma[1] * (1+ param$theta[1] * tcrossprod(A0,A))^param$p
    
    kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + 
      param$gamma[2] * (1+ param$theta[2] * tcrossprod(A1))^param$p
    kernel_cross1 <- param$gamma[2] * (1+ param$theta[2] * tcrossprod(A1,A))^param$p
  } else if (param$kernel == "RBF") {
    kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + 
      param$gamma[1] * exp(-0.5 * param$theta[1] * 
                             cost_calc_lp(A0,A0,ground_p = 2, direction = "rowwise" )^2)
    kernel_cross0 <- param$gamma[1] * exp(-0.5 *  param$theta[1] * 
                                            cost_calc_lp(A0,A,ground_p = 2, direction = "rowwise" )^2)
    
    kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + 
      param$gamma[2] * exp(-0.5 * param$theta[2] * 
                             cost_calc_lp(A1,A1,ground_p = 2, direction = "rowwise" )^2)
    kernel_cross1 <- param$gamma[2] * exp(-0.5 *  param$theta[2] * 
                                            cost_calc_lp(A1,A,ground_p = 2, direction = "rowwise" )^2)
  }
  
  pred0 <- crossprod(kernel_cross0, solve(kernel_cov0, y[z==0]))
  pred1 <- crossprod(kernel_cross1, solve(kernel_cov1, y[z==1]))
  
  tau <- if(estimand == "ATE") {
    mean(pred1 - pred0)
  } else if (estimand == "ATT") {
    mean(y[z==1] - pred0[z==1])
  } else if (estimand == "ATC") {
    mean(pred1[z==0] - y[z==0])
  }
  
  return(
    list(tau = tau,
         estimand = estimand,
         pred0 = pred0,
         pred1 = pred1,
         kernels = list(control = list(cov = kernel_cov0,
                                       cross = kernel_cross0),
                        treated = list(cov = kernel_cov1,
                                       cross = kernel_cross1)))
  )
}


testthat::test_that("gaussian process prediction works polynomial, mahalanobis", {
  
  n <- 100
  d <- 10
  x <- matrix(rnorm(n*d),n,d)
  y <- rnorm(n)
  z <- rbinom(n,1,0.5)
  n1 <- sum(z)
  n0 <- n - n1
  
  param <- list(theta = c(1.1,1.22),
                gamma = c(0.5,0.423),
                p = 2,
                sigma_2 = c(1.01, 1.033),
                is.dose = FALSE,
                kernel = "polynomial",
                metric = "mahalanobis",
                is.dose = FALSE,
                is.standardized = FALSE
                
                )
  
  rout <- gauss_pred_r(x,y,z, param)
  
  df <- data.frame(x, y = y, z = z)
  tau <- gp_pred(data = df, weights = NULL,param = param, estimand = "ATE", 
                  balance.covariates = colnames(df)[1:d],
                  treatment.indicator = "z",
                  outcome = "y")
  
  Kernel_full <- kernel_calc_pred_(X_ = as.matrix(x), 
                                   X_test_ = as.matrix(x),
                                   z = as.integer(z), p = param$p, 
                                   theta_ = param$theta, 
                                   gamma_ = param$gamma, 
                                   sigma_2_ = param$sigma_2,
                                   calc_covariance = isTRUE(param$metric == "mahalanobis"), 
                                   kernel_ = as.character(param$kernel), 
                                   estimand = as.character("ATE"))
  testthat::expect_equal(tau, rout$tau)
  testthat::expect_equal(Kernel_full[[1]]$cov, rout$kernels$control$cov)
  testthat::expect_equal(Kernel_full[[2]]$cov, rout$kernels$treated$cov)
  testthat::expect_equal(Kernel_full[[1]]$cross, rout$kernels$control$cross)
  testthat::expect_equal(Kernel_full[[2]]$cross, rout$kernels$treated$cross)
})
 
testthat::test_that("gaussian process prediction works rbf, mahalanobis", {
  
  n <- 100
  d <- 10
  x <- matrix(rnorm(n*d),n,d)
  y <- rnorm(n)
  z <- rbinom(n,1,0.5)
  n1 <- sum(z)
  n0 <- n - n1
  
  param <- list(theta = c(1.1,1.22),
                gamma = c(0.5,0.423),
                p = 2,
                sigma_2 = c(1.01, 1.033),
                is.dose = FALSE,
                kernel = "RBF",
                metric = "mahalanobis",
                is.dose = FALSE,                
                is.standardized = FALSE

  )
  
  rout <- gauss_pred_r(x,y,z, param)
  
  df <- data.frame(x, y = y, z = z)
  tau <- gp_pred(data = df, weights = NULL,param = param, estimand = "ATE", 
                 balance.covariates = colnames(df)[1:d],
                 treatment.indicator = "z",
                 outcome = "y")
  
  Kernel_full <- kernel_calc_pred_(X_=as.matrix(x), 
                                   X_test_ = as.matrix(x),
                                   z = as.integer(z), p = param$p, 
                                   theta_ = param$theta, 
                                   gamma_ = param$gamma, 
                                   sigma_2_ = param$sigma_2,
                                   calc_covariance = isTRUE(param$metric == "mahalanobis"), 
                                   kernel = as.character(param$kernel), 
                                   estimand = as.character("ATE"))
  testthat::expect_equal(tau, rout$tau)
  testthat::expect_equal(Kernel_full[[1]]$cov, rout$kernels$control$cov)
  testthat::expect_equal(Kernel_full[[2]]$cov, rout$kernels$treated$cov)
  testthat::expect_equal(Kernel_full[[1]]$cross, rout$kernels$control$cross)
  testthat::expect_equal(Kernel_full[[2]]$cross, rout$kernels$treated$cross)
})

testthat::test_that("gaussian process prediction works polynomial, Lp", {
  
  n <- 100
  d <- 10
  x <- matrix(rnorm(n*d),n,d)
  y <- rnorm(n)
  z <- rbinom(n,1,0.5)
  n1 <- sum(z)
  n0 <- n - n1
  
  param <- list(theta = c(1.1,1.22),
                gamma = c(0.5,0.423),
                p = 2,
                sigma_2 = c(1.01, 1.033),
                is.dose = FALSE,
                kernel = "polynomial",
                metric = "Lp",
                is.dose = FALSE,
                is.standardized = FALSE
                
  )
  
  rout <- gauss_pred_r(x,y,z, param)
  
  df <- data.frame(x, y = y, z = z)
  tau <- gp_pred(data = df, weights = NULL,param = param, estimand = "ATE", 
                 balance.covariates = colnames(df)[1:d],
                 treatment.indicator = "z",
                 outcome = "y")
  
  Kernel_full <- kernel_calc_pred_(X_=as.matrix(x), 
                                   X_test_ = as.matrix(x),
                                   z = as.integer(z), p = param$p, 
                                   theta_ = param$theta, 
                                   gamma_ = param$gamma, 
                                   sigma_2_ = param$sigma_2,
                                   calc_covariance = isTRUE(param$metric == "mahalanobis"), 
                                   kernel = as.character(param$kernel), 
                                   estimand = as.character("ATE"))
  testthat::expect_equal(tau, rout$tau)
  testthat::expect_equal(Kernel_full[[1]]$cov, rout$kernels$control$cov)
  testthat::expect_equal(Kernel_full[[2]]$cov, rout$kernels$treated$cov)
  testthat::expect_equal(Kernel_full[[1]]$cross, rout$kernels$control$cross)
  testthat::expect_equal(Kernel_full[[2]]$cross, rout$kernels$treated$cross)
})

testthat::test_that("gaussian process prediction works rbf, Lp", {
  
  n <- 100
  d <- 10
  x <- matrix(rnorm(n*d),n,d)
  y <- rnorm(n)
  z <- rbinom(n,1,0.5)
  n1 <- sum(z)
  n0 <- n - n1
  
  param <- list(theta = c(1.1,1.22),
                gamma = c(0.5,0.423),
                p = 2,
                sigma_2 = c(1.01, 1.033),
                is.dose = FALSE,
                kernel = "RBF",
                metric = "Lp",
                is.dose = FALSE,
                is.standardized = FALSE
  )
  
  rout <- gauss_pred_r(x,y,z, param)
  
  df <- data.frame(x, y = y, z = z)
  tau <- gp_pred(data = df, weights = NULL,param = param, estimand = "ATE", 
                 balance.covariates = colnames(df)[1:d],
                 treatment.indicator = "z",
                 outcome = "y")
  
  Kernel_full <- kernel_calc_pred_(X_=as.matrix(x), 
                                   X_test_ = as.matrix(x),
                                   z = as.integer(z), p = param$p, 
                                   theta_ = param$theta, 
                                   gamma_ = param$gamma, 
                                   sigma_2_ = param$sigma_2,
                                   calc_covariance = isTRUE(param$metric == "mahalanobis"), 
                                   kernel = as.character(param$kernel), 
                                   estimand = as.character("ATE"))
  testthat::expect_equal(tau, rout$tau)
  testthat::expect_equal(Kernel_full[[1]]$cov, rout$kernels$control$cov)
  testthat::expect_equal(Kernel_full[[2]]$cov, rout$kernels$treated$cov)
  testthat::expect_equal(Kernel_full[[1]]$cross, rout$kernels$control$cross)
  testthat::expect_equal(Kernel_full[[2]]$cross, rout$kernels$treated$cross)
})

testthat::test_that("gaussian process corrects non-positive def matrix rbf, Lp", {
  
  n <- 100
  d <- 10
  x <- matrix(rnorm(n*d),n,d)
  y <- rnorm(n)
  z <- rbinom(n,1,0.5)
  n1 <- sum(z)
  n0 <- n - n1
  
  param <- list(theta = c(1.1,1.22),
                gamma = c(0.5,0.423),
                p = 2,
                sigma_2 = -10*c(1.01, 1.033),
                is.dose = FALSE,
                kernel = "RBF",
                metric = "Lp",
                is.dose = FALSE,
                is.standardized = TRUE
  )
  
  df <- data.frame(x, y = y, z = z)
  testthat::expect_silent(tau <- gp_pred(data = df, weights = NULL,param = param, estimand = "ATE", 
                 balance.covariates = colnames(df)[1:d],
                 treatment.indicator = "z",
                 outcome = "y"))
})

# testthat::test_that("timing is better for C code, rbf lp", {
#   testthat::skip("Interactive only")
#   n <- 2^9
#   d <- 10
#   x <- matrix(rnorm(n*d),n,d)
#   y <- rnorm(n)
#   z <- rbinom(n,1,0.5)
#   n1 <- sum(z)
#   n0 <- n - n1
#   
#   param <- list(theta = c(1.1,1.22),
#                 gamma = c(0.5,0.423),
#                 p = 2,
#                 sigma_2 = c(1.01, 1.033),
#                 is.dose = FALSE,
#                 kernel = "RBF",
#                 metric = "Lp",
#                 is.dose = FALSE
#   )
#   df <- data.frame(x, y = y, z = z)
#   mb <- microbenchmark::microbenchmark(rcode = {
#     prep.data <-prep_data(df, balance.covariates = colnames(df)[1:d],
#               treatment.indicator = "z",
#               outcome = "y")
#     xx <- prep.data$df[,-c(which(colnames(prep.data$df)=="y"))]
#     yy <- prep.data$df$y
#     zz <- prep.data$z
#     gauss_pred_r(xx,yy,zz, param)
#     },
#                                  ccode = gp_pred(data = df, param = param, estimand = "ATE", 
#                                                  balance.covariates = colnames(df)[1:d],
#                                                  treatment.indicator = "z",
#                                                  outcome = "y"),
#                                  times=100)
#   print(mb)
#   microbenchmark:::autoplot.microbenchmark(mb)
#   testthat::expect_lt(median(mb$time[mb$expr=="ccode"]),
#                       median(mb$time[mb$expr=="rcode"])
#                       )
#   
# })
# testthat::test_that("timing is better for C code, rbf mahalanobis", {
#   testthat::skip("Interactive only")
#   n <- 2^9
#   d <- 10
#   x <- matrix(rnorm(n*d),n,d)
#   y <- rnorm(n)
#   z <- rbinom(n,1,0.5)
#   n1 <- sum(z)
#   n0 <- n - n1
#   
#   param <- list(theta = c(1.1,1.22),
#                 gamma = c(0.5,0.423),
#                 p = 2,
#                 sigma_2 = c(1.01, 1.033),
#                 is.dose = FALSE,
#                 kernel = "RBF",
#                 metric = "mahalanobis",
#                 is.dose = FALSE
#   )
#   df <- data.frame(x, y = y, z = z)
#   mb <- microbenchmark::microbenchmark(rcode = {
#     prep.data <-prep_data(df, balance.covariates = colnames(df)[1:d],
#                           treatment.indicator = "z",
#                           outcome = "y")
#     xx <- prep.data$df[,-c(which(colnames(prep.data$df)=="y"))]
#     yy <- prep.data$df$y
#     zz <- prep.data$z
#     gauss_pred_r(xx,yy,zz, param)
#   },
#   ccode = gp_pred(data = df, param = param, estimand = "ATE", 
#                   balance.covariates = colnames(df)[1:d],
#                   treatment.indicator = "z",
#                   outcome = "y"),
#   times=100)
#   print(mb)
#   microbenchmark:::autoplot.microbenchmark(mb)
#   testthat::expect_lt(median(mb$time[mb$expr=="ccode"]),
#                       median(mb$time[mb$expr=="rcode"])
#   )
#   
# })
# testthat::test_that("timing is better for C code, polynomial lp", {
#   testthat::skip("Interactive only")
#   n <- 2^9
#   d <- 10
#   x <- matrix(rnorm(n*d),n,d)
#   y <- rnorm(n)
#   z <- rbinom(n,1,0.5)
#   n1 <- sum(z)
#   n0 <- n - n1
#   
#   param <- list(theta = c(1.1,1.22),
#                 gamma = c(0.5,0.423),
#                 p = 2,
#                 sigma_2 = c(1.01, 1.033),
#                 is.dose = FALSE,
#                 kernel = "polynomial",
#                 metric = "Lp",
#                 is.dose = FALSE
#   )
#   df <- data.frame(x, y = y, z = z)
#   mb <- microbenchmark::microbenchmark(rcode = {
#     prep.data <-prep_data(df, balance.covariates = colnames(df)[1:d],
#                           treatment.indicator = "z",
#                           outcome = "y")
#     xx <- prep.data$df[,-c(which(colnames(prep.data$df)=="y"))]
#     yy <- prep.data$df$y
#     zz <- prep.data$z
#     gauss_pred_r(as.matrix(xx),yy,zz, param)
#   },
#   ccode = gp_pred(data = df, param = param, estimand = "ATE", 
#                   balance.covariates = colnames(df)[1:d],
#                   treatment.indicator = "z",
#                   outcome = "y"),
#   times=100)
#   print(mb)
#   microbenchmark:::autoplot.microbenchmark(mb)
#   testthat::expect_lt(median(mb$time[mb$expr=="ccode"]),
#                       median(mb$time[mb$expr=="rcode"])
#   )
#   
# })
# testthat::test_that("timing is better for C code, polynomial mahalanobis", {
#   testthat::skip("Interactive only")
#   n <- 2^9
#   d <- 10
#   x <- matrix(rnorm(n*d),n,d)
#   y <- rnorm(n)
#   z <- rbinom(n,1,0.5)
#   n1 <- sum(z)
#   n0 <- n - n1
#   
#   param <- list(theta = c(1.1,1.22),
#                 gamma = c(0.5,0.423),
#                 p = 2,
#                 sigma_2 = c(1.01, 1.033),
#                 is.dose = FALSE,
#                 kernel = "polynomial",
#                 metric = "mahalanobis",
#                 is.dose = FALSE
#   )
#   df <- data.frame(x, y = y, z = z)
#   mb <- microbenchmark::microbenchmark(rcode = {
#     prep.data <-prep_data(df, balance.covariates = colnames(df)[1:d],
#                           treatment.indicator = "z",
#                           outcome = "y")
#     xx <- prep.data$df[,-c(which(colnames(prep.data$df)=="y"))]
#     yy <- prep.data$df$y
#     zz <- prep.data$z
#     gauss_pred_r(xx,yy,zz, param)
#   },
#   ccode = gp_pred(data = df, param = param, estimand = "ATE", 
#                   balance.covariates = colnames(df)[1:d],
#                   treatment.indicator = "z",
#                   outcome = "y"),
#   times=100)
#   print(mb)
#   microbenchmark:::autoplot.microbenchmark(mb)
#   testthat::expect_lt(median(mb$time[mb$expr=="ccode"]),
#                       median(mb$time[mb$expr=="rcode"])
#   )
#   
# })

testthat::test_that("gaussian process prediction works rbf, lp", {
  testthat::skip("Interactive only")
  set.seed(3208)
  n <- 512
  d <- 10
  x <- matrix(rnorm(n*d),n,d)
  y <- 5*cos(x[,1]* 10) - 2*x[,1]^2+ 2 + rnorm(n, sd = 1)
  z <- rbinom(n,1,0.5)
  n1 <- sum(z)
  n0 <- n - n1
  
  
  m_y0   <- mean(y[z==0])
  m_y1   <- mean(y[z==1])
  sd_y0   <- sd(y[z==0])
  sd_y1   <- sd(y[z==1])
  y0 <- c(scale(y[z==0]))
  y1 <- c(scale(y[z==1]))
  
  plot(x[z==0,1], y[z==0], col = "black")
  points(x[z==1,1], y[z==1], col = "gray")
  # debugonce(RKHS_param_opt)
  # debugonce(rkhs_stan_binary_helper)
  kernel <- "RBF"
  param <- RKHS_param_opt(x[,1,drop = FALSE], y, z, power = 2:3, metric = "mahalanobis",
                          is.dose = FALSE, opt.method = "stan", kernel = kernel,
                          estimand = "ATE", verbose = TRUE, algorithm = "LBFGS")
  # cov <- cov(x)
  # A <- scale(x, center = TRUE, scale = FALSE) %*% solve(chol(cov))
  A <- scale(x[,1,drop=FALSE])
  # A <- x[,1,drop=FALSE] #scale(x[,1,drop=FALSE])
  A0 <- A[z==0,,drop = FALSE]
  A1 <- A[z==1,,drop = FALSE]

  if(kernel == "RBF") {
    kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + 
      param$gamma[1] * exp(-0.5 * param$theta[1] * 
                             cost_calc_lp(A0,A0,ground_p = 2, direction = "rowwise" )^2)
    # if(is.infinite(param$theta[2])) kernel_cov1 <- diag(param$sigma_2[2],n1,n1) * param$gamma[2]
    kernel_cross0 <- param$gamma[1] * exp(-0.5 *  param$theta[1] * 
                                            cost_calc_lp(A0,A,ground_p = 2, direction = "rowwise" )^2)
    
    kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + 
      param$gamma[2] * exp(-0.5 * param$theta[2] * 
                             cost_calc_lp(A1,A1,ground_p = 2, direction = "rowwise" )^2)
    if(is.infinite(param$theta[2])) kernel_cov1 <- diag(param$sigma_2[2],n1,n1) * param$gamma[2]
    kernel_cross1 <- param$gamma[2] * exp(-0.5 *  param$theta[2] * 
                                            cost_calc_lp(A1,A,ground_p = 2, direction = "rowwise" )^2)
    
  } else if (kernel == "polynomial") {
    kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + 
      param$gamma[1] * (1+ param$theta[1] * tcrossprod(A0))^param$p
    kernel_cross0 <- param$gamma[1] * (1+ param$theta[1] * tcrossprod(A0,A))^param$p
    
    kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + 
      param$gamma[2] * (1+ param$theta[2] * tcrossprod(A1))^param$p
    kernel_cross1 <- param$gamma[2] * (1+ param$theta[2] * tcrossprod(A1,A))^param$p
  } else if (kernel == "linear") {
      kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + tcrossprod(A0)
      kernel_cross0 <- tcrossprod(A0,A)
      
      kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + tcrossprod(A1)
      kernel_cross1 <- tcrossprod(A1,A)
  }
  
  pred0 <- crossprod(kernel_cross0, solve(kernel_cov0, y0)) * sd_y0 + m_y0
  pred1 <- crossprod(kernel_cross1, solve(kernel_cov1, y1)) * sd_y1 + m_y1
  points(x[,1], pred0, col = "blue")
  points(x[,1], pred1, col = "red")
  
  # cbind(pred0, pred1)
  # rout <- gauss_pred_r(x,y,z, param)
  
  df <- data.frame(x[,1,drop=FALSE], y = y, z = z)
  tau <- gp_pred(data = df, weights = NULL,param = param, estimand = "ATE", 
                 balance.covariates = colnames(df)[1],
                 treatment.indicator = "z",
                 outcome = "y")
  print(tau)
  w <- calc_weight(df, estimand = "ATE", method = "RKHS",
                   kernel = kernel,
                   balance.covariates = colnames(df)[1],
                   treatment.indicator = "z",
                   outcome = "y",
                   opt.hyperparam = FALSE,
                   theta = param$theta,
                   gamma = param$gamma,
                   sigma_2 = param$sigma_2,
                   p = param$p,
                   metric = "mahalanobis",
                   is.standardized = TRUE)
  
  estimate_effect(df, weights = w, target = c("ATE"),
                  balance.covariates = colnames(df)[1],
                  treatment.indicator = "z",
                  outcome = "y")
  
})
