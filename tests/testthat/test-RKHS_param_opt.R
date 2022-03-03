gp_pred_fun <- function(data, param) {
  test_pos_def_inv <- function(x,y) {
    e <- eigen(x)
    if(any(e$values <=0)) {
      min.e <- min(e$values)
      e$values <- e$values - min.e
    }
    return(e$vectors %*% diag(1/e$values) %*% crossprod(e$vectors, y))
  }
  # w0 <- weights$w0
  # w1 <- weights$w1
  prep.data <- prep_data(data)
  
  z <- prep.data$z
  y <- prep.data$df$y
  x <- as.matrix(prep.data$df[,-which(colnames(prep.data$df)=="y")])
  
  n <- length(z)
  n1 <- sum(z)
  n0 <- n - n1
  
  
  if(param$is.standardized) {
    m_y0   <- mean(y[z==0])
    m_y1   <- mean(y[z==1])
    sd_y0   <- sd(y[z==0])
    sd_y1   <- sd(y[z==1])
    y[z==0] <- c(scale(y[z==0]))
    y[z==1] <- c(scale(y[z==1]))
  } else {
    m_y0   <- 0
    m_y1   <- 0
    sd_y0   <- 1
    sd_y1   <- 1
  }
  
  
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
  } else if (param$kernel == "linear") {
    kernel_cov0 <- diag(param$sigma_2[1],n0,n0) + tcrossprod(A0)
    kernel_cross0 <- tcrossprod(A0,A)
    
    kernel_cov1 <- diag(param$sigma_2[2],n1,n1) + tcrossprod(A1)
    kernel_cross1 <- tcrossprod(A1,A)
  }
  
  
    pred0 <- crossprod(kernel_cross0, test_pos_def_inv(kernel_cov0, y[z==0])) * sd_y0 + m_y0
    pred1 <- crossprod(kernel_cross1, test_pos_def_inv(kernel_cov1, y[z==1])) * sd_y1 + m_y1
  
  
  return(list(pred0, pred1))
}

testthat::test_that("parameter optimization for RKHS", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  set.seed(290384)
  library(causalOT)
  
  n <- 2^9
  p <- 6
  overlap <- "low"
  design <- "A"
  distance <- c("mahalanobis")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  x <- data$get_x()
  y <- data$get_y()
  z <- data$get_z()
  
  # debugonce("RKHS_param_opt")
  
  # opt1 <- RKHS_param_opt(x, y, z, power = 2:3,
  #                metric = c("mahalanobis", "Lp"), is.dose = FALSE, 
  #                        opt.method = c("bayesian.optimization"),
  #                bounds = list(theta_0 = c(.Machine$double.xmin,100),
  #                     theta_1 = c(.Machine$double.xmin,100),
  #                     gamma_0 = c(.Machine$double.xmin,1000),
  #                     gamma_1 = c(.Machine$double.xmin,1000),
  #                     sigma2 = c(.Machine$double.xmin,1)),
  #                kernel = "polynomial",
  #                initPoints = 6,
  #                iters.k = 1,
  #                iters.n = 1
  #                ) 
  
  # debugonce("RKHS_param_opt")
  
  testthat::expect_warning(opt2 <- RKHS_param_opt(x, y, z, power = 2:3,
                        metric = c("mahalanobis"), is.dose = FALSE, 
                        kernel = "polynomial",
                        opt.method = c("optim"), control = list(maxit = 10)))
  testthat::expect_warning(opt2 <- RKHS_param_opt(x, y, z, power = 2:3,
                                                  metric = c("mahalanobis"), is.dose = FALSE, 
                                                  kernel = "RBF",
                                                  opt.method = c("optim"), control = list(maxit = 10)))
  
  # debugonce("RKHS_param_opt")
  testthat::expect_silent(opt3 <- causalOT:::RKHS_param_opt(x, y, z, p = 2:3,
                        metric = c("mahalanobis"), is.dose = FALSE, 
                        kernel = "polynomial",
                        opt.method = c("stan"), iter = 10))
  
  testthat::expect_message(opt3 <- causalOT:::RKHS_param_opt(x, y, z, p = 2:3,
                                                 metric = c("mahalanobis"), is.dose = FALSE, 
                                                 kernel = "RBF",
                                                 opt.method = c("stan"), iter = 10))
  
})

testthat::test_that("parameter optimization for RKHS", {
  testthat::skip_on_cran()
  set.seed(290384)
  library(causalOT)
  
  n <- 2^9
  p <- 6
  overlap <- "low"
  design <- "A"
  distance <- c("mahalanobis")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  x <- data$get_x()
  y <- data$get_y()
  z <- data$get_z()
  
  # debugonce("RKHS_param_opt")
  
  # opt1 <- RKHS_param_opt(x, y, z, power = 2:3,
  #                metric = c("mahalanobis", "Lp"), is.dose = FALSE, 
  #                        opt.method = c("bayesian.optimization"),
  #                bounds = list(theta_0 = c(.Machine$double.xmin,100),
  #                     theta_1 = c(.Machine$double.xmin,100),
  #                     gamma_0 = c(.Machine$double.xmin,1000),
  #                     gamma_1 = c(.Machine$double.xmin,1000),
  #                     sigma2 = c(.Machine$double.xmin,1)),
  #                kernel = "polynomial",
  #                initPoints = 6,
  #                iters.k = 1,
  #                iters.n = 1
  #                ) 
  
  # debugonce("RKHS_param_opt")
  
  testthat::expect_warning(opt2 <- RKHS_param_opt(x, y, z, power = 2:3,
                                                  metric = c("mahalanobis"), is.dose = FALSE, 
                                                  kernel = "polynomial",
                                                  opt.method = c("optim"), control = list(maxit = 10)))
  testthat::expect_warning(opt2 <- RKHS_param_opt(x, y, z, power = 2:3,
                                                  metric = c("mahalanobis"), is.dose = FALSE, 
                                                  kernel = "RBF",
                                                  opt.method = c("optim"), control = list(maxit = 10)))
  testthat::expect_warning(opt2 <- RKHS_param_opt(x, y, z, power = 2:3,
                                                  metric = c("mahalanobis"), is.dose = FALSE, 
                                                  kernel = "linear",
                                                  opt.method = c("optim"), control = list(maxit = 10)))
  
  # debugonce("RKHS_param_opt")
  testthat::expect_silent(opt3 <- RKHS_param_opt(x, y, z, p = 2:3,
                                                 metric = c("mahalanobis"), is.dose = FALSE, 
                                                 kernel = "polynomial",
                                                 opt.method = c("stan"), iter = 10))
  
  testthat::expect_silent(opt3 <- RKHS_param_opt(x, y, z, p = 2:3,
                                                 metric = c("mahalanobis"), is.dose = FALSE, 
                                                 kernel = "RBF",
                                                 opt.method = c("stan"), iter = 10))
  
  testthat::expect_silent(opt3 <- RKHS_param_opt(x, y, z, p = 2:3,
                                                 metric = c("mahalanobis"), is.dose = FALSE, 
                                                 kernel = "linear",
                                                 opt.method = c("stan"), iter = 10))
  
  
  testthat::expect_silent(opt3 <- RKHS_param_opt(x, y, z, p = 2:3,
                                                 metric = c("Lp"), is.dose = FALSE, 
                                                 kernel = "polynomial",
                                                 opt.method = c("stan"), iter = 10))
  
  testthat::expect_silent(opt3 <- RKHS_param_opt(x, y, z, p = 2:3,
                                                 metric = c("Lp"), is.dose = FALSE, 
                                                 kernel = "RBF",
                                                 opt.method = c("stan"), iter = 10))
  
  testthat::expect_silent(opt3 <- RKHS_param_opt(x, y, z, p = 2:3,
                                                 metric = c("Lp"), is.dose = FALSE, 
                                                 kernel = "linear",
                                                 opt.method = c("stan"), iter = 10))
  
})

testthat::test_that("parameter optimization for RKHS poly correct", {
  testthat::skip("Interactive only")
  set.seed(290384)
  library(causalOT)
  
  n <- 2^9
  p <- 6
  overlap <- "high"
  design <- "B"
  metric <- c("mahalanobis")
  kernel = "polynomial"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  x <- data$get_x()
  y <- data$get_y()
  z <- data$get_z()
  
  # debugonce("RKHS_param_opt")
  
  # opt1 <- RKHS_param_opt(x, y, z, power = 2:3,
  #                metric = c("mahalanobis", "Lp"), is.dose = FALSE, 
  #                        opt.method = c("bayesian.optimization"),
  #                bounds = list(theta_0 = c(.Machine$double.xmin,100),
  #                     theta_1 = c(.Machine$double.xmin,100),
  #                     gamma_0 = c(.Machine$double.xmin,1000),
  #                     gamma_1 = c(.Machine$double.xmin,1000),
  #                     sigma2 = c(.Machine$double.xmin,1)),
  #                kernel = "polynomial",
  #                initPoints = 6,
  #                iters.k = 1,
  #                iters.n = 1
  #                ) 
  
  # debugonce("RKHS_param_opt")
  
  #verify similarity mats
  zmat <- scale(x, scale = FALSE, center = TRUE) %*% solve(sqrt_mat(cov(x)))
  testthat::expect_equal(calc_similarity(x, z, metric = metric, kernel = kernel, is.dose = FALSE, estimand),
                         tcrossprod(zmat))
  
  # debugonce("RKHS_param_opt")
  polres <- RKHS_param_opt(x, y, z, power = 2:3,
                                                  metric = c("mahalanobis"), is.dose = FALSE, 
                                                  kernel = "polynomial",
                                                  opt.method = c("stan"))
  
  rbfres <- RKHS_param_opt(x, y, z, power = 2:3,
                           metric = c("mahalanobis"), is.dose = FALSE, 
                           kernel = "RBF",
                           opt.method = c("stan"))
  
  
  preds <- gp_pred_fun(data = data, param = polres)
  mean((preds[[1]] - y)^2)
  mean((preds[[2]] - y)^2)
  mean((preds[[1]][z==0] - y[z==0])^2)
  mean((preds[[2]][z==1] - y[z==1])^2)
  
  preds <- gp_pred_fun(data = data, param = rbfres)
  mean((preds[[1]] - y)^2)
  mean((preds[[2]] - y)^2)
  mean((preds[[1]][z==0] - y[z==0])^2)
  mean((preds[[2]][z==1] - y[z==1])^2)
  
  plot(y, preds[[1]], col = "red")
  points(y, preds[[2]], col = "blue")
  abline(0,1)
})

testthat::test_that("gp pred are ok with rbf", {
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
  
})

