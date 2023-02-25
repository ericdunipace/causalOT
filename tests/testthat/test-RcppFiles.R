
testthat::test_that("LogSumExp works", {
  # log sum exp function
  log_sum_exp <- function(x) {
    # if(is.vector(x)) {
    if(all(is.infinite(x))) return(x[1])
    mx <- max(x)
    x_temp <- x - mx
    return(log(sum(exp(x_temp)))+ mx)
    # } else if (is.matrix(x)) {
    #   mx <- apply(x, 1, max)
    #   x_temp <- x - mx
    #   return(log(rowSums(exp(x_temp)))+ mx)
    # }
  }
  x <- stats::rnorm(100)
  cpp <- causalOT:::logSumExp(x)
  R <- log_sum_exp(x)
  
  testthat::expect_equal(cpp, R)
})

testthat::test_that("rowLogSumExp works", {
  # log sum exp function by row
  row_log_sum_exp <- function(x) {
    # if(is.vector(x)) {
    if(all(is.infinite(x))) return(x[1])
    mx <- apply(x,1,max)
    mx_mat <- matrix(mx,nrow(x),ncol(x))
    x_temp <- x - mx_mat
    return(log(rowSums(exp(x_temp))) + mx)
    # } else if (is.matrix(x)) {
    #   mx <- apply(x, 1, max)
    #   x_temp <- x - mx
    #   return(log(rowSums(exp(x_temp)))+ mx)
    # }
  }
  x <- matrix(stats::rnorm(100*1001), 100, 1001)
  cpp <- causalOT:::rowLogSumExp(x)
  R <- row_log_sum_exp(x)
  
  testthat::expect_equal(cpp, R)
})

testthat::test_that("colLogSumExp works", {
  # log sum exp function by row
  col_log_sum_exp <- function(x) {
    # if(is.vector(x)) {
    if(all(is.infinite(x))) return(x[1])
    mx <- apply(x,2,max)
    mx_mat <- matrix(mx,nrow(x),ncol(x), byrow = TRUE)
    x_temp <- x - mx_mat
    return(log(colSums(exp(x_temp))) + c(mx) )
    # } else if (is.matrix(x)) {
    #   mx <- apply(x, 1, max)
    #   x_temp <- x - mx
    #   return(log(rowSums(exp(x_temp)))+ mx)
    # }
  }
  x <- matrix(stats::rnorm(100*1001), 100, 1001)
  cpp <- causalOT:::colLogSumExp(x)
  R <- col_log_sum_exp(x)
  
  testthat::expect_equal(cpp, R)
})

testthat::test_that("entropy bal weights cpp works", {
  ebal_obj <- function(par, A, delta) {
    k <- ncol(A)
    par_pos <- par[1:k]
    par_neg <- par[-c(1:k)]
    
    beta <- par_pos - par_neg
    
    lp <- A %*% beta
    obj <- sum(abs(beta)) * delta + causalOT:::logSumExp(lp)
    return(obj)
  } 
  ebal_grad <- function(par, A, delta) {
    k <- ncol(A)
    par_pos <- par[1:k]
    par_neg <- par[-c(1:k)]
    
    beta <- par_pos - par_neg
    
    lp <- A %*% beta
    s <- causalOT:::logSumExp(lp)
    cross_term <- crossprod(A, exp(lp - s))
    grad <- sign(par) * delta + c(cross_term, -cross_term)
    return(grad)
  }
  n <- 100
  d <- 5
  A <- matrix(stats::rnorm(n * d), n , d)
  delta <- 0.1
  var <- rep(0, 2 * d)
  
  cpp_obj <- causalOT:::entBW_obj_(var, A, delta)
  R_obj <- ebal_obj(var, A, delta)
  
  cpp_grad <- causalOT:::entBW_grad_(var, A, delta)
  R_grad <- ebal_grad(var, A, delta)
  
  testthat::expect_equal(cpp_obj, R_obj)
  testthat::expect_equal(cpp_grad, R_grad)
})

# testthat::test_that("cot biased entropy works", {
#   cot_obj <- function(par, source, target, cost, b, delta, lambda) {
#     n <- nrow(cost)
#     m <- ncol(cost)
#     k <- ncol(source)
#     g <- par[1:m]
#     par_A <- par[-c(1:m)]
#     par_pos <- par_A[1:k]
#     par_neg <- par_A[-c(1:k)]
#     
#     beta <- par_pos - par_neg
#     obj <- 0
#     
#     A_lp <- c(source %*% beta)/lambda
#     lp <- (matrix(g , n, m, byrow = TRUE) - cost)/lambda 
#     if ( length(A_lp) >0 ) {
#       lp <- lp - A_lp
#       obj <- -sum(par_A) * delta - sum(target * beta)
#     }
#     obj <- obj -
#       lambda * causalOT:::logSumExp(lp) +
#       sum(b * g)
#     # print(sum(b * g))
#     # print(-sum(par_A) * delta )
#     # print(lambda * causalOT:::logSumExp(lp))
#     # print(obj)
#     return(-obj)
#   } 
#   cot_grad <- function(par, source, target, cost, b, delta, lambda) {
#     n <- nrow(cost)
#     m <- ncol(cost)
#     k <- ncol(source)
#     g <- par[1:m]
#     par_A <- par[-c(1:m)]
#     par_pos <- par_A[1:k]
#     par_neg <- par_A[-c(1:k)]
#     
#     beta <- par_pos - par_neg
#     
#     A_lp <- c(source %*% beta)/lambda
#     lp <- (matrix(g , n, m, byrow = TRUE) - cost)/lambda 
#     if ( length(A_lp) >0 ) lp <- lp - A_lp
#     s <- causalOT:::logSumExp(lp)
#     lp <- lp - s
#     g_grad <- b - c(colSums(exp(lp )))
#     if(length(A_lp) > 0) {
#       cross <- crossprod(source, c(exp(causalOT:::rowLogSumExp(lp ))))
#       grad_A <- -sign(par_A) * delta + c(-target,  target) + c(cross, -cross)
#       # print(-sign(par_A) * delta)
#       # print(cross)
#       grad <- c(g_grad, 
#               grad_A)
#     } else {
#       grad <- g_grad
#     }
#     return(-grad)
#   }
#   set.seed(23048)
#   n <- 100
#   m <- 111
#   d <- 5
#   b <- rep(1/m,m)
#   A <- matrix(stats::rnorm(n * d), n , d)
#   cost <- matrix(stats::rexp(n *m), n, m)
#   delta <- 0.1
#   lambda <- 10
#   var <- as.double(1:(m + 2 * d))
#   
#   cpp_obj <- causalOT:::cotEntropy_obj_(vars_ =var, source_ = A, 
#                                         target_ = colMeans(A),
#                                         cost_ = cost, 
#                                         b_ = b,
#                                         delta = delta, lambda = lambda)
#   R_obj <- cot_obj(var, A, colMeans(A), cost, b,
#                    delta, lambda)
#   
#   cpp_grad <- causalOT:::cotEntropy_grad_(vars_ =var, source_ = A, 
#                                           target_ = colMeans(A),
#                                           cost_ = cost, 
#                                           b_ = b,
#                                           delta = delta, lambda = lambda)
#   R_grad <- cot_grad(var, A, colMeans(A), cost, b,
#                      delta, lambda)
#   
#   testthat::expect_equal(cpp_obj, R_obj)
#   testthat::expect_equal(cpp_grad, R_grad)
#   
#   # crossprod_A <- Matrix::crossprod(causalOT:::vec_to_row_constraints(n,m),A)
#   # 
#   # QQ <- cbind(Matrix::t(causalOT:::vec_to_col_constraints(n,m)),
#   #             -crossprod_A, crossprod_A
#   #             )
#   # pf <- c(b, rep(-delta, 2 * d))
#   # 
#   # testthat::expect_equal(causalOT:::cotDual_obj_(var, QQ, c(cost), pf, lambda, "entropy"), cpp_obj)
#   
#   # without bf
#   set.seed(23048)
#   n <- 100
#   m <- 111
#   d <- 5
#   b <- rep(1/m,m)
#   A <- matrix(0,0,0)
#   cost <- matrix(stats::rexp(n *m), n, m)
#   delta <- 0.1
#   lambda <- 10
#   var <- rep(0, m )
#   
#   cpp_obj <- causalOT:::cotEntropy_obj_(vars_ =var, source_ = A, 
#                                         target_ = colMeans(A),
#                                         cost_ = cost, 
#                                         b_ = b,
#                                         delta = delta, lambda = lambda)
#   R_obj <- cot_obj(var, A, colMeans(A), cost, b,
#                     delta, lambda)
#   
#   cpp_grad <- causalOT:::cotEntropy_grad_(vars_ =var, source_ = A, 
#                                           target_ = colMeans(A),
#                                           cost_ = cost, 
#                                           b_ = b,
#                                           delta = delta, lambda = lambda)
#   R_grad <- cot_grad(var, A, colMeans(A), cost, b,
#                      delta, lambda)
#   
#   testthat::expect_equal(cpp_obj, R_obj)
#   testthat::expect_equal(cpp_grad, R_grad)
#   
#   # QQ <- Matrix::t(causalOT:::vec_to_col_constraints(n,m))
#   # pf <- b
#   # 
#   # testthat::expect_equal(causalOT:::cotDual_obj_(var, QQ, c(cost), pf, lambda, "entropy"), cpp_obj)
#   
# })
