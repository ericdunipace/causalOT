testthat::test_that("costParent class forms", {
  testthat::expect_silent ( cost <- causalOT:::costParent$new())
})

testthat::test_that("costTensor class forms", {
  causalOT:::torch_check()
  set.seed(124123)
  n <- 10
  m <- 11
  d <- 4
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  
  # given just p
  testthat::expect_silent(cost <- causalOT:::costTensor$new(x = x, y = y, p = 2L))
  testthat::expect_equal(cost$data$dim(), 2)
  testthat::expect_equal(cost$data$shape, c(10L,11L))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costTensor"))
  testthat::expect_true(inherits(cost, "cost"))
  testthat::expect_equivalent(object = as.matrix(cost$data), 
                         expected = (as.matrix(stats::dist(rbind(x,y)))^2)[1:n,(n+1):(m+n)]/2,
                         ignore_attr = TRUE,
                         tolerance = 1e-5)
  
  
  # given cost_function
  cost_fun <- function(x,y,p) {
    n <- nrow(x)
    m <- nrow(y)
    as.matrix(stats::dist(rbind(x,y), method = "manhattan"))[1:n,(n+1):(m+n)]
  }
  testthat::expect_silent(cost <- causalOT:::costTensor$new(x = x, y = y, cost_function = cost_fun))
  testthat::expect_equal(cost$data$dim(), 2)
  testthat::expect_equal(cost$data$shape, c(10L,11L))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costTensor"))
  testthat::expect_true(inherits(cost, "cost"))
  testthat::expect_equivalent(object = as.matrix(cost$data), 
                              expected = as.matrix(stats::dist(rbind(x,y), method = "manhattan"))[1:n,(n+1):(m+n)],
                              ignore_attr = TRUE,
                              tolerance = 1e-5)
  testthat::expect_equal(object = cost$data, 
                              expected = causalOT:::costTensor$new(x = x, y = y, p = 1L)$data)
  
  
})

testthat::test_that("costOnline class forms", {
  testthat::skip_on_cran()
  causalOT:::rkeops_check()
  testthat::skip_on_ci()
  set.seed(124123)
  n <- 10
  m <- 11
  d <- 4
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  
  # given just p
  testthat::expect_silent(cost <- causalOT:::costOnline$new(x = x, y = y, p = 2L))
  testthat::expect_equal(dim(cost$data$x), c(n,d))
  testthat::expect_equal(dim(cost$data$y), c(m,d))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costOnline"))
  testthat::expect_true(inherits(cost, "cost"))
  causalOT:::rkeops_check()
  keops_sum <- rkeops::keops_kernel(
    formula = paste0("Sum_Reduction(", cost$fun, ", 0)"),
    args = c(
      paste0("X = Vi(",d,")"),
      paste0("Y = Vj(",d,")"))
  )
  testthat::expect_equivalent(object = sum(keops_sum(list(X = cost$data$x,
                                                 Y = cost$data$y) )), 
                              expected = sum((as.matrix(stats::dist(rbind(x,y)))[1:n,(n+1):(m+n)]^2) / 2),
                              ignore_attr = TRUE,
                              tolerance = 1e-5)
  
  
  # given cost_function
  cost_fun <- "Abs(X - Y)"
  testthat::expect_silent(cost <- causalOT:::costOnline$new(x = x, y = y, cost_function = cost_fun))
  testthat::expect_equal(dim(cost$data$x), c(n,d))
  testthat::expect_equal(dim(cost$data$y), c(m,d))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costOnline"))
  testthat::expect_true(inherits(cost, "cost"))
  causalOT:::rkeops_check()
  keops_sum2 <- rkeops::keops_kernel(
    formula = paste0("Sum_Reduction(", cost$fun, ", 0)"),
    args = c(
      paste0("X = Vi(",d,")"),
      paste0("Y = Vj(",d,")"))
  )
  testthat::expect_equivalent(object = sum(keops_sum2(list(X = cost$data$x,
                                                          Y = cost$data$y) )), 
                              expected = sum(as.matrix(stats::dist(rbind(x,y), method = "manhattan"))[1:n,(n+1):(m+n)]),
                              ignore_attr = TRUE,
                              tolerance = 1e-5)
  
  
})

testthat::test_that("cost function forms appropriate classes", {
  causalOT:::torch_check()
  set.seed(124123)
  n <- 10
  m <- 11
  d <- 4
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  
  # given just p
  testthat::expect_silent(cost <- causalOT:::cost(x = x, y = y, p = 2L))
  testthat::expect_equal(cost$data$dim(), 2)
  testthat::expect_equal(cost$data$shape, c(10L,11L))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costTensor"))
  testthat::expect_true(inherits(cost, "cost"))
  testthat::expect_equivalent(object = as.matrix(cost$data), 
                              expected = (as.matrix(stats::dist(rbind(x,y)))^2)[1:n,(n+1):(m+n)]/2,
                              ignore_attr = TRUE,
                              tolerance = 1e-5)
  
  
  # given cost_function
  cost_fun <- function(x,y,p) {
    n <- nrow(x)
    m <- nrow(y)
    as.matrix(stats::dist(rbind(x,y), method = "manhattan"))[1:n,(n+1):(m+n)]
  }
  testthat::expect_silent(cost <- causalOT:::cost(x = x, y = y, cost_function = cost_fun))
  testthat::expect_equal(cost$data$dim(), 2)
  testthat::expect_equal(cost$data$shape, c(10L,11L))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costTensor"))
  testthat::expect_true(inherits(cost, "cost"))
  testthat::expect_equivalent(object = as.matrix(cost$data), 
                              expected = as.matrix(stats::dist(rbind(x,y), method = "manhattan"))[1:n,(n+1):(m+n)],
                              ignore_attr = TRUE,
                              tolerance = 1e-5)
  testthat::expect_equal(object = cost$data, 
                         expected = causalOT:::costTensor$new(x = x, y = y, p = 1L)$data)
  
  # online
  testthat::skip_on_cran()
  causalOT:::rkeops_check()
  testthat::skip_on_ci()
  
  # given just p
  testthat::expect_silent(cost <- causalOT:::cost(x = x, y = y, p = 2L, tensorized = FALSE))
  testthat::expect_equal(dim(cost$data$x), c(n,d))
  testthat::expect_equal(dim(cost$data$y), c(m,d))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costOnline"))
  testthat::expect_true(inherits(cost, "cost"))
  causalOT:::rkeops_check()
  keops_sum <- rkeops::keops_kernel(
    formula = paste0("Sum_Reduction(", cost$fun, ", 0)"),
    args = c(
      paste0("X = Vi(",d,")"),
      paste0("Y = Vj(",d,")"))
  )
  testthat::expect_equivalent(object = sum(keops_sum(list(X = cost$data$x,
                                                          Y = cost$data$y) )), 
                              expected = sum((as.matrix(stats::dist(rbind(x,y)))[1:n,(n+1):(m+n)]^2) / 2),
                              ignore_attr = TRUE,
                              tolerance = 1e-5)
  
  
  # given cost_function
  cost_fun <- function(x,y) {
    n <- nrow(x)
    m <- nrow(y)
    as.matrix(stats::dist(rbind(x,y), method = "manhattan"))[1:n,(n+1):(m+n)]
  }
  testthat::expect_error(cost <- causalOT:::cost(x = x, y = y, cost_function = cost_fun, tensorized = FALSE))
  testthat::expect_silent(cost <- causalOT:::cost(x = x, y = y, cost_function = "Sum(Abs(X-Y))", tensorized = FALSE))
  testthat::expect_equal(dim(cost$data$x), c(n,d))
  testthat::expect_equal(dim(cost$data$y), c(m,d))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costOnline"))
  testthat::expect_true(inherits(cost, "cost"))
  causalOT:::rkeops_check()
  keops_sum2 <- rkeops::keops_kernel(
    formula = paste0("Sum_Reduction(", cost$fun, ", 0)"),
    args = c(
      paste0("X = Vi(",d,")"),
      paste0("Y = Vj(",d,")"))
  )
  testthat::expect_equivalent(object = sum(keops_sum2(list(X = cost$data$x,
                                                           Y = cost$data$y) )), 
                              expected = sum(as.matrix(stats::dist(rbind(x,y), method = "manhattan"))[1:n,(n+1):(m+n)]),
                              ignore_attr = TRUE,
                              tolerance = 1e-5)
  
  
})