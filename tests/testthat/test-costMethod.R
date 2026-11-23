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
  causalOT:::torch_check()
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
  ones <- torch::torch_ones(m, dtype = cost$data$y$dtype)
  testthat::expect_equivalent(object = sum(as.numeric(causalOT:::online_cost_mv(cost, cost$data$x, cost$data$y, ones))),
                              expected = sum((as.matrix(stats::dist(rbind(x,y)))[1:n,(n+1):(m+n)]^2) / 2),
                              ignore_attr = TRUE,
                              tolerance = 1e-5)


  # given cost_function
  cost_fun <- function(x, y, p) torch::torch_cdist(x, y, p = 1)
  testthat::expect_silent(cost <- causalOT:::costOnline$new(x = x, y = y, cost_function = cost_fun))
  testthat::expect_equal(dim(cost$data$x), c(n,d))
  testthat::expect_equal(dim(cost$data$y), c(m,d))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costOnline"))
  testthat::expect_true(inherits(cost, "cost"))
  testthat::expect_equivalent(object = sum(as.numeric(causalOT:::online_cost_mv(cost, cost$data$x, cost$data$y, ones))),
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

  # given just p
  testthat::expect_silent(cost <- causalOT:::cost(x = x, y = y, p = 2L, tensorized = FALSE))
  testthat::expect_equal(dim(cost$data$x), c(n,d))
  testthat::expect_equal(dim(cost$data$y), c(m,d))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costOnline"))
  testthat::expect_true(inherits(cost, "cost"))
  ones <- torch::torch_ones(m, dtype = cost$data$y$dtype)
  testthat::expect_equivalent(object = sum(as.numeric(causalOT:::online_cost_mv(cost, cost$data$x, cost$data$y, ones))),
                              expected = sum((as.matrix(stats::dist(rbind(x,y)))[1:n,(n+1):(m+n)]^2) / 2),
                              ignore_attr = TRUE,
                              tolerance = 1e-5)


  # given cost_function
  cost_fun <- function(x, y, p) torch::torch_cdist(x, y, p = 1)
  testthat::expect_error(cost <- causalOT:::cost(x = x, y = y, cost_function = "Sum(Abs(X-Y))", tensorized = FALSE))
  testthat::expect_silent(cost <- causalOT:::cost(x = x, y = y, cost_function = cost_fun, tensorized = FALSE))
  testthat::expect_equal(dim(cost$data$x), c(n,d))
  testthat::expect_equal(dim(cost$data$y), c(m,d))
  testthat::expect_true(inherits(cost, "R6"))
  testthat::expect_true(inherits(cost, "costOnline"))
  testthat::expect_true(inherits(cost, "cost"))
  testthat::expect_equivalent(object = sum(as.numeric(causalOT:::online_cost_mv(cost, cost$data$x, cost$data$y, ones))),
                              expected = sum(as.matrix(stats::dist(rbind(x,y), method = "manhattan"))[1:n,(n+1):(m+n)]),
                              ignore_attr = TRUE,
                              tolerance = 1e-5)
  
  
})