testthat::test_that("SCM object forms", {
  
  testthat::skip_if_not_installed("osqp")
  
  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  lambda <- 10
  wlist <- list(a,a)
  
  testthat::expect_output(scm <- causalOT:::SCM$new(source = x, target = y,
                      options = list(NULL)))
  
  testthat::expect_output(scm <- causalOT:::SCM$new(source = x, target = y,
                                                    a = a, b = b,
                                                    options = list(NULL)))
  
  testthat::expect_output(scm <- causalOT:::SCM$new(source = x, target = y,
                                                    options = list(polish = TRUE)))
  
  testthat::expect_true(scm$.__enclos_env__$private$solver$GetParams()$polish)
  
})

testthat::test_that("SCM object optimizes", {
  
  testthat::skip_if_not_installed("osqp")
  
  set.seed(234808)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  lambda <- 10
  wlist <- list(a,a)
  
  testthat::expect_output(scm <- causalOT:::SCM$new(source = x, target = y,
                                                    options = list(NULL)))
  
  testthat::expect_output(w <- scm$solve())
  
  
})
