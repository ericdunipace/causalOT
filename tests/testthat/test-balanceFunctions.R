testthat::test_that("balanceFunction class forms",{
  
  set.seed(20348)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  
  testthat::expect_silent(bf <- causalOT:::balanceFunction$new(source = x, 
                            target = y, 
                            a = a, b = b, 
                            delta = delta))
  
  testthat::expect_silent(ebw <- causalOT:::EntropyBW$new(source = x, 
                                                          target = y, 
                                                          a = a, b = b))
  
  testthat::skip_if_not_installed("osqp")
  testthat::expect_output(sbw <- causalOT:::SBW$new(source = x, 
                 target = y, 
                 a = a, b = b))
  
  
})

testthat::test_that("balanceFunction optimizes",{
  set.seed(20348)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  
  #L2
  testthat::skip_if_not_installed("osqp")
  testthat::expect_output(sbw <- causalOT:::SBW$new(source = x, 
                                                    target = y, 
                                                    a = a, b = b))
  
  delta <- sbw$gridInit(NULL, 7L)
  testthat::expect_output(w0 <- sbw$solve(delta[1], NULL))
  testthat::expect_equal(w0, a, tol = 1e-2)
  
  testthat::expect_output(w1 <- sbw$solve(delta[7], NULL))
  testthat::expect_true(is.character(all.equal(w1,a)))
  
  # entropy
  testthat::skip_if_not_installed("lbfgsb3c")
  testthat::expect_silent(ebw <- causalOT:::EntropyBW$new(source = x, 
                                                          target = y, 
                                                          a = a, b = b))
  delta <- ebw$gridInit(NULL, 7L)
  testthat::expect_warning(w0_e <- ebw$solve(delta[1], NULL))
  testthat::expect_equal(c(w0_e), a)
  testthat::expect_equal(c(w0_e), w0)
  
  testthat::expect_silent(w1_e <- ebw$solve(delta[7], NULL))
  testthat::expect_true(all(c(crossprod(ebw$A, w1_e)) < 1e-3))
  
  
})

testthat::test_that("functions catch errors", {
  set.seed(20348)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(0, n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  delta <- 0.5
  
  testthat::skip_if_not_installed("lbfgsb3c")
  testthat::expect_error(ebw <- causalOT:::EntropyBW$new(source = x, 
                                                         target = y, 
                                                         a = a, b = b))
  testthat::expect_silent(ebw <- causalOT:::EntropyBW$new(source = y, 
                                                          target = x, 
                                                          a = b, b = a))
  
  testthat::skip_if_not_installed("osqp")
  testthat::expect_error(sbw <- causalOT:::SBW$new(source = x, 
                                                    target = y, 
                                                    a = a, b = b))
  testthat::expect_output(sbw <- causalOT:::SBW$new(source = y, 
                                                   target = x, 
                                                   a = b, b = a))
  
  
})
