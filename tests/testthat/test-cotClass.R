testthat::test_that("COT objects form", {
  causalOT:::torch_check()
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
  
  #cot defaults
  testthat::expect_silent(causalOT:::COT$new(source = x, target = y))
  
  #L2 setup
  testthat::expect_output(cot <- causalOT:::COT$new(source = x, target = y,
                            options = list(penalty = "L2")))
  testthat::expect_silent(cot$weight <- c(1, rep(0,n-1)))
  cot$penalty <- 100
  testthat::expect_equal(100, cot$ot$penalty)
  testthat::expect_true(is.null(cot$.__enclos_env__$private$torch_optim$fun))
  cot$penalty <- c("lambda" = 10)
  testthat::expect_equivalent(10, cot$ot$penalty)
  testthat::expect_null(cot$bf$delta)
  testthat::expect_silent(cot$penalty <- c("lambda" = 10, "delta" = 0.5))
  testthat::expect_null(cot$bf$delta)
  testthat::expect_error(cot$penalty <- c("lambda" = 10, "carbon" = 0.5))
  testthat::expect_warning(cot$penalty <- c(0, 0.4))
  testthat::expect_null(cot$bf$delta)
  testthat::expect_equal(as.numeric(cot$weight), c(1, rep(0,n-1)))
  testthat::expect_output(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "L2",
                                                                   balance.formula = "~X1+X2")))
  testthat::expect_true(inherits(cot$bf, "balanceFunction"))
  
  # L2 with debias, expect warning
  testthat::expect_output(testthat::expect_warning(causalOT:::COT$new(source = x, target = y,
                                            options = list(penalty = "L2",
                                                           debias = TRUE))))
  
  # entropy debias
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                            options = list(penalty = "entropy",
                                                           debias = TRUE)))
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   debias = TRUE,
                                                                   balance.formula = "~X1 + X2")))
  testthat::expect_true(inherits(cot$bf, "balanceFunction"))
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   lambda = Inf,
                                                                   debias = TRUE)))
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   lambda = 0,
                                                                   debias = TRUE)))
  
  # entropy not debiased
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   debias = FALSE)))
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   debias = FALSE,
                                                                   balance.formula = "~X1*X3")))
  testthat::expect_true(inherits(cot$bf, "balanceFunction"))
  
  # debias without torch optim should throw error
  testthat::expect_error(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   debias = TRUE,
                                                                   torch.optimizer = NULL)
                                                    ))
  
  
})

testthat::test_that("L2 works", {
  causalOT:::torch_check()
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
  
  testthat::expect_output(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "L2")))
  
  testthat::expect_output(w <- cot$solve(10, wlist[1]))
  testthat::expect_output(w <- cot$solve(1e7, wlist[1]))
  testthat::expect_equal(cot$ot$penalty, 1e7)
  testthat::expect_equal(a, as.numeric(w), tol = 1e-4)
  
  testthat::expect_output(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "L2",
                                                                   balance.formula = "~X1+X2")))
  testthat::expect_output(w <- cot$solve(list(lambda = 1, delta = 0.1), wlist[1]))
  testthat::expect_true(max(crossprod(cot$bf$source, as.numeric(w))) < 0.1)
  
})

testthat::test_that("ent works", {
  causalOT:::torch_check()
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
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   debias = FALSE)))
  
  testthat::expect_silent(w <- cot$solve(1, wlist[1]))
  testthat::expect_equal(cot$ot$penalty, 1)
  testthat::expect_silent(w <- cot$solve(10, wlist[1])) 
  testthat::expect_equal(cot$ot$penalty, 1e1)
  testthat::expect_silent(w <- cot$solve(1e7, wlist[1]))
  testthat::expect_equal(cot$ot$penalty, 1e7)
  testthat::expect_equal(a, as.numeric(w), tol = 1e-4)
  testthat::expect_silent(w <- cot$solve(Inf, wlist[1]))
  testthat::expect_equal(cot$ot$penalty, Inf)
  testthat::expect_equal(a, as.numeric(w))
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   balance.formula = "~X1",
                                                                   debias = FALSE)))
  
  testthat::expect_silent(w <- cot$solve(c(lambda = 10, delta = .05), 
                                         wlist[1])) 
  testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target_vector)) < 0.051)
  testthat::expect_silent(w <- cot$solve(c(lambda = 1, delta = .05), 
                                         wlist[1])) # left off checking this
  testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target_vector)) < 0.05)
  testthat::expect_equal(cot$ot$penalty, 1)
  testthat::expect_equal(cot$bf$delta, 0.05)
  testthat::expect_silent(w <- cot$solve(c(lambda = 1e7, 
                                           delta = .025), wlist[1]))
  testthat::expect_equal(cot$ot$penalty, 1e7)
  testthat::expect_equal(cot$bf$delta, .025)
  testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target_vector)) < 0.025)
  testthat::expect_silent(w <- cot$solve(c(lambda = .1, 
                                           delta = .025), wlist[1]))
  testthat::expect_equal(cot$ot$penalty, 0.1)
  testthat::expect_equal(cot$bf$delta, .025)
  testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target_vector)) < 0.025)
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy", debias = FALSE,
                                                                   balance.formula = "~X1+X2")))
  testthat::expect_silent(w <- cot$solve(list(lambda = 1, delta = 0.05), wlist[1]))
  testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target_vector)) < 0.06)
  
  testthat::expect_silent(w <- cot$solve(list(lambda = 1e7, delta = 0.05), wlist[1]))
  testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target_vector)) < 0.06)
  testthat::expect_warning(w <- cot$solve(list(lambda = 1e-3, delta = 0.05), wlist[1]))
  testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target_vector)) < 0.4)
  
})

testthat::test_that("ent debiased works, tensor", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
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
  
  testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   debias = TRUE,
                                                                   niter = 2L,
                                                                   torch.optimizer = torch::optim_lbfgs)))
  
  
  testthat::expect_silent(w <- cot$solve(1, NULL))
  testthat::expect_silent(w2 <- cot$solve(1, wlist[1]))
  testthat::expect_equal(w,w2)
  testthat::expect_equal(cot$ot$penalty, 1)
  testthat::expect_silent(w <- cot$solve(10, wlist[1])) 
  testthat::expect_equal(cot$ot$penalty, 1e1)
  testthat::expect_silent(w <- cot$solve(1e7, wlist[1]))
  testthat::expect_equal(cot$ot$penalty, 1e7)
  testthat::expect_silent(w <- cot$solve(Inf, wlist[1]))
  testthat::expect_equal(cot$ot$penalty, Inf)
  testthat::expect_silent(w <- cot$solve(0, wlist[1]))
  testthat::expect_equal(cot$ot$penalty, 0)
  cot$weight <- wlist[[1]]
  testthat::expect_silent(w <- cot$solve(1, wlist[1]))
  w <- cot$solve(.001, wlist[1])
  
  
  # currently no BF
  # testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
  #                                                   options = list(penalty = "entropy",
  #                                                                  balance.formula = "~X1",
  #                                                                  debias = TRUE)))
  # 
  # testthat::expect_silent(w <- cot$solve(c(lambda = 10, delta = .05), 
  #                                        wlist[1])) # left off checking this
  # testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target)) < 0.051)
  # testthat::expect_silent(w <- cot$solve(c(lambda = 1, delta = .05), 
  #                                        wlist[1])) # left off checking this
  # testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target)) < 0.05)
  # testthat::expect_equal(cot$ot$penalty, 1)
  # testthat::expect_equal(cot$bf$delta, 0.05)
  # testthat::expect_silent(w <- cot$solve(c(lambda = 1e7, 
  #                                          delta = .025), wlist[1]))
  # testthat::expect_equal(cot$ot$penalty, 1e7)
  # testthat::expect_equal(cot$bf$delta, .025)
  # testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target)) < 0.025)
  # testthat::expect_silent(w <- cot$solve(c(lambda = .1, 
  #                                          delta = .025), wlist[1]))
  # testthat::expect_equal(cot$ot$penalty, 0.1)
  # testthat::expect_equal(cot$bf$delta, .025)
  # testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target)) < 0.025)
  # 
  # testthat::expect_silent(cot <- causalOT:::COT$new(source = x, target = y,
  #                                                   options = list(penalty = "entropy", debias = FALSE,
  #                                                                  balance.formula = "~X1+X2")))
  # testthat::expect_silent(w <- cot$solve(list(lambda = 1, delta = 0.05), wlist[1]))
  # testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target)) < 0.06)
  # 
  # testthat::expect_silent(w <- cot$solve(list(lambda = 1e7, delta = 0.05), wlist[1]))
  # testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target)) < 0.06)
  # testthat::expect_silent(w <- cot$solve(list(lambda = 1e-3, delta = 0.05), wlist[1]))
  # testthat::expect_true(max(abs(crossprod(cot$bf$source, as.numeric(w)) - cot$bf$target)) < 0.06)
  
})

testthat::test_that("ent debiased works, online", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
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
  
  mess<- testthat::capture_message(cot <- causalOT:::COT$new(source = x, target = y,
                      options = list(penalty = "entropy",
                                     cost.online = "online",
                                     debias = TRUE,
                                     niter = 1L,
                                     torch.optimizer = torch::optim_lbfgs)))
  
  testthat::expect_silent(cot2 <- causalOT:::COT$new(source = x, target = y,
                                                    options = list(penalty = "entropy",
                                                                   debias = TRUE,
                                                                   niter = 1L,
                                                                   torch.optimizer = torch::optim_lbfgs)))
  
  mess <- testthat::capture_output(w <- cot$solve(Inf, wlist[1]))
  mess <- testthat::capture_output(w2 <- cot2$solve(Inf, wlist[1]))
  testthat::expect_equal(cot$ot$penalty, Inf)
  testthat::expect_equal(w,w2, tol = 1e-4)
  mess <- testthat::capture_output(w <- cot$solve(0, wlist[1]))
  mess <- testthat::capture_output(w2 <- cot2$solve(0, wlist[1]))
  testthat::expect_equal(cot$ot$penalty, 0)
  testthat::expect_equal(w,w2)
  
})
