testthat::test_that("Likelihood objects form", {
  
  set.seed(20348)
  n <- 15
  d <- 3
  x <- matrix(stats::rnorm(n*d), n, d)
  z <- rbinom(n, 1, prob = 0.5)
  weights <- rep(1/n,n)
  
  data <-causalOT::dataHolder(x = x, z = z)
  
  
  testthat::expect_silent(likm <- causalOT:::likelihoodMethods(data,
                               estimand = "ATT",
                               method = "Logistic",
                               options = list(NULL)))
  
  testthat::expect_silent(likm <- causalOT:::likelihoodMethods(data,
                 estimand = "ATT",
                 method = "Logistic"))
  
  testthat::expect_true(all(slotNames(likm) %in% c("data", "estimand",
                         "method", "solver", "solver.options")))
  
  
  testthat::expect_silent(likm <- causalOT:::likelihoodMethods(data,
                               estimand = "ATT",
                               method = "CBPS",
                               options = list(NULL)))
  
  testthat::expect_silent(likm <- causalOT:::likelihoodMethods(data,
                               estimand = "ATT",
                               method = "CBPS"))
  
  testthat::expect_true(all(slotNames(likm) %in% c("data", "estimand",
                                                   "method", "solver", "solver.options")))
  
  
})

testthat::test_that("CBPS works", {
  testthat::skip_if_not_installed("CBPS")
  
  set.seed(20348)
  n <- 15
  d <- 3
  x <- matrix(stats::rnorm(n*d), n, d)
  z <- rbinom(n, 1, prob = 0.5)
  weights <- rep(1/n,n)
  
  data <-causalOT::dataHolder(x = x, z = z)
  likm <- causalOT:::likelihoodMethods(data,
                                       estimand = "ATT",
                                       method = "CBPS",
                                       options = list(NULL))
  
  mess <- testthat::capture_output(test <- likm@solver(likm))
  testthat::expect_equal(likm@solver, causalOT:::cbps_optimizer)
  
})

testthat::test_that("logistic works", {
  
  set.seed(20348)
  n <- 15
  d <- 3
  x <- matrix(stats::rnorm(n*d), n, d)
  z <- rbinom(n, 1, prob = 0.5)
  weights <- rep(1/n,n)
  
  data <-causalOT::dataHolder(x = x, z = z)
  
  testthat::expect_silent({
  likm <- causalOT:::likelihoodMethods(data,
                                       estimand = "ATT",
                                       method = "Logistic",
                                       options = list(NULL))
  
  test <- likm@solver(likm)
  
  
  likm <- causalOT:::likelihoodMethods(data,
                                       estimand = "ATC",
                                       method = "Logistic",
                                       options = list(NULL))
  
  test <- likm@solver(likm)
  
  likm <- causalOT:::likelihoodMethods(data,
                                       estimand = "ATE",
                                       method = "Logistic",
                                       options = list(NULL))
  
  test <- likm@solver(likm)
  })
})

testthat::test_that("probit works", {
  
  set.seed(20348)
  n <- 15
  d <- 3
  x <- matrix(stats::rnorm(n*d), n, d)
  z <- rbinom(n, 1, prob = 0.5)
  weights <- rep(1/n,n)
  
  data <-causalOT::dataHolder(x = x, z = z)
  
  testthat::expect_silent({
    likm <- causalOT:::likelihoodMethods(data,
                                         estimand = "ATT",
                                         method = "Probit",
                                         options = list(NULL))
    
    test <- likm@solver(likm)
    
    
    likm <- causalOT:::likelihoodMethods(data,
                                         estimand = "ATC",
                                         method = "Probit",
                                         options = list(NULL))
    
    test <- likm@solver(likm)
    
    likm <- causalOT:::likelihoodMethods(data,
                                         estimand = "ATE",
                                         method = "Probit",
                                         options = list(NULL))
    
    test <- likm@solver(likm)
  })
  
  fit <- glm(data@z ~ data@x, family = binomial(link = "probit"))
  
  probs <- predict(fit, type = "response")
  
  p1 <- probs[data@z==1]
  p0 <- 1-probs[data@z==0]
  testthat::expect_equal(test$w1, 1/p1 * 1/n)
  testthat::expect_equal(test$w0, 1/p0 * 1/n)
  
})

