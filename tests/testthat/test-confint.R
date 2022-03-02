testthat::test_that("confint works for logistic causaleffect", {
  set.seed(2323420)
  sim <- LaLonde$new(design = "NSW")
  sim$gen_data()
  
  weight <- calc_weight(sim, constraint = 0, estimand = "ATT",
                        method = "Logistic")
  
  est <- estimate_effect(data = sim, weights = weight,
                         estimand = "ATT", model = "lm", doubly.robust = FALSE,
                         matched = FALSE)
  
  # debugonce(ci_boot_ce)
  ci <- confint.causalEffect(est, level = 0.95, method = "bootstrap", n.boot = 10)
  testthat::expect_equivalent(ci$CI, c(915.7525, 3692.7361 ))
  testthat::expect_equivalent(ci$SD, 945.1653, tol = 1e-3)
})

testthat::test_that("confint works for sbw causaleffect", {
  testthat::skip_on_cran()
  set.seed(2323420)
  sim <- LaLonde$new(design = "NSW")
  sim$gen_data()
  
  weight <- calc_weight(sim, constraint = 0.05, estimand = "ATT",
                        method = "SBW", solver = "osqp",
                        adaptive_rho_interval = 10)
  
  est <- estimate_effect(data = sim, weights = weight,
                         estimand = "ATT", model = "lm", doubly.robust = FALSE,
                         matched = FALSE)
  
  # debugonce(ci_boot_ce)
  
  ci <- causalOT:::confint.causalEffect(est, level = 0.95, method = "bootstrap", n.boot = 10)
  testthat::expect_equivalent(ci$CI, c(1253.585, 3672.608    ), tol = 1e-3)
  testthat::expect_equivalent(ci$SD, 858.3696, tol = 1e-1)
})

# testthat::test_that("confint works for c wass causaleffect", {
#   testthat::skip_on_cran()
#   set.seed(2323420)
#   sim <- LaLonde$new(design = "NSW")
#   sim$gen_data()
#   
#   weight <- calc_weight(sim, constraint = 2.1, estimand = "ATT",
#                         metric = "mahalanobis", power = 2,
#                         method = "Constrained Wasserstein", solver = "osqp")
#   
#   est <- estimate_effect(data = sim, weights = weight, 
#                          estimand = "ATT", model = "lm", doubly.robust = FALSE,
#                          matched = FALSE)
#   
#   # debugonce(ci_boot_ce)
#   ci <- confint.causalEffect(est, level = 0.95, method = "bootstrap", n.boot = 10)
#   testthat::expect_equivalent(ci$CI, c(832, 2580  ), tol = 1)
#   testthat::expect_equivalent(ci$SD, 621.9692, tol = 1e-3)
# })

testthat::test_that("confint works for NNM causaleffect", {
  set.seed(2323420)
  sim <- LaLonde$new(design = "NSW")
  sim$gen_data()
  
  weight <- calc_weight(sim, estimand = "ATT",
                        metric = "mahalanobis", power = 2,
                        method = "NNM", solver = "osqp")
  
  est <- estimate_effect(data = sim, weights = weight, 
                         estimand = "ATT", model = "lm", doubly.robust = FALSE,
                         matched = FALSE)
  
  # debugonce(ci_boot_ce)
  ci <- confint.causalEffect(est, level = 0.95, method = "bootstrap", n.boot = 10)
  testthat::expect_equivalent(ci$CI, c(910, 2512), tol = 1)
  testthat::expect_equivalent(ci$SD, 554, tol = 1)
})

testthat::test_that("confint calibrated for NNM causaleffect", {
  testthat::skip("Interactive only")
  set.seed(2323420)
  sim <- Hainmueller$new(overlap = "low", design = "B",
                         n = 2^9, p = 6)
  sim$gen_data()
  
  weight <- calc_weight(sim, estimand = "ATE",
                        metric = "Lp", power = 2,
                        method = "NNM", solver = "osqp")
  
  est <- estimate_effect(data = sim, weights = weight, 
                         estimand = "ATE", model = "lm", doubly.robust = FALSE,
                         matched = FALSE)
  
  # debugonce(ci_boot_ce)
  ci <- confint.causalEffect(est, level = 0.95, method = "bootstrap", n.boot = 1000,
                             verbose = TRUE)
  print(ci)
  # ATE: boot est 0.8552763 sd = 0.3904376 CI (0.08545122, 1.60667010 )
  # ATE: true est 0.859 true sd ~ 0.549, true CI (-0.166, 1.96)
  
  # ATE DR boot: est -0.1637187 sd = 0.67,  (-1.455522,  1.079845)
  #ATE DR: est 0.572 sd 0.681 (-0.69, 1.9 )
  
  set.seed(2323420)
  sim <- Hainmueller$new(overlap = "low", design = "B",
                         n = 2^9, p = 6)
  sim$gen_data()
  
  weight <- calc_weight(sim, estimand = "ATT",
                        metric = "Lp", power = 2,
                        method = "NNM", solver = "osqp")
  
  est <- estimate_effect(data = sim, weights = weight, 
                         estimand = "ATT", model = "lm", doubly.robust = FALSE,
                         matched = FALSE)
  
  # debugonce(ci_boot_ce)
  ci <- confint.causalEffect(est, level = 0.95, method = "bootstrap", n.boot = 1000,
                             verbose = TRUE)
  # ATT boot: 0.5561051 2.565081 4.750104
  # ATT: true est 2.84 true sd 0.770, CI 1.43, 4.34
})

testthat::test_that("confint calibrated for logistic causaleffect", {
  testthat::skip("Interactive only")
  set.seed(2323420)
  sim <- Hainmueller$new(overlap = "low", design = "B",
                         n = 2^9, p = 6)
  sim$gen_data()
  
  weight <- calc_weight(sim, estimand = "ATE",
                        method = "Logistic")
  
  est <- estimate_effect(data = sim, weights = weight, 
                         estimand = "ATE", model = "lm", doubly.robust = FALSE,
                         matched = FALSE)
  
  # debugonce(ci_boot_ce)
  ci <- confint.causalEffect(est, level = 0.95, method = "bootstrap", n.boot = 1000,
                             verbose = FALSE)
  testthat::expect_equivalent(ci$CI, c(-6.538542,  1.690153 ), tol = 1e-3)
  testthat::expect_equivalent(ci$SD, 2.302784, tol = 1e-3)
  # true est 0.115, sd ~ 1.72, true CI (-3.54, 2.97)
  # boot est 0.9151212 sd 0.9915984 c(-0.9962013  2.7655312)
})

testthat::test_that("confint works for NNM causaleffect", {
  testthat::skip("Interactive only")
  set.seed(2323420)
  sim <- LaLonde$new(design = "Full")
  sim$gen_data()
  
  weight <- calc_weight(sim, estimand = "ATT",
                        metric = "mahalanobis", power = 1,
                        method = "NNM", solver = "osqp")
  
  est <- estimate_effect(data = sim, weights = weight, 
                         estimand = "ATT", model = "lm", doubly.robust = FALSE,
                         matched = FALSE)
  
  # debugonce(ci_boot_ce)
  ci <- confint.causalEffect(est, level = 0.95, method = "bootstrap", n.boot = 1000,
                             verbose = TRUE)
  testthat::expect_equivalent(ci$CI, c(396.3799, 3185.4339  ), tol = 1e-3)
  testthat::expect_equivalent(ci$SD, 706.7932, tol = 1e-3)
  #true sd = 634.4388, ci = [551, 3038].
})
