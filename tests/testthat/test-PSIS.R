testthat::test_that("PSIS diagnostics work", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  # testthat::skip_if_not_installed("gurobi")
  library(causalOT)
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC","ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p,
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  
  cost  <- cost_mahalanobis(data$get_x0(), data$get_x1(), 1)
  cost0 <- cost_mahalanobis(data$get_x(), data$get_x0(), 1)
  cost1 <- cost_mahalanobis(data$get_x(), data$get_x1(), 1)
  
  # debugonce(calc_weight_bal)
  # debugonce(qp_wass)
  test1 <- calc_weight(data = data,
                       dist = "mahalanobis",
                       p = 4,
                       constraint = 0.5,
                       estimand = "ATE",
                       method = "Wasserstein",
                       solver = "mosek")
  
  # debugonce(mosek_solver)
  # debugonce(calc_weight_bal)
  # test2 <- calc_weight(data = data,
  #                      p =1,
  #                      constraint = c(1),
  #                      estimand = "ATE",
  #                      method = "Wasserstein",
  #                      dist = "mahalanobis",
  #                      solver = "mosek")
  
  test3 <- calc_weight(data, constraint = 0.1, estimand = "ATE", method = "SBW", solver = "mosek")
  test4 <- calc_weight(data, constraint = 0.1, estimand =  "ATE", method = "Logistic", solver = "mosek")
  test4b <- calc_weight(data, constraint = 0.1, estimand =  "cATE", method = "Logistic", solver = "mosek")
  
  # test5 <- calc_weight(data, constraint = 0.1, estimand = "ATE", method = "RKHS", solver = "gurobi")
  # test5b <- calc_weight(data, constraint = 0.1, estimand = "cATE", method = "RKHS", solver = "gurobi")
  
  # test6 <- calc_weight(data = data,
  #                      p =1,
  #                      constraint = c(2.2),
  #                      estimand = "cATE",
  #                      method = "Constrained Wasserstein",
  #                      dist = "mahalanobis",
  #                      solver = "mosek")
  test7 <- calc_weight(data = data,
                       p = 4,
                       constraint = c(2.2),
                       estimand = "cATE",
                       method = "Wasserstein",
                       dist = "mahalanobis",
                       solver = "mosek")
  
  # test8 <- calc_weight(data = data,
  #                      p = 1,
  #                      constraint = 100,
  #                      estimand = "ATE",
  #                      method = "Constrained Wasserstein",
  #                      dist = "RKHS",
  #                      solver="mosek",
  #                      rkhs.args = list(p = test5$args$p,
  #                                       theta = test5$args$theta,
  #                                       gamma = test5$args$gamma))
  # 
  
  weights <- list(CW = test1,
                  # W = test2, 
                  SBW=test3, IPW = test4, 
                  # RKHS = test5,
                  # cRKHS = test5b, 
                  # cCW = test6, 
                  cW = test7
                  # , 
                  # rCW = test8
                  )
  testthat::expect_warning(ps <- lapply(weights, PSIS))
  testthat::expect_warning(ps.check <- PSIS(weights))
  testthat::expect_equivalent(ps, ps.check)
  
  diag.check <- lapply(ps, PSIS_diag)
  diag <- PSIS_diag(ps)
  testthat::expect_warning(diag.check2 <- PSIS_diag(weights))
  testthat::expect_equal(diag, diag.check)
  testthat::expect_equal(diag, diag.check2)
})

testthat::test_that("PSIS diagnostics work, feasible", {
  testthat::skip_on_cran()
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  # testthat::skip_if_not_installed("gurobi")
  library(causalOT)
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  # estimates <- c("ATT", "ATC","ATE", "feasible")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p,
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  
  cost  <- cost_mahalanobis(data$get_x0(), data$get_x1(), 1)
  cost0 <- cost_mahalanobis(data$get_x(), data$get_x0(), 1)
  cost1 <- cost_mahalanobis(data$get_x(), data$get_x1(), 1)
  
  # debugonce(calc_weight_bal)
  # debugonce(qp_wass)
  test1 <- calc_weight(data = data,
                       dist = "mahalanobis",
                       p = 1,
                       constraint = 0.01,
                       estimand = "feasible",
                       method = "SBW",
                       solver = "mosek")
  
  
  testthat::expect_silent(diag.check <- PSIS_diag(test1))
})

