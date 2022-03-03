testthat::test_that("function works as would expect", {
  set.seed(2340802)
  
  data <- causalOT::LaLonde$new(design = "NSW")
  data$gen_data()
  testthat::expect_equal(class(data$get_x()), c("matrix","array"))
  testthat::expect_equal(class(data$get_z()), "integer")
  testthat::expect_equal(class(data$get_y()), "numeric")
  testthat::expect_silent(ns <- causalOT:::get_n.DataSim(data))
  testthat::expect_silent(ps <- causalOT:::get_p.DataSim(data))
  
  testthat::expect_equal(ps, 10)
  testthat::expect_equal(ns, c(n0 = 260, n1 = 185))
  
  data <- causalOT::LaLonde$new(design = "Full")
  data$gen_data()
  testthat::expect_equal(class(data$get_x()), c("matrix","array"))
  testthat::expect_equal(class(data$get_z()), "integer")
  testthat::expect_equal(class(data$get_y()), "numeric")
  testthat::expect_silent(ns <- causalOT:::get_n.DataSim(data))
  testthat::expect_silent(ps <- causalOT:::get_p.DataSim(data))
  
  testthat::expect_equal(ps, 10)
  testthat::expect_equal(ns, c(n0 = 15992, n1 = 185))
  
  testthat::expect_equal(data$get_tau(), 1794)
})

testthat::test_that("optimal weighting works, no augmentation", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  set.seed(6464546)
  n <- 2^7
  p <- 6
  nsims <- 1
  design <- "NSW"
  distance <- c("Lp")
  power <- c(2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment = FALSE
  
  #### get simulation functions ####
  data <- causalOT::LaLonde$new()
  data$gen_data()
  
  # debugonce( data$opt_weight)
  opt_weights_mosek  <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "mosek"))
  # opt_weights_gurobi <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "gurobi"))
  # opt_weights_cplex  <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "cplex"))
  names(opt_weights_mosek) <- estimates
    # names(opt_weights_gurobi) <- estimates
    # names(opt_weights_cplex) <- estimates
  
  testthat::expect_equivalent(estimate_effect(data, weights = opt_weights_mosek[[4]], 
                                              doubly.robust = FALSE, estimand = "ATE")$estimate,
                              mean(data$get_tau()), tol = 1e-3)
  # testthat::expect_equivalent(estimate_effect(data, weights = opt_weights_gurobi[[4]], 
  #                                             doubly.robust = FALSE, estimand = "ATE")$estimate,
  #                             mean(data$get_tau()), tol = 1e-3)
  # testthat::expect_equivalent(estimate_effect(data, weights = opt_weights_cplex[[4]], 
  #                                             doubly.robust = FALSE, estimand = "ATE")$estimate,
  #                             mean(data$get_tau()), tol = 1e-3)
})

testthat::test_that("optimal weighting comparison works, no augmentation", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  set.seed(9847)
  design <- "NSW"
  distance <- c("Lp")
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment <- FALSE
  power <- 2
  
  #### get simulation functions ####
  data <- causalOT::LaLonde$new(design = design)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  testthat::expect_warning(weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 8, 
                                                       estimand = e, 
                                                       p = power,
                                                       method = "NNM",
                                                       solver = "mosek")))
  
  
  
  # opt_weights_mosek <- data$opt_weight(estimand = "ATT", augment = augment, solver = "mosek")
  # opt_weights_gurobi<- data$opt_weight(estimand = "ATE", augment = augment, solver = "gurobi")
  # opt_weights_cplex <- data$opt_weight(estimand = "ATE", augment = augment, solver = "cplex")
  
  # debugonce(data$opt_weight_dist)
  compare_mosek <- mapply(data$opt_weight_dist, weight = weights, estimand = estimates,
                          augment = augment, solver = "mosek")
  for (cc in compare_mosek) testthat::expect_lt(cc, 0.25)
})

testthat::test_that("optimal weighting comparison works. augmentation", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  set.seed(9847)
  design <- "NSW"
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment <- TRUE
  power <- 2
  
  #### get simulation functions ####
  data <- causalOT::LaLonde$new(design = design)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  testthat::expect_warning(weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 8, 
                                                       estimand = e, 
                                                       p = power,
                                                       method = "NNM",
                                                       solver = "mosek")))
  
  
  
  # opt_weights_mosek <- data$opt_weight(estimand = "ATT", augment = augment, solver = "mosek")
  # opt_weights_gurobi<- data$opt_weight(estimand = "ATE", augment = augment, solver = "gurobi")
  # opt_weights_cplex <- data$opt_weight(estimand = "ATE", augment = augment, solver = "cplex")
  
  # debugonce(data$opt_weight_dist)
  compare_mosek <- mapply(data$opt_weight_dist, weight = weights, estimand = estimates,
                          augment = augment, solver = "mosek")
  for (cc in compare_mosek) testthat::expect_lt(cc, 0.55)
})
