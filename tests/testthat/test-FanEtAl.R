testthat::test_that("function works as would expect", {
  set.seed(2340802)
  
  n <- 500
  p1 <- 4
  p <- 100
  design <- "A"

  data <- causalOT:::FanEtAl$new(n = n, p = 100, numActive = p1,
                                    design = design)
  data$gen_data()
  testthat::expect_equal(class(data$get_x()), c("matrix","array"))
  testthat::expect_equal(class(data$get_z()), "numeric")
  testthat::expect_equal(class(data$get_y()), "numeric")
  testthat::expect_silent(ns <- causalOT:::get_n.DataSim(data))
  testthat::expect_silent(ps <- causalOT:::get_p.DataSim(data))
  
  testthat::expect_equal(ps, 100)
  testthat::expect_equal(ns, c(n0 = 108,n1 = 392))
  
  design <- "A"
  p1 <- 50
  data <- causalOT:::FanEtAl$new(n = n, p = 100, numActive = p1,
                                    design = design)
  data$gen_data()
  testthat::expect_equal(class(data$get_x()), c("matrix","array"))
  testthat::expect_equal(class(data$get_z()), "numeric")
  testthat::expect_equal(class(data$get_y()), "numeric")
  testthat::expect_silent(ns <- causalOT:::get_n.DataSim(data))
  testthat::expect_silent(ps <- causalOT:::get_p.DataSim(data))
  
  
  design <- "B"
  p1 <- 4
  data <- causalOT:::FanEtAl$new(n = n, p = 100, numActive = p1,
                                    design = design)
  data$gen_data()
  testthat::expect_equal(class(data$get_x()), c("matrix","array"))
  testthat::expect_equal(class(data$get_z()), "numeric")
  testthat::expect_equal(class(data$get_y()), "numeric")
  testthat::expect_silent(ns <- causalOT:::get_n.DataSim(data))
  testthat::expect_silent(ps <- causalOT:::get_p.DataSim(data))
  
})

testthat::test_that("optimal weighting works, no augmentation", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  set.seed(6464546)
  nsims <- 1
  distance <- c("Lp")
  power <- c(2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment = FALSE
  
  #### get simulation functions ####
  testthat::expect_warning(data <- causalOT:::FanEtAl$new())
  data$gen_data()
  
  # debugonce( data$opt_weight)
  opt_weights_mosek <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "mosek"))
  # opt_weights_gurobi<- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "gurobi"))
  # opt_weights_cplex <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "cplex"))
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

testthat::test_that("optimal weighting works, augmentation", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  set.seed(6464546)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment = TRUE
  
  #### get simulation functions ####
  data <- causalOT:::FanEtAl$new(design = design, p = 100, n = 500,
                                numActive = 4)
  data$gen_data()
  
  # debugonce( data$opt_weight)
  opt_weights_mosek <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "mosek"))
  # opt_weights_gurobi<- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "gurobi"))
  # opt_weights_cplex <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "cplex"))
  names(opt_weights_mosek) <- estimates
    # names(opt_weights_gurobi) <- estimates
    # names(opt_weights_cplex) <- estimates
  
  testthat::expect_equivalent(estimate_effect(data, weights = opt_weights_mosek[[4]], 
                                              doubly.robust = augment, estimand = "ATE")$estimate,
                              mean(data$get_tau()), tol = 1e-3)
  # testthat::expect_equivalent(estimate_effect(data, weights = opt_weights_gurobi[[4]], 
  #                                             doubly.robust = augment, estimand = "ATE")$estimate,
  #                             mean(data$get_tau()), tol = 1e-1)
  # testthat::expect_equivalent(estimate_effect(data, weights = opt_weights_cplex[[4]], 
  #                                             doubly.robust = augment, estimand = "ATE")$estimate,
  #                             mean(data$get_tau()), tol = 1e-3)
})

testthat::test_that("optimal weighting comparison works, no augmentation", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  set.seed(9847)
  design <- "B"
  distance <- c("Lp")
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment <- FALSE
  power = 2
  
  #### get simulation functions ####
  data <- causalOT:::FanEtAl$new(design = design, n = 500, p = 100,
                                numActive = 4)
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
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  set.seed(9847)
  design <- "B"
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment <- TRUE
  power <- 2
  
  #### get simulation functions ####
  data <- causalOT:::FanEtAl$new(design = design, n = 500, p = 100,
                                numActive = 4)
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
