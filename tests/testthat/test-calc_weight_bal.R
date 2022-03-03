arg.names <- c("w0",  "w1",   "gamma","args", "estimand", "method")

testthat::test_that("works for Wass", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  # testthat::skip_if_not_installed(pkg="gurobi")
  testthat::skip_if_not_installed(pkg="Rmosek")
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = 2.5, 
  #                                                          estimand = e, 
  #                                                          method = "Wasserstein",
  #                                                          solver = "gurobi"))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0),
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1),
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
  #                                                          constraint = 3, 
  #                                                          estimand = e, 
  #                                                          method = "Wasserstein",
  #                                                          solver = "cplex"))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0),
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1),
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0),
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for Wass RKHS, opt", {
  testthat::skip("RKHS method deprecated")
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("RKHS")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  rkhs.argz <- NULL
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = 1e6, 
  #                                                          estimand = e, 
  #                                                          rkhs.args = rkhs.argz,
  #                                                          metric = metric,
  #                                                          method = "Wasserstein",
  #                                                          solver = "gurobi"))
  # for (w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE,
  #                        tol = 1e-3)
  # testthat::expect_equal(rep(1/n0,n0), 
  #                        weights[[3]]$w0, check.attributes = FALSE,
  #                        tol = 1e-3)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # # testthat::expect_match(all.equal(rep(1/n1,n1), 
  # #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = 1e6, 
  #                                                          estimand = e, 
  #                                                          method = "Wasserstein",
  #                                                          rkhs.args = rkhs.argz,
  #                                                          metric = metric,
  #                                                          solver = "cplex"))
  # for (w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # # testthat::expect_match(all.equal(rep(1/n0,n0), 
  # #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # # testthat::expect_match(all.equal(rep(1/n1,n1), 
  # #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 1e6, 
                                                           estimand = e, 
                                                           rkhs.args = rkhs.argz,
                                                           metric = metric,
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  for (w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for SBW", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^8
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                                      constraint = 3, 
                                                                      estimate = e, 
                                                                      method = "SBW",
                                                                      solver = "osqp"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  # testthat::skip_if_not_installed(pkg="gurobi")
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = 3, 
  #                                                          estimate = e, 
  #                                                          method = "SBW",
  #                                                          solver = "gurobi"))
  # sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "mosek"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = 3, 
  #                                                          estimate = e, 
  #                                                          method = "SBW",
  #                                                          solver = "cplex"))
  # sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  
})

testthat::test_that("works for Wass, sample weight", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  w0 <- causalOT:::renormalize(runif(n0))
  w1 <- causalOT:::renormalize(runif(n1))
  w0[seq(1,n0,2)] <- 0
  w1[seq(1,n1,2)] <- 0
  sample_weights <- list(w0 = causalOT:::renormalize(w0),
                         w1 = causalOT:::renormalize(w1))
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                                      constraint = 3, 
                                                                      estimand = e, 
                                                                      method = "Wasserstein",
                                                                      solver = "osqp",
                                                                      sample_weight = sample_weights))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # testthat::skip_if_not_installed(pkg = "gurobi")
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = 2.5, 
  #                                                          estimand = e, 
  #                                                          method = "Wasserstein",
  #                                                          solver = "gurobi",
  #                                                          sample_weight = sample_weights))
  # for (w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = 3, 
  #                                                          estimand = e, 
  #                                                          method = "Wasserstein",
  #                                                          solver = "cplex",
  #                                                          sample_weight = sample_weights))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = "mosek",
                                                           sample_weight = sample_weights))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for Wass, variance", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 100, 
                                                           estimand = e, 
                                                           penalty = "variance",
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = 30, 
  #                                                          estimand = e, 
  #                                                          penalty = "variance",
  #                                                          method = "Wasserstein",
  #                                                          solver = "cplex"))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # # testthat::expect_match(all.equal(rep(1/n0,n0), 
  # #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # # testthat::expect_match(all.equal(rep(1/n1,n1), 
  # #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 10, 
                                                           estimand = e, 
                                                           penalty = "variance",
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for Wass, entropy", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 100, 
                                                           estimand = e, 
                                                           penalty = "entropy",
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 30, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           penalty = "entropy",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 10, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           penalty = "entropy",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for Wass divergence", {
  testthat::skip_on_cran()
  skip_if_no_geomloss()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "lbfgs"
  estimates <- c("ATT", "ATC",  "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = list(penalty = 1000), 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = solver,
                                                           add.divergence = TRUE))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0),
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  
})

testthat::test_that("works for SCM", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                                      constraint = NULL, 
                                                                      estimand = e, 
                                                                      method = "SCM",
                                                                      penalty = "none",
                                                                      solver = "osqp"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0),
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1),
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # testthat::skip_if_not_installed(pkg = "gurobi")
  # testthat::expect_warning(weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = NULL, 
  #                                                          estimand = e, 
  #                                                          method = "SCM",
  #                                                          penalty = "none",
  #                                                          solver = "gurobi")))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = 3, 
  #                                                          estimand = e, 
  #                                                          method = "Wasserstein",
  #                                                          solver = "cplex"))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # # testthat::expect_match(all.equal(rep(1/n0,n0), 
  # #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # # testthat::expect_match(all.equal(rep(1/n1,n1), 
  # #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for SCM, penalties", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  testthat::expect_silent(weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = list(penalty = 10), 
                                                           estimand = e, 
                                                           method = "SCM",
                                                           penalty = "entropy",
                                                           solver = "mosek")))
  for (w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0),
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1),
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = list(penalty = 0), 
  #                                                          estimand = e, 
  #                                                          method = "SCM",
  #                                                          penalty = "L2",
  #                                                          solver = "cplex"))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0),
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1),
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::skip_if_not_installed("gurobi")
  # testthat::expect_warning(weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = list(penalty = 30), 
  #                                                          estimand = e, 
  #                                                          method = "SCM",
  #                                                          penalty = "variance",
  #                                                          solver = "gurobi")))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  
})

testthat::test_that("works for wass, joint mapping", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  # testthat::skip_if_not_installed("gurobi")
  # testthat::expect_warning(weights <- lapply(estimates, function(e) causalOT:::calc_weight_bal(data = data, 
  #                                                          constraint = list(joint = 1), 
  #                                                          estimand = e, 
  #                                                          method = "Wasserstein",
  #                                                          joint.mapping = TRUE,
  #                                                          penalty = "none",
  #                                                          solver = "gurobi")))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0),
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0),
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1),
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  weights <- lapply(estimates[1:3], function(e) causalOT:::calc_weight_bal(data = data, 
                                                           constraint = list(joint = 0.5,
                                                                             penalty = 10), 
                                                           estimand = e, 
                                                           joint.mapping = TRUE,
                                                           penalty = "L2",
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  testthat::expect_silent(weights[4] <- lapply(estimates[4], function(e) causalOT:::calc_weight_bal(data = data, 
                                                                constraint = list(joint = 0.5,
                                                                                  penalty = 10), 
                                                                estimand = e, 
                                                                joint.mapping = TRUE,
                                                                penalty = "L2",
                                                                method = "Wasserstein",
                                                                solver = "mosek")))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0),
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
})
