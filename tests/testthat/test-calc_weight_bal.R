arg.names <- c("w0",  "w1",   "gamma","estimand", "method",  "args")

testthat::test_that("works for Const Wass", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           method = "Constrained Wasserstein",
                  solver = "gurobi"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))

  
  weights2 <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 1.85, 
                                                           estimand = e, 
                                                           method = "Constrained Wasserstein",
                                                           solver = "gurobi"))
  sapply(weights2, function(w) testthat::expect_equal(names(w), arg.names))
  
  
  test_fun <- function(w1,w2){testthat::expect_match(all.equal(w1$w1, w2$w1), "Mean relative difference")}
  mapply(test_fun, w1 = weights[c(2,4)], w2 = weights2[c(2,4)])
  
  test_fun2 <- function(w1,w2){testthat::expect_match(all.equal(w1$w0, w2$w0), "Mean relative difference")}
  # test_fun2 <- function(w1,w2){all.equal(w1$w0, w2$w0)}
  
  mapply(test_fun2, w1 = weights[c(1,4)], w2 = weights2[c(1,4)])
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights2[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights2[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
})

testthat::test_that("works for Const Wass RKHS", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("RKHS")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  rkhs.argz <- list(p = 2,
                    theta = c(1,2),
                    gamma = c(10,5))
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 30, 
                                                           estimand = e, 
                                                           metric = metric,
                                                           method = "Constrained Wasserstein",
                                                           solver = "gurobi",
                                                           rkhs.args = rkhs.argz)
                                                            )
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  
  weights2 <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                            constraint = 0.001, 
                                                            estimand = e, 
                                                            metric = metric,
                                                            rkhs.args = rkhs.argz,
                                                            method = "Constrained Wasserstein",
                                                            solver = "gurobi"))
  sapply(weights2, function(w) testthat::expect_equal(names(w), arg.names))
  
  
  test_fun <- function(w1,w2){testthat::expect_match(all.equal(w1$w1, w2$w1), "Mean relative difference")}
  mapply(test_fun, w1 = weights[c(2,4)], w2 = weights2[c(2,4)])
  
  test_fun2 <- function(w1,w2){testthat::expect_match(all.equal(w1$w0, w2$w0), "Mean relative difference")}
  # test_fun2 <- function(w1,w2){all.equal(w1$w0, w2$w0)}
  
  mapply(test_fun2, w1 = weights[c(1,4)], w2 = weights2[c(1,4)])
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0),
  #                                  weights2[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights2[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
})

testthat::test_that("works for Const Wass RKHS, opt", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("RKHS")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  rkhs.argz <- NULL
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 300, 
                                                           estimand = e, 
                                                           metric = metric,
                                                           method = "Constrained Wasserstein",
                                                           solver = "gurobi",
                                                           rkhs.args = rkhs.argz)
  )
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  
  weights2 <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                            constraint = 40, 
                                                            estimand = e, 
                                                            metric = metric,
                                                            rkhs.args = rkhs.argz,
                                                            method = "Constrained Wasserstein",
                                                            solver = "gurobi"))
  sapply(weights2, function(w) testthat::expect_equal(names(w), arg.names))
  
  
  test_fun <- function(w1,w2){testthat::expect_match(all.equal(w1$w1, w2$w1, check.attributes = FALSE), "Mean relative difference")}
  test_fun1 <- function(w1,w2){testthat::expect_true(all.equal(w1$w1, w2$w1, check.attributes = FALSE))}
  
  mapply(test_fun1, w1 = weights[c(1,3)], w2 = weights2[c(1,3)])
  
  test_fun2 <- function(w1,w2){testthat::expect_match(all.equal(w1$w0, w2$w0, check.attributes = FALSE), "Mean relative difference")}
  # test_fun2 <- function(w1,w2){all.equal(w1$w0, w2$w0)}
  test_fun2 <- function(w1,w2){all.equal(w1$w0, w2$w0, check.attributes = FALSE)}
  mapply(test_fun2, w1 = weights[c(4)], w2 = weights2[c(4)])
  
  mapply(test_fun2, w1 = weights[c(1,3)], w2 = weights2[c(1,3)])
  
  test_fun2(list(w0 = rep(1/n0,n0)), 
                                   weights2[[1]])
  test_fun2(list(w0 = rep(1/n0,n0)), 
            weights2[[3]])
  test_fun2(list(w0 = rep(1/n0,n0)), 
            weights2[[4]])
  
  test_fun1(list(w1 = rep(1/n1,n1)),
            weights2[[1]])
  test_fun1(list(w1 = rep(1/n1,n1)), 
           weights2[[3]])
  # test_fun1(list(w1 = rep(1/n1,n1)), 
  #          weights2[[4]])
  
})

testthat::test_that("works for Wass", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 2.5, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = "gurobi"))
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = "cplex"))
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("RKHS")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  rkhs.argz <- NULL
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 1e6, 
                                                           estimand = e, 
                                                           rkhs.args = rkhs.argz,
                                                           metric = metric,
                                                           method = "Wasserstein",
                                                           solver = "gurobi"))
  for (w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE,
                         tol = 1e-3)
  testthat::expect_equal(rep(1/n0,n0), 
                         weights[[3]]$w0, check.attributes = FALSE,
                         tol = 1e-3)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 1e6, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           rkhs.args = rkhs.argz,
                                                           metric = metric,
                                                           solver = "cplex"))
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "gurobi"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "mosek"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "cplex"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  
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
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  w0 <- renormalize(runif(n0))
  w1 <- renormalize(runif(n1))
  w0[seq(1,n0,2)] <- 0
  w1[seq(1,n1,2)] <- 0
  sample_weights <- list(w0 = renormalize(w0),
                         w1 = renormalize(w1))
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 2.5, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = "gurobi",
                                                           sample_weight = sample_weights))
  for (w in weights) testthat::expect_equal(names(w), arg.names)
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = "cplex",
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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

testthat::test_that("works for Const Wass, sample weight", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  w0 <- renormalize(runif(n0))
  w1 <- renormalize(runif(n1))
  w0[seq(1,n0,2)] <- 0
  w1[seq(1,n1,2)] <- 0
  sample_weights <- list(w0 = renormalize(w0),
                         w1 = renormalize(w1))
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           method = "Constrained Wasserstein",
                                                           solver = "gurobi",
                                                           sample_weight = sample_weights))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  
  weights2 <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                            constraint = 1.85, 
                                                            estimand = e, 
                                                            method = "Constrained Wasserstein",
                                                            solver = "gurobi",
                                                            sample_weight = sample_weights))
  sapply(weights2, function(w) testthat::expect_equal(names(w), arg.names))
  
  
  test_fun <- function(w1,w2){testthat::expect_match(all.equal(w1$w1, w2$w1), "Mean relative difference")}
  mapply(test_fun, w1 = weights[c(2,4)], w2 = weights2[c(2,4)])
  
  test_fun2 <- function(w1,w2){testthat::expect_match(all.equal(w1$w0, w2$w0), "Mean relative difference")}
  # test_fun2 <- function(w1,w2){all.equal(w1$w0, w2$w0)}
  
  mapply(test_fun2, w1 = weights[c(1,4)], w2 = weights2[c(1,4)])
  
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights2[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights2[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights2[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights2[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights2[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights2[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
})

testthat::test_that("works for Wass, variance", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 30, 
                                                           estimand = e, 
                                                           penalty = "variance",
                                                           method = "Wasserstein",
                                                           solver = "cplex"))
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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

testthat::test_that("works for const Wass, variance", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 100, 
                                                           estimand = e, 
                                                           penalty = "variance",
                                                           method = "Constrained Wasserstein",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0),
  #                                  weights[[2]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_warning(weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 1, 
                                                           estimand = e, 
                                                           penalty = "variance",
                                                           method = "Constrained Wasserstein",
                                                           solver = "cplex")))
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 10, 
                                                           estimand = e, 
                                                           penalty = "variance",
                                                           method = "Constrained Wasserstein",
                                                           solver = "mosek"))
  for (w in weights) testthat::expect_equal(names(w), arg.names)
  
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for Const Wass, entropy", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 100, 
                                                           estimand = e, 
                                                           penalty = "entropy",
                                                           method = "Constrained Wasserstein",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 30, 
                                                           estimand = e, 
                                                           method = "Constrained Wasserstein",
                                                           penalty = "entropy",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 10, 
                                                           estimand = e, 
                                                           method = "Constrained Wasserstein",
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
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = NULL, 
                                                           estimand = e, 
                                                           method = "SCM",
                                                           penalty = "none",
                                                           solver = "gurobi"))
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = "cplex"))
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
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
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  testthat::expect_warning(weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = list(penalty = 3), 
                                                           estimand = e, 
                                                           method = "SCM",
                                                           penalty = "L2",
                                                           solver = "cplex"))
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = list(penalty = 30), 
                                                           estimand = e, 
                                                           method = "SCM",
                                                           penalty = "variance",
                                                           solver = "gurobi"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
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
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = list(joint = 1), 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           joint.mapping = TRUE,
                                                           penalty = "none",
                                                           solver = "gurobi"))
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
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = list(joint = 1, penalty = 1), 
                                                           estimand = e, 
                                                           penalty = "variance",
                                                           method = "Wasserstein",
                                                           joint.mapping = TRUE,
                                                           solver = "cplex"))
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
  
  weights <- lapply(estimates[1:3], function(e) calc_weight_bal(data = data, 
                                                           constraint = list(joint = 0.5,
                                                                             penalty = 10), 
                                                           estimand = e, 
                                                           joint.mapping = TRUE,
                                                           penalty = "L2",
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  testthat::expect_warning(weights[4] <- lapply(estimates[4], function(e) calc_weight_bal(data = data, 
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

testthat::test_that("works for Wass, neg.weights", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 100, 
                                                           estimand = e, 
                                                           penalty = "L2",
                                                           method = "Wasserstein",
                                                           solver = "mosek",
                                                           neg.weights = TRUE))
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
  testthat::expect_true(any(weights[[4]]$w1 < 0))
  testthat::expect_true(any(weights[[4]]$w0 < 0))
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 30, 
                                                           estimand = e, 
                                                           penalty = "L2",
                                                           method = "Wasserstein",
                                                           solver = "cplex",
                                                           neg.weights = TRUE))
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
  testthat::expect_true(any(weights[[4]]$w1 < 0))
  testthat::expect_true(any(weights[[4]]$w0 < 0))
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 1e-6, 
                                                           estimand = e, 
                                                           penalty = "variance",
                                                           method = "Wasserstein",
                                                           solver = "mosek",
                                                           neg.weights = TRUE))
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
  testthat::expect_true(any(weights[[4]]$w1 < 0))
  testthat::expect_true(any(weights[[4]]$w0 < 0))
})

testthat::test_that("works for CWass, neg.weights", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 1, 
                                                           estimand = e, 
                                                           penalty = "L2",
                                                           method = "Constrained Wasserstein",
                                                           solver = "mosek",
                                                           neg.weights = TRUE))
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
  testthat::expect_true(any(weights[[4]]$w1 < 0))
  testthat::expect_true(any(weights[[4]]$w0 < 0))
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 1, 
                                                           estimand = e, 
                                                           penalty = "L2",
                                                           method = "Constrained Wasserstein",
                                                           solver = "cplex",
                                                           neg.weights = TRUE))
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
  testthat::expect_true(any(weights[[4]]$w1 < 0))
  testthat::expect_true(any(weights[[4]]$w0 < 0))
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = list(joint = 0, penalty = 1e-100), 
                                                           estimand = e, 
                                                           penalty = "variance",
                                                           method = "Constrained Wasserstein",
                                                           solver = "mosek",
                                                           neg.weights = TRUE))
  for (w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0),
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   # weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1),
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1),
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_true(any(weights[[2]]$w1 < 0))
  testthat::expect_true(any(weights[[1]]$w0 < 0))
})

testthat::test_that("works for SCM, neg.weights", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = list(penalty = 100), 
                                                           estimand = e, 
                                                           penalty = "L2",
                                                           method = "SCM",
                                                           solver = "mosek",
                                                           neg.weights = TRUE))
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
  # testthat::expect_true(any(weights[[4]]$w1 < 0))
  # testthat::expect_true(any(weights[[4]]$w0 < 0))
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = list(penalty = 30), 
                                                           estimand = e, 
                                                           penalty = "L2",
                                                           method = "SCM",
                                                           solver = "cplex",
                                                           neg.weights = TRUE))
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
  testthat::expect_true(any(weights[[4]]$w1 < 0))
  testthat::expect_true(any(weights[[4]]$w0 < 0))
  
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = list(penalty = 1), 
                                                           estimand = e, 
                                                           penalty = "SCM",
                                                           method = "SCM",
                                                           solver = "gurobi",
                                                           neg.weights = TRUE))
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
  # testthat::expect_true(any(weights[[4]]$w1 < 0))
  # testthat::expect_true(any(weights[[4]]$w0 < 0))
})