testthat::test_that("works for Const Wass", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
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
  sapply(weights, function(w) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma")))

  
  weights2 <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 1.85, 
                                                           estimand = e, 
                                                           method = "Constrained Wasserstein",
                                                           solver = "gurobi"))
  sapply(weights2, function(w) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma")))
  
  
  test_fun <- function(w1,w2){testthat::expect_match(all.equal(w1$w1, w2$w1), "Mean relative difference")}
  mapply(test_fun, w1 = weights[2:4], w2 = weights2[2:4])
  
  test_fun2 <- function(w1,w2){testthat::expect_match(all.equal(w1$w0, w2$w0), "Mean relative difference")}
  # test_fun2 <- function(w1,w2){all.equal(w1$w0, w2$w0)}
  
  mapply(test_fun2, w1 = weights[c(1,3:4)], w2 = weights2[c(1,3:4)])
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
})

testthat::test_that("works for Const Wass RKHS", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("RKHS")
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
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           distance = distance,
                                                           method = "Constrained Wasserstein",
                                                           solver = "gurobi",
                                                           rkhs.args = rkhs.argz)
                                                            )
  sapply(weights, function(w) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma")))
  
  
  weights2 <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                            constraint = .2, 
                                                            estimand = e, 
                                                            distance = distance,
                                                            rkhs.args = rkhs.argz,
                                                            method = "Constrained Wasserstein",
                                                            solver = "gurobi"))
  sapply(weights2, function(w) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma")))
  
  
  test_fun <- function(w1,w2){testthat::expect_match(all.equal(w1$w1, w2$w1), "Mean relative difference")}
  mapply(test_fun, w1 = weights[2:4], w2 = weights2[2:4])
  
  test_fun2 <- function(w1,w2){testthat::expect_match(all.equal(w1$w0, w2$w0), "Mean relative difference")}
  # test_fun2 <- function(w1,w2){all.equal(w1$w0, w2$w0)}
  
  mapply(test_fun2, w1 = weights[c(1,3:4)], w2 = weights2[c(1,3:4)])
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights2[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights2[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
})

testthat::test_that("works for Const Wass RKHS, opt", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("RKHS")
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
                                                           distance = distance,
                                                           method = "Constrained Wasserstein",
                                                           solver = "gurobi",
                                                           rkhs.args = rkhs.argz)
  )
  sapply(weights, function(w) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma")))
  
  
  weights2 <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                            constraint = 40, 
                                                            estimand = e, 
                                                            distance = distance,
                                                            rkhs.args = rkhs.argz,
                                                            method = "Constrained Wasserstein",
                                                            solver = "gurobi"))
  sapply(weights2, function(w) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma")))
  
  
  test_fun <- function(w1,w2){testthat::expect_match(all.equal(w1$w1, w2$w1, check.attributes = FALSE), "Mean relative difference")}
  mapply(test_fun, w1 = weights[2:4], w2 = weights2[2:4])
  
  test_fun2 <- function(w1,w2){testthat::expect_match(all.equal(w1$w0, w2$w0, check.attributes = FALSE), "Mean relative difference")}
  # test_fun2 <- function(w1,w2){all.equal(w1$w0, w2$w0)}
  mapply(test_fun2, w1 = weights[c(4)], w2 = weights2[c(4)])
  test_fun2 <- function(w1,w2){all.equal(w1$w0, w2$w0)}
  mapply(test_fun2, w1 = weights[c(1,3)], w2 = weights2[c(1,3)])
  
  test_fun2(list(w0 = rep(1/n0,n0)), 
                                   weights2[[1]])
  test_fun2(list(w0 = rep(1/n0,n0)), 
            weights2[[3]])
  test_fun2(list(w0 = rep(1/n0,n0)), 
            weights2[[4]])
  
  test_fun(list(w1 = rep(1/n1,n1)), 
            weights2[[2]])
  test_fun(list(w1 = rep(1/n1,n1)), 
           weights2[[3]])
  test_fun(list(w1 = rep(1/n1,n1)), 
           weights2[[4]])
  
})

testthat::test_that("works for Wass", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
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
  for(w in weights) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma"))
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
  for(w in weights) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma"))
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
  for(w in weights) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma"))
  
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

# testthat::test_that("works for Wass RKHS", {
#   set.seed(23483)
#   n <- 2^7
#   p <- 6
#   nsims <- 1
#   overlap <- "low"
#   design <- "A"
#   distance <- c("RKHS")
#   power <- c(1,2)
#   solver <- "gurobi"
#   estimates <- c("ATT", "ATC", "cATE", "ATE")
#   rkhs.argz <- list(p = 2,
#                     theta = c(1,2),
#                     gamma = c(10,5))
#   
#   #### get simulation functions ####
#   data <- causalOT::Hainmueller$new(n = n, p = p, 
#                                     design = design, overlap = overlap)
#   data$gen_data()
#   ns <- data$get_n()
#   n0 <- ns["n0"]
#   n1 <- ns["n1"]
#   
#   weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
#                                                            constraint = 30, 
#                                                            estimand = e, 
#                                                            rkhs.args = rkhs.argz,
#                                                            distance = distance,
#                                                            method = "Wasserstein",
#                                                            solver = "gurobi"))
#   for(w in weights) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma"))
#   testthat::expect_match(all.equal(rep(1/n0,n0), 
#                                    weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n0,n0), 
#                                    weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n0,n0), 
#                                    weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
#   
#   testthat::expect_match(all.equal(rep(1/n1,n1), 
#                                    weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n1,n1), 
#                                    weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n1,n1), 
#                                    weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
#   
#   weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
#                                                            constraint = 30, 
#                                                            estimand = e, 
#                                                            method = "Wasserstein",
#                                                            rkhs.args = rkhs.argz,
#                                                            distance = distance,
#                                                            solver = "cplex"))
#   for(w in weights) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma"))
#   testthat::expect_match(all.equal(rep(1/n0,n0), 
#                                    weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n0,n0), 
#                                    weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n0,n0), 
#                                    weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
#   
#   testthat::expect_match(all.equal(rep(1/n1,n1), 
#                                    weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n1,n1), 
#                                    weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n1,n1), 
#                                    weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
#   
#   weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
#                                                            constraint = 30, 
#                                                            estimand = e, 
#                                                            rkhs.args = rkhs.argz,
#                                                            distance = distance,
#                                                            method = "Wasserstein",
#                                                            solver = "mosek"))
#   for(w in weights) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma"))
#   
#   testthat::expect_match(all.equal(rep(1/n0,n0), 
#                                    weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n0,n0), 
#                                    weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n0,n0), 
#                                    weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
#   
#   testthat::expect_match(all.equal(rep(1/n1,n1), 
#                                    weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n1,n1), 
#                                    weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
#   testthat::expect_match(all.equal(rep(1/n1,n1), 
#                                    weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
# })

testthat::test_that("works for Wass RKHS, opt", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("RKHS")
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
                                                           constraint = 10, 
                                                           estimand = e, 
                                                           rkhs.args = rkhs.argz,
                                                           distance = distance,
                                                           method = "Wasserstein",
                                                           solver = "gurobi"))
  for(w in weights) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma"))
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
                                                           constraint = 10, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           rkhs.args = rkhs.argz,
                                                           distance = distance,
                                                           solver = "cplex"))
  for(w in weights) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma"))
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
                                                           constraint = 10, 
                                                           estimand = e, 
                                                           rkhs.args = rkhs.argz,
                                                           distance = distance,
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma"))
  
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

testthat::test_that("works for SBW", {
  set.seed(23483)
  n <- 2^8
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
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
  sapply(weights, function(w) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma")))
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "mosek"))
  sapply(weights, function(w) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma")))
  weights <- lapply(estimates, function(e) calc_weight_bal(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "cplex"))
  sapply(weights, function(w) testthat::expect_equal(names(w), c("w0",    "w1",    "gamma")))
  
  
})

# 
# library(causalOT)
# set.seed(23483)
# n <- 2^9
# p <- 6
# nsims <- 1
# overlap <- "low"
# design <- "B"
# distance <- c("Lp")
# power <- c(1,2)
# solver <- "gurobi"
# estimates <- c("ATT", "ATC","ATE")
# 
# #### get simulation functions ####
# data <- causalOT::Hainmueller$new(n = n, p = p,
#                                   design = design, overlap = overlap)
# data$gen_data()
# ns <- data$get_n()
# n0 <- ns["n0"]
# n1 <- ns["n1"]
# n <- sum(ns)
# #
# #
# cost  <- cost_mahalanobis(data$get_x0(), data$get_x1(), 2)
# cost0 <- cost_mahalanobis(data$get_x0(), data$get_x(), 1)
# cost1 <- cost_mahalanobis(data$get_x1(), data$get_x(), 1)
# #
# # # debugonce(calc_weight_bal)
# # debugonce(qp_wass)
# # debugonce(mosek_solver)
# test1 <- calc_weight(data = data,
#                          dist = "mahalanobis",
#                          p = 2,
#                         constraint = 0.5,
#                         estimand = "ATE",
#                         method = "Wasserstein",
#                         solver = "mosek")
# 
# # debugonce(mosek_solver)
# # debugonce(calc_weight_bal)
# test2 <- calc_weight(data = data,
#                          p =1,
#                         constraint = c(1),
#                         estimand = "ATE",
#                         method = "Constrained Wasserstein",
#                         dist = "mahalanobis",
#                         solver = "mosek")
# 
# # all.equal(test2$w0, test1$w0)
# # a <- rep(1/n,n)
# # b0 <- test2$w0
# # b1 <- test2$w1
# #
# # transport::wasserstein(a,b0,p = 1,costm= cost0)
# # transport::wasserstein(a,b1,p = 1,costm= cost1)
# # transport::wasserstein(a, rep(1/length(test1$w0), length(test1$w0)), p = 1,costm = cost0)
# # transport::wasserstein(a, rep(1/length(test1$w1), length(test1$w1)), p = 1,costm = cost1)
# #
# #
# # b0 <- test1$w0
# # b1 <- test1$w1
# #
# # transport::wasserstein(a,b0,p = 1,costm= cost0)
# # transport::wasserstein(a,b1, p = 1,costm= cost1)
# # transport::wasserstein(a, rep(1/length(test1$w0), length(test1$w0)), p = 1, costm = cost0)
# # transport::wasserstein(a, rep(1/length(test1$w1), length(test1$w1)), p = 1,costm = cost1)
# 
# # debugonce(mosek_solver)
# # debugonce(calc_weight_bal)
# test3 <- calc_weight(data, constraint = 0.1, estimand = "ATE", method = "SBW", solver = "mosek")
# test4 <- calc_weight(data, constraint = 0.1, estimand =  "ATE", method = "Logistic", solver = "mosek")
# test4b <- calc_weight(data, constraint = 0.1, estimand =  "cATE", method = "Logistic", solver = "mosek")
# 
# test5 <- calc_weight(data, constraint = 0.1, estimand = "ATE", method = "RKHS", solver = "gurobi")
# test5b <- calc_weight(data, constraint = 0.1, estimand = "cATE", method = "RKHS", solver = "gurobi")
# 
# 
# 
# # test6_temp <- list(calc_weight_bal(data = data,
# #                          p =1,
# #                          constraint = c(2.2),
# #                          estimand = "ATC",
# #                          method = "Constrained Wasserstein",
# #                          dist = "mahalanobis",
# #                          solver = "mosek"),
# #                    calc_weight_bal(data = data,
# #                                    p =1,
# #                                    constraint = c(2.2),
# #                                    estimand = "ATT",
# #                                    method = "Constrained Wasserstein",
# #                                    dist = "mahalanobis",
# #                                    solver = "mosek")
# #                     )
# test6 <- calc_weight(data = data,
#                          p =1,
#                          constraint = c(2.2),
#                          estimand = "cATE",
#                          method = "Constrained Wasserstein",
#                          dist = "mahalanobis",
#                          solver = "mosek")
# test7 <- calc_weight(data = data,
#                          p =1,
#                          constraint = c(1.2),
#                          estimand = "cATE",
#                          method = "Wasserstein",
#                          dist = "mahalanobis",
#                          solver = "mosek")
# 
# test8 <- calc_weight(data = data,
#                          p = 1,
#                          constraint = .1,
#                          estimand = "ATE",
#                          method = "Constrained Wasserstein",
#                          dist = "RKHS",
#                          solver="mosek",
#                          rkhs.args = list(p = test5$addl.args$p,
#                                           theta = test5$addl.args$theta,
#                                           gamma = test5$addl.args$gamma))
# 
# outcome_model(data, weights = test1, target = "ATE", matched = TRUE, doubly.robust = TRUE)
# outcome_model(data, weights = test2, target = "ATE", matched = TRUE, doubly.robust = TRUE)
# outcome_model(data, weights = test3, target = "ATE", matched = TRUE, doubly.robust = TRUE)
# outcome_model(data, weights = test4, target = "ATE", matched = TRUE, doubly.robust = TRUE)
# outcome_model(data, weights = test5, target = "ATE", matched = TRUE, doubly.robust = TRUE)
# outcome_model(data, weights = test5b, target = "ATE", matched = TRUE, doubly.robust = TRUE)
# outcome_model(data, weights = test6, target = "ATE", matched = TRUE, doubly.robust = TRUE)
# outcome_model(data, weights = test7, target = "ATE", matched = TRUE, doubly.robust = TRUE)
# outcome_model(data, weights = test8, target = "ATE", matched = TRUE, doubly.robust = TRUE)
# 
# library(causalOT)
# options(mc.cores = 3)
# ess <- lapply(weights, function(w) ESS(w)/ns)
# weights <- list(W = test1, CW = test2, SBW=test3, IPW = test4, RKHS = test5,
#                 cRKHS = test5b, cCW = test6, cW = test7, rCW = test8)
# ps <- mapply(function(w, r) {PSIS(w, r_eff = r)}, w = weights, r = ess, SIMPLIFY = FALSE)
# ps.check <- PSIS(weights, r_eff = ess)
# testthat::expect_equivalent(ps, ps.check)
# 
# diag.check <- lapply(ps, PSIS_diag)
# diag <- PSIS_diag(ps.check)
# diag.check2 <- PSIS_diag(weights)
# 
# 
# sapply(diag.check2, function(d) c(d$w0$n_eff, d$w1$n_eff))
# sapply(diag.check2, function(d) c(d$w0$pareto_k, d$w1$pareto_k))
# 
# 
# sapply(diag, function(d) c(d$w0$n_eff, d$w1$n_eff))
# sapply(diag, function(d) c(d$w0$pareto_k, d$w1$pareto_k))
# sapply(weights, function(w) ESS(w))
