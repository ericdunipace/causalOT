arg.names <- c("w0",  "w1",   "gamma","estimand", "method",  "args")

testthat::test_that("works for Const Wass", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 3, 
                                                           estimand = e, 
                                                           p = power,
                                                           method = "Constrained Wasserstein",
                                                           solver = "gurobi"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  
  weights2 <- lapply(estimates, function(e) calc_weight(data = data, 
                                                            constraint = 1.85, 
                                                            estimand = e, 
                                                            p = power,
                                                            method = "Constrained Wasserstein",
                                                            solver = "gurobi"))
  sapply(weights2, function(w) testthat::expect_equal(names(w), arg.names))
  
  
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
  power <- c(2)
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
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 3, 
                                                       estimand = e, 
                                                       p = power,
                                                       method = "Constrained Wasserstein",
                                                       solver = "gurobi",
                                                       metric = distance,
                                                       rkhs.args = rkhs.argz))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  
  weights2 <- lapply(estimates, function(e) calc_weight(data = data, 
                                                        constraint = .001, 
                                                        estimand = e, 
                                                        p = power,
                                                        method = "Constrained Wasserstein",
                                                        solver = "gurobi",
                                                        metric = distance,
                                                        rkhs.args = rkhs.argz))
  sapply(weights2, function(w) testthat::expect_equal(names(w), arg.names))
  
  
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
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 2, 
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
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 2, 
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
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 2, 
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
#   weights <- lapply(estimates, function(e) calc_weight(data = data, 
#                                                        constraint = 3, 
#                                                        estimand = e, 
#                                                        method = "Wasserstein",
#                                                        solver = "gurobi",
#                                                        metric = distance,
#                                                        rkhs.args = rkhs.argz))
#   for(w in weights) testthat::expect_equal(names(w), arg.names)
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
#   weights <- lapply(estimates, function(e) calc_weight(data = data, 
#                                                        constraint = 3, 
#                                                        estimand = e, 
#                                                        method = "Wasserstein",
#                                                        solver = "cplex",
#                                                        metric = distance,
#                                                        rkhs.args = rkhs.argz))
#   for(w in weights) testthat::expect_equal(names(w), arg.names)
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
#   weights <- lapply(estimates, function(e) calc_weight(data = data, 
#                                                        constraint = 3, 
#                                                        estimand = e, 
#                                                        method = "Wasserstein",
#                                                        solver = "mosek",
#                                                        metric = distance,
#                                                        rkhs.args = rkhs.argz))
#   for(w in weights) testthat::expect_equal(names(w), arg.names)
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
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "gurobi"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "mosek"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "cplex"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  
})

testthat::test_that("works for SBW grid", {
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
  constraint <- seq(0,1, 10)
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = constraint, 
                                                       grid.search = TRUE,
                                                       grid = 10,
                                                       estimate = e, 
                                                       method = "SBW",
                                                       solver = "gurobi"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  sapply(weights, function(w) testthat::expect_lte(w$args$standardized.mean.difference , 10))
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = constraint,
                                                       grid.search = TRUE,
                                                       grid = 10,
                                                       estimate = e, 
                                                       method = "SBW",
                                                       solver = "mosek"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  sapply(weights, function(w) testthat::expect_lte(w$args$standardized.mean.difference , 10))
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = constraint,
                                                       grid.search = TRUE,
                                                       grid = 10,
                                                       estimate = e, 
                                                       method = "SBW",
                                                       solver = "cplex"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  sapply(weights, function(w) testthat::expect_lte(w$args$standardized.mean.difference , 10))
  
  
  
})