test_that("RKHS grid works", {
  set.seed(9870)
  library(causalOT)
  
  n <- 2^9
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  ground_power <- 2
  solver <- "gurobi"
  estimates <- c("ATT", "ATC","feasible")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  
  testthat::expect_silent(cplex.check <- RKHS_grid_search(data = data, grid = NULL, estimate = "ATE", n.boot = 100, solver = "cplex"))
  testthat::expect_silent(gurobi.check <- RKHS_grid_search(data = data, grid = NULL, estimate = "ATE", n.boot = 100, solver = solver))
})

test_that("RKHS grid gives expected value", {
  set.seed(9870)
  library(causalOT)
  
  n <- 2^9
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  solver <- "gurobi"
  
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  weight <- RKHS_grid_search(data, grid = seq(0,100, length.out = 11), 
                            n.boot = 100,
                            solver = solver,
                            theta = c(1,1),
                            gamma = c(1,1),
                            p = 1,
                            dist = distance)
  testthat::expect_equal(weight[1:3], do.call("calc_weight_RKHS", 
                                              list(data = data, 
                                                   method = "RKHS",
                                                   estimate = "ATE",
                                                   solver = solver,
                                                   lambda = weight$lambda,
                                                   theta = c(1,1),
                                                   gamma = c(1,1),
                                                   p = 1,
                                                   dist = distance)))
})